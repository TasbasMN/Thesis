#!/usr/bin/env python

"""
Last update: 02/15/2017
author: jun ding
usage: python TarPmiR.py -a<miRNA> -b<mRNA> -m<model> -p <probability_cut>
"""
import pdb,sys,os,re,math,getopt
import urllib,urllib2
import subprocess
import pickle 
import uuid
import pymysql
import ssl
import numpy as np
try:
	ssl._create_default_https_context = ssl._create_unverified_context
except:
	pass 

#------------------------------------------------
try:
	db1=pymysql.connect("genome-mysql.cse.ucsc.edu",'genome')
	cur1=db1.cursor()
	cur1.execute("use hg19")
except:
	pass

#---------------------------------------------------------------------

#***********************************************************************
#=======================================================================
# global functions#

def filterBS(A):
	def mergeBS(a,b):
		a_bs=a[2].split(',')
		b_bs=b[2].split(',')
		if float(a[-1])>float(b[-1]):
			return a
		elif float(a[-1])<float(b[-1]):
			return b
		else:
			if int(a_bs[1])-int(a_bs[0])>int(b_bs[1])-int(b_bs[0]):
				return a
			else:
				return b
			
	def overlap(a,b):
		a=a.split(',')
		b=b.split(',')
		flag=(int(a[1])-int(b[0]))*(int(a[0])-int(b[1]))
		if flag<=0:
			return True
		else:
			return False

	for i in range(len(A)-1):
		for j in range(i+1,len(A)):
			if A[i]!=[] and A[j]!=[]:
				if overlap(A[i][2],A[j][2]):
					A[i]=mergeBS(A[i],A[j])
					A[j]=[]
	A=[item for item in A if item!=[]]
	return A	

def get_conservation(pos):
	url="https://genome.ucsc.edu/cgi-bin/hgTables?"
	form_data={'hgsid':'385332177_fMwuaEoYIaAfIxVs6VEa6hgHN1kA',
				'clade':'mammal',
				'org':'human',
				'db':'hg19',
				'hgta_group':'compGeno',
				'hgta_track':'cons46way',
				'hgta_table':'phyloP46wayAll',
				'hgta_regionType':'range',
				'position':pos,
				'hgta_outputType':'wigData',
				'hgta_outFileName':'1',
				'hgta_doTopSubmit':'get output'}
	params=urllib.urlencode(form_data)
	response=urllib2.urlopen(url,params)
	data=response.read()
	data=data.split('\n')
	data=data[1:]
	out=[]
	for i in data:
		if (len(i)>0) and (i[0]!="#")and (i[0]!='v'):
			ii=i.split()
			ii=[int(ii[0]),float(ii[1])]
			out.append(ii)
	return out



def scan(seq,m,cut):
	dN={"A":0,"C":1,"G":2,"T":3}
	maxS=maxScore(m)
	loc=[]
	for i in range(len(seq)-len(m)+1):
		sc_i=0
		for j in range(len(m)):
			seqij=seq[i+j]
			if seqij!="N":
				sc=m[j][dN[seqij]]
			else:
				sc=0.0
			sc_i+=sc
		if sc_i/maxS>=cut:
			loc.append(i)
	return loc

	def maxScore(m):
		s=0
		for i in m:
			s+=max(i)
		return s

	
def get_ensGene(name):
	pn=re.compile('N._')
	row=None
	if name[0:4]=='ENST':
		cur1.execute('select * from hg19.ensGene where name="%s"'%(name))
		row=cur1.fetchone()
	elif pn.search(name):
		row=cur1.execute('select * from hg19.refGene where name="%s"'%(name))
	if row!=None:
		row=[row[2],row[3],row[4],row[5],row[9],row[10]]
	else:
		row=[]
	return row
	
def revcom(A):
	dN={"A":"T","C":"G","G":"C","T":"A","N":"N"}
	A=[dN[item] for item in A]
	A=A[::-1]
	A="".join(A)
	return A
	
def get_seq(pos,strand):
	link="http://genome.ucsc.edu/cgi-bin/das/hg19/dna?segment="+pos
	f=urllib.urlopen(link)
	lf=f.read()
	lf=lf.split("<DNA")[1].split("</DNA>")[0].split('>')[1]
	lf=lf.split("\n")
	lf=[item for item in lf if item!='']
	lf="".join(lf)
	lf=lf.upper()
	if strand=="+":
		seq=lf
	else:
		seq=revcom(lf)
	return seq

def check_seed(seedx,x):
	dN={'A':'T','C':'G','G':'CT','U':'AG','T':'AG','N':'N'}
	seedx=[dN[item] for item in seedx]
	rc_seedx=seedx[::-1]
	flag=1
	for i in range(len(x)):
		if x[i] not in rc_seedx[i]:
			flag=0
			break 
	return flag
		
def progressbar(count, total, suffix=""):
    bar_len = 60
    filled_len = int(round(bar_len * count / float(total)))

    percents = round(100.0 * count / float(total), 1)
    bar = '=' * filled_len + '-' * (bar_len - filled_len)

    sys.stdout.write('[%s] %s%s ...%s\r' % (bar, percents, '%', suffix))
    sys.stdout.flush()
    
def longestStem(bb,bm):
	#bb is the pairing output of miRNA from RNAduplex
	ct=0
	A=[]
	B=[]
	pc=0
	key='('
	for i in bb:
		if i==key:
			ct+=1
		else:
			A.append([ct,pc-ct])
			ct=0
		pc+=1
	A.append([ct,pc-ct])	
	A=max(A)
	
	
	pa=[] # pairing location in miRNA
	pb=[] # pairing location in mRNA
	for i in range(len(bb)):
		if bb[i]=='(':
			pa.append(i)

	for i in range(len(bm)-1,-1,-1):
		if bm[i]==')':
			pb.append(i)
	
	ss=bb[0:A[1]].count('(') # stem start
	se=bb[0:A[1]+A[0]-1].count('(')# stem end
	#pdb.set_trace()	
	return [pb[se],pb[ss]]
		
#=======================================================================
#***********************************************************************



#=======================================================================
# File processing

class File:
	def __init__(self,file_path):
		self.path=file_path
		
	def readInfo(self,sep):
		f=open(self.path,'r')
		lf=f.readlines()
		f.close()
		lf=[item.strip().split(sep) for item in lf]
		return lf
	
	def readSeq(self):
		f=open(self.path,'r')
		lf=f.readlines()
		f.close()
		lf="".join(lf)
		lf=lf.split(">")[1:]
		for i in range(len(lf)):
			lf[i]=lf[i].split('\n')
			lf[i]=["".join(lf[i][0]),''.join(lf[i][1:])]
		return lf


	
class miRNA:
	def __init__(self,mi,mi_seq):
		self.name=mi
		self.seq=mi_seq
	
class mRNA:		
	def __init__(self,m,seq):
		self.name=m
		self.seq=seq
		self.acc=self.accessibility_plfold() 
		self.phy=self.phyloP()
		
	def accessibility_plfold(self):
			W=80
			L=40
			sessionID='fn_'+str(uuid.uuid1())
			input_mseq='>'+sessionID+'\n'+self.seq
			pp=subprocess.Popen(["Vienna\RNAplfold","-W", "80", "-L", "40","-u", "16"],stdout=subprocess.PIPE, stdin=subprocess.PIPE,stderr=subprocess.PIPE)
			#pdb.set_trace()
			out,err=pp.communicate(input=input_mseq)
			pp.wait()
			#pp.terminate()
			f=open("%s_lunp"%(sessionID),'r')
			lf=f.readlines()
			f.close()
			lf=[item.strip() for item in lf]
			lf=[item for item in lf if len(item)!=0]
			lf=lf[2:]
			lf=[item.split('\t') for item in lf]
			lf=[item[-1] for item in lf]
			
			#===
			try:
				os.remove("%s_dp.ps"%(sessionID))
				os.remove("%s_lunp"%(sessionID))
			except:
				pass 
			
			return lf
				
	def phyloP(self):
		try:
			row=get_ensGene(self.name)
			span=80000
			currentS=row[2]
			currentE=currentS+span
			cvm=[]
			while(currentE<row[3]):
				gpos=str(row[0])+':'+str(currentS)+'-'+str(currentE)
				cvm+=get_conservation(gpos)
				currentS=currentE
				currentE=currentS+span
			gpos=str(row[0])+':'+str(currentS)+'-'+str(row[3])
			cvm+=get_conservation(gpos)
			dc={item[0]:item[1] for item in cvm}
			cv=[]
			if len(row)>0:
				exons=row[4].split(',')
				exons=[item for item in exons if item!='']
				exone=row[5].split(',')
				exone=[item for item in exone if item!='']
				for i in range(len(exons)):
					for j in range(int(exons[i]),int(exone[i])):
						if j in dc:
							cv.append(dc[j])
						else:
							cv.append(0)
			return cv	
		except:
			return []
		

class Interactions:
	def __init__(self,mir,m):
		self.mir=mir
		self.m=m
		
	def cand_site(self):
		ss=self.seed_site()
		es=self.energy_site()
		return ss+es
		
	def seed_site(self):
		seed=self.mir.seq[1:16]
		ss=[]
		for i in range(len(self.m.seq)-len(self.mir.seq)+1):
			x=self.m.seq[i:i+len(seed)]
			flag=check_seed(seed,x)
			if flag!=0:
				ss.append(i+len(seed))
		out=[]
		x=24
		for i in ss:
			out.append(str(max(0,i-x))+","+str(i+1))
		return out
				
	def energy_site(self):
		S=4
		C=-15
		ip_seq=">"+self.mir.name+"\n"+self.mir.seq+"\n"+">"+self.m.name+"\n"+self.m.seq
		pp1=subprocess.Popen(["Vienna\RNAduplex","-e", "5", "-s"],stdout=subprocess.PIPE, stdin=subprocess.PIPE,stderr=subprocess.PIPE)
		out,err=pp1.communicate(input=ip_seq)
		pp1.wait()
		#pp1.terminate()
		lf=[item for item in out.split('\n') if len(item)>0]
		lf=[item.strip() for item in lf if item[0]!=">"]
		lf=[item.split() for item in lf]
		lf=[item for item in lf if float(item[-1][1:-1])<C]
		
		lf=lf[:S]
		lf=[item[-2] for item in lf]
		return lf



class binding:
	def __init__(self,mir,m,pos):
		self.mir=mir
		self.m=m
		self.p=pos
		# features of each miRNA-mRNA binding site
		[self.en,self.sp,self.bb,self.bm]=self.mfe()
		self.seed=self.hasSeed()
		self.au=self.AU_content()
		[self.nc,self.pnc]=self.consecutive_pairs()
		self.me=self.me_motif()
		self.ac=self.acc()
		self.np=self.nbp()
		self.bl=self.bl()
		[self.pe,self.dpse]=self.prThreeEnd()
		[self.phys,self.phyf]=self.cv()
		
	def getFeatures(self):
		F=[self.en,self.seed,self.ac,self.au,self.phys,self.phyf,self.me,self.np,self.bl,self.nc,self.pnc,self.pe,self.dpse]
		return F
		
	def mfe(self):
		# RNAduplex
		t=self.p.split(",")
		t=[int(item) for item in t]
		seqt=self.m.seq[t[0]:t[1]]
		ipp_seq=">"+self.mir.name+"\n"+self.mir.seq+"\n"+">"+self.p+"\n"+seqt
		pp2=subprocess.Popen(["Vienna\RNAduplex","-e", "5", "-s"],stdout=subprocess.PIPE, stdin=subprocess.PIPE,stderr=subprocess.PIPE)
		out,err=pp2.communicate(input=ipp_seq)
		pp2.wait()
		#pp2.terminate()
		lf=[item for item in out.split('\n') if len(item)>0]
		
		lf=[item.strip() for item in lf if item[0]!=">"]
		#pdb.set_trace()
		lf=lf[0].split()
		pm=lf[3].split(',')
		pm=[int(item) for item in pm]
		pmi=lf[1].split(',')
		pmi=[int(item) for item in pmi]
		
		en=float(lf[-1][1:-1])
		bb=lf[0].split('&')[0]
		bm=lf[0].split('&')[1]
		###########################
		sd=pmi[0]-2
		msd=pm[1]+sd
		#================================================
		bb="."*(pmi[0]-1)+bb+'.'*(len(self.mir.seq)-pmi[1])
		sp=bb[1:7].count('(')
		sp1=bb[2:8].count('(')
		sp=max(sp,sp1)
		#================================================
		pm=[item-1+t[0] for item in pm]
	
		#en: energy
		#sp: # of pairings in the seed region
		#bb: pairing details in miRNA
		#bm
		
		return [en,sp,bb,bm]
	
	def hasSeed(self):
		return self.sp/6
		
	def AU_content(self):
		t=self.p.split(',')
		offset=30 # 30 nt
		x=8 # seed length
		t=[int(item) for item in t]
		seed_pos=[max(0,t[1]-1-x-offset),t[1]-1+offset]
		#t=[max(0,t[0]-offset),t[1]+offset]
		A=[0,0,0,0]
		seqt=self.m.seq[seed_pos[0]:seed_pos[1]]
		dN={"A":0,"C":1,"G":2,"T":3}
		for i in seqt:
			if i in dN:
				A[dN[i]]+=1
		au=round(float(A[0]+A[3])/sum(A),3)
		return au
		
	def consecutive_pairs(self):
		#bm is the pairing of mRNA output from RNAduplex
		bm=self.bm
		ct=0
		qt=0
		am=2
		A=[]
		B=[]
		pc=0
		for i in bm:
			if i==")":
				ct+=1
				qt+=1
			else:
				if am<0:
					#pdb.set_trace()
					B.append([qt,len(bm)-pc])
					qt=0
				am=am-1
				A.append([ct,len(bm)-pc])
				ct=0
			pc+=1
		A.append([ct,len(bm)-pc])
		A=max(A)
		return A
		
	def me_motif(self):
		# bb : pairing of miRNA
		bb=self.bb
		dBP={'(':'m','.':'e'}
		bp=[dBP[item] for item in bb]
		bp=''.join(bp)
		# Averge probability of pairing in each position learned from CLASH
		A=[0.464,0.762,0.826,0.954,0.984,0.87,0.824,0.91,0.742,0.776,0.944,0.966,0.784,0.702,0.676,0.558,0.554,0.456,0.562,0.552,0.422,0.096]
		dN={'m':0,'e':1}
		x=16
		p0=0
		for j in range(x):
			ind=dN[bp[j]]
			if ind==0:
				p=A[j]
			else:
				p=1-A[j]
			p0+=math.log(p)
		return p0	
		
	def acc(self):
		t=self.p.split(',')
		t=[int(item) for item in t]
		return float(self.m.acc[min(max(15,t[1]-1),len(self.m.seq)-1)])
		
	def nbp(self):
		return self.bb.count('(')
	
	def bl(self):
		return(len(self.bm))
	
	def prThreeEnd(self):
		pe=self.bb[-8:].count('(')
		ps=self.bb[0:8].count('(')
		dpse=abs(pe-ps)
		return [pe,dpse]
		
	def cv(self):
		t=self.p.split(',')
		t=[int(item) for item in t]
		lstem=longestStem(self.bb,self.bm)
		lstem=[item+t[0] for item in lstem]
		try:
			phys=self.m.phy[lstem[0]:lstem[1]+1]
			phyf1=self.m.phy[max(0,lstem[0]-50):max(0,lstem[0]-10)]
			phyf2=self.m.phy[lstem[1]+10:lstem[1]+50]
			Avg_phys=sum(phys)/len(phys)
			Avg_phyf1=sum(phyf1)/len(phyf1)
			Avg_phyf2=sum(phyf2)/len(phyf2)
			Avg_phyf=(Avg_phyf1+Avg_phyf2)/2
			return [Avg_phys, Avg_phyf]
		except:
			return [0,0]
			

		
#=======================================================================


#=======================================================================
def main():
	
	# parsing parameters
	opts, args = getopt.getopt(sys.argv[1:], 'a:b:m:p:')
	
	if len(opts)!=4:
		print("please check your inputs!")
		print("python TarPmiR.py -a <mir> -b <mRNA> -m <model> -p <probability_cut>")
		sys.exit(0)
		
	for i in opts:
		if i[0]=='-a':
			mir_path=i[1]
		if i[0]=='-b':
			mrna_path=i[1]
		if i[0]=='-m':
			mode_path=i[1]
		if i[0]=='-p':
			pb_cut=float(i[1])
			
	global sessionID 
	sessionID= 'fn_'+str(uuid.uuid1())
	#=======================================================================
	
	# read in microRNAs.txt
	mir=File(mir_path)
	mir_info=mir.readSeq()

	# read in mRNA
	mRNA_file=File(mrna_path)
	mRNA_seq=mRNA_file.readSeq()

	# read in model
	f=open(mode_path,'r')
	RR=pickle.load(f)
	#=======================================================================

	#pb_cut=0.5 # confidential cut of binding probability

	out=[]
	ct=0

	N=len(mRNA_seq)
	mrna_name=mrna_path.split('/')[-1]
	mrna_name1=mrna_path.split('\\')[-1]
	if len(mrna_name1)<len(mrna_name):
		mrna_name=mrna_name1
	
	if os.path.exists(mir_path+'_'+mrna_name+'.bp'):
		os.remove(mir_path+'_'+mrna_name+'.bp')
	
	f=open(mir_path+'_'+mrna_name+'.bp','a')
	
	for i in mRNA_seq:
		# instance of mRNA class 
		m=mRNA(i[0],i[-1])
		out=[]
		for j in mir_info:
			mir=miRNA(j[0],j[-1])
			ii=Interactions(mir,m)
			bs=ii.cand_site()
			FS=[]
			out_bs=[]
			for k in bs:
				bijk=binding(mir,m,k)
				F=bijk.getFeatures()
				FF=np.array(F).reshape(1,-1)
				pb=RR.predict_proba(FF)[0][1]
				if pb>pb_cut:
					out_bs.append([bijk.mir.name,bijk.m.name,k,str(pb)]+[str(item) for item in F])
			out_bs=filterBS(out_bs)
			out+=out_bs
		out=['\t'.join(item) for item in out]
		out='\n'.join(out)
		f.write(out+'\n')
		ct+=1
		progressbar(ct,N)  
	
	f.close()
	print('\r\n')
	#=======================================================================
	# removing temp files
	try:
		os.remove("%s_dp.ps"%(sessionID))
		os.remove("%s_lunp"%(sessionID))
		os.remove("%s.input"%(sessionID))
		os.remove("%s.out"%(sessionID))
	except:
		pass 
	#=============================			
if __name__=='__main__':
	main()
