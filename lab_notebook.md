## TODO:


- [x] target_abundance kolonu eklenecek (hedef UTR'da o miRNA'nın kaç tane non-overlapping hedefi var?)
- [x] figure out a way to implement target_conservation (TargetScan uses Pct)

# 2- 9 August 2023




## todo
rename midpoint (mrna percentage?)
rename close proximity
fetch annotations
convert grch38 miRNAs into grch37 (cli tool liftover
https://github.com/sritchie73/liftOverPlink/blob/master/README.md)

parallelize find_matches_for_vcf function (abstractten sonra)

xgb neden floating point output veriyor?
mirna_conservation na'leri sekansla en yakın olan mirna'nın değerini alacak şekilde doldur

check au_content function

convert 1-0 predictions back into floating point and find most different ones & showcase them with their dot bracket scheme

2 figures, feature importance & showcase


## this week:

midpoint ve close proximity kolonları kaldırıldı, çünkü transkript temelli değil DNA temelli gidiyoruz.

  örn: chr1 50m-50m+8 gibi pozisyonlarımız var. Bu pozisyonun midpointini, yani transkript üzerindeki görece yerini (0-1 arası) elde edemiyorum.  

- removed pred_energy scaling
- extended extension window for vcf positions from 11 to 30
- added few more hyperparams to tune
- in generate_positive_data, dropped all the results that have at least 2 difference between predicted & true positions
- removed _mutated suffix from mutated vcf entries as they have is_mutated column
- implemented sana's pipeline and predicted 2 rows with it
- case 1 almost completed
- case 2 file preprocess script completed and added into its respective folder in /data
- case 2 mutation processing script completed

  

case 1: miRNA'ya düşen mutasyon
case 2: mRNA'ya düşen mutasyon






## 21 - 26 july 2023:

- updated preprocess_clash.ipynb
  - refactored code
  - added full sequences and chr:start:end values using pyensembl
  - converted all biological indices to 0-based

- added accessions to the mirbase csv


## 15 June - 6 July 2023

- deprecated nucleotide_toolkit.py
- deprecated utilsv2_for_jupyter.py 
- created utils_latest.py that contains all the helper functions


todo:
 dot bracket şemasını printleyecek metod. ortaya git & işaretini bul, mismatchlere karşı satıra space ekle
     df = find_CLASH_V_sites(df) flag columnu 1 yapmıyor

düzeldikten sonra CLASH type trueleri bizim clash type predictionlarla eşleştir, karşılığına bak
yeni kolon fikri: RNAduplex'ten gelen satır sayısı


3 algoritmadan gelen sonuçları yan yana koy karşılaştırr
yeni kolonlar ekle

 
## 8-15 June 2023

- merge_pos_and_neg_data.ipynb shows how the previous seed match heuristic isn't working
- generate_data shows the steps of the data processing pipeline
- merge_pos_and_neg_data


## 20-27 April 2023

- revamped data folder
- new data, grosswendt 2014


opsiyonel:

- https://www.nature.com/articles/nsmb.2230 
- https://compgen.bio.ub.edu/datasets/454/docs/Parameters.html


breast cancer vcf'lerini indir, vcfteki adresten sekansı çek

çok sonraki işler:
- signaturelerle ilişkilendirilmesi
- sig based miRNA binding


13 Nisan

(N ve - mismatch olacak)


## 13-20 April 2023:

- wrote main_v2.0.ipynb and all its functions

- moved tools.py into nucleotide_toolkit.py to better represent its contents
  - added docstrings 

- added features.py that contains feature generating functions

- changed usages of "accession" column into "name" column for uniformity purposes between mb and ts datas

- added driver.py
- added supplementary.py for creating supplementary tsv files
- added utils_v2.py, removed unnecessary functions from utils.py

- changed create_supplementary_files.py into create_supplementary_files() in driver.py 
- deprecated generate_avg_position_column()


## 6-13 April 2023:

  - started hybrid.ipynb
  - added align_sequences()
  - added find_matches()
  - find_k_consecutive_bps()







## weeks

## 30 march - 6 april 2023

- Implement different canonical sites from CLASH
- blast.ipynb
- tried DP alignment algo
- CLASH_V_matches
- deprecated unpack_results_df()


## 21-30 march 2023

todo:

- create negative dataset

done:

- opened issue on mirbind github
- modified & ran find matches on clash

future:

- generate negative set
- find noncanonical matches (don't use dict matching)


## 9-21 march 2023

done:

- get_pair_probabilities_of_sequence() now appends "NaN" instead of 0
- added parse_clash_data() into create_supplementary_files.py
- added preprocess_clash_data() to utils.py
- added get_ensembl60_data() to utils.py
- added tools.py
- examined pyensembl use cases in pyensembl tutorial.ipynb 
- moved script to main.ipynb
- moved some functions to tools.py
- removed irrelevant functions from utils.py

## 2-9 march 2023

### done:

- moved new functions into utils.py
- tidied up main script

- added invoke_rnaduplex()
- added invoke_rnacofold()

- added generate_rnacofold_energy_column()
- added generate_rnaduplex_energy_column()

- renamed complement() into rna_complement()

- found & implemented a new feature, generate_stable_duplex_column()

### questions:

- 3 basamağa yuvarladım log10'u. Az mı? (accessibility functionında)
- append 0 mantıklı mı yoksa append math.nan mı gerekir (accessibility functionında)




## 23 february - 2 march 2023

### done

- fixed inconsistencies column selectors
  - all df["column"].tolist() uses changed into df["column"].values.tolist(), as the method is faster
- added one_hot_encode_match_types() function
- added get_pair_probabilities_of_sequence() function
- added get_accessibility_column() function with optional log10-scaled output

todo for next days:

- add rnaduplex wrapper functions for thermodynamical calcs
- start exploring xgb

## 16-23 february 2023

### General Stuff

- configured linux environment
- learned & explored docker containers for packaging up the project
  - could work inside a container if ViennaRNA gets problematic

Why containerizing my work is important?

- Reproducibility
- Easy collaboration
- Portability (runs on any machine that can run docker)
- Scalability (just deploy additional containers & manage with container orchestration tools such as kubernetes)

example containers for bioinformatics (from <https://biocontainers.pro/>)

- <https://biocontainers-edu.readthedocs.io/en/latest/introduction.html>

### Thesis Stuff

- extracted now useless db comparison methods to db_comparison.py
- generate_flanking_dinucleotides_columns() is fixed, removed arbitrary scoring of nucleotides

### ViennaRNA Stuff

- installed & ran ViennaRNA python plugin
- viennarna_examples.ipynb is used to figure out how ViennaRNA python package works
- Even asked ChatGPT for example uses, it delivered 👍🏻

### Tarpmir

- Thoroughly examined tarpmir.py to see how they implemented ViennaRNA calculations
  - they used subprocess.Popen() to run ViennaRNA instead of implementing into python
- Extracted useful snippets from tarpmir into a notebook
