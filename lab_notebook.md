## TODO:

- küçük harf olan fastalarda script patlar
- çizgi olan fastalarda script patlar
- N olan fastalarda script patlar
- discarded featureleri biz de deneyip sonucu etkilemediğini görsek mi?


- [x] target_abundance kolonu eklenecek (hedef UTR'da o miRNA'nın kaç tane non-overlapping hedefi var?)
- [x] figure out a way to implement target_conservation (TargetScan uses Pct)




opsiyonel:

- https://www.nature.com/articles/nsmb.2230 
- https://compgen.bio.ub.edu/datasets/454/docs/Parameters.html


fonksiyonları birleştirip sonuç df'ini hazırla

negatif dataseti hallet

build xgb

breast cancer vcf'lerini indir, vcfteki adresten sekansı çek

çok sonraki işler:
- signaturelerle ilişkilendirilmesi
- sig based miRNA binding




## Problems:

- ta_sps tsv can't be created because it is an excel file. find a way

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
