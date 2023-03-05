
# 2-9 march 2023

- moved new functions into utils.py
- added calculate_free_energy()
- removed fasta header addition & removal from calculate_free_energy()
- added generate_hybridization_energy_column()


## todo

unpack_results_df unpacks more than needed. deprecate.
- name of complement() function is not great. rename it.
- küçük harf olan fastalarda script patlar
- çizgi olan fastalarda script patlar
- N olan fastalarda script patlar



# 23 february - 2 march 2023


## done

- fixed inconsistencies column selectors
  - all df["column"].tolist() uses changed into df["column"].values.tolist(), as the method is faster
- added one_hot_encode_match_types() function
- added get_pair_probabilities_of_sequence() function
- added get_accessibility_column() function with optional log10-scaled output

todo for next days:

- add rnaduplex wrapper functions for thermodynamical calcs
- start exploring xgb

# 16-23 february 2023

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
