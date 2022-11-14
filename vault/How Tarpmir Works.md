# How Tarpmir Works

TarPmiR predicts miRNA target sites in three steps with the input of a set of miRNAs and a set of mRNAs.

1. First, TarPmiR generates candidate target sites based on seed match or minimal folding energy ([Enright](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5018371/#btw318-B13) _[et al.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5018371/#btw318-B13)_[, 2004](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5018371/#btw318-B13); [Grimson](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5018371/#btw318-B16) _[et al.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5018371/#btw318-B16)_[, 2007](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5018371/#btw318-B16); [Yousef](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5018371/#btw318-B52) _[et al.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5018371/#btw318-B52)_[, 2007](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5018371/#btw318-B52)). For a given miRNA, TarPmiR scans an mRNA sequence with the seed region of the miRNA (positions 2–7) to find perfect seed-matching sites.
2. In addition, TarPmiR applies RNA-duplex from the Vienna RNA package ([Hofacker, 2003](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5018371/#btw318-B19)) to obtain the top target sites with the lowest folding energy.
3. For each candidate target sites, TarPmiR calculates the values of the 13 selected features ([Supplementary File S1](http://bioinformatics.oxfordjournals.org/lookup/suppl/doi:10.1093/bioinformatics/btw318/-/DC1)).
4. Finally, TarPmiR applies the trained random-forest based predictor to predict target sites.

## Dataset

_(Mapping the human miRNA interactome by CLASH reveals frequent noncanonical binding, 2013, https://doi.org/10.1016/j.cell.2013.03.043)_

## Methods

To select important features, we applied the following four machine learning methods: step-wise logistic regression ([Ralston and Wilf, 1960](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5018371/#btw318-B38)), least absolute shrinkage and selection operator (LASSO) ([Tibshirani, 1996](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5018371/#btw318-B45)), randomized logistic regression ([Meinshausen and Bühlmann, 2010](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5018371/#btw318-B32)) and random forests ([Svetnik et al., 2003](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5018371/#btw318-B43)).

### Training

We downloaded 18 514 miRNA target sites of 399 miRNAs from CLASH experiments (Helwak et al., 2013). These target sites were considered as positive target sites.

We also generated 18 514 corresponding negative or ‘false’ target sites in a manner similar to a previous study (Li et al., 2014), with the following criteria:
- A positive site and its corresponding negative site are on the same mRNA
- The positive and its corresponding negative site has similar CG dinucleotide frequency
- The positive and its corresponding negative site has similar number of the nucleotide G; (iv) A negative site does not overlap with any positive site
- With multiple candidate negative sites in an mRNA, select the one with the lowest folding energy.

### Testing

To determine which method to be used, we randomly chose 10 000 positive sites and 10 000 negative sites for training and the remaining positive and negative sites for testing. We repeated this process five times and selected the method with the F2 scores.

We also collected two independent PAR-CLIP datasets from the human HEK293 cell line for testing. PAR-CLIP datasets were used because a large number of potential miRNA target regions called crosslink-centered regions (CCRs) could be obtained from PAR-CLIP. CCRs were considered as positive target sites.

To test TarPmiR on general datasets, we compared the TarPmiR predictions with the experimentally validated miRNA targets by general methods in TarBase 7.0 (Vlachos et al., 2014).

The output of the random-forest model is the predicted probability that a candidate target site is a true target site. We have compared nine probability cutoffs to define target sites using the F2 score, since we put more emphasis on the recall than the precision. The cutoffs 0.5 and 0.6 have almost the similar F2 scores, while the cutoff 0.5 has the largest recall ([Supplementary File S2](http://bioinformatics.oxfordjournals.org/lookup/suppl/doi:10.1093/bioinformatics/btw318/-/DC1)). Therefore, we used 0.5 for the following analyses. We provide a parameter –_p_ in TarPmiR, users can choose other cutoffs based on their own needs.

---

Tried Methods:
- Stepwise SVM
- LASSO
- Randomized logistic regression
- Random Forest
