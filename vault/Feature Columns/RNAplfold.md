# RNAplfold

*(Predicting effective microRNA target sites in mammalian mRNAs, 2015, https://doi.org/10.7554/eLife.05005)*

## RNA Structure Prediction

[Request a detailed protocol](https://bio-protocol.org/eLIFErap05005?item=s4-5)

3′ UTRs were folded locally using RNAplfold ([Bernhart et al., 2006](https://elifesciences.org/articles/05005#bib13)), allowing the maximal span of a base pair to be 40 nucleotides, and averaging pair probabilities over an 80 nt window (parameters -L 40 -W 80), parameters found to be optimal when evaluating siRNA efficacy ([Tafer et al., 2008](https://elifesciences.org/articles/05005#bib235)). For each position 15 nt upstream and downstream of a target site, and for 1–15 nt windows beginning at each position, the partial correlation of the log10(unpaired probability) to the log2(mRNA fold change) associated with the site was plotted, controlling for known determinants of targeting used in the context+ model, which include min_dist, local_AU, 3P_score, SPS, and TA ([Garcia et al., 2011](https://elifesciences.org/articles/05005#bib60)). For the final predicted SA score used as a feature, we computed the log10 of the probability that a 14-nt segment centered on the match to sRNA positions 7 and 8 was unpaired.
