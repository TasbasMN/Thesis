# RNAplfold

*(Predicting effective microRNA target sites in mammalian mRNAs, 2015, https://doi.org/10.7554/eLife.05005)*

## RNA Structure Prediction

[Request a detailed protocol](https://bio-protocol.org/eLIFErap05005?item=s4-5)

3′ UTRs were folded locally using RNAplfold ([Bernhart et al., 2006](https://elifesciences.org/articles/05005#bib13)), allowing the maximal span of a base pair to be 40 nucleotides, and averaging pair probabilities over an 80 nt window (parameters -L 40 -W 80), parameters found to be optimal when evaluating siRNA efficacy ([Tafer et al., 2008](https://elifesciences.org/articles/05005#bib235)). For each position 15 nt upstream and downstream of a target site, and for 1–15 nt windows beginning at each position, the partial correlation of the log10(unpaired probability) to the log2(mRNA fold change) associated with the site was plotted, controlling for known determinants of targeting used in the context+ model, which include min_dist, local_AU, 3P_score, SPS, and TA ([Garcia et al., 2011](https://elifesciences.org/articles/05005#bib60)). For the final predicted SA score used as a feature, we computed the log10 of the probability that a 14-nt segment centered on the match to sRNA positions 7 and 8 was unpaired.

---

A second feature that we re-evaluated was the predicted structural accessibility of the site. As scored previously, the degree to which the site nucleotides were predicted to be free of pairing to flanking 3′-UTR regions was not informative after controlling for the contribution of local AU content ([Grimson et al., 2007](https://elifesciences.org/articles/05005#bib78)). However, analysis inspired by work on siRNA site accessibility ([Tafer et al., 2008](https://elifesciences.org/articles/05005#bib235)) suggested an improved scoring scheme for this feature. For this analysis we used RNAplfold ([Bernhart et al., 2006](https://elifesciences.org/articles/05005#bib13)) to predict the unpaired probabilities for variable-sized windows in the proximity of the site and then examined the relationship between these probabilities and the repression associated with sites in our compendium of normalized datasets, while controlling for local AU content and other features of the context+ model ([Figure 4A](https://elifesciences.org/articles/05005#fig4)). Based on these results, which resembled those reported previously ([Tafer et al., 2008](https://elifesciences.org/articles/05005#bib235)), we scored predicted structural accessibility (SA) as proportional to the log10 value of the unpaired probability for a 14-nt region centered on the match to miRNA nucleotides 7 and 8.


