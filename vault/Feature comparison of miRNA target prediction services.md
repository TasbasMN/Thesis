
## Table
Notes: sRNA means small RNAs. In this context, sRNA is interchangable with miRNA.


| x                                                                  | Description                                                                                                                                                         | TargetScan | Tarpmir |
| ------------------------------------------------------------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ---------- | ------- |
| Seed match                                                         | Perfect match between seed & MRE sequences                                                                                                                          | x          | x       |
| Folding Energy                                                     | minimal folding energy calculated by RNAduplex from ViennaRNA package                                                                                               | x          | x       |
| Site Accessibility                                                 | The accessibility is a measurement of whether the target site region in the mRNA sequence is open for miRNA to binding.                                             | x          | x       |
| Local AU content                                                   | AU content near the site (Grimson et al., 2007; Nielsen et al., 2007)                                                                                               | x          | x       |
| Stem conservation                                                  | The Stem conservation was calculated as the Average PhyloP sore in the miRNA-mRNA binding stem region.                                                              |            | x       |
| Flanking conservation                                              | calculated as the Average PhyloP sore in the 40nt upstream and 40nt downstream of the binding site.                                                                 |            | x       |
| m/e motif                                                          | pairing probabilities at different positions of miRNA                                                                                                               |            | x       |
| Length of target mRNA region                                       | the length of miRNA binding target site region                                                                                                                      |            | x       |
| Total no of pairs                                                  | Total number of paired positions for each miRNA-mRNA binding site.                                                                                                  |            | x       |
| Length of largest consecutive pairing                              | Length of largest consecutive pairs between sRNA & mRNA                                                                                                             |            | x       |
| Pos. of largest consecutive pairing related to the 5' end of miRNA | the relative position of largest consecutive pairs to to 5’ end of miRNA                                                                                            |            | x       |
| No. of pairs at 3' end of miRNA (3'end=last 7 bases)               | the number of paired positions in the miRNA 3’ end.                                                                                                                 |            | x       |
| Numerical diff. btw. pairs of seed region and 3' end               | difference in # of paired positions between the seed region and the miRNA 3’ end region.                                                                            |            | x       |
| 3′-UTR target-site abundance                                       | Number of sites in all annotated 3′ UTRs (Arvey et al., 2010; Garcia et al., 2011)                                                                                  | x          |         |
| Predicted seed-pairing stability SPS                               | Predicted thermodynamic stability of seed pairing [(Garcia et al., 2011)](https://elifesciences.org/articles/05005#bib60)                                           | x          |         |
| sRNA position 1                                                    | Identity of nucleotide at position 1 of the sRNA                                                                                                                    | x          |         |
| sRNA position 8                                                    | sRNA position 8                                                                                                                                                     | x          |         |
| Site position 8                                                    | Identity of nucleotide at position 8 of the site                                                                                                                    | x          |         |
| 3′ supplementary pairing                                           | Supplementary pairing at the miRNA 3′ end [(Grimson et al., 2007)](https://elifesciences.org/articles/05005#bib78)                                                  | x          |         |
| Minimum distance                                                   | log10(Minimum distance of site from stop codon or polyadenylation site) (Gaidatzis et al., 2007; Grimson et al., 2007; Majoros and Ohler, 2007)                     | x          |         |
| Probability of conserved targeting PCT                             | Probability of site conservation, controlling for dinucleotide evolution and site context [(Friedman et al., 2009)](https://elifesciences.org/articles/05005#bib58) | x          |         |
| ORF length                                                         | log10(Length of the ORF)                                                                                                                                            | x          |         |
| 3′-UTR length                                                      | log10(Length of the 3′ UTR) [(Hausser et al., 2009)](https://elifesciences.org/articles/05005#bib90)                                                                | x          |         |
| 3′-UTR offset-6mer sites                                           | Number of offset-6mer sites in the 3′ UTR [(Friedman et al., 2009)](https://elifesciences.org/articles/05005#bib58)                                                 | x          |         |
| ORF 8mer sites                                                     | Number of 8mer sites in the ORF (Lewis et al., 2005; Reczko et al., 2012)                                                                                           | x          |         |

## Features of Tarpmir

Common features column indicates features that have been selected **significant** by all of the 4 ML models.
Common among other tools column indicates features that had been used in existing tools [(Peterson et al., 2014)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5018371/#btw318-B36).

| All Considered Features                                                                                                 | Final Features | Common features | Common among other tools |
| ----------------------------------------------------------------------------------------------------------------------- | -------------- | --------------- | ------------------------ |
| folding energy                                                                                                          | x              | x               | x                        |
| seed match                                                                                                              | x              | x               | x                        |
| site accessibility                                                                                                      | x              | x               | x                        |
| AU content                                                                                                              | x              | x               | x                        |
| stem conservation                                                                                                       | x              | x               | x                        |
| flanking conservation                                                                                                   | x              | x               | x                        |
| difference between stem and flanking conservation                                                                       |                |                 | x                        |
| m/e motif                                                                                                               | x              | x               |                          |
| the total number of paired positions                                                                                    | x              |                 |                          |
| the length of the target mRNA region                                                                                    | x              | x               |                          |
| the length of the largest consecutive pairs                                                                             | x              | x               |                          |
| the position of the largest consecutive pairs relative to the miRNA 5′ end                                              | x              |                 |                          |
| the length of the largest consecutive pairs allowing 2 mismatches                                                       |                |                 |                          |
| the position of the largest consecutive pairs allowing 2 mismatches                                                     |                |                 |                          |
| the number of paired positions at the miRNA 3′ end, where 3′ miRNA end was defined as the last 7 positions of the miRNA | x              |                 |                          |
| the total number of paired positions in the seed region and the miRNA 3′ end                                            |                |                 |                          |
| the difference between the number of paired positions in the seed region and that in the miRNA 3′ end                   | x              | x                |                          |
| exon preference                                                                                                         |                |                 |                          |

### Biological interpretation of m/e motif:

>The inclusion of the m/e motif implied that there existed preferred matching positions shared by all miRNAs. The length of the target site was selected, showing the importance of the binding preference of miRNAs to mRNA regions with specific lengths. The length of the largest consecutive pairing positions mattered, which extended the concept of seed match, as seed match was just a simple case with a long consecutive pairing positions. The difference between the number of paired positions in the seed region and that in the miRNA 3′ end also suggested that the seed match may be unimportant, given a high-quality 3′ end region matching. This also supported the idea that a long consecutive matching region is critical for functional miRNA target sites.

[dingTarPmiRNewApproach2016, https://doi.org/10.1093/bioinformatics/btw318]

### What are the effects of new columns of Tarpmir?
>The selected new features significantly improved the prediction accuracy of TarPmiR. To show the contribution of the new features to the accuracy of TarPmiR, we removed the seven new features and retrained random forests in TarPmiR. Compared with the original TarPmiR with 13 features, the recall and precision of the modified TarPmiR dropped 8.6% and 9.7%, respectively.

[dingTarPmiRNewApproach2016, https://doi.org/10.1093/bioinformatics/btw318]

## Extra Notes:

- **The accessibility** was proposed in Kertesz, M., Iovino, N., Unnerstall, U., Gaul, U. and Segal, E. (2007) The role of site accessibility in microRNA target recognition. Nature genetics, 39, 1278-1284.
-  The **local AU content** reflects the transcript AU content 30nt upstream and downstream of predicted site.
- **Length of target mRNA region:** For example, if miRNA x binds to mRNA y and the binding site between x and y are 28 nts region on mRNA y, this feature is 28.
![](images/me_motif.png)