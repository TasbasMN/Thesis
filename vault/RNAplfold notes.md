# My Notes

RNAplfold works by sliding a short RNA folding window of a certain length over the longer
mRNA sequence.

w = window size
L = constraint that controls base pairs between nucleotides separated by, at most L nucleotides
u = the region of mRNA that's being assessed

RNAplfold assesses accessibility of a region (of length u) by computing the p (prob. of that region u is unpaired) and averages this p over all windows W (W > L) that contain the stretch u

The RNA is folded locally in a window w
Within w, base pairing is restricted to at most L distance (w > L should be satisfied)
u is the stretch of nucleotides that are being assessed (miRNA seed in our case)

![](images/Pasted%20image%2020230226142409.png)
![](images/Pasted%20image%2020230226142443.png)

Accessibility is measured as the probability that a region of predefined length u is free of base pairing in thermodynamic equilibrium.


## RNAplfold output

```
#unpaired probabilities

#i$ l=1 2 3 4 5 6 7 8 9 10 11 12 13 14

1 0.8963457 NA NA NA NA NA NA NA NA NA NA NA NA NA

2 0.8866609 0.8319123 NA NA NA NA NA NA NA NA NA NA NA NA
```


The following block is taken from RNAplfold manpage, which I couldn't understand


> The output is a plain text matrix containing on each line a position i followed by the probability that i is unpaired, [i-1..i] is unpaired [i-2..i] is unpaired and so on to the probability that [i-x+1..i] is unpaired.

Another text from same page says that

> compute the mean probability that regions of length 1 to a given length are unpaired.

---

-o makes the RNAplfold generate a basic output with 3 columns, 2 nucleotide positions and a probability column which contains the base pairing prob. of said 2 nucleotides.

```
1 21 0.0573305

1 23 0.0114782

1 33 0.024513

2 20 0.0362148

2 31 0.0462913

2 32 0.0149115
```


# Literature Notes

*(Local RNA base pairing probabilities in large sequences, 2006, https://doi.org/10.1093/bioinformatics/btk014)*

The parameter L should be a bit larger than the structures of interest.

---

*(The impact of target site accessibility on the design of effective siRNAs, 2008, https://doi.org/10.1038/nbt1404)*

RNAplfold allows fast computation of local base-pair probabilities and accessibilities within mRNA transcripts by sliding a short RNA folding window of a certain length over the longer
mRNA sequence

More precisely, in our method, a window of size W is moved along the mRNA sequence, and the partition function for all local structures within the window is computed under the constraint that base pairing is allowed only between positions separated by, at the most, L nucleotides.

RNAplfold assesses the accessibility of a region of length u by computing the probability that this stretch is unpaired in thermodynamic equilibrium and averages this accessibility over all windows containing the region.

![](images/Pasted%20image%2020230226142409.png)
![](images/Pasted%20image%2020230226142443.png)

## The Optimal Parameters

with target site accessibility over a wide range of parameters analyzed, with the most significant separation resulting from 80 nucleotides (nts) and 40 nts for W and L,
respectively (Fig. 1b and Supplementary Fig. 1 online). These values may seem small as they exclude any structure with base-pairs that span more than 40 nts. However, it is clear that actively translated coding regions are largely devoid of long range structures, because these
structures are destroyed by the passing ribosome and are slow to reform .

## Optimal -u Parameters

When varying the length u of the unpaired region, we observed two parameter ranges with especially good separation (Fig. 1c). The first region with a P-value peak measures the accessibility of the 6–8 nucleotides starting at the 3¢ end of the target site and, therefore,
corresponds to the so-called seed region.This is in agreement with previous observations that the 5'-seed region of both siRNAs and microRNAs is the major determinant for RISC-mediated target recognition.

Furthermore, a second region with a P-value peak was observed for u values of 12–16, reminiscent of biochemical data showing that accessibility of the first 16 nts within the target site is required for highly efficient RISC-mediated cleavage 3.

## Methods

### Local Folding Using RNAplfold and Parameter Selection.

Accessibility parameters for mRNAs were determined using the RNAplfold program, which
implements a local folding algorithm to compute base-pairing probabilities and accessibility in a scanning window approach (window size W). The locality of the structure is determined by L, which specifies the maximal distance between two pairing positions. Accessibility is measured as the probability that a region of predefined length u is free of base pairing in thermodynamic equilibrium. Base-pairing probabilities (restricted by L) and accessibility (of the region u) are averaged over all windows W (W > L) that contain the stretch u. To assess
optimal folding parameters for the partition of functional and nonfunctional siRNAs, we performed a grid search by varying the window size W from 20 to 200 nts in a step size of 20 nts, the maximal binding lengths to L 1⁄4 1/4W, L 1⁄4 1/2W, L 1⁄4 3/4W and W-5 and the size of the region for which the accessibility was computed from 1 to 19 nts in steps of 1 nt. A Wilcoxon test (a nonparametric test for assessing whether two samples of observations come from the same distribution) was applied on the two sets of functional and nonfunctional siRNAs for all W, L, u, triplets (Supplementary Fig. 1). The parameter set (W 1⁄4 80 nts, L 1⁄4 40 nts, u 1⁄4 16 nts for data set 1 or 8 nts for data set 2) gave the overall best P-values as determined by a Wilcoxon test and were therefore considered as the optimal parameters for siRNA prediction and subsequently used in RNAxs. Note, that for an siRNA to be selected by RNAxs it is required that both accessibility scores (u 1⁄4 8 and u 1⁄4 16) lie
above thresholds.

---
---

*(Predicting effective microRNA target sites in mammalian mRNAs, 2015, https://doi.org/10.7554/eLife.05005)*

 we scored predicted structural accessibility (SA) as proportional to the log10 value of the unpaired probability for a 14-nt region centered on the match to miRNA nucleotides 7 and 8.

parameters -L 40 -W 80

parameters found to be optimal when evaluating siRNA efficacy ([Tafer et al., 2008](https://elifesciences.org/articles/05005#bib235)).

	 For the final predicted SA score used as a feature, we computed the log10 of the probability that a 14-nt segment centered on the match to sRNA positions 7 and 8 was unpaired.
