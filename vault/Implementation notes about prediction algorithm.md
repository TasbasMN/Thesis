# Implementation Notes About Prediction Algorithm

This is a random collection about stuff that will be useful during the implementation.

---

TargetScan doesn't use the noncanonical interaction "centered sites" as it doesn't appear to be functional.

>Centered sites, which do not appear to be functional for most miRNAs are no longer considered or displayed in TargetScan.

*(The biochemical basis of microRNA targeting efficacy, 2019, https://doi.org/10.1126/science.aav1741)*

---

>3'UTR sites within 15 nt of the stop codon were less effective on the array, compared to sites elsewhere in the 3'UTR. These results confirmed that the segment immediately following the stop codon was inhospitable for targeting.

*(MicroRNA Targeting Specificity in Mammals: Determinants beyond Seed Pairing, 2007, https://doi.org/10.1016/j.molcel.2007.06.017)*

---

There are four commonly used features for miRNA target prediction tools:
1. seed match
2. conservation
3. free energy
4. site accessibility

*(Common features of microRNA target prediction tools, 2014, https://doi.org/10.3389/fgene.2014.00023)*

---

Noncanonical interactions are as important as canonical ones.

>The binding of most miRNAs includes the 5′ seed region, but around 60% of seed interactions are noncanonical, containing bulged or mismatched nucleotides. Moreover, seed interactions are generally accompanied by specific, nonseed base pairing. 18% of miRNA-mRNA interactions involve the miRNA 3′ end, with little evidence for 5′ contacts, and some of these were functionally validated.

*(Mapping the human miRNA interactome by CLASH reveals frequent noncanonical binding, 2013, https://doi.org/10.1016/j.cell.2013.03.043)*

---

miRNA - miRNA complexes exist, but they are very rare.

>Interactions between pairs of distinct miRNAs were not very frequent (∼3%), but some were highly reproducible and apparently isoform specific

*(Mapping the human miRNA interactome by CLASH reveals frequent noncanonical binding, 2013, https://doi.org/10.1016/j.cell.2013.03.043)*

---

TargetScan discards other mRNA regions than 3' UTR. TargetScan also discards first 15 nt of 3' UTR, because of steric hindrance caused by ribosomes.

>Sites in 5'UTRs, in ORFs, or within 15 nt of the stop codon were excluded, in accord with our results showing that such sites were generally not effective.

*(Predicting effective microRNA target sites in mammalian mRNAs, 2015, https://doi.org/10.7554/eLife.05005)*

---

>For instance, although matching seed is not always sufficient for a functional miRNA–mRNA interaction (Brennecke et al., 2005; Didiano and Hobert, 2006), it has been thought to be necessary for most animal miRNA–mRNA binding. **However, studies have shown non-canonical pairings that allow G:U wobbles and even mismatches can be functional (Brennecke et al., 2005; Didiano and Hobert, 2006).**

*(TarPmiR: A new approach for microRNA target site prediction, 2016, https://doi.org/10.1093/bioinformatics/btw318)*

---

There are some implementation notes found in [miRmap](Literature%20Notes/miRmap.md).
