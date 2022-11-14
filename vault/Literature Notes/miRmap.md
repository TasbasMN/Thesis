# miRmap

*(miRmap: Comprehensive prediction of microRNA target repression strength, 2012, https://doi.org/10.1093/nar/gks901)*

This note contains implementation notes of features:
- Thermodynamics
- Motif occurence
- Conservation

---

## Target Prediction Features

### Thermodynamics of miRNA–mRNA Interactions

>The miRNA–mRNA pair forms an RNA duplex. Using the Vienna RNA Secondary Structure library ([29](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3526310/#gks901-B29)), we computed the minimum free folding energy (MFE) of this duplex (with the ‘cofold’ function), and named it ‘ΔG duplex’.

>The RISC is much larger than the miRNA ([30](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3526310/#gks901-B30)) and must bind to the mRNA in an extended single-stranded form. We computed the energy required to unfold the 3′-UTR region of the target site (this area can be optionally extended), named ‘ΔG open’, similarly to PITA ([12](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3526310/#gks901-B12)), with the ‘pf_fold’ function from the Vienna Library ([29](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3526310/#gks901-B29)). The computation of ‘ΔG open’ requires two energy calculations; the free energy of the mRNA constrained to maintain the target site single stranded is subtracted from the free energy of the same unconstrained mRNA. The single-strand constraint was placed on a segment of 70 nucleotides centred on the target site. Finally, ‘ΔG open’ summed with ‘ΔG duplex’ or ‘ΔG binding’ gives the total system energy: we named it ‘ΔG total’ (named ΔΔG in PITA ([12](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3526310/#gks901-B12))).

### Probability of the Motif Occurrence

> We modelled the 3′-UTR sequence as a Markov process (order 1, as 3′-UTR sequences are too short to parameterize higher orders) and determined the expected probability of finding at least *n* occurrences of the motif defined as either an exact seed match or the full miRNA binding site, using two different methods. In the first method, the probability distribution was approximated with a binomial distribution, as in Marín and Vanícek ([18](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3526310/#gks901-B18)), while in the second method, we computed the exact probability distribution based on the theoretical work of Nuel *et al.* ([31](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3526310/#gks901-B31)).

### Conservation of the Target Site

> Using the UCSC ([27](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3526310/#gks901-B27)) MultiZ multiple genome sequence alignments (hg19, MultiZ 46-way; mm9, MultiZ 30-way), we searched for conserved miRNA target sites in the alignment blocks defined by the 3′-UTRs of the reference species (human or mouse for the HITS-CLIP data). From a mammalian species tree (UCSC ([27](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3526310/#gks901-B27))), we first pruned all the species that did not contain the target site. We then summed the lengths of the remaining branches (as in ([32](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3526310/#gks901-B32))) to obtain the branch length score (BLS). As implemented by Friedman *et al.* ([19](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3526310/#gks901-B19)), we summed the branch lengths of the species topology fitted for each 3′-UTR alignment with the REV model using the PhyloFit program from the PHAST suite ([33](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3526310/#gks901-B33)). Tree manipulations were done with the DendroPy ([34](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3526310/#gks901-B34)) library.
