
This feature assigns m (match) or e (no match) characters to each of the bases of seed sequence.

**Examples:**

8mer: MMMMMMMM
7mer-A1: MMMMMMME

---

The m/e motif describes how different positions in miRNAs match the corresponding positions in target sites. Here two positions match means that nucleotides at the two positions are complement to each other. For instance, nucleotides at positions in miRNA seed regions tend to match the nucleotides at the corresponding positions in target sites and nucleotides at positions in other miRNA regions tend to form mismatches or bulges with the corresponding positions in target sites. We thus have a sequential pattern composed of two letters ‘m’ and ‘e’ to describe preferred matching and non-matching positions, respectively. To calculate the m/e scores, for each position in miRNAs, we calculate a probability _pi_ that this position matches the corresponding position in target sites by using all positive target sites in the training dataset.

  ![](../images/me_motif.png)

_(TarPmiR: A new approach for microRNA target site prediction, 2016, https://doi.org/10.1093/bioinformatics/btw318)_

# Biological Interpretation of m/e Motif:

>The inclusion of the m/e motif implied that there existed preferred matching positions shared by all miRNAs. The length of the target site was selected, showing the importance of the binding preference of miRNAs to mRNA regions with specific lengths. The length of the largest consecutive pairing positions mattered, which extended the concept of seed match, as seed match was just a simple case with a long consecutive pairing positions. The difference between the number of paired positions in the seed region and that in the miRNA 3′ end also suggested that the seed match may be unimportant, given a high-quality 3′ end region matching. This also supported the idea that a long consecutive matching region is critical for functional miRNA target sites.

_(TarPmiR: A new approach for microRNA target site prediction, 2016, https://doi.org/10.1093/bioinformatics/btw318)_
