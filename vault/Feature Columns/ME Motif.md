# Me Score

The m/e motif describes how different positions in miRNAs match the corresponding positions in target sites. Here two positions match means that nucleotides at the two positions are complement to each other. For instance, nucleotides at positions in miRNA seed regions tend to match the nucleotides at the corresponding positions in target sites and nucleotides at positions in other miRNA regions tend to form mismatches or bulges with the corresponding positions in target sites. We thus have a sequential pattern composed of two letters ‘m’ and ‘e’ to describe preferred matching and non-matching positions, respectively. To calculate the m/e scores, for each position in miRNAs, we calculate a probability _pi_ that this position matches the corresponding position in target sites by using all positive target sites in the training dataset.

  ![](../images/me_motif.png)

_(TarPmiR: A new approach for microRNA target site prediction, 2016, https://doi.org/10.1093/bioinformatics/btw318)_
