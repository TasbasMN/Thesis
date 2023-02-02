# deepTarget, End-to-end Learning Framework for microRNA Target Prediction Using Deep Recurrent Neural Networks

According to Menor et al., as many as 151 kinds of features appear in the literature, which can
be broadly grouped into four common types:
- the degree of Watson-Crick matches of a seed sequence
- the degree of sequence conservation across species
- Gibbs free energy, which measures the stability of the binding of a miRNA-mRNA pair
- the site accessibility, which measures the hybridization possibility of a pair from their secondary structures

---

The miRNA-mRNA binding sites can be classified into three types :
- 5’-dominant canonical
- 5’-dominant seed only
- 3’compensatory

![](../images/3%20types%20of%20miRNA-mRNA%20interactions.jpg)

> Figure 2. Approximate secondary structures of the three main types of target site duplex. **(a)** Canonical sites have good or perfect complementarity at both the 5′ and 3′ ends of the miRNA, with a characteristic bulge in the middle. **(b)** Dominant seed sites have perfect seed 5′ complementarity to the miRNA but poor 3′ complementarity. **(c)** Compensatory sites have a mismatch or wobble in the 5′ seed region but compensate through excellent complementarity at the 3′ end.

*(Prediction of microRNA targets, 2007, https://doi.org/10.1016/j.drudis.2007.04.002)*

---

Biological sequences (such as DNA, RNA, and protein sequences) naturally fit the recurrent neural networks that are capable of temporal modeling. Nonetheless, prior work on applying deep learning to bioinformatics utilized only convolutional and fully connected neural networks. The biggest novelty of our work lies in its use of recurrent neural networks to model RNA sequences and further learn their sequence-to-sequence interactions, without laborious feature engineering (e.g., more than 151 features of miRNA-target pairs have been proposed in the literature). As shown in our experimental results, even without any of the known features, deepTarget delivered
substantial performance boosts (over 25% increase in F-measure) over existing miRNA target detectors, demonstrating the effectiveness of recent advances in end-to-end learning methodologies.

---

Notably, deepTarget does not depend on any sequence alignment operation, which has been used in many bioinformatics pipelines as a holy grail to reveal similarity/interactions between sequences.
Although effective in general, sequence alignment is susceptible to changes in parameters (e.g., gap/mismatch penalty and match premium) and often fails to reveal the true interactions between sequences, as is often observed in most of the alignment-based miRNA target detectors. By processing miRNA and RNA sequences with RNN-based auto-encoders without alignment, DeepTarget successfully discover the inherent sequence representations, which are effectively used in the next step of deepTarget for interaction learning.
