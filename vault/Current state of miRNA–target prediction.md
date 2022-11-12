# Current State of miRNA–target Prediction

*(RPmirDIP: Reciprocal Perspective improves miRNA targeting prediction, 2020, https://doi.org/10.1038/s41598-020-68251-4)*

---

Notable examples of classical machine learning include:
- TarPMir[20](https://www.nature.com/articles/s41598-020-68251-4#ref-CR20 "Ding, J., Li, X. & Hu, H. Tarpmir: a new approach for microrna target site prediction. Bioinformatics 32, 2768–2775 (2016).")
- RFMirTarget[21](https://www.nature.com/articles/s41598-020-68251-4#ref-CR21 "Mendoza, M. R. et al. RFMirTarget: predicting human microRNA target genes with a random forest classifier. PLoS One https://doi.org/10.1371/journal.pone.0070153 (2013).")
- MirTarget[22](https://www.nature.com/articles/s41598-020-68251-4#ref-CR22 "Liu, W. & Wang, X. Prediction of functional microRNA targets by integrative modeling of microRNA binding and target expression data. Genome Biol. 20, 18 (2019).")

with recent models leveraging deep learning, including
- MiRTDL[23](https://www.nature.com/articles/s41598-020-68251-4#ref-CR23 "Cheng, S. et al. MiRTDL: a deep learning approach for miRNA target prediction. IEEE/ACM Trans. Comput. Biol. Bioinform. 13, 1161–1169 (2015).")
- DeepMirTar[24](https://www.nature.com/articles/s41598-020-68251-4#ref-CR24 "Wen, M., Cong, P., Zhang, Z., Lu, H. & Li, T. Deepmirtar: a deep-learning approach for predicting human miRNA targets. Bioinformatics 34, 3781–3787 (2018).")
- miRAW[25](https://www.nature.com/articles/s41598-020-68251-4#ref-CR25 "Pla, A., Zhong, X. & Rayner, S. miRAW: a deep learning-based approach to predict microRNA targets by analyzing whole microRNA transcripts. PLoS Comput. Biol. 14, e1006185 (2018).").

MirTarget is a Support Vector Machine (SVM) trained on CLIP experimentally validated interactions and miRNA overexpression data. The miRNA overexpression data provides a complimentary view to understanding functional targets as the elucidation of target interaction does not necessarily result in gene down-regulation[22](https://www.nature.com/articles/s41598-020-68251-4#ref-CR22 "Liu, W. & Wang, X. Prediction of functional microRNA targets by integrative modeling of microRNA binding and target expression data. Genome Biol. 20, 18 (2019).").

The miRTDL method implemented a Convolutional Neural Network (CNN) with selected features obtained from the convolved feature maps[23](https://www.nature.com/articles/s41598-020-68251-4#ref-CR23 "Cheng, S. et al. MiRTDL: a deep learning approach for miRNA target prediction. IEEE/ACM Trans. Comput. Biol. Bioinform. 13, 1161–1169 (2015).").

DeepMirTar used Stacked denoising Autoencoders (SdA) to learn a lower-dimensional representation of latent features[24](https://www.nature.com/articles/s41598-020-68251-4#ref-CR24 "Wen, M., Cong, P., Zhang, Z., Lu, H. & Li, T. Deepmirtar: a deep-learning approach for predicting human miRNA targets. Bioinformatics 34, 3781–3787 (2018).")

miRAW leveraged autoencoders without the denoising step[25](https://www.nature.com/articles/s41598-020-68251-4#ref-CR25 "Pla, A., Zhong, X. & Rayner, S. miRAW: a deep learning-based approach to predict microRNA targets by analyzing whole microRNA transcripts. PLoS Comput. Biol. 14, e1006185 (2018).").

TarPMir used a Random Forest (RF) classifier trained on an experimentally validated dataset[20](https://www.nature.com/articles/s41598-020-68251-4#ref-CR20 "Ding, J., Li, X. & Hu, H. Tarpmir: a new approach for microrna target site prediction. Bioinformatics 32, 2768–2775 (2016).").

RFMirTarget also used a Random Forest classifier, however it was trained on data originally pre-computed by miRanda, thus acting as a cascaded refinement of *ab initio* predictions[21](https://www.nature.com/articles/s41598-020-68251-4#ref-CR21 "Mendoza, M. R. et al. RFMirTarget: predicting human microRNA target genes with a random forest classifier. PLoS One
https://doi.org/10.1371/journal.pone.0070153
(2013).").
