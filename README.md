
# dAsthma source codes for "Latent-space embedding of expression data identifies gene signatures from sputum samples of asthmatic patients"

**Backgrounds**: The pathogenesis of asthma is a complex process involving multiple genes and pathways. Identifying biomarkers from asthma datasets, especially those that include heterogeneous subpopulations, is challenging. In this work, we developed a framework that incorporates a denoising autoencoder and a supervised learning approach to identify gene signatures related to asthma severity. The autoencoder embeds high-dimensional gene expression data into a lower-dimensional latent space in an unsupervised fashion, enabling us to extract distinguishing features from gene expression data.

**Results**: Using the trained autoencoder model, we found that the weights on hidden units in the latent space correlate well with previously defined and clinically relevant clusters of patients. Moreover, pathway analysis based on each gene's contribution to the hidden units showed significant enrichment in known asthma-related pathways. We then used genes that contribute most to the hidden units to develop a secondary supervised classifier (based on random forest) for directly predicting asthma severity. The random-forest importance metric from this classifier identified a signature based on 50 key genes, which can predict severity with an AUROC of 0.81 and thus have potential as diagnostic biomarkers. Furthermore, the key genes could also be used for successfully estimating, via support-vector-machine regression, the FEV1/FVC ratios across patients, achieving pre- and post-treatment correlations of 0.56 and 0.65, respectively (between predicted and observed values). 

**Conclusions**: The denoising autoencoder framework could extract meaningful functional genes and patient groups from the gene expression profile of asthma patients. These patterns may provide potential sources for biomarkers for asthma severity.


