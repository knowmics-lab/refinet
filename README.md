# REFINET: Reference Tissue Identification using Neural Transformers

## Source code usage

To use the REFINET source code, you first need to download the dataset from the 
following link:
[https://refinet.s3.eu-west-1.amazonaws.com/refinet_dataset.tar.xz](https://refinet.s3.eu-west-1.amazonaws.com/refinet_dataset.tar.xz). 

After downloading the dataset, you can extract the files under the `src/datasets` folder
with the following command:

```bash
tar -Jxvf <path_to_downloaded_file>/refinet_dataset.tar.xz -C <path_to_your_project>/src/datasets
```

## Background

Studying complex diseases often requires RNA-seq analysis to understand the dynamic mechanisms behind pathogenesis, drug resistance, and other disease-related phenomena. However, analyzing such data involves not only collecting samples from patients but also obtaining reference samples. These reference samples are crucial for determining whether changes have occurred in normal expression patterns. Therefore, careful attention must be paid to the selection of controls.

Unfortunately, obtaining reference samples can be challenging, especially when involving patients, as healthy individuals may not be able to undergo invasive or potentially dangerous procedures for various reasons.

Fortunately, databases such as GEO, SRA, and EBI, along with projects like TCGA, TARGET, and GTEx, have assembled numerous samples from healthy individuals. If utilized effectively, these resources could be invaluable in addressing this issue.

In this context, Zeng W. Z. D. et al. (2019), further expanded upon by Zeng B. et al. (2021), recently proposed a methodology for identifying source tissues and selecting reference samples using autoencoders.

## Methods

In this study, we developed REFINET (Reference Tissue Identification using Neural Transformers), a novel AI-based tool that utilizes the generalization capabilities of transformer neural networks for reference tissue identification. It also embeds the data in a lower-dimensional space to facilitate the selection of the most suitable reference samples.

First, we created an internal sample database containing expression data from normal and tumor samples sourced from TCGA, TARGET, GTEx, and GEO. This database comprises 19,253 samples across 60,498 genes. The data were normalized using a TPM transformation.

Given user-provided input samples in the form of raw counts, the algorithm first transforms these data into TPM-normalized values. It then applies a neural transformer autoencoder (the architecture details are shown in panel A of the figure), which is trained to encode the data into a lower-dimensional latent space of 64 dimensions while predicting the reference tissue.
Next, we select all control samples from our internal database that match the identified reference tissue. When the number of available control samples is fewer than the number of input samples, we return the entire set of control samples. However, if there are enough control samples, we implement a greedy selection procedure that ranks all control samples based on Spearman’s correlation applied to the latent space. We then collect control samples with the highest correlation values and repeat this process until we have gathered sufficient control samples.

## Results

To train our network, we divided our sample data into three sets: a 10% holdout set, a 10% test set, and a 90% training set. We utilized the ADAM optimizer for 400 epochs to train the network using the training and test sets. The training converged after approximately 250 epochs. We then evaluated the results on the holdout set for both the autoencoder and classifier components of the network. The results showed a mean squared error (MSE) of 0.1514 for the autoencoder and an accuracy of 98.42% for the classifier (detailed training results are available on GitHub).

Subsequently, we downloaded 34 independent datasets from the GEO database, which included a total of 2,148 samples from 16 different tissues. We used this dataset to further test the classifier, yielding an accuracy of 92.52% classifier (detailed training results are available on GitHub).

Finally, we employed Limma to identify Differentially Expressed Genes (DEGs) with a FDR < 0.05. We compared these DEGs with those obtained using controls predicted by our methodology, as well as those identified with Zeng et al. (2021) and a random selection. The results, presented in panels B and C of the figure, demonstrate that our methodology identifies a set of controls that shows greater consistency than the other approaches.
