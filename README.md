# microbiome-network-analysis-research

## Overview

This research project aims to understand how microbial compositions shift in response to carbon source perturbation. We use a Compositional Data Analysis (CoDA) approach instead of the standard microbial pipeline to gain better insight into the gradient shifts and respect the inherent compositional nature of high-throughput sequencing data.

## Theory

The steps in a standard microbiome analysis toolkit have compositional replacements.

![Compositional Data Analysis Pipeline](img/compositional_data_analysis_pipeline.jpg)
Image Source: [Microbiome Datasets Are Compositional: And This Is Not Optional (Gloor 2017)](https://www.frontiersin.org/articles/10.3389/fmicb.2017.02224/full)

In this analysis, we use:
1. Isometric-Log Ratio (ILR) Transform
2. Aitchison Distance
3. PCA (Variance)
4. ANCOM Differential Abundance

### Theory References

**Compositional Data Analysis:**  
* Introduction to CoDA: [Microbiome Datasets Are Compositional: And This Is Not Optional (Gloor 2017)](https://www.frontiersin.org/articles/10.3389/fmicb.2017.02224/full)
* Case Study: [Compositional Data Analysis Approaches to Improve Microbiome Studies: from Collection to Conclusions](https://www.youtube.com/watch?v=j1IbfQrT2Cs)

**Balance Trees: Isometric-Log Ratio Transform**  
In this analysis, we use ILR. In microbial ecology, PCA on various distance metrics are used to explore sources of variation in the data. The ILR transform is the only compositional transform that is suitable as a feature space for PCA exploratory data analysis. [(Source 1)](https://stats.stackexchange.com/questions/305965/can-i-use-the-clr-centered-log-ratio-transformation-to-prepare-data-for-pca) [(Source 2)](http://www.statsathome.com/2017/08/09/we-can-do-better-than-the-alr-or-softmax-transform/)

## Code Overview

### Process Data Files

**Dependencies:**
* Python skbio: [set-up](http://scikit-bio.org/)
* Python pandas, numpy, re

**Run:** `jupyter notebook scripts/process_data.ipynb` to create the following files:

* Sample Name x Meta Information
    * `data/processed/sample_metadata.tsv`
* OTU ID x Sample ID x Absolute Counts
    * `data/processed/OTU_counts.tsv`
* OTU ID x Taxonomy, Functions
    * `data/processed/feature_metadata.tsv`

### Make PCA on LR Data

**Dependencies:**
* R robCompositions: [issue: does not work with microbiome data](https://github.com/matthias-da/robCompositions/issues/10)

Actually going to follow analysis from 


Steps:
* Compositional Data --> ILR Coordinates
* PCA on ILR Coordinates
* Transform PCA_ILR back to CLR

---

## Use QIIME2

**Dependencies:**
* BIOM: [set-up](https://biom-format.org/index.html)
* QIIME2: [set-up](https://docs.qiime2.org/2021.11/install/)
* GNEISS: `pip install git+https://github.com/biocore/gneiss.git`

Activate QIIME Environment: `conda activate qiime2-2021.11`

### Creating Balances
The selection of balance tree controls for variation and allows us to identify interesting differentially abundant partitions of taxa. Through gneiss, there are three approaches:

* **Correlation Clustering:** If we don’t have relevant prior information about how to cluster together organisms, we can group together organisms based on how often they co-occur with each other. This is available in the `correlation-clustering` command and creates tree input for `ilr-hierarchical`.
* **Gradient Clustering:** Use a metadata category to cluster taxa found in similar sample types. For example, if we want to evaluate if pH is a driving factor, we can cluster according to the pH that the taxa are observed in, and observe whether the ratios of low-pH organisms to high-pH organisms change as the pH changes. This is available in the `gradient-clustering` command and creates tree input for `ilr-hierarchical`.
* **Phylogenetic Analysis:** A phylogenetic tree (e.g. `rooted-tree.qza`) created outside of gneiss can also be used. In this case you can use your phylogenetic tree as input for `ilr-phylogenetic`.

**Resources**
* [QIIME Docs Tutorial: Differential abundance analysis with gneiss](https://docs.qiime2.org/2021.11/tutorials/gneiss/)
* [GNEISS Docs Tutorial: Linear regression on balances in the 88 soils](https://biocore.github.io/gneiss/docs/v0.4.0/tutorials/qiime2/88soils-qiime2-tutorial.html)
* [GNEISS Docs Tutorial: Linear mixed effects models on balances in a CF study](https://biocore.github.io/gneiss/docs/v0.4.0/tutorials/qiime2/cfstudy-qiime2-tutorial.html)
* [GNEISS Docs Tutorial: Linear regression on balances in the Chronic Fatigue Syndrome](https://biocore.github.io/gneiss/docs/v0.4.0/tutorials/qiime2/cfs-qiime2-tutorial.html)

Activate the environment: `conda activate qiime2-2021.11`
Then run `jupyter notebook scripts/qiime.ipynb` to do the following analysis:

## Use gneiss

**Dependencies:**
* gneiss: [set-up](https://github.com/biocore/gneiss)

**Resources**
* [GNEISS Docs Tutorial: Linear regression on balances in the 88 soils](https://biocore.github.io/gneiss/docs/v0.4.0/tutorials/python/88soils-python-tutorial.html)
* [GNEISS Docs Tutorial: Linear mixed effects models on balances in a CF study](https://biocore.github.io/gneiss/docs/v0.4.0/tutorials/python/cfstudy-python-tutorial.html)
* [GNEISS Docs Tutorial: Linear regression on balances in the Chronic Fatigue Syndrome](https://biocore.github.io/gneiss/docs/v0.4.0/tutorials/python/cfs-python-tutorial.html)