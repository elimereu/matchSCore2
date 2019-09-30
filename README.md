
## Quick introduction

MatchSCore2 is a R package for the comparison of single cell RNA-seq data across tools and experiments. The package allows a gene marker-based projection of single cells onto a reference sample and, thus, the identification of cell types in unknown cells.  A more detailed version of the method is at the bioRxiv paper:  http://dx.doi.org/10.1101/630087. The code to reproduce the downstream analyses of the paper can be found at the folder Benchmarking_of_scRNA-seq_protocols in this repository. The snakemake workflow and additional code used for the manuscript is also at the following github: https://github.com/ati-lz/HCA_Benchmarking.

If you want to use our PBMC+HEK293T QCed data (UMI counts + matchSCore2 annotations) from all 13 protocols, you can download the data object at the following dropbox link: https://www.dropbox.com/s/indoehhesx6f457/sce.all_technologies.RData?dl=0. 
It’s a unique SingleCellExperiment object from which you can extract the counts from each protocol.
After you load it into R, if you check the colData(object) there are three metadata, which are: nnet2, ident, batch .

1. nnet2 is the annotation from the matchSCore2 classifier. 
2. ident is the Seurat clustering result. The clusters are manually annotated by looking are known gene markers.
3. batch is the protocol. 

To subset each technology you can use this command: 

sce_technology <- sce[,which(sce@colData$batch==“technology.name”)]

or with the function “subset”.. (depends on the version you have installed).

where sce is the name of the SingleCellExperiment object. 

## Installation

matchSCore2 is a new version of the older matchSCore package. This is a development version running on R (>= 3.5.1). The package makes use of nnet, Matrix, splatter, MixSim and ggplot2 libraries. Other libraries are also used for clustering and identification of markers (e.g. Seurat).

matchSCore2 can be installed by devtools:

```{r,eval=FALSE}

library("devtools")
install_github('elimereu/matchSCore2')


```


## Functionalities of matchSCore2

MatchSCore2 has the following main functions:

1. Annotation of cells throught a supervised classification model based on a reference dataset. (**Cell annotations**)
2. Weighted GO pathways enrichment of cluster-specific markers based on the Kolmogorov-Smirnov test of genes p-values. (**GOannotation function**)
3. Measure and visualize the matching between two clustering. (**Cluster matching**) 
4. Track the accuracy trend of a tool in clustering and marker identification compared with the optimal solution provided by a simulated data set. In this case, matchSCore2 works in combination with the Splatter package - https://github.com/Oshlack/splatter/blob/master/vignettes/splatter.Rmd - (**Benchmarking**)


## Usage

To predict cell identities matchSCore2 requires:

1. ``` scale.data ```: A matrix of log-normalized and scaled gene expression values from the reference dataset (like the matrix scale.data in Seurat object).
2. ``` clus ```: A named factor with reference identities (like in the @ident slot in Seurat object).
3. ``` gene_cl.ref ```: A named list of markers. Each element of the list contains cell type specific gene markers (Usually top100 ranked markers of each cell type). An example of gene_cl.ref can be found at the data folder of this repository. If you have the output of FindAllMarkers from Seurat, you could use the cut_markers function to get gene_cl.ref by the command ``` cut_markers(levels(output.seu$cluster),output.seu,ntop=100) ```.


![Scheme](matchSCore2_Overview.png)


```{r,eval=FALSE}

library(matchSCore2)
library(nnet)
library(Matrix)

### Training of the model  
mod <- train_model(scale.data = ref@scale.data,clus = ref@ident,gene_cl.ref = gene_cl.ref,prop = 0.75)

## Cell projection
out <- identity_map(scale.data = obs@scale.data,model = mod,gene_cl.ref)

### cell identities
ids <- out$ids 


## To each cell probabilities are assigned for any possible identity class
probabilities <- out$fit.prob


```
MatchSCore2 can be used to improve the clustering annotation, detecting subtle differences between cell types. 
![Annotations](Clustering_vs_matchSCore_annotations.png)


## Clustering Annotation

We can also measure the grade of matching across clusters from two considered datasets.
For example, you could use the top 100 ranked markers we got (by using Seurat) from the Chromium and MARS-Seq dataset.
You can load files directly from the data folder in this repository. 

```{r,eval=TRUE}

## We used the Chromium data as reference
load(file="data/gene_cl.ref_Chromium.RData")

## And the MARS-Seq as test
load(file="data/gene_cl_MARS-Seq.RData")


## The matchSCore2 function computes the clustering comparison and produce the heatmap table with Jaccard Indexes for each group combination

ms <- matchSCore2(gene_cl.ref = gene_cl.ref,gene_cl.obs = gene_cl,ylab = "Chromium",xlab = "MARS-Seq")

## The matchSCore heatmap is stored in the ggplot slot of ms. 
ms$ggplot



```
![Heatmap](Heatmap.png)





