
## Quick introduction

MatchSCore2 is R package for the comparison of single cell RNA-seq data across tools and experiments. The package allows a gene marker-based projection of single cells onto a reference sample and, thus, the identification of cell types in unknown cells.  

## Installation

matchSCore2 is a new version of the older matchSCore package. 
This is a development version running on R (>= 3.5.1). The package makes use of nnet, Matrix, splatter, MixSim and ggplot2 libraries. Other libraries are also used for clustering and identification of markers (e.g. Seurat).

matchSCore2 can be installed by devtools:

```{r,eval=FALSE}

library("devtools")
install_github('elimereu/matchSCore2')


```

The purpose of this vignette is to guide you to the use of the matchSCore2 with some examples.


## Functionalities of matchSCore2

MatchSCore2 has tree main functions:

1. Annotation of cells throught a supervised classification model based on a reference dataset. (**Cell annotation**) 
2. Measure and visualize the matching between your clusters and clusters from a reference dataset. (**Cluster annotation**) 
3. Track the accuracy trend of a tool in clustering and marker identification compared with the optimal solution provided by a simulated data set. In this case, matchSCore2 works in combination with the Splatter package - https://github.com/Oshlack/splatter/blob/master/vignettes/splatter.Rmd - (**Benchmarking**).

```{r,eval=FALSE}

library(matchSCore2)

### Training of the model by using a reference dataset (ref) where cell types (ref@ident) and their markers are known 

mod <- train_model(scale.data = ref@scale.data,clus = ref@ident,gene_cl.ref = gene_cl.ref,prop = 0.75)

## Predict cell identities in a test data (obs)

out <- identity_map(scale.data = obs@scale.data,model = mod,gene_cl.ref)

### cell identities
ids <- out$ids 


## To each cell probabilities are assigned for any possible identity class
probabilities <- out$fit.prob


```



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

ms <- matchSCore2(gene_cl.ref = gene_cl.ref,gene_cl.obs = gene_cl.obs,ylab = "Chromium",xlab = "MARS-Seq")

## The matchSCore heatmap is stored in the ggplot slot of ms. 
ms$ggplot



```
![Heatmap](Heatmap.png)


## Benchmarking by using simulated data


```{r,eval=FALSE}

library(matchSCore)
library(splatter)
library(scater)
library(MixSim)
## To plot
library(ggplot2) 
library(ggpubr)


## For example we create a simulated data set with 4 groups of equal proportions. 
sim <- splatSimulate(batchCells=rep(500,2),group.prob=rep(0.25,4),method = "groups")
n_groups <- 4

## Compute the ranking of genes per group
rank_df <- rank_sim(sim)

## Set the proportion of top ranked genes (specificity) we want to output
specificity=0.2
## feature info data frame
fd <- rowData(sim)
##List of markers per group
markers_pos <- markers_by_specificity(rank_df,specificity,n_groups) ## by position in fd$Gene

## List with the top 10% ranked markers per group
sim.markers <- lapply(markers_pos,function(x) fd$Gene[x])
sim.markers

```

For example, we could use Seurat to analyze this data set by running the function seurat_run

```{r,eval=FALSE}

##Run Seurat
out_seu <-seurat_run(sim,ntop=100,out_seu = NULL,res = NULL,dims.use = NULL,test.de = "wilcox")

pd=colData(sim)
lab.sim <- as.integer(factor(pd$Group))
idc <- out_seu$clusters ## Seurat clusters
gene_cl <- out_seu$gene_cl ## Seurat cluster genes

## identify the group labels (defined by the simulated groups) for the clusters 
lab <- compute_labels(lab.sim,idc)


RandIndex(lab.sim,idc) ## we can compute the Rand Index (RI), adjusted Rand Index (ARI) or F score

## Now, we can compute the matchSCore value 
matchSCore(markers_pos,gene_cl,lab)


## If you want to increase the number of ntop to see if the matchSCore is higher, you can run the seurat_run function, but passing its old output to be faster
gene_cl=seurat_run(sim,ntop=500,out_seu,res = NULL,dims.use = NULL,test.de = "wilcox")

## And then re-compute the matchSCore 
matchSCore(markers_pos,gene_cl,lab)


## Let's plot the matchSCore trend by using different values of top ranked markers (ntop) and specificity (proportion of true group markers). 

ntop=seq(250,2000,250)
## The score function compute a matchSCore at all different values of ntop.
ms=scores(sim,out_seu,ntop,seurat_run,lab)

## Now we can plot the matchSCore curve

scores_df=data.frame(ms)
tscores_plots(scores_df)


```



