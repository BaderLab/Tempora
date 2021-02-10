Tempora: cell trajectory inference using time-series single-cell RNA
sequencing data
================

  - [Introduction](#introduction)
  - [Usage](#usage)
      - [Installation](#installation)
      - [Sample data](#sample-data)
      - [Input data](#input-data)
      - [Calculate clusters’ pathway enrichment
        profiles](#calculate-clusters-pathway-enrichment-profiles)
      - [Build and visualize
        trajectory](#build-and-visualize-trajectory)
      - [Identify temporally dependent
        pathways](#identify-temporally-dependent-pathways)

## Introduction

Tempora is a novel cell trajectory inference method that orders cells
using time information from time-series scRNAseq data. Tempora uses
biological pathway information to help identify cell type relationships
and can identify important time-dependent pathways to help interpret the
inferred trajectory.

## Usage

### Installation

You can install Tempora using devtools:

``` r
# install devtools
install.packages("devtools")

# install Tempora
devtools::install_github("BaderLab/Tempora")

library(Tempora)
```

### Sample data

The Tempora package was validated using three datasets: an *in vitro*
differentiation of human skeletal muscle myoblasts, *in vivo* early
development of murine cerebral cortex and *in vivo* embryonic and
postnatal development of murine cerebellum. These processed datasets can
be accessed on the Bader Lab website at
<https://www.baderlab.org/Software/Tempora>.

The MouseCortex dataset will be used in this vignette as an example.

### Input data

Tempora takes processed scRNAseq data as input, either as a gene
expression matrix with separate time and cluster labels for all cells,
or a Seurat or SingleCellExperiment object containing gene expression
data and a clustering result. Tempora does not implement clustering or
batch effect correction as part of its pipeline and assumes that the
user has input a well-annotated cluster solution free of batch effect
into the method.

``` r
#install Seurat package when using the MouseCortex data.
if (!require('Seurat')) {
  install.packages('Seurat')
} 

#Load MouseCortex sample data
load("MouseCortex.RData")
```

We can the import the Seurat object containing the murine cerebral
cortex development data into a Tempora object to start the analysis.
Here, as the clusters have been manually annotated prior to running
Tempora, a vector of cluster label is given to the function. If you have
yet to annotate your clusters but have a list of marker genes for
expected cell types in the data, you can input the list of marker genes
to this function to run automated cluster labeling with GSVA.

``` r
#Import MouseCortex data 
#As this is a Seurat v2 object, set assayType to ""
#See ?Tempora::ImportSeuratObject for additional arguments to import Seurat v3 or SingleCellExperiment obbjects
cortex_tempora <- ImportSeuratObject(MouseCortex, clusters = "res.0.6",
                                     timepoints = "Time_points", 
                                     assayType = "",
                                     cluster_labels = c("Neurons","Young neurons","APs/RPs",
                                                        "IPs","APs/RPs", "Young neurons", "IPs"),
                                     timepoint_order = c("e11", "e13", "e15", "e17"))
```

From the specified clustering result, Tempora will automatically
calculate the temporal score of each cluster, which is based on its
composition of cells from each timepoint. This information will be
stored in the *cluster.metadata* slot of the Tempora object.

### Calculate clusters’ pathway enrichment profiles

Next, the pathway enrichment profiles of the clusters are calculated
using
[GSVA](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-14-7)
and stored in the *cluster.pathways* slot of the Tempora object. The
default pathway gene set database Tempora uses is the Bader Lab pathway
gene set database without electronic annotation Gene Ontology terms,
which can be accessed on the [Bader Lab
website](http://download.baderlab.org/EM_Genesets/current_release/).

To automatically pull the latest version of the gmt file:
```{r download baderlab gmt file, message=FALSE, warning=FALSE}
if (!require('RCurl')) {
  install.package('RCurl')
} 
gmt_url = "http://download.baderlab.org/EM_Genesets/current_release/Mouse/symbol/"

#list all the files on the server
filenames = getURL(gmt_url)
tc = textConnection(filenames)
contents = readLines(tc)
close(tc)

#get the gmt that has all the pathways and does not include terms inferred from electronic annotations(IEA)
#start with gmt file that has pathways only
rx = gregexpr("(?<=<a href=\")(.*.GOBP_AllPathways_no_GO_iea.*.)(.gmt)(?=\">)",
  contents, perl = TRUE)

gmt_file = unlist(regmatches(contents, rx))
dest_gmt_file <- file.path(getwd(),gmt_file )
download.file(
    paste(gmt_url,gmt_file,sep=""),
    destfile=dest_gmt_file
)
```


This function also performs principal component analysis (PCA) on the
clusters pathway enrichment profiles to remove redundancy due to
overrepresentation of certain pathways in the database. The PCA result
is stored in the *cluster.pathways.dr* slot. Tempora also outputs a
scree plot to help users identify the number of principal components
(PCs) to be used in downstream trajectory construction.

``` r
#Estimate pathway enrichment profiles of clusters
cortex_tempora <- CalculatePWProfiles(cortex_tempora, 
                gmt_path = gmt_file,
                method="gsva", min.sz = 5, max.sz = 200, parallel.sz = 1)
```

### Build and visualize trajectory

We can now build the trajectory based on the clusters’ pathway
enrichment profiles. Tempora employs the mutual information (MI) rank
and data processing inequality approach implemented in
[ARACNE](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-7-S1-S7)
to calculate MI between all cluster pairs present in the data as well as
remove edges with weak MIs. The trajectory is stored as a dataframe of
edge lists in the *trajectory* slot. Tempora then assigns directions to
all edges in the network so that edges point from clusters with low
temporal scores to clusters with high temporal scores.

``` r
#Build trajectory with 6 PCs 
cortex_tempora <- BuildTrajectory(cortex_tempora, n_pcs = 6, difference_threshold = 0.01)
```

After building the trajectory, we can visualize it as a network, with
the piechart at each node representing the composition of cells
collected at different time points in the experiment and the arrow
connecting each pair of nodes representing lineage relationship between
them.

``` r
#Visualize the trajectory
cortex_tempora <- PlotTrajectory(cortex_tempora)
```

This function will add a slot *layouts* containing the x and y
coordinates of all nodes, determined using the Sugiyama layered graph
drawing algorithm.

### Identify temporally dependent pathways

Finally, we can use Tempora to investigate time-dependent pathways.
Tempora fits a generalized additive model to the data to identify
pathways whose expressions change over the temporal axis. The results of
this analysis is stored in the *varying.pws* slot of the Tempora object.

``` r
#Fit GAMs on pathway enrichment profile
cortex_tempora <- IdentifyVaryingPWs(cortex_tempora, pval_threshold = 0.05)

#Plot expression trends of significant time-varying pathways
PlotVaryingPWs(cortex_tempora)
```
