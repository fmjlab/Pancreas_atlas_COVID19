# The following code allows for the analysis of 6 single cell RNAseq datasets of the human pancreas
# Information on these datasets can be found in the following locations:
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE81076
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE85241
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE86469
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84133
# https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-5061/
# This code was written by Fahd Qadir PhD. on 06/03/2020 email: mqadir@tulane.edu

# 1. installation and loading of packages
# Devtools
install.packages('devtools')
library(devtools)

# Seuratdata
devtools::install_github('satijalab/seurat-data')

# Seurat wrappers
devtools::install_github('satijalab/seurat-wrappers')

# Load packages
library(Seurat)
library(ggplot2)
library(patchwork)
library(SeuratData)
library(SeuratWrappers)
library(future)

# Set RAM to 50GB
# options(future.globals.maxSize = 40 * 1024^3)

# check the current active plan
# plan()

# change the current plan to access parallelization
# future::availableCores()
# future::availableWorkers()
# plan("multiprocess", workers = 15)
# plan()

# Loading of refrence datasets
#GSE81076 <- read.csv("C:/Users/mqadir/Box/Lab 2301/sCell Analysis Project/Refrence Human Pancreas/scRNAseq datasets/pancreas/GSE81076.csv", header = TRUE, sep = ",", row.names = 1)
GSE85241 <- read.csv("C:/Users/mqadir/Box/Lab 2301/sCell Analysis Project/Refrence Human Pancreas/scRNAseq datasets/pancreas/GSE85241.csv", header = TRUE, sep = ",", row.names = 1)
GSE86469 <- read.csv("C:/Users/mqadir/Box/Lab 2301/sCell Analysis Project/Refrence Human Pancreas/scRNAseq datasets/pancreas/GSE86469.csv", header = TRUE, sep = ",", row.names = 1)
GSE84133 <- read.csv("C:/Users/mqadir/Box/Lab 2301/sCell Analysis Project/Refrence Human Pancreas/scRNAseq datasets/pancreas/GSE84133.csv", header = TRUE, sep = ",", row.names = 1)
EMTAB5061 <- read.csv("C:/Users/mqadir/Box/Lab 2301/sCell Analysis Project/Refrence Human Pancreas/scRNAseq datasets/pancreas/EMTAB5061.csv", header = TRUE, sep = ",", row.names = 1)
GSE131886 <- read.csv("C:/Users/mqadir/Box/Lab 2301/sCell Analysis Project/Refrence Human Pancreas/scRNAseq datasets/pancreas/GSE131886.csv", header = TRUE, sep = ",", row.names = 1)

# Create Seurat objects
#GSE81076 <- CreateSeuratObject(counts = GSE81076, project = "SeuratProject", assay = "RNA", min.cells = 3, min.features = 200)
GSE85241 <- CreateSeuratObject(counts = GSE85241, project = "SeuratProject", assay = "RNA", min.cells = 3, min.features = 200)
GSE86469 <- CreateSeuratObject(counts = GSE86469, project = "SeuratProject", assay = "RNA", min.cells = 3, min.features = 200)
GSE84133 <- CreateSeuratObject(counts = GSE84133, project = "SeuratProject", assay = "RNA", min.cells = 3, min.features = 200)
EMTAB5061 <- CreateSeuratObject(counts = EMTAB5061, project = "SeuratProject", assay = "RNA", min.cells = 3, min.features = 200)
GSE131886 <- CreateSeuratObject(counts = GSE131886, project = "SeuratProject", assay = "RNA", min.cells = 3, min.features = 200)
# Load in Luca's data
adult_pancreas <- readRDS("C:/Users/mqadir/Box/Lab 2301/sCell Analysis Project/Refrence Human Pancreas/scRNAseq datasets/pancreas/adult_pancreas.rds")
chronic_pancreatitis <- readRDS("C:/Users/mqadir/Box/Lab 2301/sCell Analysis Project/Refrence Human Pancreas/scRNAseq datasets/pancreas/chronic_pancreatitis.rds")
neonatal_pancreas <- readRDS("C:/Users/mqadir/Box/Lab 2301/sCell Analysis Project/Refrence Human Pancreas/scRNAseq datasets/pancreas/neonatal_pancreas.rds")

# Sample specific Metadata addition
#GSE81076$sample <- "GSE81076"
GSE85241$sample <- "GSE85241"
GSE86469$sample <- "GSE86469"
GSE84133$sample <- "GSE84133"
EMTAB5061$sample <- "EMTAB5061"
GSE131886$sample <- "GSE131886"
adult_pancreas$sample <- "EGAS00001004653_adult"
chronic_pancreatitis$sample <- "EGAS00001004653_CP"
neonatal_pancreas$sample <- "EGAS00001004653_NP"

# Sex segregation specific Metadata addition
# For GSE85241
levels(GSE85241)
male <- c("D28.1", "D28.2", "D28.3", "D28.4", "D28.5", "D28.6", "D28.7", "D28.8",
          "D29.1", "D29.2", "D29.3", "D29.4", "D29.5", "D29.6", "D29.7", "D29.8",
          "D31.1", "D31.2", "D31.3", "D31.4", "D31.5", "D31.6", "D31.7", "D31.8")
female <- c("D30.1", "D30.2", "D30.3", "D30.4", "D30.5", "D30.6", "D30.7", "D30.8")
GSE85241@meta.data$sex[GSE85241@meta.data$orig.ident %in% male] <- "male"
GSE85241@meta.data$sex[GSE85241@meta.data$orig.ident %in% female] <- "female"

# For EMTAB5061
levels(EMTAB5061)
male <- c("AZ", "HP1502401", "HP1504101T2D", "HP1504901", "HP1507101", "HP1509101", "HP152301T2D")
female <- c("HP1506401", "HP1508501T2D", "HP1526901T2D")
EMTAB5061@meta.data$sex[EMTAB5061@meta.data$orig.ident %in% male] <- "male"
EMTAB5061@meta.data$sex[EMTAB5061@meta.data$orig.ident %in% female] <- "female"

# fOR GSE131886
levels(GSE131886)
male <- c("HPD3")
female <- c("HPD1", "HPD2")
GSE131886@meta.data$sex[GSE131886@meta.data$orig.ident %in% male] <- "male"
GSE131886@meta.data$sex[GSE131886@meta.data$orig.ident %in% female] <- "female"

# fOR GSE84133
levels(GSE84133)
male <- c("m1", "m3")
female <- c("f2", "f4")
GSE84133@meta.data$sex[GSE84133@meta.data$orig.ident %in% male] <- "male"
GSE84133@meta.data$sex[GSE84133@meta.data$orig.ident %in% female] <- "female"

# fOR GSE86469
levels(GSE86469)
male <- c("H1", "H2", "H3", "H4", "H6", "H7", "H8")
female <- c("H5", "H9", "H10", "H11", "H12", "H13")
GSE86469@meta.data$sex[GSE86469@meta.data$orig.ident %in% male] <- "male"
GSE86469@meta.data$sex[GSE86469@meta.data$orig.ident %in% female] <- "female"

# Ref-dataset specific Metadata addition
#GSE81076$ref <- "ref"
GSE85241$ref <- "ref"
GSE86469$ref <- "ref"
GSE84133$ref <- "ref"
EMTAB5061$ref <- "ref"
GSE131886$ref <- "ref"
adult_pancreas$ref <- "ref"
chronic_pancreatitis$ref <- "ref"
neonatal_pancreas$ref <- "ref"

#Subset out to only save male and female
Idents(pancreas.integrated) <- "sex"
pancreas.integrated <- subset(pancreas.integrated, idents = c("male", "female"))

# Create a list of datasets containing seurat objects
pancreas.list <- list(#"GSE81076" = GSE81076, 
  "GSE85241" =GSE85241, "GSE86469" = GSE86469, 
  "GSE84133" = GSE84133, "EMTAB5061" = EMTAB5061, "GSE131886" = GSE131886, "EGAS00001004653_adults" = adult_pancreas,
  "EGAS00001004653_CP" = chronic_pancreatitis, "EGAS00001004653_NP" = neonatal_pancreas)
#,"panc_sex_cau_m1" = panc_sex_cau_m1, "panc_sex_cau_f1" = panc_sex_cau_f1)

pancreas.list
pancreas.list <- lapply(X = pancreas.list, FUN = function(x) {
  x <- NormalizeData(x, verbose = TRUE)
  x <- FindVariableFeatures(x, verbose = TRUE)
})

features <- SelectIntegrationFeatures(object.list = pancreas.list)
pancreas.list <- lapply(X = pancreas.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

anchors <- FindIntegrationAnchors(object.list = pancreas.list, reference = c(6, 7), reduction = "rpca", 
                                  dims = 1:50)
pancreas.integrated <- IntegrateData(anchorset = anchors, dims = 1:50)

pancreas.integrated <- ScaleData(pancreas.integrated, verbose = TRUE)
pancreas.integrated <- RunPCA(pancreas.integrated, verbose = TRUE)
pancreas.integrated <- RunUMAP(pancreas.integrated, dims = 1:50)

DimPlot(pancreas.integrated, group.by = "sample")
DimPlot(pancreas.integratedx, group.by = "sex")

# Remove NAs
pancreas.integratedx <- subset(pancreas.integrated, subset = sex != "NA")
pancreas.integrated <- pancreas.integratedx

# Normalize based on RNA
pancreas.integrated <- NormalizeData(pancreas.integrated, normalization.method = "LogNormalize", assay = "RNA", scale.factor = 1e4, 
                                     verbose = TRUE)

#Clustering
pancreas.integrated <- FindNeighbors(pancreas.integrated, dims = 1:30)
pancreas.integrated <- FindClusters(pancreas.integrated, resolution = 1.2)

# For UMAP visualization
DefaultAssay(object = pancreas.integrated) <- "RNA"
FeaturePlot(object = pancreas.integrated, 
            features = c("ADRB1"),
            pt.size = 1,
            cols = c("darkgrey", "red"),
            min.cutoff = 0,
            max.cutoff = 20,
            order = TRUE)

# Visualization Clustering
plots <- DimPlot(pancreas.integrated, group.by = c("ref", "sample"))
plots & theme(legend.position = "right") & guides(color = guide_legend(nrow = 14, byrow = TRUE,
                                                                       override.aes = list(size = 5)))
Idents(pancreas.integrated) <- "CellType"
DimPlot(pancreas.integrated, label = TRUE)

# Organize clusters
Idents(pancreas.integrated) <- "seurat_clusters"
plot <- DimPlot(pancreas.integrated, reduction = "umap")
pancreas.integrated <- CellSelector(plot = plot, object = pancreas.integrated, ident = "Beta")
pancreas.integrated <- CellSelector(plot = plot, object = pancreas.integrated, ident = "Alpha")
pancreas.integrated <- CellSelector(plot = plot, object = pancreas.integrated, ident = "Delta")
pancreas.integrated <- CellSelector(plot = plot, object = pancreas.integrated, ident = "Epsilon")
pancreas.integrated <- CellSelector(plot = plot, object = pancreas.integrated, ident = "Gamma")
pancreas.integrated <- CellSelector(plot = plot, object = pancreas.integrated, ident = "Ductal")
pancreas.integrated <- CellSelector(plot = plot, object = pancreas.integrated, ident = "Acinar")
pancreas.integrated <- CellSelector(plot = plot, object = pancreas.integrated, ident = "Ducto-Acinar")
pancreas.integrated <- CellSelector(plot = plot, object = pancreas.integrated, ident = "Ducto-Endocrine")
pancreas.integrated <- CellSelector(plot = plot, object = pancreas.integrated, ident = "Unclassified-Endocrine")
pancreas.integrated <- CellSelector(plot = plot, object = pancreas.integrated, ident = "Bcells")
pancreas.integrated <- CellSelector(plot = plot, object = pancreas.integrated, ident = "Macrophage")
pancreas.integrated <- CellSelector(plot = plot, object = pancreas.integrated, ident = "Tcells")
pancreas.integrated <- CellSelector(plot = plot, object = pancreas.integrated, ident = "Tuftcells")
pancreas.integrated <- CellSelector(plot = plot, object = pancreas.integrated, ident = "Endothelial")
pancreas.integrated <- CellSelector(plot = plot, object = pancreas.integrated, ident = "Quiescent stellate")
pancreas.integrated <- CellSelector(plot = plot, object = pancreas.integrated, ident = "Activated stellate")
pancreas.integrated <- CellSelector(plot = plot, object = pancreas.integrated, ident = "Schwann")
pancreas.integrated <- CellSelector(plot = plot, object = pancreas.integrated, ident = "Mast")
levels(pancreas.integrated)

# Saving this information in the metadata slot
head(Idents(pancreas.integrated))
pancreas.integrated$CellType <- Idents(pancreas.integrated)
head(pancreas.integrated@meta.data)

# Run find variable features again running this is questionable, as only the var features from integrated data is useful
# But Seurat recommends re-running this
DefaultAssay(object = pancreas.integrated) <- "RNA"
pancreas.integrated <- FindVariableFeatures(pancreas.integrated, selection.method = "vst", nfeatures = 3000)

# Define an order of cluster identities remember after this step-
# cluster re-assignment occurs, which re-assigns clustering in my_levels
my_levels <- c("Beta", "Alpha", "Delta", "Gamma", "Epsilon", 
               "Ductal", "Acinar", "Quiescent stellate", "Activated stellate", 
               "Schwann", "Endothelial", "Macrophage", "Mast", "Tcells", "Bcells",
               "Tuftcells")
head(pancreas.integrated@meta.data$CellType)

# Re-level object@meta.data this just orders the actual metadata slot, so when you pull its already ordered
pancreas.integrated@meta.data$CellType <- factor(x = pancreas.integrated@meta.data$CellType, levels = my_levels)
DimPlot(pancreas.integrated)

#Save Object
saveRDS(pancreas.integrated, "C:/Users/mqadir/Box/Lab 2301/sCell Analysis Project/Refrence Human Pancreas/scRNAseq datasets/Workspace/pancreas.integrated.rds")
pancreas.integrated <- readRDS("C:/Users/mqadir/Box/Lab 2301/sCell Analysis Project/Refrence Human Pancreas/scRNAseq datasets/Workspace/pancreas.integrated.rds")

# Subsetting Our cells out
sex <- subset(pancreas.integrated, subset = ref == "panc_sex")
DimPlot(sex)

# Check metadata
head(pancreas.integrated@meta.data)
table(pancreas.integrated$sample)
table(Idents(pancreas.integrated))

# Check activeidents
head(Idents(pancreas.integrated))

# Change active idents to CellType
Idents(pancreas.integrated) <- "sex"

# For UMAP visualization
DefaultAssay(object = pancreas.integrated) <- "RNA"
FeaturePlot(object = pancreas.integrated, 
            features = c("PGR"),
            pt.size = 1,
            cols = c("darkgrey", "red"),
            min.cutoff = 0,
            max.cutoff = 20,
            order = FALSE)

# Visualize information
table(pancreas.integrated$sample)
DefaultAssay(object = pancreas.integrated) <- "RNA"
VlnPlot(pancreas.integrated, c("PGR"), group.by = "CellType", split.by = "sex", assay = "RNA", slot = "data", ncol = 1, pt.size = 1)

# Average expression of all cells within a cluster
males <- subset(pancreas.integrated, subset = (sex == "male"))
females <- subset(pancreas.integrated, subset = (sex == "female"))
Idents(female) <- "CellType"
Idents(males) <- "CellType"
cluster.averages.males <- AverageExpression(males)
cluster.averages.females <- AverageExpression(females)
head(cluster.averages.males[["RNA"]])
head(cluster.averages.females[["RNA"]])
cluster.averages.males[["RNA"]][c("PGR"),]
cluster.averages.females[["RNA"]][c("PGR"),]

# Issue 371
# Subset your cluster of interest for as an example I am subsetting a cluster called 'beta'
# The following creates a seurat object of only the cluster 'beta'
betacells <- subset(pancreas.integrated, subset = (CellType == c("Beta")) & (sex == "female") & (sample == "EGAS00001004653_CP"))
#betacells <- subset(pancreas.integrated, subset = (CellType == c("Beta")) & (sex == "female"))
betacells <- subset(pancreas.integrated, subset = (CellType == c("Alpha")) & (sample == "EGAS00001004653_CP"))

# Point your new cluster towards the object you will use to perform calculations.
# I like doing this because otherwise, you have to write lengths of redundant code
# Also I'm really lazy
ThisWayIsTotallyMentalButItWorks <- betacells
GOI1 <- 'ACE2' #you will have to name your first gene here, im choosing PDX1 as an example
GOI2 <- 'TMPRSS2' #you will have to name your first gene here, im choosing INS as an example
GOI1.cutoff <- .1
GOI2.cutoff <- .1

# Enjoy!
GOI1.cells <- length(which(FetchData(ThisWayIsTotallyMentalButItWorks, vars = GOI1) > GOI1.cutoff))
GOI2.cells <- length(which(FetchData(ThisWayIsTotallyMentalButItWorks, vars = GOI2) > GOI2.cutoff))
GOI1_GOI2.cells <- length(which(FetchData(ThisWayIsTotallyMentalButItWorks, vars = GOI2) > GOI2.cutoff & FetchData(ThisWayIsTotallyMentalButItWorks, vars = GOI1) > GOI1.cutoff))
all.cells.incluster <- table(ThisWayIsTotallyMentalButItWorks@active.ident)
GOI1.cells/all.cells.incluster*100 # Percentage of cells in Beta that express GOI1
GOI2.cells/all.cells.incluster*100 #Percentage of cells in Beta that express GOI2
GOI1_GOI2.cells/all.cells.incluster*100 #Percentage of cells in Beta that co-express GOI1 + GOI2

# Some cool code for total percentage (need to x100)
betacells <- subset(pancreas.integrated, subset = (sample == "EGAS00001004653_CP"))

PrctCellExpringGene <- function(object, genes, group.by = "all"){
  if(group.by == "all"){
    prct = unlist(lapply(genes,calc_helper, object=object))
    result = data.frame(Markers = genes, Cell_proportion = prct)
    return(result)
  }

  else{
    list = SplitObject(object, group.by)
    factors = names(list)

    results = lapply(list, PrctCellExpringGene, genes=genes)
    for(i in 1:length(factors)){
      results[[i]]$Feature = factors[i]
    }
    combined = do.call("rbind", results)
    return(combined)
  }
}

calc_helper <- function(object,genes){
  counts = object[['RNA']]@counts
  ncells = ncol(counts)
  if(genes %in% row.names(counts)){
    sum(counts[genes,]>0)/ncells
  }else{return(NA)}
}

PrctCellExpringGene(betacells, c("ACE2", "TMPRSS2"), group.by = "CellType")
calc_helper(pancreas.integrated, c("ACE2", "TMPRSS2"))

# Plotting one gene on a dimplot
betacells <- subset(pancreas.integrated, subset = (sex == "female"))
betacells <- subset(pancreas.integrated, subset = (sex == "female"))
FeaturePlot(object = betacells, 
            features = c("ACE2"),
            pt.size = 1,
            cols = c("darkgrey", "red"),
            min.cutoff = 0,
            max.cutoff = 3,
            order = TRUE)

# Set cell identity to sample identity so that you can extraxt cell type information for plotting
Idents(object = pancreas.integrated) <- pancreas.integrated@meta.data$celltype

# How can I extract expression matrix for all beta cells
betacells <- subset(pancreas.integrated, idents = c("Beta"))

# Violin plot
DefaultAssay(object = betacells) <- "RNA"
VlnPlot(object = betacells, features = c("ACE2", "TMPRSS2"), group.by = "sample", slot = "data")

# How can I extract expression matrix for all beta cells
alphacells <- subset(pancreas.integrated, idents = c("alpha"))

# Violin plot
DefaultAssay(object = alphacells) <- "RNA"
Idents(pancreas.integrated) <- "sex"
VlnPlot(object = pancreas.integrated, features = c("XIST"), group.by = "sample", split.by = "sex", slot = "data")

# Set cell identity to sample identity
Idents(object = pancreas.integrated) <- pancreas.integrated@meta.data$celltype

# Find if SRD genes are differentially expressed
beta.integrated.markers <- FindAllMarkers(object = pancreas.integrated, slot = 'data', test.use = 'wilcox')

# How can I calculate the average expression of all cells within a cluster?
cluster.averages <- AverageExpression(pancreas.integrated, assay= "RNA", slot = "data")
head(cluster.averages[["RNA"]][c("ACE2", "TMPRSS2"), 1:14])
