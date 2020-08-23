#######################################################################################################################
### INSTALL & LOAD REQUIRED PACKAGES
#######################################################################################################################

#install the package maSigPro (for gene expression differences in time course)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("maSigPro")
#load the packge maSigPro
library("maSigPro")
library("mclust")

######################################################################################################################
### LOAD THE GENE EXPRESSION DATA FROM CSV FILE
######################################################################################################################
gene_exp_data <- read.csv("out_v3.csv") #read the csv file 
gene_exp_data_2 <- gene_exp_data[,-1] #remove column "X"
gene_exp_data_4 <- data.frame(gene_exp_data_2[,-1], row.names=gene_exp_data_2$Id) #read the file as a dataframe and make the gene name as a row name

dd <- data.frame(t(gene_exp_pt$influ.info))#extract influential genes and remove those rows
dd2 <- tibble::rownames_to_column(dd, "id")

gene_exp_data_5 <- gene_exp_data_4[!(as.character(gene_exp_data_2$Id) %in% as.character(dd2$id)), ]

######################################################################################################################
### MAKE AN EXPERIMENTAL DESIGN FILE (showing conditions A & B, replicates A_1 to A_12 & B_1 to B_12)
######################################################################################################################
design_gene_expression <- read.csv("design_gene_expression.csv") #make an experimental design file with time series and conditionA vs. conditionB
design_gene_expression_1 <- data.frame(design_gene_expression[,-1], row.names=design_gene_expression$X) #read the file as a dataframe and make the row names from A_1 to B_12
design_gene_expression_2 <- as.matrix(data.frame(design_gene_expression_1))


#####################################################################################################################
### 1) MAKE A DESIGN MATRIX FOR DEFINING THE MODEL
#####################################################################################################################
design_exp <- make.design.matrix(design_gene_expression_2, degree=2, time.col=1, repl.col=2, group.cols=c(3,4)) #generate the design matrix
#degree of polynomial is 2, sincere there are 3 time points (quintic regression model)
a <- design_exp$groups.vector
plot <- plot(a) #plot the design matrix (shows a regression model)

############################################1st regression step
##############################################################################################################
### 2) FIND STATISTICALLY SIGNIFICANT GENE MODELS + THEIR REGRESSION COEFFICIENTS
##############################################################################################################
gene_exp_p <- p.vector(gene_exp_data_5, design_exp, counts=FALSE, MT.adjust="BH", Q=0.05, min.obs=23, family=NULL) #compute a regression fit & p-value associated to the F-Statistic for each gene & adjust the p-value for multiple comparison Benjamini-Hochberg
help(p.vector)
#counts are FALSE (since we have normalised the data)
gene_exp_p$i #no of significant genes
gene_exp_p$p.adjusted #give adjusted p-values
gene_exp_p$SELEC #list of significant genes with their p-values

############################################2nd regression step
#####################################################################################################################
### 3) FIND BEST FIT OF THE VARIABLES INTO THE GENE MODELS
#####################################################################################################################
gene_exp_pt <- T.fit(gene_exp_p)


####################################################################################################################
### 4) GET A LIST OF GENES (for which variables best fit the model)
####################################################################################################################
get <- get.siggenes(gene_exp_pt, rsq=0.95, vars="groups")
get$summary #see the genes in clusters
dfs <- data.frame(get$summary)
hhs <- dfs$ConditionvsControl
write.csv(hhs, "outputkkk.csv")
#cut-off value for R^2 is set to 0.7 (i.e. 70% of differential gene expression could be explained by condition A vs. condition B)
#vars (variables) are set of be grouped by condition A vs. Condition B


##############################################################################################################
### 5) VISUALISE SIGNIFICANT GENES AS CLUSTERS
##############################################################################################################

h <- see.genes(get$sig.genes$ConditionvsControl, cluster.method="hclust")#cluster analysis of significant genes
gg <- names(which(h$cut==10)) #extract the genes from a specific cluster
gg
