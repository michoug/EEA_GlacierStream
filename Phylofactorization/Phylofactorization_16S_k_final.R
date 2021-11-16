library(phylofactor)
library(ggplot2)
library(ggpubr)


X <- read.csv('my_metadata_table.csv') %>% as.data.table
Data <- read.table('my_SV_table.txt',row.names = 1, header=T) %>% as.matrix()
tree <- read.tree('my_tree.nwk')
taxonomy <- read.csv('my_taxonomy_table.csv') %>% as.data.table()
taxonomy[,taxonomy:=paste(Kingdom,Phylum,Class,Order,Family,Genus,Species,sep='; ')]
colnames(taxonomy)[1] <- 'id'
taxonomy <- taxonomy[,c('id','taxonomy')]


###### Remove samples with 10K or fewer reads and asvs with 6 or fewer sequences 
Data <- Data[,colSums(Data)>1e4]
Data <- Data[,grepl('_sed',colnames(Data))]
Data <- Data[rowSums(Data)>0,]
pa <- rowSums(Data>0)

tb <- cumsum(rev(table(rowSums(Data>0))))
Data <- Data[pa>6,]


##### Align rows of Data with tree tip labels and taxonomy
tree <- keep.tip(tree,rownames(Data))
Data <- Data[tree$tip.label,]
taxonomy <- taxonomy[id %in% tree$tip.label]
all(rownames(Data) %in% tree$tip.label) #should return TRUE, similarly for taxonomy matches

##### Align columns of data with Metadata codes
X <- X[match(colnames(Data),code),]
all(colnames(Data) %in% X$code) #check, should return TRUE

## data + meta-data + tree aligned. Let's phylofactor!

############## Phylofactorization

#X[,logK:=log10(K)] #in our specific example we log-transformed decomposition rates (k) before analysis for normality reasons
#X[,logChla:=log10(chla)] #the same for chl-a concentrations

#The objective funtion for phylofactorization, with expedition considered
K3_objective <- function(y,X,PF.output=FALSE,...){
  dataset <- cbind('Data'=y,X)
  gg <- mgcv::gam(logK~ expedition + Data,data=dataset,...)
  if (PF.output){
    return(gg)
    break
  } else {
    output <- NULL
    s <- summary(gg)
    output$objective <- abs(s$p.table['Data','t value'])
    output$stopStatistics <- s$p.table['Data','Pr(>|t|)']
    return(output)
  }
}

load.mgcv <- 'library(mgcv)'

#The actual phylofactorization happens here
PF.K3 <- PhyloFactor(Data,tree,X,choice.fcn=K3_objective,
                     cluster.depends = load.mgcv,nfactors=500,ncores=31)

save(list=ls(),file='data/K_phylofactorization_16S.Rd')


######### Refining information on phylofactor object, extracting p-values and correcting
get_Pvals <- function(n=1,PF.K.=PF.K3) summary(PF.K3$custom.output[[n]])$p.table['Data','Pr(>|t|)']
PF.K3$factors$Pval <- sapply(1:PF.K3$nfactors,get_Pvals)

holm_adjust <- function(p,PF.K.=PF.K3){
  n_tips <- length(PF.K3$tree$tip.label)
  n_edges_tested <- 2*n_tips-1 - 2*(1:PF.K3$nfactors)
  p.adjust(p,method = 'holm',n = n_edges_tested) %>% return
}

PF.K3$factors$Holm_FWER <- holm_adjust(PF.K3$factors$Pval)
sum(PF.K$factors$Holm_FWER<0.05) #the number of significant phylogenetic factors after Holm correction

#To refine further, we pulled out a list of factors with > 50 ASVs each and FWER < 0.01. 
##View(PF.K3$factors) ##Example of visual examination, we are interested in Factors where Group 1 contains ~50 or more SVs and Holm_FWER column is less than 0.01
#We then explored their taxonomies with pf.taxa()
##pf.taxa(PF.K3, taxonomy = taxonomy, factor=7) #example where we pull the consensus taxonomies of factor 7. Again, we are interested in Group 1 (Group 2 is the outgroup)

#Then we plotted the selected non-nested factors. We first pruned the extra long (length>0.5) edges of the tree for visualization purposes
#PF.K3_tree <- pf.tree(PF.K3, tree=tree, factors = c(2,7,11,13,30,31,40,49,51,64,69,77,98,100,135,150,152,168,434)) ##Our specific example
#PF.K3_tree_outgroup <- pf.tree(PF.K3, tree=tree, factors = 7, groups = 2) ##This additionally plots the big negatively correlated clade
#pp_K3 <- PF.K3_tree$ggplot ##The overall plot object
#pp_K3_outgroup <- PF.K3_tree_outgroup$ggplot #A similar object that just marks the big outgroup

#Then we just copied the outgroup graph to the overall tree and we manually painted and annotated the groups' taxonomies
