oe.lof<-read.csv("/Users/Joaco/Desktop/Joaco/Medicina/UOC/TFM/TFM/gnomad.v2.1.1.lof_metrics.by_gene.txt",sep = "")
PVS1.lof<-read.csv("/Users/Joaco/Desktop/Joaco/Medicina/UOC/TFM/TFM/PVS1.LOF.genes.hg19.csv", header = F)
miss.aa<-read.csv("/Users/Joaco/Desktop/Joaco/Medicina/UOC/TFM/TFM/missense_aa.csv",sep = ",")
uniprot<-read.csv("/Users/Joaco/Desktop/Joaco/Medicina/UOC/TFM/TFM/protein.domain.csv",sep = "\t")
#rmsk<-read.csv("/Users/Joaco/Desktop/Joaco/Medicina/UOC/TFM/TFM/rmsk.csv", sep="\t")
BP1.gene<-read.csv("/Users/Joaco/Desktop/Joaco/Medicina/UOC/TFM/TFM/BP1.genes.hg19.csv", sep="", header=F)
gwas.clean<-read.csv("/Users/Joaco/Desktop/Joaco/Medicina/UOC/TFM/TFM/gwas.clean.csv")
omim.clean<-read.csv("/Users/Joaco/Desktop/Joaco/Medicina/UOC/TFM/TFM/omim.clean.csv")

ensembl = useEnsembl(biomart='ensembl', 
                     dataset="hsapiens_gene_ensembl",GRCh=37) 
cadd.pvs<-c("CANONICAL_SPLICE","STOP_GAINED","STOP_LOST")
review<-c("criteria provided, multiple submitters, no conflicts",
          "criteria provided, single submitter","reviewed by expert panel")