
# Load all the requeried files
oe.lof<-read.csv("gnomad.v2.1.1.lof_metrics.by_gene.txt",sep = "")
PVS1.lof<-read.csv("PVS1.LOF.genes.hg19.csv", header = F)
miss.aa<-read.csv("missense_aa.csv",sep = ",")
uniprot<-read.csv("protein.domain.csv",sep = "\t")
rmsk<-read.csv("rmsk_clean.csv", sep=",")
BP1.gene<-read.csv("BP1.genes.hg19.csv", sep="", header=F)
gwas.clean<-read.csv("gwas.clean.csv")
omim.clean<-read.csv("omim.clean.csv")

ensembl = useEnsembl(biomart='ensembl', 
                     dataset="hsapiens_gene_ensembl",GRCh=37) 
cadd.pvs<-c("CANONICAL_SPLICE","STOP_GAINED","STOP_LOST")
review<-c("criteria provided, multiple submitters, no conflicts",
          "criteria provided, single submitter","reviewed by expert panel")
