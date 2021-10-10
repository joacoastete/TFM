library(myvariant)
library(biomaRt)
library(knitr)

# Carga los dataset necesarios
oe.lof<-read.csv("gnomad.v2.1.1.lof_metrics.by_gene.txt",sep = "")
PVS1.lof<-read.csv("PVS1.LOF.genes.hg19.csv", header = F)
miss.aa<-read.csv("missense_aa.csv",sep = ",")
uniprot<-read.csv("protein.domain.csv",sep = "\t")
rmsk<-read.csv("rmsk.csv", sep="\t")
BP1.gene<-read.csv("BP1.genes.hg19.csv", sep="", header=F)
gwas.clean<-read.csv("gwas.clean.csv")
omim.clean<-read.csv("omim.clean.csv")

ensembl = useEnsembl(biomart='ensembl', 
                     dataset="hsapiens_gene_ensembl",GRCh=37) 
cadd.pvs<-c("CANONICAL_SPLICE","STOP_GAINED","STOP_LOST")
review<-c("criteria provided, multiple submitters, no conflicts",
          "criteria provided, single submitter","reviewed by expert panel")

# Valores de corte
oe.cutoff=0.35
mis.cutoff=3.09
mis.z<-oe.lof[oe.lof$mis_z>mis.cutoff,] 
af.cutoff=0.01
prior.prob<-0.1
odds.path<-350

# Carga los el archivo con los ids de las variantes
validation<-read.csv("validation.csv", sep=";")
valid<-as.vector(validation$id)

#Genera el objeto con los datos de las variantes
variants<-getVariants(valid, fields = "all", return.as = "records")

# Genera el data.frame vacío
names<-c("ID","Gene", "Clasification","Post_P","PVS1",
         "PS1","PS2","PS3", "PS4",
         "PM1","PM2", "PM3", "PM4","PM5","PM6", 
         "PP1", "PP2", "PP3","PP4","PP5",
         "BA1","BS1","BS2", "BS3","BS4",
         "BP1", "BP2","BP3","BP4","BP5","BP6", "BP7",
         "PS","PM","PP","BS","BP", 
         "odd.PVS", "odd.PS","odd.PM","odd.PP", "odd.BS","odd.BP", "comb.odd")

df<-as.data.frame(matrix(c(rep(NA,4),rep(0,33),rep(NA,7)), length(variants), length(names), byrow = T))
colnames(df)<-names

for(var in 1:length(variants)){
  # Completa los datos de la variante
  df$ID[var]<-variants[[var]]$`_id`
  if(is.null(variants[[var]]$hg19)==F){
    positions <- data.frame(chromosome =variants[[var]]$chrom,
                            start = variants[[var]]$hg19$start,
                            end = variants[[var]]$hg19$end)
    results <- getBM(attributes = c("hgnc_symbol"), 
                     filters = c("chromosome_name", "start", "end"),
                     values = list(positions[,1], positions[,2], positions[,3]),
                     mart = ensembl)
    df$Gene[var]<-results[1]
  }else if(is.null(variants[[var]]$cadd$gene$genename)==F){
    df$Gene[var]<-variants[[var]]$cadd$gene$genename
  }else if(is.null(variants[[var]]$snpeff$ann[[1]]$gene_id)==F){
    df$Gene[var]<-variants[[var]]$snpeff$ann[[1]]$gene_id
  }else if(is.null(variants[[var]]$dbnsfp$genename)==F){
    df$Gene[var]<-variants[[var]]$dbnsfp$genename
  }
  
  #BA1 busca en las distintas base de datos si la frecuencia alélica es >0.05
  
  
  if((is.null(variants[[var]]$cadd$esp$af)|is.null(variants[[var]]$gnomad_genome$af$af)|is.null(variants[[var]]$gnomad_exome$af$af)|
      is.null(variants[[var]]$dbnsfp$`1000gp3`$af)|is.null(variants[[var]]$cadd$`1000g`$af))==F){
    if(variants[[var]]$cadd$esp$af|variants[[var]]$gnomad_genome$af$af|variants[[var]]$gnomad_exome$af$af|
       variants[[var]]$dbnsfp$`1000gp3`$af|variants[[var]]$cadd$`1000g`$af>0.05){
      df$BA1[var]=1
    }
  }
  
  
  
  # PVS1. Si el gen está en la lista de LoF y si es una variante nula
  
  
  if (is.null(variants[[var]]$cadd$consequence)==F){
    if(is.na(df$Gene[var])==F){
      if(length(grep(df$Gene[var],PVS1.lof,value=F))>0){
        if(length(grep(variants[[var]]$cadd$consequence,cadd.pvs,value=F))>0){
          df$PVS1[var]=1
        }
      }
    }
  } else if(is.null(variants[[var]]$snpeff$ann)==F){
    if(length(variants[[var]]$snpeff$ann$effect)==0){
      if(is.na(df$Gene[var])==F){
        if(length(grep(df$Gene[var],PVS1.lof[,1],value=F))>0){
          if(variants[[var]]$snpeff$ann[[1]]$effect=="frameshift_variant"){
            df$PVS1[var]=1
          }
        }
      }
    }else if(length(grep(df$Gene[var],PVS1.lof[,1],value=F))>0){
      if(variants[[var]]$snpeff$ann$effect=="frameshift_variant"){
        df$PVS1[var]=1
      }
    }
  }
  
  
  # PS1 & PM5. Si la variante es missense, mira el tipo de cambio de aminoacido
  
  if(is.null(variants[[var]]$cadd$consdetail)==F){
    if(variants[[var]]$cadd$consdetail=="missense"){
      for(p in 1:nrow(miss.aa)){
        if(is.null(variants[[var]]$cgi)==F){
          if(length(variants[[var]]$cgi$protein_change)==0){
            if(variants[[var]]$cgi[[1]]$protein_change==miss.aa$join[p]){
              df$PS1[var]=1
              df$PM5[var]=0
              if (df$PS1[var]==1){
                break
              }
            }
          }else if(variants[[var]]$cgi$protein_change==miss.aa$join[p]){
            df$PS1[var]=1
            df$PM5[var]=0
            if (df$PS1[var]==1){
              break
            }
          }
        }
        if(is.null(variants[[var]]$dbnsfp$hgvsp)==F & is.na(df$Gene[var])==F){
          if((df$Gene[var]==miss.aa$gene[p])==T& 
             (variants[[var]]$dbnsfp$hgvsp[1]==miss.aa$hgvs.p[p])==F&
             (variants[[var]]$dbnsfp$hgvsp[1]==miss.aa$aa.one[p])==F){
            df$PS1[var]=0
            df$PM5[var]=1
          }
          if((df$Gene[var]==miss.aa$gene[p])==T& 
             (variants[[var]]$dbnsfp$hgvsp[1]==miss.aa$hgvs.p[p]|
              variants[[var]]$dbnsfp$hgvsp[1]==miss.aa$aa.one[p])==T){
            df$PS1[var]=1
            df$PM5[var]=0
            if (df$PS1[var]==1){
              break
            }
            if(is.null(variants[[var]]$dbnsfp$mutpred$aa_change)==F & is.na(df$Gene[var])==F){
              if((df$Gene[var]==miss.aa$gene[p])==T& 
                 (variants[[var]]$dbnsfp$mutpred$aa_change[1]==miss.aa$one.clean[p])==F)
                df$PS1[var]=0
              df$PM5[var]=1
            }
            if((df$Gene[var]==miss.aa$gene[o])==T& 
               (variants[[var]]$dbnsfp$mutpred$aa_change[1]==miss.aa$one.clean[o])==T){
              df$PS1[var]=1
              df$PM5[var]=0
              if (df$PS1[var]==1){
                break
              }
              if(is.null(variants[[var]]$snpeff$ann$hgvs_p)==F & is.na(df$Gene[var])==F){
                if((variants[[var]]$df$Gene[var]==miss.aa$gene[o])==T& 
                   (variants[[var]]$snpeff$ann$hgvs_p==miss.aa$one.clean[o])==F)
                  df$PS1[var]=0
                df$PM5[var]=1
              }
              if((df$Gene[var]==miss.aa$gene[o])==T& 
                 (variants[[var]]$snpeff$ann$hgvs_p==miss.aa$one.clean[o])==T){
                df$PS1[var]=1
                df$PM5[var]=0
                if (df$PS1[var]==1){
                  break
                }
              }
            }
          }
          
        } 
      }
    }
  }
  
  # BS1 Si la frecuencia alelica es mayor a la frec PopMax de Gnomad
  
  
  if(is.null(variants[[var]]$gnomad_exome$af$af_popmax)==F&is.null(variants[[var]]$gnomad_exome$af$af)==F){
    if(variants[[var]]$gnomad_exome$af$af>variants[[var]]$gnomad_exome$af$af_popmax){
      df$BS1[var]=1
    }
  }else if(is.null(variants[[var]]$gnomad_exome$af$af_popmax)==T&is.null(variants[[var]]$gnomad_exome$af$af)==F){
    if(variants[[var]]$gnomad_exome$af$af>af.cutoff){
      df$BS1[var]=1
    }
  }
  
  
  
  # PM2. Si la variante no aparece en ninguna de las basese poblacionales
  
  
  if((is.null(variants[[var]]$cadd$esp$af)&is.null(variants[[var]]$gnomad_genome$af$af)&is.null(variants[[var]]$gnomad_exome$af$af)|
      is.null(variants[[var]]$dbnsfp$`1000gp3`$af)&is.null(variants[[var]]$cadd$`1000g`$af))==T){
    df$PM2[var]=1
  }
  
  #BS2. Si se observa la variante en algun control sano según el tipo de patologia
  
  if(df$PM2[var]==0){
    if (is.null(variants[[var]]$gnomad_exome)==F){
      if(is.na(df$Gene[var])==F){
        for(o in 1:nrow(omim.clean)){
          if(df$Gene[var]==omim.clean$gene[o]){
            if(omim.clean$inheritance[o]=="AR"&variants[[var]]$gnomad_exome$hom$hom>3){
              df$BS2[var]=1
            }else if(omim.clean$inheritance[o]=="AD"&variants[[var]]$gnomad_exome$ac$ac>5){
              df$BS2[var]=1
            }else if(omim.clean$inheritance[o]=="XL"&variants[[var]]$gnomad_exome$ac$ac_male>3){
              df$BS2[var]=1
            }
          }
        }
      } 
    } else if (is.null(variants[[var]]$gnomad_genome)==F){
      if(is.na(df$Gene[var])==F){
        for(o in 1:nrow(omim.clean)){
          if(df$Gene[var]==omim.clean$gene[o]){
            if(omim.clean$inheritance[o]=="AR"&variants[[var]]$gnomad_genome$hom$hom>3){
              df$BS2[var]=1
            }else if(omim.clean$inheritance[o]=="AD"&variants[[var]]$gnomad_genome$ac$ac>5){
              df$BS2[var]=1
            }else if(omim.clean$inheritance[o]=="XL"&variants[[var]]$gnomad_genome$ac$ac_male>3){
              df$BS2[var]=1
            }
          }
        }
      }
    }
  }
  
  
  
  # PS4. Si la variante se encuentra en el listado de OR>0.5 de la base de datos GWAS
  
  
  
  if(is.null(variants[[var]]$dbsnp$rsid)==F){
    if(length(grep(variants[[var]]$dbsnp$rsid, gwas.clean, value=F))>0){
      df$PS4[var]=1
    }
  }else if(is.null(variants[[var]]$dbnsfp$rsid)==F){
    if(length(grep(variants[[var]]$dbnsfp$rsid, gwas.clean, value=F))>0){
      df$PS4[var]=1
    }
  }else if(is.null(variants[[var]]$clinvar$rsid)==F){
    if(length(grep(variants[[var]]$clinvar$rsid, gwas.clean, value=F))>0){
      df$PS4[var]=1
    }
  }
  
  
  # PM1. Para variantes missense que alteran el dominio proteico ya conocido de patogenicidad
  
  
  if(is.null(variants[[var]]$cadd$consdetail)==F&is.null(variants[[var]]$dbnsfp$interpro_domain)==F&
     is.na(df$Gene[var])==F){
    if(variants[[var]]$cadd$consdetail=="missense"){
      for(dom in 1:nrow(uniprot)){
        if(df$Gene[var]==uniprot$Gene[dom]&(length(grep(variants[[var]]$dbnsfp$interpro_domain,uniprot$Interpro_domain[dom], value=F))>0)){
          df$PM1[var]=1
        }
      }
    }
  }
  
  
  # PM4 BP3. Segun si la variante se afecta un región de repetición
  
  
  if(is.null(variants[[var]]$cadd$consequence)==F){
    if(variants[[var]]$cadd$consequence=="STOP_LOST"){
      df$PM1[var]=1
    }
  }
  
  if(is.null(variants[[var]]$dbsnp$vartype)==F& is.null(variants[[var]]$snpeff$ann)==F){
    if(length(variants[[var]]$snpeff$ann$effect)==0){
      if(variants[[var]]$dbsnp$vartype=="delins"&
         isFALSE(variants[[var]]$snpeff$ann[[1]]$effect=="frameshift_variant")==T){
        name<-paste0("chr", variants[[var]]$chrom, sep = "")
        for(r in 1:nrow(rmsk[rmsk$genoName==name,])){
          if(variants[[var]]$hg19$start>=rmsk[rmsk$genoName==name,]$genoStart&
             variants[[var]]$hg19$start<rmsk[rmsk$genoName==name,]$genoEnd){
            df$BP3[var]=1
          }else{
            df$PM1[var]=1
          }
        }
      }
    }else if(variants[[var]]$dbsnp$vartype=="delins"&isFALSE(variants[[var]]$snpeff$ann$effect=="frameshift_variant")==T){
      name<-paste0("chr", variants[[var]]$chrom, sep = "")
      for(r in 1:nrow(rmsk[rmsk$genoName==name,])){
        if(variants[[var]]$hg19$start>=rmsk[rmsk$genoName==name,]$genoStart&
           variants[[var]]$hg19$start<rmsk[rmsk$genoName==name,]$genoEnd){
          df$BP3[var]=1
        }else{
          df$PM1[var]=1
        }
      }
    }
  }
  
  # PP3 BP4. Según la predicciones de distintos algoritmos
  
  pred.s<-c()
  if(is.null(variants[[var]]$dbnsfp$provean$pred[1])==F){
    if (variants[[var]]$dbnsfp$provean$pred=="N"){
      pred.s<-c(pred.s,0)
    }else{pred.s<-c(pred.s,1)}
  }
  if(is.null(variants[[var]]$cadd$sift$cat)==F){
    if (variants[[var]]$cadd$sift$cat=="tolerated"){
      pred.s<-c(pred.s,0)
    }else{pred.s<-c(pred.s,1)}
  }
  if(is.null(variants[[var]]$dbnsfp$fathmm$pred)==F){
    if (variants[[var]]$dbnsfp$fathmm$pred[1]=="T"){
      pred.s<-c(pred.s,0)
    }else{pred.s<-c(pred.s,1)}
  }
  if(is.null(variants[[var]]$cadd$polyphen$cat)==F){
    if (variants[[var]]$cadd$polyphen$cat=="benign"){
      pred.s<-c(pred.s,0)
    }else{pred.s<-c(pred.s,1)}
  }
  if(is.null(variants[[var]]$cadd$consscore)==F){
    if (variants[[var]]$cadd$consscore<20){
      pred.s<-c(pred.s,0)
    }else{pred.s<-c(pred.s,1)}
  } 
  if(is.null(variants[[var]]$dbnsfp$vest4$rankscore)==F){
    if (variants[[var]]$dbnsfp$vest4$rankscore<0.5){
      pred.s<-c(pred.s,0)
    }else{pred.s<-c(pred.s,1)}
  }
  
  if(is.null(variants[[var]]$cadd$gerp$rs)==F){
    if (variants[[var]]$cadd$gerp$rs>2){
      pred.s<-c(pred.s,0)
    }else{pred.s<-c(pred.s,1)}
  }
  
  
  
  if(is.null(pred.s)==F){
    pred.tot<-sum(pred.s)/length(pred.s)
    if(pred.tot>=0.5){
      df$PP3[var]=1
    }else if(pred.tot<0.5){
      df$BP4[var]=1}
  }
  
  pred.s<-c()
  
  
  #PP2. Si la variante missense afecta un gen suceptible de missense
  
  
  if(is.null(variants[[var]]$cadd$consdetail)==F){
    if(variants[[var]]$cadd$consdetail=="missense"){
      if(is.na(df$Gene[var])==F){
        if(length(grep(df$Gene[var],mis.z[,1],value=F))>0){
          df$PP2[var]=1
        }
      }
    }
  }
  
  
  #BP1. Si la variante missense afecta un gen no susceptible de missense
  
  if(df$PP2[var]==0){
  if(is.null(variants[[var]]$cadd$consdetail)==F){
    if(variants[[var]]$cadd$consdetail=="missense"){
      if(is.na(df$Gene[var])==F){
        if(length(grep(df$Gene[var],BP1.gene[,1],value=F))>0){
          df$BP1[var]=1
        }
      }
    }
  }
  }
  
  
  #PP5 BP6. Si la variante ya ha sido clasificada con evidencia
  
  rcv.clasi<-c()
  rcv.review<-c()
  
  if(is.null(variants[[var]]$clinvar$rcv)==F){
    if(length(variants[[var]]$clinvar$rcv$clinical_significance)==0){
      for(i in 1:length(variants[[var]]$clinvar$rcv)){
        rcv.clasi<-c(rcv.clasi,variants[[var]]$clinvar$rcv[[i]]$clinical_significance)
      }
    }else {rcv.clasi<-variants[[var]]$clinvar$rcv$clinical_significance}
    if(length(grep("Pathogenic",rcv.clasi, value = F))>0){
      if(length(variants[[var]]$clinvar$rcv$review_status)==0){
        for(i in 1:length(variants[[var]]$clinvar$rcv)){
          rcv.review<-c(rcv.review,variants[[var]]$clinvar$rcv[[i]]$review_status)
        }
        if(length(grep(review[1], rcv.review))>0){
          df$PP5[var]=1
        }else if(length(grep(review[2], rcv.review))>0){
          df$PP5[var]=1
        }else if(length(grep(review[3], rcv.review))>0){
          df$PP5[var]=1
        }
      }else if(length(grep(review[1], variants[[var]]$clinvar$rcv$review_status))>0){
        df$PP5[var]=1
      }else if(length(grep(review[2], variants[[var]]$clinvar$rcv$review_status))>0){
        df$PP5[var]=1
      }else if(length(grep(review[3], variants[[var]]$clinvar$rcv$review_status))>0){
        df$PP5[var]=1
      }
    }else if(length(grep("Benign",rcv.clasi, value = F))>0){
      if(length(variants[[var]]$clinvar$rcv$review_status)==0){
        for(i in 1:length(variants[[var]]$clinvar$rcv)){
          rcv.review<-c(rcv.review,variants[[var]]$clinvar$rcv[[i]]$review_status)
        }
        if(length(grep(review[1], rcv.review))>0){
          df$BP6[var]=1
        }else if(length(grep(review[2], rcv.review))>0){
          df$BP6[var]=1
        }else if(length(grep(review[3], rcv.review))>0){
          df$BP6[var]=1
        }
      }else if(length(grep(review[1], variants[[var]]$clinvar$rcv$review_status))>0){
        df$BP6[var]=1
      }else if(length(grep(review[2], variants[[var]]$clinvar$rcv$review_status))>0){
        df$BP6[var]=1
      }else if(length(grep(review[3], variants[[var]]$clinvar$rcv$review_status))>0){
        df$BP6[var]=1
      }
    }
  }
  
  # Calcula la Probabilidad Post
  
  df$PS[var]<-sum(df[var,c("PS1","PS2","PS3","PS4")])
  df$PM[var]<-sum(df[var,c("PM1","PM2","PM3","PM4", "PM5", "PM6")])
  df$PP[var]<-sum(df[var,c("PP1","PP2","PP3","PP4", "PP5")])
  df$BS[var]<-sum(df[var,c("BS1","BS2","BS3","BS4")])
  df$BP[var]<-sum(df[var,c("BP1","BP2","BP3","BP4", "BP5", "BP6", "BP7")])
  df$odd.PVS[var]<-350^df$PVS1[var]
  df$odd.PS[var]<-18.7^df$PS[var]
  df$odd.PM[var]<-4.33^df$PM[var]
  df$odd.PP[var]<-2.08^df$PP[var]
  df$odd.BS[var]<-18.7^((-1)*df$BS[var])
  df$odd.BP[var]<-2.08^((-1)*df$BP[var])
  df$comb.odd[var]<-df$odd.PVS[var]*df$odd.PS[var]*df$PM[var]*df$odd.PP[var]*
    df$odd.BS[var]*df$odd.BP[var]
  df$Post_P[var]<-round((df$comb.odd[var]*prior.prob)/((df$comb.odd[var]-1)*prior.prob+1),3)
  
  # Clseifica según la probabilidad post
  if(df$BA1[var]==1){
    df$Clasification[var]<-"Benign"
  }else if(df$Post_P[var]>0.99){
    df$Clasification[var]<-"Pathogenic"
  }else if(df$Post_P[var]>0.90&df$Post_P[var]<=0.99){
    df$Clasification[var]<-"Likely Pathogenic"
  }else if(df$Post_P[var]>=0.10&df$Post_P[var]<=0.90){
    df$Clasification[var]<-"Uncertain significance"
  }else if(df$Post_P[var]>=0.001&df$Post_P[var]<=0.10){
    df$Clasification[var]<-"Likely Benign"
  }else if(df$Post_P[var]<0.001){
    df$Clasification[var]<-"Benign"
  }
}


df
comp<-cbind(df$ID,df$Clasification,df$Post_P,validation$clasification)
colnames(comp)<-c("ID", "Clasification", "Post_P", "Clasification Clinvar")
kable(comp)

