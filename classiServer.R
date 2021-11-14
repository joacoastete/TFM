first.table<-function(inpt,mis.cut,af.cut,PS2,PS3,PM3,BP2,PP1,PP4,BP5, BP7){

variants<-getVariants(inpt, fields = "all", return.as = "records")
oe.cutoff=0.35
mis.cutoff=mis.cut
mis.z<-oe.lof[oe.lof$mis_z>mis.cutoff,] 
af.cutoff=af.cut
prior.prob<-0.1
odds.path<-350

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
    df
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
df$PS2[var]<-ifelse(PS2=="Paternity confirmed",1,0)
df$PM6[var]<-ifelse(PS2=="Paternity non confirmed",1,0)
df$PS3[var]<-ifelse(PS3=="Studies supportive of a damaging effect",1,0)
df$BS3[var]<-ifelse(PS3=="Studies shows no damaging effect",1,0)
df$PM3[var]<-ifelse(PM3=="Yes",1,0)
df$BP2[var]<-ifelse(BP2=="Yes",1,0)
df$PP1[var]<-ifelse(PP1=="Co-segregation with disease in multiple affected family members",1,0)
df$BS4[var]<-ifelse(PP1=="Lack of segregation in affected members of a family",1,0)
df$PP4[var]<-ifelse(PP4=="Yes",1,0)
df$BP5[var]<-ifelse(BP5=="Yes",1,0)
df$BP7[var]<-ifelse(BP7=="Yes",1,0)
}
df
}  

  # Calcula la Probabilidad Post
classification.func<-function(tab1,pat.pred,ben.pred,pat.new,ben.new){
  pat.pred<-as.vector(pat.pred)
  ben.pred<-as.vector(ben.pred)
  pat.new<-as.vector(pat.new)
  ben.new<-as.vector(ben.new)
  if(isTRUE(pat.pred==pat.new & ben.pred==ben.new)){
    df2<-tab1
  }else{
  names<-c("ID","Gene", "Clasification","Post_P","PVS1",
           "PS1","PS2","PS3", "PS4",
           "PM1","PM2", "PM3", "PM4","PM5","PM6", 
           "PP1", "PP2", "PP3","PP4","PP5",
           "BA1","BS1","BS2", "BS3","BS4",
           "BP1", "BP2","BP3","BP4","BP5","BP6", "BP7",
           "PS","PM","PP","BS","BP", 
           "odd.PVS", "odd.PS","odd.PM","odd.PP", "odd.BS","odd.BP", "comb.odd")
  
  df2<-as.data.frame(matrix(c(rep(NA,4),rep(0,33),rep(NA,7)), nrow(tab1), length(names), byrow = T))
  colnames(df2)<-names
  df2$ID<-tab1$ID
  df2$Gene<-tab1$Gene
  df2$PVS1<-ifelse(length(grep("PVS1",pat.new,value=F))>0,1,0)
  df2$PS1<-ifelse(length(grep("PS1",pat.new,value=F))>0,1,0)
  df2$PS2<-ifelse(length(grep("PS2",pat.new,value=F))>0,1,0)
  df2$PS3<-ifelse(length(grep("PS3",pat.new,value=F))>0,1,0)
  df2$PS4<-ifelse(length(grep("PS4",pat.new,value=F))>0,1,0)
  df2$PM1<-ifelse(length(grep("PM1",pat.new,value=F))>0,1,0)
  df2$PM2<-ifelse(length(grep("PM2",pat.new,value=F))>0,1,0)
  df2$PM3<-ifelse(length(grep("PM3",pat.new,value=F))>0,1,0)
  df2$PM4<-ifelse(length(grep("PM4",pat.new,value=F))>0,1,0)
  df2$PM5<-ifelse(length(grep("PM5",pat.new,value=F))>0,1,0)
  df2$PM6<-ifelse(length(grep("PM6",pat.new,value=F))>0,1,0)
  df2$PP1<-ifelse(length(grep("PP1",pat.new,value=F))>0,1,0)
  df2$PP2<-ifelse(length(grep("PP2",pat.new,value=F))>0,1,0)
  df2$PP3<-ifelse(length(grep("PP3",pat.new,value=F))>0,1,0)
  df2$PP4<-ifelse(length(grep("PP4",pat.new,value=F))>0,1,0)
  df2$PP5<-ifelse(length(grep("PP5",pat.new,value=F))>0,1,0)
  df2$PP6<-ifelse(length(grep("PP6",pat.new,value=F))>0,1,0)
  df2$BA1<-ifelse(length(grep("BA1",ben.new,value=F))>0,1,0)
  df2$BS1<-ifelse(length(grep("BS1",ben.new,value=F))>0,1,0)
  df2$BS2<-ifelse(length(grep("BS2",ben.new,value=F))>0,1,0)
  df2$BS3<-ifelse(length(grep("BS3",ben.new,value=F))>0,1,0)
  df2$BS4<-ifelse(length(grep("BS4",ben.new,value=F))>0,1,0)
  df2$BP1<-ifelse(length(grep("BP1",ben.new,value=F))>0,1,0)
  df2$BP2<-ifelse(length(grep("BP2",ben.new,value=F))>0,1,0)
  df2$BP3<-ifelse(length(grep("BP3",ben.new,value=F))>0,1,0)
  df2$BP4<-ifelse(length(grep("BP4",ben.new,value=F))>0,1,0)
  df2$BP5<-ifelse(length(grep("BP5",ben.new,value=F))>0,1,0)
  df2$BP6<-ifelse(length(grep("BP6",ben.new,value=F))>0,1,0)
  df2$BP7<-ifelse(length(grep("BP7",ben.new,value=F))>0,1,0)
  }
  for(var in 1:nrow(tab1)){
    
  df2$PS[var]<-sum(df2[var,c("PS1","PS2","PS3","PS4")])
  df2$PM[var]<-sum(df2[var,c("PM1","PM2","PM3","PM4", "PM5", "PM6")])
  df2$PP[var]<-sum(df2[var,c("PP1","PP2","PP3","PP4", "PP5")])
  df2$BS[var]<-sum(df2[var,c("BS1","BS2","BS3","BS4")])
  df2$BP[var]<-sum(df2[var,c("BP1","BP2","BP3","BP4", "BP5", "BP6", "BP7")])
  df2$odd.PVS[var]<-350^df2$PVS1[var]
  df2$odd.PS[var]<-18.7^df2$PS[var]
  df2$odd.PM[var]<-4.33^df2$PM[var]
  df2$odd.PP[var]<-2.08^df2$PP[var]
  df2$odd.BS[var]<-18.7^((-1)*df2$BS[var])
  df2$odd.BP[var]<-2.08^((-1)*df2$BP[var])
  df2$comb.odd[var]<-df2$odd.PVS[var]*df2$odd.PS[var]*df2$PM[var]*df2$odd.PP[var]*
    df2$odd.BS[var]*df2$odd.BP[var]
  df2$Post_P[var]<-round((df2$comb.odd[var]*prior.prob)/((df2$comb.odd[var]-1)*prior.prob+1),3)
  
  # Clseifica según la probabilidad post
  if(df2$BA1[var]==1){
    df2$Clasification[var]<-"Benign"
  }else if(df2$Post_P[var]>0.99){
    df2$Clasification[var]<-"Pathogenic"
  }else if(df2$Post_P[var]>0.90&df2$Post_P[var]<=0.99){
    df2$Clasification[var]<-"Likely Pathogenic"
  }else if(df2$Post_P[var]>=0.10&df2$Post_P[var]<=0.90){
    df2$Clasification[var]<-"Uncertain significance"
  }else if(df2$Post_P[var]>=0.001&df2$Post_P[var]<=0.10){
    df2$Clasification[var]<-"Likely Benign"
  }else if(df2$Post_P[var]<0.001){
    df2$Clasification[var]<-"Benign"
  }
  }
  df2
}

group.pat<-c("PVS1",
             "PS1","PS2","PS3", "PS4",
             "PM1","PM2", "PM3", "PM4","PM5","PM6", 
             "PP1", "PP2", "PP3","PP4","PP5")
check.pat<-function(tab1){
  choice.pat<-c()
  if(tab1$PVS1==1){choice.pat<-c(choice.pat,"PVS1")}
  if(tab1$PS1==1){choice.pat<-c(choice.pat,"PS1")}
  if(tab1$PS2==1){choice.pat<-c(choice.pat,"PS2")}
  if(tab1$PS3==1){choice.pat<-c(choice.pat,"PS3")}
  if(tab1$PS4==1){choice.pat<-c(choice.pat,"PS4")}
  if(tab1$PM1==1){choice.pat<-c(choice.pat,"PM1")}
  if(tab1$PM2==1){choice.pat<-c(choice.pat,"PM2")}
  if(tab1$PM3==1){choice.pat<-c(choice.pat,"PM3")}
  if(tab1$PM4==1){choice.pat<-c(choice.pat,"PM4")}
  if(tab1$PM5==1){choice.pat<-c(choice.pat,"PM5")}
  if(tab1$PM6==1){choice.pat<-c(choice.pat,"PM6")}
  if(tab1$PP1==1){choice.pat<-c(choice.pat,"PP1")}
  if(tab1$PP2==1){choice.pat<-c(choice.pat,"PP2")}
  if(tab1$PP3==1){choice.pat<-c(choice.pat,"PP3")}
  if(tab1$PP4==1){choice.pat<-c(choice.pat,"PP4")}
  if(tab1$PP5==1){choice.pat<-c(choice.pat,"PP5")}
  choice.pat
}
group.ben<-c("BA1","BS1","BS2", "BS3","BS4",
  "BP1", "BP2","BP3","BP4","BP5","BP6", "BP7")
check.ben<-function(tab1){
  choice.ben<-c()
  if(tab1$BA1==1){choice.ben<-c(choice.ben,"BA1")}
  if(tab1$BS1==1){choice.ben<-c(choice.ben,"BS1")}
  if(tab1$BS2==1){choice.ben<-c(choice.ben,"BS2")}
  if(tab1$BS3==1){choice.ben<-c(choice.ben,"BS3")}
  if(tab1$BS4==1){choice.ben<-c(choice.ben,"BS4")}
  if(tab1$BP1==1){choice.ben<-c(choice.ben,"BP1")}
  if(tab1$BP2==1){choice.ben<-c(choice.ben,"BP2")}
  if(tab1$BP3==1){choice.ben<-c(choice.ben,"BP3")}
  if(tab1$BP4==1){choice.ben<-c(choice.ben,"BP4")}
  if(tab1$BP5==1){choice.ben<-c(choice.ben,"BP5")}
  if(tab1$BP6==1){choice.ben<-c(choice.ben,"BP6")}
  if(tab1$BP7==1){choice.ben<-c(choice.ben,"BP7")}
  choice.ben
}
