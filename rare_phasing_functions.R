
add_autism_trios_ancestry = function(data){
  fam = read.delim("ASC_v13_JAK8_sexFixed_probands_noDup_removeMissingFams.fam", stringsAsFactors = F,header=F)
  colnames(fam) = c("sa.fam.famID", "Sample", "PatID", "MatID","Sex","Pheno")
  ancestry = fam %>% 
    left_join(pca_complete %>%
                select(Sample,ancPCA,ancAdmix)) %>% 
    select(sa.fam.famID,ancPCA,ancAdmix) %>%
    filter(!duplicated(sa.fam.famID))
  return(data %>% left_join(ancestry))
}


getSummary <-function(x){
  return(x %>% 
           summarise(NtrioHap=sum(sameTrioHap)+sum(diffTrioHap),
                     nSameHapTrio= sum(sameTrioHap),
                     nDiffHapTrio = sum(diffTrioHap),
                     percSameTrioHap = sum(sameTrioHap) / NtrioHap,
                     #                    NcoInExAC= sum(coInExAC),
                     #                     NnotCoInExAC = sum(notCoInExAC),
                     #                     nSameHapTrioAndExAC = sum(sameHapTrioAndExAC),
                     #                    nDiffHapTrioAndExAC = sum(diffHapTrioAndExAC),
                     #                    percCorrectSegExAC = (nSameHapTrioAndExAC+nDiffHapTrioAndExAC) / (NcoInExAC+NnotCoInExAC),
                     #                     percMissSeg = 1-((NcoInExAC+NnotCoInExAC)/NtrioHap),
                     nSameHapExACEM = sum(sameHapExACEM),
                     nDiffHapExACEM = sum(diffHapExACEM),
                     nSameHapTrioAndExACEM = sum(sameHapTrioAndExACEM),
                     nDiffHapTrioAndExACEM = sum(diffHapTrioAndExACEM),
                     percCorrectEM = (nSameHapTrioAndExACEM+nDiffHapTrioAndExACEM)/(nSameHapExACEM+nDiffHapExACEM),
                     percCorrectSameHapEM = nSameHapTrioAndExACEM / nSameHapExACEM,
                     percCorrectDiffHapEM = nDiffHapTrioAndExACEM / nDiffHapExACEM,
                     percMissEM = 1-((nSameHapExACEM+nDiffHapExACEM)/NtrioHap)
                     
           ))
}

getSummary2 <-function(x){
  return(x %>% 
           #mutate(single_het_ratio = (AaBB + AABb) / (AaBB + AABb + AaBb)) %>%
           summarise(NtrioHap=n(),
                     nSameHapTrio= sum(same_trio_haplotype),
                     nDiffHapTrio = NtrioHap - nSameHapTrio,
                     percSameTrioHap = nSameHapTrio / NtrioHap,
                     #                    NcoInExAC= sum(coInExAC),
                     #                     NnotCoInExAC = sum(notCoInExAC),
                     #                     nSameHapTrioAndExAC = sum(sameHapTrioAndExAC),
                     #                    nDiffHapTrioAndExAC = sum(diffHapTrioAndExAC),
                     #                    percCorrectSegExAC = (nSameHapTrioAndExAC+nDiffHapTrioAndExAC) / (NcoInExAC+NnotCoInExAC),
                     #                     percMissSeg = 1-((NcoInExAC+NnotCoInExAC)/NtrioHap),
                     nSameHapEM = sum(same_em_haplotype, na.rm = T),
                     nDiffHapEM = sum(!same_em_haplotype, na.rm = T),
                     nCorrectSameHap = sum(same_em_haplotype & same_trio_haplotype, na.rm = T),
                     nCorrectDiffHap = sum(!same_em_haplotype & !same_trio_haplotype, na.rm = T),
                     percCorrectEM = (nCorrectSameHap+nCorrectDiffHap)/(nSameHapEM+nDiffHapEM),
                     percCorrectSameHap = nCorrectSameHap / nSameHapEM,
                     percCorrectDiffHap = nCorrectDiffHap / nDiffHapEM,
                     percMissEM = 1-((nSameHapEM+nDiffHapEM)/NtrioHap),
                     nSameHapR = sum(single_het_ratio < 0.75, na.rm=T),
                     nDiffHapR = sum(single_het_ratio > 0.75, na.rm=T),
                     nCorrectSameHapR = sum(single_het_ratio < 0.75 & same_trio_haplotype, na.rm = T),
                     nCorrectDiffHapR = sum(single_het_ratio > 0.75 & !same_trio_haplotype, na.rm = T),
                     percCorrectR = (nCorrectSameHapR+nCorrectDiffHapR)/(nSameHapR+nDiffHapR),
                     percCorrectSameHapR = nCorrectSameHapR / nSameHapR,
                     percCorrectDiffHapR = nCorrectDiffHapR / nDiffHapR
                     
           ))
}

getSummary3 <-function(x, compute_EM_stats=T, compute_ratio_stats=T, compute_overlap_stats=F, compute_like_stats=F){
  
  res = x %>%
    summarise(nVariantPairs=length(unique(paste(chrom1,pos1,ref1,alt1,chrom2,pos2,ref2,alt2))),
              nHapTrio=n(),
              nSameHapTrio= sum(same_trio_haplotype),
              nDiffHapTrio = nHapTrio - nSameHapTrio,
              percSameTrioHap = nSameHapTrio / nHapTrio
    )
  
  if(compute_EM_stats){
    res %<>%
      inner_join(
        x %>%
          summarise(
            nHapTrio=n(),
            nSameHapTrio= sum(same_trio_haplotype),
            nDiffHapTrio = nHapTrio - nSameHapTrio,
            nCorrectSameHapEM = sum(same_em_haplotype & same_trio_haplotype, na.rm = T),
            nCorrectDiffHapEM = sum(!same_em_haplotype & !same_trio_haplotype, na.rm = T),
            nWrongSameHapEM = sum(same_em_haplotype & !same_trio_haplotype, na.rm = T),
            nWrongDiffHapEM = sum(!same_em_haplotype & same_trio_haplotype, na.rm = T),
            nMissingSameHapEM = sum(is.na(same_em_haplotype) & same_trio_haplotype, na.rm=T),
            nMissingDiffHapEM = sum(is.na(same_em_haplotype) & !same_trio_haplotype, na.rm=T),
            precisionSameHapEM = nCorrectSameHapEM / (nCorrectSameHapEM + nWrongSameHapEM),
            recallSameHapEM = nCorrectSameHapEM / nSameHapTrio,
            fSameHapEM = 2*(precisionSameHapEM * recallSameHapEM) / (precisionSameHapEM + recallSameHapEM),
            precisionDiffHapEM = nCorrectDiffHapEM / (nCorrectDiffHapEM + nWrongDiffHapEM),
            recallDiffHapEM = nCorrectDiffHapEM / nDiffHapTrio,
            fDiffHapEM = 2*(precisionDiffHapEM * recallDiffHapEM) / (precisionDiffHapEM + recallDiffHapEM),
            precisionEM = (nCorrectSameHapEM + nCorrectDiffHapEM) / (nCorrectSameHapEM + nCorrectDiffHapEM + nWrongSameHapEM + nWrongDiffHapEM),
            recallEM = (nCorrectSameHapEM + nCorrectDiffHapEM) / nHapTrio,
            fEM = 2*(precisionEM * recallEM) / (precisionEM + recallEM)
          )
      )
  }
  
  if(compute_ratio_stats){
    res %<>%
      inner_join(
        x %>%
          summarise(
            nHapTrio=n(),
            nSameHapTrio= sum(same_trio_haplotype),
            nDiffHapTrio = nHapTrio - nSameHapTrio,
            nCorrectSameHapR = sum(same_ratio_haplotype & same_trio_haplotype, na.rm = T),
            nCorrectDiffHapR = sum(!same_ratio_haplotype & !same_trio_haplotype, na.rm = T),
            nWrongSameHapR = sum(same_ratio_haplotype & !same_trio_haplotype, na.rm = T),
            nWrongDiffHapR = sum(!same_ratio_haplotype & same_trio_haplotype, na.rm = T),
            nMissingSameHapR = sum(is.na(same_ratio_haplotype) & same_trio_haplotype, na.rm=T),
            nMissingDiffHapR = sum(is.na(same_ratio_haplotype) & !same_trio_haplotype, na.rm=T),
            precisionSameHapR = nCorrectSameHapR / (nCorrectSameHapR + nWrongSameHapR),
            recallSameHapR = nCorrectSameHapR / nSameHapTrio,
            fSameHapR = 2*(precisionSameHapR * recallSameHapR) / (precisionSameHapR + recallSameHapR),
            precisionDiffHapR = nCorrectDiffHapR / (nCorrectDiffHapR + nWrongDiffHapR),
            recallDiffHapR = nCorrectDiffHapR / nDiffHapTrio,
            fDiffHapR = 2*(precisionDiffHapR * recallDiffHapR) / (precisionDiffHapR + recallDiffHapR),
            precisionR = (nCorrectSameHapR + nCorrectDiffHapR) / (nCorrectSameHapR + nCorrectDiffHapR + nWrongSameHapR + nWrongDiffHapR),
            recallR = (nCorrectSameHapR + nCorrectDiffHapR) / nHapTrio,
            fR = 2*(precisionR * recallR) / (precisionR + recallR)
          )
      )
  }
  
  if(compute_overlap_stats){
    res %<>%
      inner_join(
        x %>%
          summarise(
            nHapTrio=n(),
            nSameHapTrio= sum(same_trio_haplotype),
            nDiffHapTrio = nHapTrio - nSameHapTrio,
            nCorrectSameHapBoth = sum(same_em_haplotype & same_ratio_haplotype & same_trio_haplotype, na.rm = T),
            nCorrectDiffHapBoth = sum(!same_em_haplotype & !same_ratio_haplotype & !same_trio_haplotype, na.rm = T),
            nWrongSameHapBoth = sum((same_em_haplotype & same_ratio_haplotype) & !same_trio_haplotype, na.rm = T),
            nWrongDiffHapBoth = sum((!same_em_haplotype & !same_ratio_haplotype) & same_trio_haplotype, na.rm = T),
            nMissingSameHapBoth = sum(( (same_em_haplotype != same_ratio_haplotype) | is.na(same_em_haplotype) | is.na(same_ratio_haplotype)) & same_trio_haplotype, na.rm=T),
            nMissingDiffHapBoth = sum(( (same_em_haplotype != same_ratio_haplotype) | is.na(same_em_haplotype) | is.na(same_ratio_haplotype)) & !same_trio_haplotype, na.rm=T)
          )
      )
  }
  
  if(compute_like_stats){
    res %<>%
      inner_join(
        x %>%
          summarise(
            nHapTrio=n(),
            nSameHapTrio= sum(same_trio_haplotype),
            nDiffHapTrio = nHapTrio - nSameHapTrio,
            nCorrectSameHapLike = sum(same_like_haplotype & same_trio_haplotype, na.rm = T),
            nCorrectDiffHapLike = sum(!same_like_haplotype & !same_trio_haplotype, na.rm = T),
            nWrongSameHapLike = sum(same_like_haplotype & !same_trio_haplotype, na.rm = T),
            nWrongDiffHapLike = sum(!same_like_haplotype & same_trio_haplotype, na.rm = T),
            nMissingSameHapLike = sum(is.na(same_like_haplotype) & same_trio_haplotype, na.rm=T),
            nMissingDiffHapLike = sum(is.na(same_like_haplotype) & !same_trio_haplotype, na.rm=T),
            precisionSameHapLike = nCorrectSameHapLike / (nCorrectSameHapLike + nWrongSameHapLike),
            recallSameHapLike = nCorrectSameHapLike / nSameHapTrio,
            fSameHapLike = 2*(precisionSameHapLike * recallSameHapLike) / (precisionSameHapLike + recallSameHapLike),
            precisionDiffHapLike = nCorrectDiffHapLike / (nCorrectDiffHapLike + nWrongDiffHapLike),
            recallDiffHapLike = nCorrectDiffHapLike / nDiffHapTrio,
            fDiffHapLike = 2*(precisionDiffHapLike * recallDiffHapLike) / (precisionDiffHapLike + recallDiffHapLike),
            precisionLike = (nCorrectSameHapLike + nCorrectDiffHapLike) / (nCorrectSameHapLike + nCorrectDiffHapLike + nWrongSameHapLike + nWrongDiffHapLike),
            recallLike = (nCorrectSameHapLike + nCorrectDiffHapLike) / nHapTrio,
            fLike = 2*(precisionLike * recallLike) / (precisionLike + recallLike)
          )
      )
  }
  
  return(res)
}

get_mutation_string = function(ref, alt, cpg){
  
  if(nchar(ref) > 1 | nchar(alt) > 1){
    return(paste(ref,alt,sep=">"))
  }
  
  if(ref == "G" & alt =="A"){ 
    if(cpg){return("C>T@CpG") }
    else{return("C>T") }
  }
  if(ref == "G" & alt =="T"){ return("C>A") }
  if(ref == "T" & alt =="C"){ return("A>G") }
  if(ref == "T" & alt =="G"){ return("A>C") }
  if(ref == "T" & alt =="A"){ return("A>T") }
  if(ref == "G" & alt =="C"){ return("C>G") }
  if(ref =="C" & alt == "T" & cpg){return("C>T@CpG")}
  
  return(paste(ref,alt,sep=">"))
  
}


#Compute compound het likelihood based on the following 3 haplotypes model:
#Haplotype aB freq: p
#Haplotype Ab freq: q
#Haplotype AB freq: 1-p-q-e
#Haplotype ab freq: e
#
#       BB                    Bb                 bb
#AA (1-p-q-e)^2       2*(1-p-q-e)*q          q^2
#Aa 2*(1-p-q-e)*p   2*(p*q + (1-p-q-e)*e)    2*q*e
#aa    p^2                  2*p*e                e^2
compute_chet_likelihood = function(AABB, AaBB, aaBB, AABb, AaBb, aaBb, AAbb, Aabb, aabb, e = 1e-6, debug=F, contrast=F){
  n = 2* (AABB + AaBB + aaBB + AABb + AaBb + aaBb + AAbb + Aabb + aabb)
  p = (AaBB + AaBb + Aabb+  2*(aaBB + aaBb + aabb)) / n
  q = (AABb + AaBb + aaBb + 2*(AAbb + Aabb + aabb)) / n
  
  #Is this right, or should it be this?
  #p = (AaBB + AaBb + aaBb + 2*aaBB) / n
  #q = (AABb + AaBb + Aabb + 2*AAbb) / n
  
  if(debug){
    print("---Chet hap likelihoods---")
    print(paste("p:", p))
    print(paste("q:", q))
    print(paste("AABB",log10(1-p-q-e)*2*AABB))
    print(paste("AABb",log10(2*(1-p-q-e)*q)*AABb))
    print(paste("AAbb",log10(q)*2*AAbb))
    print(paste("AaBB",log10(2*(1-p-q-e)*p)*AaBB))
    print(paste("AaBb",log10(2*(p*q + (1-p-q-e)*e))*AaBb))
    print(paste("Aabb",log10(2*q*e)*Aabb))
    print(paste("aaBB",log10(p)*2*aaBB))
    print(paste("aaBb",log10(2*p*e)*aaBb))
    print(paste("aabb",log10(e)*2*aabb))
  }
  
  if(contrast){
    return(
      log10(2*(p*q + (1-p-q-e)*e))*AaBb +
      ifelse(p>q,
             log10(2*(1-p-q-e)*q)*AABb,
             log10(2*(1-p-q-e)*p)*AaBB)
    )
  }
  
  return(
    log10(1-p-q-e)*2*AABB +
    log10(2*(1-p-q-e)*q)*AABb +
    log10(q)*2*AAbb +
    log10(2*(1-p-q-e)*p)*AaBB +
    log10(2*(p*q + (1-p-q-e)*e))*AaBb +
    log10(2*q*e)*Aabb +
    log10(p)*2*aaBB +
    log10(2*p*e)*aaBb +
    log10(e)*2*aabb
  )
}

#Compute same haplotype likelihood based on the following 3 haplotypes model:
#where p >= q
#Haplotype aB freq: p-q ## 0 and Interchangeable with Ab when perfectly in LD
#Haplotype Ab freq: e
#Haplotype AB freq: 1-p-e
#Haplotype ab freq: q
#
#       BB             Bb                   bb
#AA   (1-p-e)^2        2*(1-p-e)*e             e^2
#Aa  2*(1-p-e)*(p-q) 2*((1-p-e)*q + e*(p-q))  2*q*e
#aa   (p-q)^2        2*(p-q)*q             q^2
compute_same_hap_likelihood = function(AABB, AaBB, aaBB, AABb, AaBb, aaBb, AAbb, Aabb, aabb, e = 1e-6, debug=F, contrast=F){
  
  compute_ordered_same_hap_likelihood = function(AABB, AaBB, aaBB, AABb, AaBb, aaBb, AAbb, Aabb, aabb,p,q,e,debug,contrast){
    
    if(debug){
      print(paste("p:",p))
      print(paste("q:",q))
      print(paste("e:",e))
      print(paste("AABB",log10(1-p-e)*2*AABB))
      print(paste("AABb",log10(2*(1-p-e)*e)*AABb))
      print(paste("AAbb",log10(e)*2*AAbb))
      print(paste("AaBB",log10(2*(1-p-e)*(p-q))*AaBB))
      print(paste("AaBb",log10(2*((1-p-e)*q + e*(p-q)))*AaBb))
      print(paste("Aabb",log10(2*q*e)*Aabb))
      print(paste("aaBB",log10(p-q)*2*aaBB))
      print(paste("aaBb",log10(2*(p-q)*q)*aaBb))
      print(paste("aabb",log10(q)*2*aabb))
    }
    
    if(contrast){
      return(log10(2*((1-p-e)*q + e*(p-q)))*AaBb + 
               log10(2*(1-p-e)*e)*AABb)
    }
    
    result = log10(1-p-e)*2*AABB +
      log10(2*(1-p-e)*e)*AABb +
      log10(e)*2*AAbb +
      log10(2*((1-p-e)*q + e*(p-q)))*AaBb
    
    if(q>0){
      result = result + 
        log10(2*q*e)*Aabb +
        log10(q)*2*aabb
      if(p!=q ){
        result = result + log10(2*(p-q)*q)*aaBb
      }
    }
    
    if(p!=q ){
      result = result +
          log10(2*(1-p-e)*(p-q))*AaBB +
          log10(p-q)*2*aaBB 
    }

    return(result)
   
  }
  
  n = 2* (AABB + AaBB + aaBB + AABb + AaBb + aaBb + AAbb + Aabb + aabb)
  f1 = (AaBB + AaBb + Aabb + 2*(aabb + aaBB + aaBb)) / n
  f2 = (AABb + AaBb + aaBb + 2*(AAbb + Aabb + aabb)) / n
  
  if(debug){
    print("---Same hap likelihoods---")
    print(paste("counts: ", AABB, AaBB, aaBB, AABb, AaBb, aaBb, AAbb, Aabb, aabb))
    print(paste("f1: ",f1))
    print(paste("f2: ",f2))
  }
  
  q = (AaBb + Aabb + aaBb + 2*aabb) / n
  if(is.na(f1) | is.na(f2)){
    return(NA)
  }
  
  if(f1>f2){
    return(compute_ordered_same_hap_likelihood(AABB, AaBB, aaBB, AABb, AaBb, aaBb, AAbb, Aabb, aabb,f1,q,e,debug,contrast))
  }
  else{
    return(compute_ordered_same_hap_likelihood(AABB, AABb, AAbb, AaBB, AaBb, Aabb, aaBB, aaBb,aabb,f2,q,e,debug,contrast))
  }
}


preprocess_gnomad = function(data, trios_meta, 
                             same_hap_em_threshold = 0.25, 
                             diff_hap_em_threshold = 0.25, 
                             same_hap_ratio_threshold = 0.5, 
                             diff_hap_ratio_threshold = 0.9, 
                             same_hap_like_threshold = 0.95, 
                             diff_hap_like_threshold = 1.05,
                             dprime_threshold = 0.25){
  
  bad_fams = unique(trios_meta %>% filter(!drop_reasons %in% c('kept','related','child','NR',"NO_PERM")) %>% .$famID)
  
  return(
    data %>%
      filter(!fam %in% bad_fams) %>%
      mutate(
        ac_notrios1 = ifelse(is.na(ac_notrios1),0,ac_notrios1),
        ac_notrios2 = ifelse(is.na(ac_notrios2),0,ac_notrios2),
        #ac1 = AaBB + AaBb + Aabb + 2*(aaBB + aaBb + aabb),
        #ac2 = AABb + AaBb + aaBb + 2*(AAbb + Aabb + aabb),
        #an = 2*(AABB + AaBB + aaBB + AABb + AaBb + aaBb + AAbb + Aabb + aabb),
        #af1 = ifelse(an > 0, ac1 / an, 0),
        #af2 = ifelse(an > 0, ac2 / an, 0), 
        n = AB + Ab + aB + ab,
        pAB = AB/n,
        pAb = Ab/n,
        paB = aB/n,
        pab = ab/n,
        pa = pab+paB,
        pb = pab+pAb,
        d = pab - (pa * pb),
        kid_adj = kid_v1_gq >= 20 & kid_v2_gq >= 20 & kid_v1_dp >= 10 & kid_v2_dp >= 10 & (kid_v1_gt != 1 | kid_v1_ad1/kid_v1_dp > 0.2)  & (kid_v2_gt != 1 | kid_v2_ad1/kid_v2_dp > 0.2),
        mom_adj = mom_v1_gq >= 20 & mom_v2_gq >= 20 & mom_v1_dp >= 10 & mom_v2_dp >= 10 & (mom_v1_gt != 1 | mom_v1_ad1/mom_v1_dp > 0.2)  & (mom_v2_gt != 1 | mom_v2_ad1/mom_v2_dp > 0.2),
        dad_adj = dad_v1_gq >= 20 & dad_v2_gq >= 20 & dad_v1_dp >= 10 & dad_v2_dp >= 10 & (dad_v1_gt != 1 | dad_v1_ad1/dad_v1_dp > 0.2)  & (dad_v2_gt != 1 | dad_v2_ad1/dad_v2_dp > 0.2),
        pass1 = pass1 == "true",
        pass2 = pass2 == "true",
        cpg1 = ifelse(!is.na(cpg1) & cpg1 == "true", T, F),
        cpg2 = ifelse(!is.na(cpg2) & cpg2 == "true", T, F),
        wasSplit1 = wasSplit1 == "true",
        wasSplit2 = wasSplit2 == "true",
        same_trio_haplotype = same_trio_haplotype == "true",
        same_em_haplotype = case_when(
          prob_same_haplotype > same_hap_em_threshold ~ T, 
          prob_same_haplotype < diff_hap_em_threshold ~ F, 
          T ~ NA),
        single_het_ratio = ifelse(AaBB > AABb,
                                  (AABb + AAbb) / (AABb + AAbb + AaBb + aaBb + Aabb + aabb),
                                  (AaBB + aaBB) / (AaBB + aaBB + AaBb + aaBb + Aabb + aabb))
        ) %>%
      rowwise() %>% 
      mutate(
        #single_het_ratio = min(AaBB,AABb) / (min(AaBB,AABb) + AaBb),
             dprime = d / ifelse( 
               d<0,
               max(-pa*pb, -(1-pa)*(1-pb)),
               min(pa*(1-pb), (1-pa)*pb)
             ),
             chet_likelihood = compute_chet_likelihood(AABB,AaBB,aaBB,AABb,AaBb,aaBb,AAbb,Aabb,aabb),
             same_hap_likelihood = compute_same_hap_likelihood(AABB,AaBB,aaBB,AABb,AaBb,aaBb,AAbb,Aabb,aabb),
             like_ratio = chet_likelihood / same_hap_likelihood,
             same_like_haplotype = case_when(
               like_ratio > same_hap_like_threshold ~ T,
               like_ratio < diff_hap_like_threshold ~ F,
               T ~ NA
             )
             ) %>%
      ungroup() %>%
      mutate(
        same_ratio_haplotype = case_when(single_het_ratio > diff_hap_ratio_threshold ~ F,
                                         single_het_ratio < same_hap_ratio_threshold ~ T,
                                         T ~ NA)
      )
    )
}

get_em_table_for_ppt = function(data){
  data %>% transmute(
    pop = pop,
    same_same = nCorrectSameHapEM,
    same_same_perc = nCorrectSameHapEM / nSameHapTrio,
    same_diff = nDiffHap - nCorrectDiffHapEM,
    same_diff_perc = same_diff / nSameHapTrioEM,
    same_missing = nSameHapTrio - same_diff - same_same,
    same_missing_perc = same_missing / nSameHapTrio,
    diff_same = nSameHapEM - nCorrectSameHapEM,
    diff_same_perc = diff_same / nDiffHapTrio,
    diff_diff = nCorrectDiffHap,
    diff_diff_perc = nCorrectDiffHap / nDiffHapTrio,
    diff_missing = nDiffHapTrio - diff_diff - diff_same,
    diff_missing_perc = diff_missing / nDiffHapTrio
    
    
  )
  
}


get_ratio_table_for_ppt = function(data){
  data %>% transmute(
    pop = pop,
    same_same = nCorrectSameHap,
    same_same_perc = nCorrectSameHap / nSameHapTrio,
    same_diff = nDiffHapEM - nCorrectDiffHap,
    same_diff_perc = same_diff / nSameHapTrio,
    same_missing = nSameHapTrio - same_diff - same_same,
    same_missing_perc = same_missing / nSameHapTrio,
    diff_same = nSameHapEM - nCorrectSameHap,
    diff_same_perc = diff_same / nDiffHapTrio,
    diff_diff = nCorrectDiffHap,
    diff_diff_perc = nCorrectDiffHap / nDiffHapTrio,
    diff_missing = nDiffHapTrio - diff_diff - diff_same,
    diff_missing_perc = diff_missing / nDiffHapTrio
    
    
  )
  
}



