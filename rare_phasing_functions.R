
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
                     NsameTrioHap= sum(sameTrioHap),
                     NdiffTrioHap = sum(diffTrioHap),
                     percSameTrioHap = sum(sameTrioHap) / NtrioHap,
                     #                    NcoInExAC= sum(coInExAC),
                     #                     NnotCoInExAC = sum(notCoInExAC),
                     #                     NsameHapTrioAndExAC = sum(sameHapTrioAndExAC),
                     #                    NdiffHapTrioAndExAC = sum(diffHapTrioAndExAC),
                     #                    percCorrectSegExAC = (NsameHapTrioAndExAC+NdiffHapTrioAndExAC) / (NcoInExAC+NnotCoInExAC),
                     #                     percMissSeg = 1-((NcoInExAC+NnotCoInExAC)/NtrioHap),
                     NsameHapExACEM = sum(sameHapExACEM),
                     NdiffHapExACEM = sum(diffHapExACEM),
                     NsameHapTrioAndExACEM = sum(sameHapTrioAndExACEM),
                     NdiffHapTrioAndExACEM = sum(diffHapTrioAndExACEM),
                     percCorrectEM = (NsameHapTrioAndExACEM+NdiffHapTrioAndExACEM)/(NsameHapExACEM+NdiffHapExACEM),
                     percCorrectSameHapEM = NsameHapTrioAndExACEM / NsameHapExACEM,
                     percCorrectDiffHapEM = NdiffHapTrioAndExACEM / NdiffHapExACEM,
                     percMissEM = 1-((NsameHapExACEM+NdiffHapExACEM)/NtrioHap)
                     
           ))
}

getSummary2 <-function(x){
  return(x %>% 
           #mutate(single_het_ratio = (AaBB + AABb) / (AaBB + AABb + AaBb)) %>%
           summarise(NtrioHap=n(),
                     NsameTrioHap= sum(same_trio_haplotype),
                     NdiffTrioHap = NtrioHap - NsameTrioHap,
                     percSameTrioHap = NsameTrioHap / NtrioHap,
                     #                    NcoInExAC= sum(coInExAC),
                     #                     NnotCoInExAC = sum(notCoInExAC),
                     #                     NsameHapTrioAndExAC = sum(sameHapTrioAndExAC),
                     #                    NdiffHapTrioAndExAC = sum(diffHapTrioAndExAC),
                     #                    percCorrectSegExAC = (NsameHapTrioAndExAC+NdiffHapTrioAndExAC) / (NcoInExAC+NnotCoInExAC),
                     #                     percMissSeg = 1-((NcoInExAC+NnotCoInExAC)/NtrioHap),
                     NsameHapEM = sum(same_em_haplotype, na.rm = T),
                     NdiffHapEM = sum(!same_em_haplotype, na.rm = T),
                     NcorrectSameHap = sum(same_em_haplotype & same_trio_haplotype, na.rm = T),
                     NcorrectDiffHap = sum(!same_em_haplotype & !same_trio_haplotype, na.rm = T),
                     percCorrectEM = (NcorrectSameHap+NcorrectDiffHap)/(NsameHapEM+NdiffHapEM),
                     percCorrectSameHap = NcorrectSameHap / NsameHapEM,
                     percCorrectDiffHap = NcorrectDiffHap / NdiffHapEM,
                     percMissEM = 1-((NsameHapEM+NdiffHapEM)/NtrioHap),
                     NsameHapR = sum(single_het_ratio < 0.75, na.rm=T),
                     NdiffHapR = sum(single_het_ratio > 0.75, na.rm=T),
                     NcorrectSameHapR = sum(single_het_ratio < 0.75 & same_trio_haplotype, na.rm = T),
                     NcorrectDiffHapR = sum(single_het_ratio > 0.75 & !same_trio_haplotype, na.rm = T),
                     percCorrectR = (NcorrectSameHapR+NcorrectDiffHapR)/(NsameHapR+NdiffHapR),
                     percCorrectSameHapR = NcorrectSameHapR / NsameHapR,
                     percCorrectDiffHapR = NcorrectDiffHapR / NdiffHapR
                     
           ))
}

getSummary3 <-function(x){
  return(x %>% 
           #mutate(single_het_ratio = (AaBB + AABb) / (AaBB + AABb + AaBb)) %>%
           summarise(NtrioHap=n(),
                     NsameTrioHap= sum(same_trio_haplotype),
                     NdiffTrioHap = NtrioHap - NsameTrioHap,
                     percSameTrioHap = NsameTrioHap / NtrioHap,
                     #                    NcoInExAC= sum(coInExAC),
                     #                     NnotCoInExAC = sum(notCoInExAC),
                     #                     NsameHapTrioAndExAC = sum(sameHapTrioAndExAC),
                     #                    NdiffHapTrioAndExAC = sum(diffHapTrioAndExAC),
                     #                    percCorrectSegExAC = (NsameHapTrioAndExAC+NdiffHapTrioAndExAC) / (NcoInExAC+NnotCoInExAC),
                     #                     percMissSeg = 1-((NcoInExAC+NnotCoInExAC)/NtrioHap),
                     NcorrectSameHapEM = sum(same_em_haplotype & same_trio_haplotype, na.rm = T),
                     NcorrectDiffHapEM = sum(!same_em_haplotype & !same_trio_haplotype, na.rm = T),
                     NwrongSameHapEM = sum(same_em_haplotype & !same_trio_haplotype, na.rm = T),
                     NwrongDiffHapEM = sum(!same_em_haplotype & same_trio_haplotype, na.rm = T),
                     NmissingSameHapEM = sum(is.na(same_em_haplotype) & same_trio_haplotype, na.rm=T),
                     NmissingDiffHapEM = sum(is.na(same_em_haplotype) & !same_trio_haplotype, na.rm=T),
                     precisionSameHapEM = NcorrectSameHapEM / (NcorrectSameHapEM + NwrongSameHapEM),
                     recallSameHapEM = NcorrectSameHapEM / NsameTrioHap,
                     precisionDiffHapEM = NcorrectDiffHapEM / (NcorrectDiffHapEM + NwrongDiffHapEM),
                     recallDiffHapEM = NcorrectDiffHapEM / NdiffTrioHap,
                     precisionEM = (NcorrectSameHapEM + NcorrectDiffHapEM) / (NcorrectSameHapEM + NcorrectDiffHapEM + NwrongSameHapEM + NwrongDiffHapEM),
                     recallEM = (NcorrectSameHapEM + NcorrectDiffHapEM) / NtrioHap,
                     NcorrectSameHapR = sum(same_ratio_haplotype & same_trio_haplotype, na.rm = T),
                     NcorrectDiffHapR = sum(!same_ratio_haplotype & !same_trio_haplotype, na.rm = T),
                     NwrongSameHapR = sum(same_ratio_haplotype & !same_trio_haplotype, na.rm = T),
                     NwrongDiffHapR = sum(!same_ratio_haplotype & same_trio_haplotype, na.rm = T),
                     NmissingSameHapR = sum(is.na(same_ratio_haplotype) & same_trio_haplotype, na.rm=T),
                     NmissingDiffHapR = sum(is.na(same_ratio_haplotype) & !same_trio_haplotype, na.rm=T),
                     precisionSameHapR = NcorrectSameHapR / (NcorrectSameHapR + NwrongSameHapR),
                     recallSameHapR = NcorrectSameHapR / NsameTrioHap,
                     precisionDiffHapR = NcorrectDiffHapR / (NcorrectDiffHapR + NwrongDiffHapR),
                     recallDiffHapR = NcorrectDiffHapR / NdiffTrioHap,
                     precisionR = (NcorrectSameHapR + NcorrectDiffHapR) / (NcorrectSameHapR + NcorrectDiffHapR + NwrongSameHapR + NwrongDiffHapR),
                     recallR = (NcorrectSameHapR + NcorrectDiffHapR) / NtrioHap,
                     NcorrectSameHapBoth = sum(same_em_haplotype & same_ratio_haplotype & same_trio_haplotype, na.rm = T),
                     NcorrectDiffHapBoth = sum(!same_em_haplotype & !same_ratio_haplotype & !same_trio_haplotype, na.rm = T),
                     NwrongSameHapBoth = sum((same_em_haplotype & same_ratio_haplotype) & !same_trio_haplotype, na.rm = T),
                     NwrongDiffHapBoth = sum((!same_em_haplotype & !same_ratio_haplotype) & same_trio_haplotype, na.rm = T),
                     NmissingSameHapBoth = sum(( (same_em_haplotype != same_ratio_haplotype) | is.na(same_em_haplotype) | is.na(same_ratio_haplotype)) & same_trio_haplotype, na.rm=T),
                     NmissingDiffHapBoth = sum(( (same_em_haplotype != same_ratio_haplotype) | is.na(same_em_haplotype) | is.na(same_ratio_haplotype)) & !same_trio_haplotype, na.rm=T)
           ))
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

preprocess_gnomad = function(data, trios_meta, prob_threshold = 0.25, ratio_threshold = 0.25){
  
  bad_fams = unique(trios_meta %>% filter(!drop_reasons %in% c('kept','related','child','NR',"NO_PERM")) %>% .$famID)
  
  return(
    data %>%
      filter(!fam %in% bad_fams) %>%
      filter(!fam %in% bad_fams) %>%
      mutate(
        #ac1 = AaBB + AaBb + Aabb + 2*(aaBB + aaBb + aabb),
        #ac2 = AABb + AaBb + aaBb + 2*(AAbb + Aabb + aabb),
        #an = 2*(AABB + AaBB + aaBB + AABb + AaBb + aaBb + AAbb + Aabb + aabb),
        #af1 = ifelse(an > 0, ac1 / an, 0),
        #af2 = ifelse(an > 0, ac2 / an, 0),
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
          prob_same_haplotype > 1-prob_threshold ~ T, 
          prob_same_haplotype < prob_threshold ~ F, 
          T ~ NA)) %>%
      rowwise() %>% 
      mutate(single_het_ratio = min(AaBB,AABb) / (min(AaBB,AABb) + AaBb)) %>%
      ungroup() %>%
      mutate(
        same_ratio_haplotype = case_when(single_het_ratio > 1-ratio_threshold ~ F,
                                         single_het_ratio < ratio_threshold ~ T,
                                         T ~ NA)
      )
    )
}

get_em_table_for_ppt = function(data){
  data %>% transmute(
    pop = pop,
    same_same = NcorrectSameHapEM,
    same_same_perc = NcorrectSameHapEM / NsameTrioHap,
    same_diff = NdiffHap - NcorrectDiffHapEM,
    same_diff_perc = same_diff / NsameTrioHapEM,
    same_missing = NsameTrioHap - same_diff - same_same,
    same_missing_perc = same_missing / NsameTrioHap,
    diff_same = NsameHapEM - NcorrectSameHapEM,
    diff_same_perc = diff_same / NdiffTrioHap,
    diff_diff = NcorrectDiffHap,
    diff_diff_perc = NcorrectDiffHap / NdiffTrioHap,
    diff_missing = NdiffTrioHap - diff_diff - diff_same,
    diff_missing_perc = diff_missing / NdiffTrioHap
    
    
  )
  
}


get_ratio_table_for_ppt = function(data){
  data %>% transmute(
    pop = pop,
    same_same = NcorrectSameHap,
    same_same_perc = NcorrectSameHap / NsameTrioHap,
    same_diff = NdiffHapEM - NcorrectDiffHap,
    same_diff_perc = same_diff / NsameTrioHap,
    same_missing = NsameTrioHap - same_diff - same_same,
    same_missing_perc = same_missing / NsameTrioHap,
    diff_same = NsameHapEM - NcorrectSameHap,
    diff_same_perc = diff_same / NdiffTrioHap,
    diff_diff = NcorrectDiffHap,
    diff_diff_perc = NcorrectDiffHap / NdiffTrioHap,
    diff_missing = NdiffTrioHap - diff_diff - diff_same,
    diff_missing_perc = diff_missing / NdiffTrioHap
    
    
  )
  
}



