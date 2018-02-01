library(cowplot)

### Data Loading ####
#Read fam file
exomes_fam = read.delim("gnomad_exomes.qctrios.fam",stringsAsFactors = F, header = F)
colnames(exomes_fam) = c("famID","sample","matID","patID","sex","pheno")
exomes_fam %<>% select(sample,famID,matID,patID)

#Read samples metadata
exomes_meta = read.delim("super_meta_april_01_2017.txt.bgz", stringsAsFactors = F)
exomes_meta_pca = read.delim("all_data.tsv.gz",stringsAsFactors = F)
exomes_meta %<>% join(exomes_meta_pca %>% select(sample, pid, pc1, pc2, pc3,pc4,pc5,pc6,pc7,pc8,pc9,pc10))
exomes_meta$trio = exomes_meta$sample %in% exomes_fam$sample
rm(exomes_meta_pca)

#Compile trio metadata
trios_meta = exomes_fam %>% inner_join(exomes_meta)
#trios_meta = exomes_fam %>% inner_join(gnomad_pca %>% mutate(sample = gsub('exome_','',sample)))

#PC plot
ggplot() + geom_point(data = exomes_meta %>% mutate(population = as.factor(population)) %>% filter(!trio), 
                      mapping = aes(pc1,pc2,col=population), alpha=0.05 ,shape = 4) + 
  geom_point(data = exomes_meta %>% mutate(population = as.factor(population)) %>% filter(trio & population == "oth"), 
             mapping = aes(pc1,pc2,fill=population), shape = 21)


#Read trio data
#Preprocess data:
#Convert true / false to T / F
#Remove bad samples / families
bad_fams = unique(trios_meta %>% filter(!drop_reasons %in% c('kept','related','child','NR',"NO_PERM")) %>% .$famID)
gnomad = read.delim("exome_trios_full_byPop.txt.bgz", stringsAsFactors = F) %>%
  mutate(prob_diff_haplotype = 1 - prob_same_haplotype) %>%
  preprocess_gnomad(trios_meta,
                    same_hap_ratio_threshold = 0.5, 
                    diff_hap_ratio_threshold = 0.75,
                    same_hap_em_threshold = 0.95, 
                    diff_hap_em_threshold = 0.5)
gnomad %<>% filter((is.na(an_release1) | an_release1 > 221670) & 
                            (is.na(an_release2) | an_release2 > 221670) & 
                            kid_adj & mom_adj & dad_adj)

gnomad %<>% mutate(
  indel1 = nchar(ref1) != nchar(alt1),
  indel2 = nchar(ref2) != nchar(alt2)
)


gnomad_no_pop = read.delim("exome_trios_full.txt.bgz", stringsAsFactors = F) %>%
  mutate(prob_diff_haplotype = 1 - prob_same_haplotype) %>%
  preprocess_gnomad(trios_meta,
                    same_hap_ratio_threshold = 0.5, 
                    diff_hap_ratio_threshold = 0.75,
                    same_hap_em_threshold = 0.95, 
                    diff_hap_em_threshold = 0.5)
gnomad_no_pop %<>% filter((is.na(an_release1) | an_release1 > 221670) & 
                            (is.na(an_release2) | an_release2 > 221670) & 
                            kid_adj & mom_adj & dad_adj)
gnomad_no_pop %<>% mutate(
  indel1 = nchar(ref1) != nchar(alt1),
  indel2 = nchar(ref2) != nchar(alt2)
)

## Old version
gnomad_long = gnomad %>%
  mutate(method = "byPop") %>%
  rbind(
    gnomad_no_pop %>%
      mutate(method = "noPop")
  )

gnomad_long2 = gnomad %>% 
  mutate(byPop = T) %>%
  rbind(
    gnomad_no_pop %>%
      mutate(byPop = F)
  ) %>%
  gather("method","score", "prob_diff_haplotype", "single_het_ratio") %>%
  mutate(
    method = ifelse(method == "prob_diff_haplotype", "EM", "ratio"),
    same_pred_haplotype = ifelse(method == "prob_diff_haplotype", same_em_haplotype, same_ratio_haplotype)
  ) %>%
  select(-c(same_em_haplotype, same_ratio_haplotype)) %>%
  rowwise() %>%
  mutate(
    min_ac_release = min(ac_release1, ac_release2)
  ) %>%
  ungroup()

#gnomad_adj = read.delim("exome_trios_adj.txt.bgz", stringsAsFactors = F) %>%
#  preprocess_gnomad(trios_meta, prob_threshold = 0.25)

#### Likelihood model ####
### TESTS 
gnomad %>%
  filter(chrom1 == "7" & pos1 == 105662766 & pos2 == 105665004) %>%
  rowwise() %>%
  mutate(
    chet_likelihood = compute_chet_likelihood(AABB,AaBB,aaBB,AABb,AaBb,aaBb,AAbb,Aabb,aabb,debug=T, contrast = T),
    same_hap_likelihood = compute_same_hap_likelihood(AABB,AaBB,aaBB,AABb,AaBb,aaBb,AAbb,Aabb,aabb,debug=T, contrast = T),
    same_like_haplotype = same_hap_likelihood > chet_likelihood
  ) %>% select(
    AABB,AaBB,aaBB,AABb,AaBb,aaBb,AAbb,Aabb,aabb,chet_likelihood,same_hap_likelihood,
    same_like_haplotype,same_em_haplotype,same_trio_haplotype
  )


gnomad %>% 
  rowwise() %>%
  mutate(
    chet_likelihood = compute_chet_likelihood(AABB,AaBB,aaBB,AABb,AaBb,aaBb,AAbb,Aabb,aabb,contrast=T),
    same_hap_likelihood = compute_same_hap_likelihood(AABB,AaBB,aaBB,AABb,AaBb,aaBb,AAbb,Aabb,aabb,contrast=T),
    same_like_haplotype = same_hap_likelihood > chet_likelihood,
    like_ratio =  same_hap_likelihood - chet_likelihood
  ) %>% select(
    #AABB,AaBB,aaBB,AABb,AaBb,aaBb,AAbb,Aabb,aabb,
    AaBB,AABb,AaBb,
    chet_likelihood,
    same_hap_likelihood,
    like_ratio,
    same_like_haplotype,
    same_em_haplotype,
    same_trio_haplotype
  ) %>%
  filter(same_like_haplotype != same_trio_haplotype) %>%
  head(100) %>%
  View()

#Plot number of true hets by likelihood ratio
gnomad %>%
  mutate(like_ratio = chet_likelihood - same_hap_likelihood,
         bin = ntile(like_ratio, 100)) %>%
  group_by(bin)  %>%
  summarise(
    prop_same_hap = sum(same_trio_haplotype) / n(),
    min_ratio = min(like_ratio),
    max_ratio = max(like_ratio)
  ) %>%
  ggplot(aes(min_ratio, prop_same_hap)) + geom_point()

gnomad %<>%
  mutate(
    same_like_haplotype = case_when(
      like_ratio >= 1 ~ T,
      like_ratio < 0.95 ~ F,
      T ~ NA
    )
  )
  
##Plot ROC curves
#Chets
p = gnomad %>%
  mutate(like_ratio = - chet_likelihood + same_hap_likelihood,
         bin = ntile(like_ratio, 100)) %>%
  arrange(like_ratio) %>%
  filter(!is.na(like_ratio)) %>%
  mutate(
    cumtp = cumsum(!same_trio_haplotype),
    cumfp = cumsum(same_trio_haplotype),
    p = sum(!same_trio_haplotype),
    n = sum(same_trio_haplotype),
    tpr = cumtp / p,
    fpr = cumfp / n
  ) %>%
  sample_n(5000, replace = F) %>%
  ggplot(aes(fpr,tpr,text=like_ratio)) + geom_point()

p = gnomad %>%
  mutate(like_ratio = chet_likelihood / same_hap_likelihood,
         bin = ntile(like_ratio, 100)) %>%
  arrange(-like_ratio) %>%
  filter(!is.na(like_ratio)) %>%
  mutate(
    cumtp = cumsum(same_trio_haplotype),
    cumfp = cumsum(!same_trio_haplotype),
    p = sum(same_trio_haplotype),
    n = sum(!same_trio_haplotype),
    tpr = cumtp / p,
    fpr = cumfp / n
  ) %>%
  sample_n(5000, replace = F) %>%
  ggplot(aes(fpr,tpr,text=like_ratio)) + geom_point()

  ggplotly(p)


#### Threshold selection ####


#Violin plot of scores by AC
gnomad_long2 %>%
  mutate(ac_bin = ntile(min_ac_release, 5)) %>%
  group_by(ac_bin) %>%
  mutate(min_ac_bin = min(min_ac_release),
         max_ac_bin = max(min_ac_release),
         label = sprintf("[%d,%d[",first(min_ac_bin), first(max_ac_bin))
         ) %>%
  filter(!byPop) %>%
  ggplot(aes(score, fill=method)) + geom_density(position="dodge")

  ggplot(aes(method, score)) + geom_violin() + facet_grid(label ~ .) + 
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) +
  xlab("Minimum allele count")
  
  filter(!byPop & min_ac_release > 1060) %>%
  ggplot(aes(score, fill=method)) + geom_histogram(position="dodge")

  
#Investigate weird bi-modality
pdf("likelihood_bimodality.pdf")
for(x in c("ac_notrios1","ac_notrios2","AaBB","AABb","AaBb")){
  p = gnomad %>%
  filter(same_trio_haplotype & like_ratio > 0.9 & ac_notrios1 < 100) %>%
    mutate(lr = like_ratio < 1.5) %>%
  ggplot(aes_string(x, fill="lr")) + geom_histogram(position="dodge") 
  print(p)
}
dev.off()
  
gnomad %>%
  #filter(same_trio_haplotype & like_ratio > 0.9 & ac_notrios1 < 100) %>%
  filter(same_trio_haplotype & abs(same_hap_likelihood - chet_likelihood) < 100) %>%
  sample_n(500) %>% 
  select(chrom1, pos1, pos2, AaBB, AABb, AaBb, same_hap_likelihood, chet_likelihood) %>%
  View()
  
  


ratio_bins = gnomad_long %>%
  group_by(method) %>%
  mutate(het_ratio_bin = ntile(single_het_ratio,500)) %>%
  group_by(method, het_ratio_bin) %>%
  summarise(
    n = n(),
    propSameHap = sum(same_trio_haplotype)/n,
    propDiffHap = sum(!same_trio_haplotype)/n,
    min_ratio = min(single_het_ratio),
    max_ratio = max(single_het_ratio)
  ) 

pdf("threshold_selection_ratio.pdf", width=6.5, height = 3)
ratio_bins %>%
  mutate(same_hap = case_when(
    max_ratio < 0.5 ~ "red",
    max_ratio > 0.85 ~ "blue"
  )) %>%
  filter(method == "byPop") %>%
  ggplot(aes(max_ratio , propDiffHap, col=same_hap)) + geom_point() +
  xlab("Proportion of individuals carrying only the rarest allele") + 
  ylab("Proportion of\ncompound heterozygotes") +
  geom_vline(xintercept = 0.5, linetype="dashed") + 
  geom_vline(xintercept = 0.85, linetype="dashed") +
  theme_bw(base_size = 14) +
  theme(legend.position="None") 
dev.off()

  em_bins = gnomad_long %>%
  group_by(method) %>%
  mutate(em_bin = ntile(prob_same_haplotype,100)) %>%
  group_by(method, em_bin) %>%
  summarise(
    n = n(),
    propSameHap = sum(same_trio_haplotype)/n,
    propDiffHap = sum(!same_trio_haplotype)/n,
    min_ratio = min(prob_same_haplotype),
    max_ratio = max(prob_same_haplotype)
  ) 

em_bins %>%
  ggplot(aes(em_bin, propDiffHap, col=method)) + geom_point()

#Merge numbers from haplotyping using global vs by population
global_pop_columns = c(4,19,34,37:49,98:100)
gnomad_all = gnomad %>% full_join(
  gnomad_no_pop %>% 
    rename_(.dots=setNames(names(.)[global_pop_columns], paste0("global_", names(.)[global_pop_columns])))
) 


#View summary
gnomad_long %>%
  filter(!cpg1 & !cpg2) %>%
  #filter(ac_notrios1 > 0 & ac_notrios2 > 0) %>%
  group_by(method,pop) %>%
  getSummary3() %>%
  View

gnomad_long %>%
  filter(ac_release1 > 1060 & ac_release2 > 1060) %>%
  #filter(ac_notrios1 > 0 & ac_notrios2 > 0) %>%
  group_by(method,pop) %>%
  getSummary3() %>%
  View

#Look at distribution of genotype qualities
gnomad_all %>%
  transmute(ratio_correct = same_ratio_haplotype == same_trio_haplotype,
            em_correct = same_em_haplotype == same_trio_haplotype,
            kid_v1_gq = kid_v1_gq,
            kid_v2_gq = kid_v2_gq,
            mom_v1_gq = mom_v1_gq,
            mom_v2_gq = mom_v2_gq,
            dad_v1_gq = dad_v1_gq,
            dad_v2_gq = dad_v2_gq
            ) %>%
  gather("sample","gq",3:8) %>% 2
ggplot(aes(gq, col =ratio_correct)) + geom_density()

#Look at callrate
gnomad_all %>%
  transmute(ratio_correct = same_ratio_haplotype == same_trio_haplotype,
            em_correct = same_em_haplotype == same_trio_haplotype,
            an_release1 = an_release1,
            an_release2 = an_release2) %>%
  rowwise() %>%
  mutate(callrate = min(an_release1,an_release2)/246300) %>%
  ggplot(aes(callrate, col=em_correct)) + geom_density()

#Mismatches by gene
gnomad_by_gene = gnomad %>%
  filter(!cpg1 & !cpg2) %>%
  group_by(gene) %>%
  getSummary3()

gnomad_by_gene %>%
  filter(nVariantPairs >= 10) %>%
  ggplot(aes(precisionDiffHapR, recallDiffHapR, col = nHapTrio)) + geom_point()

#Look at prob vs rate on same haplotype
em_by_vp = gnomad %>%
  mutate(vp = paste(chrom1,pos1,ref1, alt1, chrom2, pos2, ref2, alt2)) %>%
  filter(!cpg1 & !cpg2) %>%
  group_by(vp) %>%
    summarise(n = n(),
              ac1 = first(ac_notrios1),
              ac2 = first(ac_notrios2),
              em_prob = first(prob_same_haplotype),
              single_het_ratio= first(single_het_ratio),
              trio_prob = sum(same_trio_haplotype)/n)

em_by_vp %>% filter(n > 1 & !is.na(em_prob) & ac1 > 1000 & ac2  >1000) %>%
  ggplot(aes(single_het_ratio, trio_prob)) + geom_point()

#General properties
gnomad_long %>% 
  group_by(method) %>% 
  getSummary3() %>% 
  glimpse()

gnomad_long %>% 
  group_by(method) %>% 
  getSummary3() %>% 
  glimpse()

#### Plot PR as a function of AF across methods #####
gnomad_af = gnomad_long %>%
  #filter(ac_notrios1 > 0 & ac_notrios2 > 0) %>%
  rowwise() %>%
  mutate(
    min_ac_release = min(ac_release1, ac_release2),
    ac_group = log10(min_ac_release)) %>%
  ungroup() %>%
  mutate(ac_bin = ntile(ac_group, 30)) %>%
  group_by(ac_bin) %>%
  mutate(min_ac_bin = min(min_ac_release),
         max_ac_bin = max(min_ac_release),
         label = sprintf("[%d,%d[",first(min_ac_bin), first(max_ac_bin))) %>%
  ungroup()

ac_bins_labels = gnomad_af %>% 
  group_by(ac_bin) %>% 
  summarise(label = first(label)) %>% 
  arrange(ac_bin) %$% 
  factor(label, levels=label)

#Prepare plotting data
pre_plotting_data = gnomad_af %>%
  #filter(!cpg1 & !cpg2) %>%
  mutate(ac_bin = factor(label, levels = levels(ac_bins_labels))) %>%
  group_by(method, ac_bin) %>%
  getSummary3()

#Low ACs, one by one
pre_plotting_data = gnomad_af %>%
  filter(!cpg1 & !cpg2) %>%
  filter(ac_release1 < 30 | ac_release2 < 30) %>%
  rowwise() %>%
  mutate(ac_bin = min(ac_release1, ac_release2)) %>%
  ungroup() %>%
  group_by(method, ac_bin) %>%
  getSummary3()

plotting_data = data.frame(method=character(),
                           ac_bin=factor(levels = levels(ac_bins_labels)),
                           #ac_bin=numeric(),
                           algo=character(),
                           data=character(),
                           metric=character(),
                           value=double())
for (dataType in c("SameHap", "DiffHap", "")){
  for (metricType in c("precision","recall","f")){
    plotting_data %<>% 
      bind_rows(
          pre_plotting_data %>%
        gather("key","value", one_of(c(paste(metricType,dataType,"EM",sep=""),
                                       paste(metricType,dataType,"R",sep="")))) %>%
        mutate(
                  algo = ifelse(grepl("EM",key),"EM","ratio"),
                  data = ifelse(nchar(dataType)==0,"All",dataType),
                  metric = metricType,
                  value = value
        ) %>%
          mutate_(
            n = paste("n",ifelse(nchar(dataType)==0,"Hap",dataType),"Trio",sep="")
          ) %>%
          select(
            method, ac_bin, algo, data, metric, value, n 
          )
      )
  }
}

plotting_data %>%
  mutate(method = paste(algo, method)) %>% 
  ggplot(aes(ac_bin, value, col = method, size = n)) + 
  geom_point(alpha=0.5) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  facet_grid(data ~ metric)
  
plotting_data %>%
  filter(method == "noPop" & algo == "ratio") %>%
  mutate(method = algo) %>% 
  filter(metric != "f") %>%
  ggplot(aes(ac_bin, value, col = data)) + 
  geom_point(alpha=0.5) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  facet_grid(. ~ metric) +
  xlab("Allele Count bin")

#### Plot PR as a function of AC of both variants #####

#Test as a table 
gnomad_long %>%
  mutate(
    ac_bin1 = case_when(
      ac_notrios1 == 0 ~ "absent",
      ac_notrios1 == 1 ~ "singleton",
      ac_notrios1 == 2 ~ "doubleton",
      ac_notrios1 / an_release1 < 0.001 ~ "AF<0.1%",
      T ~ "AF<1%"
    ),
    ac_bin2 = case_when(
      ac_notrios2 == 0 ~ "absent",
      ac_notrios2 == 1 ~ "singleton",
      ac_notrios2 == 2 ~ "doubleton",
      ac_notrios2 / an_release2 < 0.001 ~ "AF<0.1%",
      T ~ "AF<1%"
    )
  ) %>%
  group_by(method,
           pop,
           ac_bin1,
           ac_bin2) %>%
  getSummary3(compute_EM_stats=F) %>% 
  View()

gnomad_long %>%
  filter(ac_notrios1==1 & ac_notrios2 ==1) %>%
  group_by(method) %>%
  getSummary3(compute_EM_stats=F) %>% 
  View()

gnomad_long %>%
  group_by(pop, method, cpg1 + cpg2) %>%
  getSummary3(compute_EM_stats=F) %>% 
  filter(method =="noPop" & pop %in% c('eas','nfe')) %>%
  View()


gnomad_af2 = gnomad_long %>%
  filter(ac_release1 > 0 & ac_release2 > 0) %>%
  mutate(ac_bin1 = ntile(ac_release1, 100),
         ac_bin2 = ntile(ac_release2, 100)) %>%
  group_by(ac_bin1, ac_bin2) %>%
  mutate(min_ac1_bin = min(ac_release1),
         max_ac1_bin = max(ac_release1),
         min_ac2_bin = min(ac_release2),
         max_ac2_bin = max(ac_release2),
         label1 = sprintf("[%d,%d[",first(min_ac1_bin), first(max_ac1_bin)),
         label2 = sprintf("[%d,%d[",first(min_ac2_bin), first(max_ac2_bin))) %>%
  ungroup()

gnomad_af2 %>% 
  filter(method=="byPop") %>%
  group_by(ac_bin1, ac_bin2) %>% 
  getSummary3() %>%
  ggplot(aes(ac_bin1, ac_bin2, fill = recallDiffHapR)) + geom_hex(stat="identity")

#merge with pLI
gnomad %<>% left_join(pli %>% select(gene,pLI,pRec))
pli = read.delim("fordist_cleaned_exac_r03_march16_z_pli_rec_null_nregions_data.txt", stringsAsFactors = F)
high_pli_genes = pli %>% filter(pLI > 0.9) %>% .$gene %>% unique

#### Phase of variants missing from gnomAD ####
gnomad %>% 

#### Distance vs phase ####
ggplot(gnomad,aes(prob_same_haplotype, col = same_trio_haplotype)) + geom_density() + scale_y_log10()
ggplot(gnomad %>% filter(is.na(ac1) & is.na(ac2)), aes(distance, col=same_trio_haplotype)) + geom_density() + scale_x_log10()


##Plot showing the relationship between distance and chet status
pdf("distance_vs_compound_het.pdf", width=6.5, height = 3)
gnomad %>%
  mutate(distance = 10^round(log10(ifelse(distance==0,1,distance))),
         variants_in_gnomad = factor((ac_notrios1 > 0) + (ac_notrios2 > 0))) %>%
  group_by(distance) %>%
  getSummary3() %>%
  ggplot(aes(distance, 1 - percSameTrioHap)) + 
  geom_point(size=2) + 
  geom_line(size=1.25) + 
  xlab("Distance between variants in pair (log scale)") +
  ylab("Proportion of \ncompound heterozygotes") + 
  scale_x_log10(breaks=c(1,10,100,1000,10000,100000, 1000000),
                labels=c("1b","10b","100b","1kb","10kb","100kb","1mb")) +
  theme_bw(base_size = 14)
dev.off()

##Plot comparing distance/chet variants in/missing from gnomAD
gnomad %>%
  mutate(distance = 10^round(log10(ifelse(distance==0,1,distance))),
         variants_in_gnomad = factor((ac_notrios1 > 0) + (ac_notrios2 > 0))) %>%
  group_by(distance, variants_in_gnomad) %>%
  getSummary3() %>%
  ggplot(aes(distance, 1 - percSameTrioHap, col = variants_in_gnomad)) + 
  geom_point() + 
  geom_line() + 
  xlab("Distance between variants in pair") +
  ylab("Proportion of compound heterozygous pairs") + 
  scale_x_log10() + 
  scale_color_discrete("Number of variants \nin the pair in gnomAD") + 
  theme_bw(base_size = 16) + 
  theme(legend.position = c(0.7,0.2))

  
  
#Plot propotion of compound hets in pairs with one variant in and one variant missing from gnomAD
pdf("chets_missing_in_gnomad.pdf", width=6.5, height = 3)
gnomad %>%
  filter(ac_notrios2 == 0, distance >= 1000) %>%
  mutate(
    ac_bin1 = factor(case_when(
      ac_notrios1 == 0 ~ "Absent",
      ac_notrios1 == 1 ~ "Singleton",
      T ~ "AF<1%"
    ), levels = c("Absent","Singleton","AF<1%"))
  ) %>%
  group_by(ac_bin1) %>%
  getSummary3() %>%
  ggplot(aes(ac_bin1, 1 - percSameTrioHap, fill=ac_bin1)) + 
  geom_bar(stat="identity") + 
  geom_line() + 
  xlab("Frequency of second variant in gnomAD") +
  ylab("Proportion of\ncompound heterozygotes") +
  theme_classic(base_size = 14) +
  scale_fill_brewer() + 
  theme(legend.position="None") +
  geom_text(aes(ac_bin1, 0.2, label=sprintf("N=%d",nHapTrio)), size=4)
  
dev.off()

pdf("chets_missing_in_gnomad_old.pdf", width=6.5, height = 3)
gnomad %>%
  filter(ac_notrios2 == 0) %>%
  mutate(
    ac_bin1 = factor(case_when(
      ac_notrios1 == 0 ~ "absent",
      ac_notrios1 == 1 ~ "singleton",
      ac_notrios1 == 2 ~ "doubleton",
      ac_notrios1 / an_release1 < 0.001 ~ "AF<0.1%",
      T ~ "AF<1%"
    ), levels = c("absent","singleton","doubleton","AF<0.1%","AF<1%")),
    distant = distance >= 1000
  ) %>%
  group_by(ac_bin1, distant) %>%
  getSummary3() %>%
  ggplot(aes(ac_bin1, 1 - percSameTrioHap, col=distant)) + 
  geom_point(size=3) + 
  geom_line() + 
  xlab("Frequency of second variant in gnomAD") +
  ylab("Proportion of compound heterozygotes") +
  scale_color_discrete(name="Inter-variants distance", 
                       labels=c("<1kb",">1kb")) +
  theme(legend.position = c(0.8,0.2))
dev.off()



gnomad%>% 
  mutate(n_variants_missing =  (ac_notrios1 ==0) + (ac_notrios2==0)) %>%
  group_by(n_variants_missing, distance < 1000) %>% getSummary3(F,F) %>% View

gnomad_dist = gnomad %>%
  mutate(n_in_gnomad = as.factor(2-(is.na(ac1) + is.na(ac2))),
    dist_bin = .bincode(distance,breaks=quantile(distance,probs = seq(0,1,0.05))),include.lowest=T) %>%
  group_by(n_in_gnomad) %>%
  mutate(tot = n()) %>%
  group_by(n_in_gnomad, dist_bin) %>%
  summarise(
    n=n(),
    prop_tot = n/first(tot),
    prop_same_hap = sum(same_trio_haplotype)/n,
    min_dist = min(distance),
    max_dist = max(distance)
  )
#x %<>% mutate(dist_bin = ifelse(is.na(dist_bin),1,dist_bin))
ggplot(gnomad_dist,aes(max_dist,prop_same_hap, size = prop_tot, col=n_in_gnomad)) + geom_point() +
  scale_x_log10()

#Look at mutational spectrum
gnomad_mut = gnomad %>%
  filter(alleleType1 == "SNP" & alleleType2 =="SNP" & same_trio_haplotype) %>%
  rowwise() %>%
  mutate(
    mut1 = get_mutation_string(ref1,alt1,cpg1),
    mut2 = get_mutation_string(ref2,alt2,cpg2)
  ) %>%
  select(distance, mut1, mut2, ac1, ac2) %>%
  gather(var, mut, mut1,mut2) %>%
  mutate(mnp = distance < 100, 
         in_gnomad = ifelse(!is.na(ac1) & !is.na(ac2),"At least one variant in gnomAD","Both variants absent from gnomAD"),
         mut = factor(mut,levels= c(
           "C>T@CpG","C>T","A>G","C>A","A>C","A>T","C>G"
         ))) %>% 
  group_by(mnp,in_gnomad) %>%
  mutate(tot = n()) %>%
  group_by(mnp, in_gnomad, mut) %>%
  summarise(n = n(),
            prop = n / first(tot))

ggplot(gnomad_mut,aes(mut, prop,fill=mnp)) + geom_bar(stat="identity", position="dodge") +
  xlab("Substitution") + ylab("Proportion") + 
  facet_grid(. ~ in_gnomad) + 
  geom_text(aes(mut, 0.025, label = sprintf("n = %d",n), angle = 90), position=position_dodge(width=0.9))

ggplot(gnomad_mut %>% filter(in_gnomad== "Both variants absent from gnomAD"),aes(mut, prop,fill=mnp)) + geom_bar(stat="identity", position="dodge") +
  xlab("") + ylab("Proportion") + 
  scale_fill_brewer(name="Distance between variants", labels=c(">1kb","<=1kb")) +
  geom_text(aes(mut, 0.025, label = sprintf("n = %d",n), angle = 90), position=position_dodge(width=0.9)) +
  theme_bw(base_size=16) + 
  theme(legend.position=c(0.7, 0.85))

#### Overall summary ####
#gnomAD -- all Genotypes
gnomad %>% getSummary3() %>% glimpse()
gnomad %>% group_by(pop) %>% getSummary3() %>% glimpse()
gnomad %>% getSummary2()
gnomad %>% mutate(hasCpG = cpg1 | cpg2) %>% group_by(hasCpG, pop) %>% getSummary2() %>% data.frame()
gnomad %>% filter(ac1==0 & ac2==0) %>% getSummary2()
gnomad %>% filter(ac1==0 & ac2==1) %>% getSummary2()
gnomad %>% filter(ac1==1 & ac2==1) %>% getSummary2()
gnomad %>% filter(ac1<100 & ac2 < 100) %>% getSummary2()

#gnomAD -- Adj Genotypes
gnomad_adj %>% getSummary2()
gnomad_adj %>% mutate(hasCpG = cpg1 | cpg2) %>% group_by(pop) %>% getSummary2() %>% data.frame()

#Summary by family
gnomad_by_fam = gnomad %>% group_by(sa.fam.famID) %>% getSummary()

#### Single vs Variant pairs impact ####
single_var = read.delim("gnomad_exomes.single_variants.txt.bgz", stringsAsFactors = F)
var_pairs = read.delim("gnomad_exomes_bySample.variant_pairs.txt.bgz",stringsAsFactors = F)

#Look at LoF per gene
single_var_per_gene = single_var %>%
  filter(pass == "true") %>%
  group_by(pop, gene, impact) %>%
  summarise(
    nUniqueVariants = length(unique(pos)),
    nSampleVariants = n(),
    minAC = min(ac),
    maxAC = max(ac),
    medianAC = median(ac)
  )

single_var_per_gene %>% 
  filter(impact == "high" & gene %in% high_pli_genes) %>% 
  group_by(pop) %>% 
  summarise(genes = n(),
            variants = sum(nUniqueVariants),
            knockouts = sum(nSampleVariants),
            ratio = knockouts / variants)


var_pairs_per_gene = var_pairs %>%
  filter(pass1 == "true" & pass2 == "true") %>%
  group_by(pop, gene, impact1, impact2) %>%
  summarise(
    nUniqueVariants = length(unique(paste(pos1,pos2))),
    nSampleVariants = n(),
    minAC = min(ac1, ac2),
    maxAC = max(ac1, ac2),
    medianAC = median(ac1, ac2)
  )

var_pairs_per_gene %>% 
  filter(impact1 == "high" & impact2 == "high") %>% 
  group_by(pop) %>% 
  summarise(genes = n(),
            variants = sum(nUniqueVariants),
            knockouts = sum(nSampleVariants),
            ratio = knockouts / variants)

var_pairs_per_gene %>% 
  filter(impact1 == "high" & impact2 == "high" & gene %in% high_pli_genes) %>% 
  group_by(pop) %>% 
  summarise(genes = n(),
            variants = sum(nUniqueVariants),
            knockouts = sum(nSampleVariants),
            ratio = knockouts / variants)
                                                                                                                                    nUniqueVariants = sum(nUniqueVariants),
                                                                                                                                    nSampleVariants = sum(nSampleVariants))
ko_per_gene = single_var_per_gene %>%
  filter(impact == "high") %>%
  full_join(
    var_pairs_per_gene %>%
      filter(impact1 == "high" & impact2 == "high"),
    by = "gene"
  ) %>%
  replace_na(list(nUniqueVariants.x = 0, nUniqueVariants.y = 0, nSampleVariants.x = 0, nSampleVariants.y = 0))

ggplot(ko_per_gene, aes(nUniqueVariants.x, nUniqueVariants.y)) + geom_point() + 
  xlab("Number of homozygous PTVs") + ylab("Number of compound het PTVs") 
  
ggplot(ko_per_gene, aes(nUniqueVariants.x - nUniqueVariants.y)) + geom_histogram(binwidth = 1) +
  xlab("Number of homozygous PTVs - Number of compound het PTVs")

#### Investigate genotype counts ####

ggplot(gnomad %>% 
         filter(ac1 > 0 & ac2> 0 & !is.na(same_em_haplotype)) %>%
         mutate(x = ifelse(same_trio_haplotype, "Same trio haplotype","Different trio haplotypes"),
                y = ifelse(same_em_haplotype, "Predicted same haplotype","Predicted different haplotypes"),
                Correct = same_trio_haplotype == same_em_haplotype), 
       aes((AaBB + AABb) / (AaBB + AABb + AaBb), col = Correct)) + 
  geom_density() + 
  facet_grid(y ~ x, scales = "free") + 
  theme_set(theme_bw(base_size = 16))

#global AC for variants misclassified as on the same haplotype
gnomad %>%
  filter(global_same_em_haplotype & !same_trio_haplotype & !same_em_haplotype) %>%
  group_by(pop) %>%
  summarise(
    n_ac0 = sum(other_ac1 * other_ac2 == 0),
    n_AaBb = sum(global_AaBb - AaBb > 0),
    n = n()
  )

gnomad %>%
  filter(global_same_em_haplotype & !same_trio_haplotype & !same_em_haplotype) %>%
  mutate(hasCpG = cpg1 | cpg2) %>%
ggplot(aes(af2, other_af2, col = global_AaBb > AaBb)) + geom_point(alpha=0.3) + geom_abline(slope = 1, intercept = 0) +
  scale_x_continuous(name="Population AF") + scale_y_continuous(name="All other populations AF")

