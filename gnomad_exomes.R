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
  preprocess_gnomad(trios_meta, prob_threshold = 0.25)
gnomad_no_pop = read.delim("exome_trios_full.txt.bgz", stringsAsFactors = F) %>%
  preprocess_gnomad(trios_meta, prob_threshold = 0.25)
gnomad_adj = read.delim("exome_trios_adj.txt.bgz", stringsAsFactors = F) %>%
  preprocess_gnomad(trios_meta, prob_threshold = 0.25)


#Merge numbers from haplotyping using global vs by population
global_pop_columns = c(4,19,34,37:49,98:100)
gnomad_all = gnomad %>% full_join(
  gnomad_no_pop %>% 
    rename_(.dots=setNames(names(.)[global_pop_columns], paste0("global_", names(.)[global_pop_columns])))
) 

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

gnomad_long = gnomad %>% 
  mutate(method = "byPop") %>%
  rbind(
    gnomad_no_pop %>%
      mutate(method = "noPop")
  )

#Look at callrate
gnomad_all %>%
  transmute(ratio_correct = same_ratio_haplotype == same_trio_haplotype,
            em_correct = same_em_haplotype == same_trio_haplotype,
            an_release1 = an_release1,
            an_release2 = an_release2) %>%
  rowwise() %>%
  mutate(callrate = min(an_release1,an_release2)/246300) %>%
  ggplot(aes(callrate, col=em_correct)) + geom_density()


#General properties
gnomad_long %>% 
  filter(an_release1 > 221670 & an_release2 > 221670 & kid_adj & mom_adj & dad_adj) %>%
  group_by(method) %>% 
  getSummary3() %>% 
  glimpse()

#Plot PR as a function of AF
gnomad_af = gnomad_long %>%
  filter(an_release1 > 221670 & an_release2 > 221670 & kid_adj & mom_adj & dad_adj) %>%
  filter(ac_release1 > 0 & ac_release2 > 0) %>%
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
  mutate(ac_bin = factor(label, levels = levels(ac_bins_labels))) %>%
  group_by(method, ac_bin) %>%
  getSummary3()

plotting_data = data.frame(method=character(),
                           ac_bin=factor(levels = levels(ac_bins_labels)),
                           algo=character(),
                           data=character(),
                           metric=character(),
                           value=double())
for (dataType in c("SameHap", "DiffHap", "")){
  for (metricType in c("precision","recall")){
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
          select(
            method, ac_bin, algo, data, metric, value
          )
      )
  }
}

plotting_data %>%
  mutate(method = paste(algo, method)) %>% 
  ggplot(aes(ac_bin, value, col = method)) + 
  geom_point() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  facet_grid(data ~ metric)
  

#merge with pLI
gnomad %<>% left_join(pli %>% select(gene,pLI,pRec))
pli = read.delim("fordist_cleaned_exac_r03_march16_z_pli_rec_null_nregions_data.txt", stringsAsFactors = F)
high_pli_genes = pli %>% filter(pLI > 0.9) %>% .$gene %>% unique

#### Distance vs phase ####
ggplot(gnomad,aes(prob_same_haplotype, col = same_trio_haplotype)) + geom_density() + scale_y_log10()
ggplot(gnomad %>% filter(is.na(ac1) & is.na(ac2)), aes(distance, col=same_trio_haplotype)) + geom_density() + scale_x_log10()

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
  xlab("Substitution") + ylab("Proportion") + 
  scale_fill_discrete(name="Distance < 100bp") +
  geom_text(aes(mut, 0.025, label = sprintf("n = %d",n), angle = 90), position=position_dodge(width=0.9))

#### Overall summary ####
#gnomAD -- all Genotypes
gnomad %>% filter(an_release1 > 221670 & an_release2 > 221670 & kid_adj & mom_adj & dad_adj) %>% getSummary3() %>% glimpse()
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

