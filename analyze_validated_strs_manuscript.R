# checking out the validated str dnms

options(stringsAsFactors = F)
setwd("~/Documents/Harris/Simons/for_michael/")

library(ggplot2)
library(tidyr)
library(dplyr)
library(ggridges)
library(cowplot)

parental_ages <- read.table("../SSC_parental_age_at_birth_in_years.txt", header = T)
ancestry_data <- read.delim("../autism_inferred_ancestry.k8.with_superpops.txt", sep = "\t", header = T)
ancestry_data <-
  ancestry_data %>% separate(eichler.id, c("family", "indiv"), "\\.", remove = F)
all_admixture_data <- read.delim("./full.dataset.for.admixture.autism.8.Q", sep = " ", header = F)
colnames(all_admixture_data) <- paste0("K", c(1:8))
k8_ancestry_mappings <- read.delim("./k8_superpop_mapping.txt", sep = "\t")
colnames(all_admixture_data) <- k8_ancestry_mappings[colnames(all_admixture_data), 1]
all_autism_ped <- read.delim("./full.dataset.for.admixture.autism.fam", sep = " ", header = F)
all_admixture_data <- cbind(id=all_autism_ped$V1, all_admixture_data)
ssc_admixture_data <- subset(all_admixture_data, all_admixture_data$id %in% ancestry_data$X2)
ancestry_data <- merge(ancestry_data, ssc_admixture_data[, c("id", "AFR")], by.x = "X2", by.y = "id", all.x = T)

ped_file <- read.delim("./ssc_copy.ped", header = F, sep = " ")
colnames(ped_file) <- c("family", "Sample", "father.id", "mother.id", "asd_state", "sex", "ssc_id")

# ses_data <- read.delim("./SSC_SES_forMG.csv", sep = ",")
# ses_data$"family" <- sub("\\..*", "", ses_data$individual)

str_data <- read.delim("./dnms_w_filter_info.csv", sep = ",", header = F)
colnames(str_data) <-
  c("chrom", "pos", "period", "prior", "family", "child", 
    "phenotype", "posterior", "newallele", "mutsize", "poocase", 
    "child_gt", "mat_gt", "pat_gt", "encl_child", "encl_mother", 
    "encl_father", "encl_parent", "long_mother", "long_father", 
    "phase", "ref_allele_len", "is_in_asd_gene_coding_region", 
    "len_repeat_unit", "repeat_unit", "reptiming", "reptiming_quartile", 
    "mom_allele_counts", "dad_allele_counts", "pro_allele_counts", 
    "sib_allele_counts", "validation", "all_dnms_discoverable")

str_data$"expansion_or_deletion" <- 
  sapply(
    str_data$mutsize,
    function(x) if(x < 0) "deletion" else "expansion"
  )
# str_data$"reptiming_early_late" <-
#   sapply(
#     str_data$reptiming_quartile,
#     function(x) if(is.na(x)) NA else if(x > 2) "early" else "late"
#   )
# str_data <- merge(str_data, ancestry_data[, c("X2", "SuperPop", "AFR")], by.x = "child", by.y = "X2", all.x = T)
# str_data$"len_repeat_unit_condensed" <- 
#   sapply(str_data$len_repeat_unit,
#          function(x) if(x >= 5) ">=5" else x)
# str_data$len_repeat_unit_condensed <-
#   factor(str_data$len_repeat_unit_condensed, levels = c("1", "2", "3", "4", ">=5"))
str_data$"is_homopolymer" <-
  (str_data$len_repeat_unit == 1)

str_data$mat_homo_het <-
  sapply(
    str_data$mat_gt,
    function(x)
      length(unique(strsplit(x, ",")[[1]]))
  )
str_data$pat_homo_het <-
  sapply(
    str_data$pat_gt,
    function(x)
      length(unique(strsplit(x, ",")[[1]]))
  )
str_data$"mutate_from_homo_or_het" <-
  apply(
    str_data,
    1,
    function(x)
      if(x[["poocase"]] == 2) x[["pat_homo_het"]]
    else if(x[["poocase"]] == 3) x[["mat_homo_het"]]
    else 0
  )
str_data$"n_parents_hets" <-
  str_data$mat_homo_het + str_data$pat_homo_het - 2
# str_data$"accurate_mutsize" <-
#   str_data$mutate_from_homo_or_het == 1 | str_data$mat_gt == str_data$pat_gt
# str_data$"accurate_expansion_or_deletion" <-
#   apply(
#     str_data,
#     1,
#     function(x)
#       if(x[["accurate_mutsize"]]) return(TRUE)
#     else{
#       possible_alleles=c()
#       if(x[["poocase"]] == 4) {
#         possible_alleles = 
#           as.integer(unlist(strsplit(c(x[["mat_gt"]], x[["pat_gt"]]), ",")))
#       }
#       else if(x[["poocase"]] == 2){
#         possible_alleles = 
#           as.integer(unlist(strsplit(x[["pat_gt"]], ",")))
#       }
#       else {
#         possible_alleles = 
#           as.integer(unlist(strsplit(x[["mat_gt"]], ",")))
#       }
#       return(x[["newallele"]] < min(possible_alleles) | 
#                x[["newallele"]] > max(possible_alleles))
#     }
#   )
# str_data$"quality_bin" <-
#   cut(str_data$posterior, breaks = quantile(str_data$posterior, probs = seq(0, 1, by=0.5)))
# str_data$"parents_share_allele" <-
#   apply(
#     str_data,
#     1,
#     function(x)
#       any(unlist(strsplit(x[["mat_gt"]], ",")) %in% 
#             unlist(strsplit(x[["pat_gt"]], ",")))
#   )
str_data$"any_gc" <-
  sapply(
    str_data$repeat_unit,
    function(x)
      nchar(sub("[CG]", "", x)) < nchar(x)
  )

# table(
#   str_data$len_repeat_unit_condensed,
#   str_data$validation
# )

str_data_by_indiv <-
  data.frame(
    ssc_id = unique(str_data$child)
  )

str_data_by_indiv <- 
  merge(
    str_data_by_indiv,
    as.data.frame.table(
      table(
        subset(str_data[, c("child", "poocase")],
               str_data$validation == "true_de_novo")
      )
    ) %>% spread(poocase, Freq),
    by.x = "ssc_id", by.y = "child", all.x = T
  ) %>% 
  rename(., "paternal" = "2", "maternal" = "3", "unphased" = "4") %>%
  mutate(., paternal = replace_na(paternal, 0),
         maternal = replace_na(maternal, 0), 
         unphased = replace_na(unphased, 0), .keep = "unused"
         )
str_data_by_indiv$"total" <- rowSums(str_data_by_indiv[, c("paternal", "maternal", "unphased")])
str_data_by_indiv <-
  merge(
    str_data_by_indiv,
    as.data.frame.table(
      table(
        subset(str_data[, c("child", "poocase")],
               str_data$validation == "true_de_novo" &
                 !str_data$is_homopolymer)
      )
    ) %>% spread(poocase, Freq),
    by.x="ssc_id", by.y="child"
  ) %>% 
  rename(., paternal_nonhomopolymer = "2", maternal_nonhomopolymer = "3", unphased_nonhomopolymer = "4") %>%
  mutate(., paternal_nonhomopolymer = replace_na(paternal_nonhomopolymer, 0),
         maternal_nonhomopolymer = replace_na(maternal_nonhomopolymer, 0),
         unphased_nonhomopolymer = replace_na(unphased_nonhomopolymer, 0),
         .keep = "unused")

str_data_by_indiv <-
  merge(str_data_by_indiv, 
        ped_file[, c("Sample", "father.id", "mother.id", "ssc_id", "sex")], 
        all.x = T,
        by="ssc_id")
str_data_by_indiv <-
  merge(
    str_data_by_indiv,
    parental_ages,
    by.x="Sample", by.y="sample",
    all.x = T
  )

# wilcox.test(subset(str_data_by_indiv$motherAgeAtBirth, str_data_by_indiv$sex == "M"),
#             subset(str_data_by_indiv$motherAgeAtBirth, str_data_by_indiv$sex == "F"))
# no detected difference in parental age at birth between sexes of kids, not going to worry about excluding chrX from analyses

# phased deletion or expansions
str_data_by_indiv <-
  merge(
    str_data_by_indiv,
    as.data.frame.table(
      table(
        subset(str_data[, c("child", "expansion_or_deletion", "poocase")], 
               str_data$poocase %in% c(2, 3) & 
                 str_data$validation == "true_de_novo" &
                 !str_data$is_homopolymer)
      )
    ) %>% mutate(expansion_or_deletion_poocase=paste(expansion_or_deletion, poocase, sep = "_"), 
                 .keep = "unused") %>% 
      spread(expansion_or_deletion_poocase, Freq),
    by.x = "ssc_id", by.y = "child", all.x = T
  ) %>% mutate(
    ., deletion_2 = replace_na(deletion_2, 0),
    deletion_3 = replace_na(deletion_3, 0),
    expansion_2 = replace_na(expansion_2, 0),
    expansion_3 = replace_na(expansion_3, 0),
    .keep = "unused"
  )
str_data_by_indiv$"motherAge_quantile" <-
  cut(str_data_by_indiv$motherAgeAtBirth, 
      breaks = quantile(str_data_by_indiv$motherAgeAtBirth, probs = seq(0, 1, by=0.2), na.rm = T))
str_data_by_indiv$"fatherAge_quantile" <-
  cut(str_data_by_indiv$fatherAgeAtBirth, 
      breaks = quantile(str_data_by_indiv$fatherAgeAtBirth, probs = seq(0, 1, by=0.2), na.rm = T))
str_data_by_indiv$"family" <-
  sub("\\..*", "", str_data_by_indiv$Sample)

str_data_by_indiv <- 
  merge(
    str_data_by_indiv,
    as.data.frame.table(
      table(
        str_data$child
      )
    ) %>% rename(., c("ssc_id" = "Var1", "total_unfiltered" = "Freq")),
    by = "ssc_id", all.x = T
  )
str_data <-
  merge(str_data,
        str_data_by_indiv[, c("ssc_id", "fatherAgeAtBirth", "motherAgeAtBirth")],
        by.x = "child", by.y = "ssc_id", all.x = T)

# validation_data <- 
#   read.delim("./mitra_et_al_2021_validations.txt")
# validation_data <- 
#   rename(validation_data, c("validation_family_number" = "Family",
#                                       "chrom" = "Chrom",
#                                       "pos" = "TR.start.position..hg38."))
# validation_data$"family" <-
#   sapply(validation_data$validation_family_number,
#          function(x)
#            if (x == 1) 14686
#          else if (x == 2) 11963
#          else if (x == 3) 12605
#          else if (x == 4) 14329
#          else if (x == 5) 11450)
# 
# str_data$"is_in_mitra_validation" <-
#   sapply(
#     c(1:nrow(str_data)),
#     function(x)
#       if(!str_data$family[x] %in% validation_data$family) FALSE
#     else
#       nrow(plyr::match_df(validation_data[, c("chrom", "pos", "family")], 
#                     str_data[x, c("chrom", "pos", "family")])) > 0
#   )
# str_data$"pass_mitra_validation" <-
#   sapply(
#     c(1:nrow(str_data)),
#     function(x)
#       if(!str_data$family[x] %in% validation_data$family) NA
#     else
#       if(!str_data$pos[x] %in% subset(validation_data$pos, validation_data$family == str_data$family[x]))
#         NA
#     else
#       str_data$pos[x] %in% subset(validation_data$pos, validation_data$Mutation.validated == "Y")
#   )

# ggsave(
#   "../figures/mitra_vs_read1_flank2_validations_20220816.pdf",
#   ggplot(
#     data = subset(str_data, str_data$is_in_mitra_validation),
#     aes(x = pass_mitra_validation)
#   ) +
#     geom_bar(aes(fill = validation), position = "fill")
# )

ggsave(
  "../figures/frac_from_homozygous_str_dnms_one_parent_het_20230725.pdf",
  ggplot(
    data = 
      as.data.frame.table(
        table(
          str_data[, c("child", "n_parents_hets", "mutate_from_homo_or_het")] %>%
            filter(., n_parents_hets == 1 & mutate_from_homo_or_het != 0)
        ) 
      ) %>% spread(mutate_from_homo_or_het, Freq) %>%
      rename(., c("homo" = "1", "het" = "2")) %>%
      mutate(., frac_homozygous = homo / (homo + het), .keep = "unused") 
    ,
    aes(x = frac_homozygous)
  ) + geom_histogram() +
    theme_minimal()
)

ggsave(
  "../figures/frac_from_homozygous_str_dnms_one_parent_het_validated_nonhomopolymer_20230725.pdf",
  ggplot(
    data = 
      as.data.frame.table(
        table(
          str_data[, c("child", "n_parents_hets", "mutate_from_homo_or_het", "validation", "is_homopolymer")] %>%
            filter(., n_parents_hets == 1 & mutate_from_homo_or_het != 0 & validation == "true_de_novo" & !is_homopolymer)
        ) 
      ) %>% spread(mutate_from_homo_or_het, Freq) %>%
      rename(., c("homo" = "1", "het" = "2")) %>%
      mutate(., frac_homozygous = homo / (homo + het), .keep = "unused") 
    ,
    aes(x = frac_homozygous)
  ) + geom_histogram() +
    theme_minimal()
)

summary(
  glm(
    maternal_nonhomopolymer ~ motherAgeAtBirth, 
    data = str_data_by_indiv, 
    family = poisson(link = "identity")
  )
)
summary(
  glm(
    paternal_nonhomopolymer ~ fatherAgeAtBirth, 
    data = str_data_by_indiv, 
    family = poisson(link = "identity")
  )
)
summary(
  glm(
    maternal_nonhomopolymer ~ motherAgeAtBirth + fatherAgeAtBirth, 
    data = str_data_by_indiv, 
    family = poisson(link = "identity")
  )
)

summary(
  glm(
    maternal_nonhomopolymer ~ motherAgeAtBirth, 
    data = str_data_by_indiv, 
    family = poisson(link = "log")
  )
)
summary(
  glm(
    paternal_nonhomopolymer ~ fatherAgeAtBirth, 
    data = str_data_by_indiv, 
    family = poisson(link = "log")
  )
)

library(AER)
dispersiontest(
  glm(
    maternal_nonhomopolymer ~ motherAgeAtBirth, 
    data = str_data_by_indiv, 
    family = poisson(link = "identity")
  )
)

ggsave(
  "../figures/mat_str_vs_mat_age_20230605.pdf",
  ggplot(
    data = 
      str_data_by_indiv,
    aes(x = motherAgeAtBirth, y = maternal_nonhomopolymer)
  ) +
    geom_jitter(height = 0.1, alpha = 0.1) +
    theme_cowplot() +
    scale_y_continuous(breaks=seq(0,9,1)) +
    stat_smooth(method = "glm", formula = y ~ x, method.args = list(family = poisson(link = "identity"))),
  width = 7, height = 5
)
ggsave(
  "../figures/pat_str_vs_pat_age_20230725.pdf",
  ggplot(
    data = 
      str_data_by_indiv,
    aes(x = fatherAgeAtBirth, y = paternal_nonhomopolymer)
  ) +
    geom_jitter(height = 0.1, alpha = 0.1) +
    theme_cowplot() +
    scale_y_continuous(breaks=seq(0,24,4)) +
    stat_smooth(method = "glm", formula = y ~ x, method.args = list(family = poisson(link = "identity"))),
  width = 7, height = 5
)

ggsave(
  "../figures/mat_str_vs_mat_age_homopolymer_vs_nonhomopolymer_20230605.pdf",
  ggplot(
    data = 
      merge(
        as.data.frame.table(
          table(
            subset(
              str_data[, c("child", "is_homopolymer")],
              str_data$validation == "true_de_novo" &
                str_data$poocase == 3
            )
          )
        ) %>%
          pivot_wider(., names_from = "is_homopolymer", values_from = "Freq"),
        str_data_by_indiv[, c("ssc_id", "motherAgeAtBirth")],
        by.x = "child", by.y = "ssc_id", all.x = T, all.y = T) %>%
      pivot_longer(., cols = c("FALSE", "TRUE"), names_to = "is_homopolymer", values_to = "Freq") %>%
      mutate(., is_homopolymer = replace_na(is_homopolymer, 0)), 
    aes(x = motherAgeAtBirth, y = Freq)
  ) +
    geom_jitter(aes(color = is_homopolymer), height = 0.1, alpha = 0.1) +
    theme_cowplot() +
    scale_y_continuous(breaks = seq(0, 9, 1)) +
    stat_smooth(aes(color = is_homopolymer), method = "glm", formula = y ~ x, method.args = list(family = poisson(link = "identity"))),
  width = 7, height = 5
)
ggsave(
  "../figures/pat_str_vs_pat_age_homopolymer_vs_nonhomopolymer_20230607.pdf",
  ggplot(
    data = 
      merge(
        as.data.frame.table(
          table(
            subset(
              str_data[, c("child", "is_homopolymer")],
              str_data$validation == "true_de_novo" &
                str_data$poocase == 2
            )
          )
        ) %>%
          pivot_wider(., names_from = "is_homopolymer", values_from = "Freq"),
        str_data_by_indiv[, c("ssc_id", "fatherAgeAtBirth")],
        by.x = "child", by.y = "ssc_id", all.x = T, all.y = T) %>%
      pivot_longer(., cols = c("FALSE", "TRUE"), names_to = "is_homopolymer", values_to = "Freq") %>%
      mutate(., is_homopolymer = replace_na(is_homopolymer, 0)), 
    aes(x = fatherAgeAtBirth, y = Freq)
  ) +
    geom_jitter(aes(color = is_homopolymer), height = 0.1, alpha = 0.1) +
    theme_cowplot() +
    scale_y_continuous(breaks = seq(0, 20, 4)) +
    stat_smooth(aes(color = is_homopolymer), method = "glm", formula = y ~ x, method.args = list(family = poisson(link = "identity"))),
  width = 7, height = 5
)

summary(
  glm(
    Freq ~ fatherAgeAtBirth, family = poisson(link = "identity"),
    data = 
      merge(
        as.data.frame.table(
          table(
            subset(
              str_data$child,
              str_data$validation == "true_de_novo" &
                str_data$poocase == 2 &
                str_data$is_homopolymer
            )
          )
        ),
        str_data_by_indiv[, c("ssc_id", "fatherAgeAtBirth")],
        by.x = "Var1", by.y = "ssc_id", all.x = T, all.y = T) %>%
      mutate(., Freq = replace_na(Freq, 0))
  )
)
summary(
  glm(
    Freq ~ motherAgeAtBirth, family = poisson(link = "identity"),
    data = 
      merge(
        as.data.frame.table(
          table(
            subset(
              str_data$child,
              str_data$validation == "true_de_novo" &
                str_data$poocase == 3 &
                str_data$is_homopolymer
            )
          )
        ),
        str_data_by_indiv[, c("ssc_id", "motherAgeAtBirth")],
        by.x = "Var1", by.y = "ssc_id", all.x = T, all.y = T) %>%
      mutate(., Freq = replace_na(Freq, 0))
  )
)
# Neither maternal nor paternal homopolymer mutation rates are associated with age

# Expansions or deletions more significantly affected by parental age?
expansion_deletion_parental_age_matrix <-
  rbind(
    rbind(
      str_data_by_indiv[, c("ssc_id", "deletion_2", "fatherAgeAtBirth")] %>%
        rename(., c("Freq" = "deletion_2")),
      str_data_by_indiv[, c("ssc_id", "expansion_2", "fatherAgeAtBirth")] %>%
        rename(., c("Freq" = "expansion_2"))
    ) %>% rename(., c("parental_age" = "fatherAgeAtBirth")),
    rbind(
      str_data_by_indiv[, c("ssc_id", "deletion_3", "motherAgeAtBirth")] %>%
        rename(., c("Freq" = "deletion_3")),
      str_data_by_indiv[, c("ssc_id", "expansion_3", "motherAgeAtBirth")] %>%
        rename(., c("Freq" = "expansion_3"))
    ) %>% rename(., c("parental_age" = "motherAgeAtBirth"))
  )
expansion_deletion_parental_age_matrix$"deletion_or_expansion" <-
  rep(rep(c("expansion", "deletion"), each = nrow(str_data_by_indiv)), times = 2)
expansion_deletion_parental_age_matrix$"parental_sex" <-
  rep(c("father", "mother"), each = nrow(str_data_by_indiv) * 2)

summary(
  glm(
    Freq ~ deletion_or_expansion * parental_age,
    data = subset(expansion_deletion_parental_age_matrix, expansion_deletion_parental_age_matrix$parental_sex == "mother"),
    family = poisson(link = "identity")
  )
)
summary(
  glm(
    Freq ~ deletion_or_expansion * parental_age,
    data = subset(expansion_deletion_parental_age_matrix, expansion_deletion_parental_age_matrix$parental_sex == "father"),
    family = poisson(link = "identity")
  )
)
ggsave(
  "../figures/expansion_deletion_parental_age_20230607.pdf",
  ggplot(
    data = expansion_deletion_parental_age_matrix,
    aes(x = parental_age, y = Freq)
  ) +
    geom_jitter(aes(color = deletion_or_expansion), height = 0.1, alpha = 0.05) +
    theme_cowplot() +
    stat_smooth(formula = y ~ x, method = "glm", aes(color = deletion_or_expansion)) +
    facet_grid(rows = vars(parental_sex), scales = "free_y")
)
# No

# Comparison to Sun et al. 2012
kong_markers <- read.delim("./kong_et_al_2002_markers_hg38.txt", header = T)

# note: looks like two of the markers are called as separate loci in the Kong list but now map to a single locus in hg38
# subsetting str data for specific analysis because the merge function is very slow
# no homopolymers in the kong panel
str_data_for_kong_analysis <-
  subset(str_data, 
         str_data$validation == "true_de_novo" &
           str_data$poocase %in% c(2, 3) & 
           !str_data$is_homopolymer)
str_data_for_kong_analysis$"chrom_pos" <-
  paste(str_data_for_kong_analysis$chrom,
        str_data_for_kong_analysis$pos, sep = "_")


# set.seed(0)
# subset_kong <-
#   kong_markers[sample(c(1:nrow(kong_markers)), size = 2477, replace = F), ]
# 
# str_data_for_kong_analysis$"overlap_kong" <-
#   apply(
#     str_data_for_kong_analysis,
#     1,
#     function(x)
#       nrow(
#         filter(kong_markers, X.chrom == x[["chrom"]] & start <= as.numeric(x[["pos"]]) & end >= as.numeric(x[["pos"]]))
#       )
#   )
# str_data_for_kong_analysis$"overlap_kong_subset" <-
#   apply(
#     str_data_for_kong_analysis,
#     1,
#     function(x)
#       nrow(
#         filter(subset_kong, X.chrom == x[["chrom"]] & start <= as.numeric(x[["pos"]]) & end >= as.numeric(x[["pos"]]))
#       )
#   )


# kong_subset_str_data_by_indiv <-
#   merge(
#     str_data_by_indiv[, c("ssc_id", "motherAgeAtBirth", "fatherAgeAtBirth")],
#     as.data.frame.table(
#       table(
#         subset(str_data_for_kong_analysis[, c("child", "poocase")], 
#                str_data_for_kong_analysis$overlap_kong >= 1 &
#                  str_data_for_kong_analysis$validation == "true_de_novo")
#       )
#     ) %>%
#       spread(poocase, Freq) %>%
#       rename(., "paternal" = "2",
#              "maternal" = "3"),
#     by.x = "ssc_id", by.y = "child", all.x = T
#   ) %>%
#   mutate(., paternal = replace_na(paternal, 0),
#          maternal = replace_na(maternal, 0), .keep = "unused")

# kong_subset_str_data_by_indiv <-
#   merge(
#     str_data_by_indiv[, c("ssc_id", "motherAgeAtBirth", "fatherAgeAtBirth")],
#     as.data.frame.table(
#       table(
#         subset(str_data_for_kong_analysis[, c("child", "poocase")], 
#                str_data_for_kong_analysis$overlap_kong_subset >= 1 &
#                  str_data_for_kong_analysis$validation == "true_de_novo")
#       )
#     ) %>%
#       spread(poocase, Freq) %>%
#       rename(., "paternal_subset" = "2",
#              "maternal_subset" = "3"),
#     by.x = "ssc_id", by.y = "child", all.x = T
#   ) %>%
#   mutate(., paternal_subset = replace_na(paternal_subset, 0),
#          maternal_subset = replace_na(maternal_subset, 0), .keep = "unused")

# str_panel_w_kong <- read.delim("./sites_annotated_reptiming_dnm_counts_ru_borders_kong.csv", header = F)
# names(str_panel_w_kong) <-
#   c("chrom", "pos", "repeat_unit", "len_ref_ru", "reptiming", "n_dnms", "X", "border_homology", "in_kong")
# 
# denom_kong <- nrow(subset(str_panel_w_kong,
#                           str_panel_w_kong$in_kong & nchar(str_panel_w_kong$repeat_unit) > 1))
# denom_kong_subset <-
#   sum(
#     apply(
#       subset(
#         str_panel_w_kong,
#         str_panel_w_kong$in_kong
#       ),
#       1,
#       function(x)
#         nrow(
#           filter(subset_kong, X.chrom == x[["chrom"]] & start <= as.numeric(x[["pos"]]) & end >= as.numeric(x[["pos"]]))
#         )
#     )
#   )
# 
# kong_subset_str_data_by_indiv$"maternal_rate" <-
#   kong_subset_str_data_by_indiv$maternal / denom_kong
# kong_subset_str_data_by_indiv$"paternal_rate" <-
#   kong_subset_str_data_by_indiv$paternal / denom_kong
# 
# kong_subset_str_data_by_indiv$"maternal_subset_rate" <-
#   kong_subset_str_data_by_indiv$maternal_subset / denom_kong_subset
# kong_subset_str_data_by_indiv$"paternal_subset_rate" <-
#   kong_subset_str_data_by_indiv$paternal_subset / denom_kong_subset

# summary(
#   glm(
#     maternal_rate ~ motherAgeAtBirth, data = kong_subset_str_data_by_indiv
#   )
# )
# # significant now!
# summary(
#   glm(
#     maternal_subset_rate ~ motherAgeAtBirth, data = kong_subset_str_data_by_indiv
#   )
# )
# # but subset is not significant; P value unaffected by changes in denom (which is not accurate above)
# summary(
#   glm(
#     paternal_rate ~ fatherAgeAtBirth, data = kong_subset_str_data_by_indiv
#   )
# )

# kong_subset_str_data_by_indiv$"maternal_age_decile" <-
#   cut(kong_subset_str_data_by_indiv$motherAgeAtBirth, 
#       breaks = quantile(kong_subset_str_data_by_indiv$motherAgeAtBirth, probs = seq(0, 1, by=0.1), na.rm = T))
# kong_subset_str_data_by_indiv$"paternal_age_decile" <-
#   cut(kong_subset_str_data_by_indiv$fatherAgeAtBirth, 
#       breaks = quantile(kong_subset_str_data_by_indiv$fatherAgeAtBirth, probs = seq(0, 1, by=0.1), na.rm = T))
# 
# kong_subset_str_data_by_indiv$"median_maternal_age" <-
#   quantile(kong_subset_str_data_by_indiv$motherAgeAtBirth,
#            probs = seq(0.05, 0.95, by=0.1), na.rm = T)[as.numeric(kong_subset_str_data_by_indiv$maternal_age_decile)]
# kong_subset_str_data_by_indiv$"median_paternal_age" <-
#   quantile(kong_subset_str_data_by_indiv$fatherAgeAtBirth,
#            probs = seq(0.05, 0.95, by=0.1), na.rm = T)[as.numeric(kong_subset_str_data_by_indiv$paternal_age_decile)]

# ggsave(
#   "../figures/kong_markers_in_ssc_families_parental_age_effect_20230510.pdf",
#   ggplot(
#     data =
#       rbind.data.frame(
#         data.frame(
#           rate = kong_subset_str_data_by_indiv$paternal_rate,
#           median_age = kong_subset_str_data_by_indiv$median_paternal_age,
#           age = kong_subset_str_data_by_indiv$fatherAgeAtBirth,
#           lineage = "paternal"
#         ),
#         data.frame(
#           rate = kong_subset_str_data_by_indiv$maternal_rate,
#           median_age = kong_subset_str_data_by_indiv$median_maternal_age,
#           age = kong_subset_str_data_by_indiv$motherAgeAtBirth,
#           lineage = "maternal"
#         )
#       ),
#     aes(x = median_age, y = rate, color = lineage, group = lineage)
#   ) + 
#     stat_summary(fun = mean,
#                  fun.min = function(x) max(mean(x) - sd(x), 0), 
#                  fun.max = function(x) mean(x) + sd(x), 
#                  geom = "pointrange") +
#     stat_smooth()
# )
# 
# ggsave(
#   "../figures/kong_markers_subset_in_ssc_families_parental_age_effect_20230605.pdf",
#   ggplot(
#     data =
#       rbind.data.frame(
#         data.frame(
#           rate = kong_subset_str_data_by_indiv$paternal_subset_rate,
#           median_age = kong_subset_str_data_by_indiv$median_paternal_age,
#           age = kong_subset_str_data_by_indiv$fatherAgeAtBirth,
#           lineage = "paternal"
#         ),
#         data.frame(
#           rate = kong_subset_str_data_by_indiv$maternal_subset_rate,
#           median_age = kong_subset_str_data_by_indiv$median_maternal_age,
#           age = kong_subset_str_data_by_indiv$motherAgeAtBirth,
#           lineage = "maternal"
#         )
#       ),
#     aes(x = median_age, y = rate, color = lineage, group = lineage)
#   ) + 
#     stat_summary(fun = mean,
#                  fun.min = function(x) max(mean(x) - sd(x), 0), 
#                  fun.max = function(x) mean(x) + sd(x), 
#                  geom = "pointrange") +
#     stat_smooth()
# )

# # Bootstrap down to 2477
# str_panel_w_kong$"chrom_pos" <-
#   paste(str_panel_w_kong$chrom,
#         str_panel_w_kong$pos, sep = "_")
# str_panel_w_kong_subset <-
#   subset(
#     str_panel_w_kong,
#     str_panel_w_kong$in_kong
#   )
# kong_markers_subset <- 
#   kong_markers[sample(c(1:nrow(kong_markers)), size = 2477, replace = F),]
# str_panel_w_kong_subset$"curr_subset" <-
#   apply(
#     str_panel_w_kong_subset,
#     1,
#     function(x)
#       nrow(
#         filter(kong_markers_subset, X.chrom == x[["chrom"]] & start <= as.numeric(x[["pos"]]) & end >= as.numeric(x[["pos"]]))
#       )
#   )
# str_data_for_kong_analysis$"curr_kong_subset" <-
#   str_data_for_kong_analysis$chrom_pos %in% str_panel_w_kong_subset$chrom_pos[str_panel_w_kong_subset$curr_subset > 0]

ac_sites_in_kong <- read.delim("./ac_dinucs_overlapping_kong.txt", header = F) %>%
  rename(., chrom = "V1", pos = "V2") %>%
  mutate(., chrom_pos = paste(chrom, pos, sep = "_"))

str_data_for_kong_analysis$"in_ac_kong_overlap" <- 
  str_data_for_kong_analysis$chrom_pos %in% ac_sites_in_kong$chrom_pos

kong_str_data_by_indiv <-
  merge(
    str_data_by_indiv[, c("ssc_id", "motherAgeAtBirth", "fatherAgeAtBirth")],
    as.data.frame.table(
      table(
        subset(str_data_for_kong_analysis[, c("child", "poocase")], 
               str_data_for_kong_analysis$in_ac_kong_overlap)
      )
    ) %>%
      spread(poocase, Freq) %>%
      rename(., "paternal" = "2",
             "maternal" = "3"),
    by.x = "ssc_id", by.y = "child", all.x = T
  ) %>%
  mutate(., paternal_subset = replace_na(paternal, 0),
         maternal_subset = replace_na(maternal, 0), .keep = "unused")

kong_str_data_by_indiv$"maternal_rate" <-
  kong_str_data_by_indiv$maternal / nrow(ac_sites_in_kong)
kong_str_data_by_indiv$"paternal_rate" <-
  kong_str_data_by_indiv$paternal / nrow(ac_sites_in_kong)

kong_str_data_by_indiv$"maternal_age_decile" <-
  cut(kong_str_data_by_indiv$motherAgeAtBirth, 
      breaks = quantile(kong_str_data_by_indiv$motherAgeAtBirth, probs = seq(0, 1, by=0.1), na.rm = T))
kong_str_data_by_indiv$"paternal_age_decile" <-
  cut(kong_str_data_by_indiv$fatherAgeAtBirth, 
      breaks = quantile(kong_str_data_by_indiv$fatherAgeAtBirth, probs = seq(0, 1, by=0.1), na.rm = T))

kong_str_data_by_indiv$"median_maternal_age" <-
  quantile(kong_str_data_by_indiv$motherAgeAtBirth,
           probs = seq(0.05, 0.95, by=0.1), na.rm = T)[as.numeric(kong_str_data_by_indiv$maternal_age_decile)]
kong_str_data_by_indiv$"median_paternal_age" <-
  quantile(kong_str_data_by_indiv$fatherAgeAtBirth,
           probs = seq(0.05, 0.95, by=0.1), na.rm = T)[as.numeric(kong_str_data_by_indiv$paternal_age_decile)]

summary(
  glm(
    maternal_rate ~ motherAgeAtBirth, data = kong_str_data_by_indiv
  )
)
summary(
  glm(
    paternal_rate ~ fatherAgeAtBirth, data = kong_str_data_by_indiv
  )
)

ggsave(
  "../figures/kong_markers_in_ssc_families_parental_age_effect_ac_overlap_20230729.pdf",
  ggplot(
    data =
      rbind.data.frame(
        data.frame(
          rate = kong_str_data_by_indiv$paternal_rate,
          median_age = kong_str_data_by_indiv$median_paternal_age,
          age = kong_str_data_by_indiv$fatherAgeAtBirth,
          lineage = "paternal"
        ),
        data.frame(
          rate = kong_str_data_by_indiv$maternal_rate,
          median_age = kong_str_data_by_indiv$median_maternal_age,
          age = kong_str_data_by_indiv$motherAgeAtBirth,
          lineage = "maternal"
        )
      ),
    aes(x = median_age, y = rate, color = lineage, group = lineage)
  ) + 
    stat_summary(fun = mean,
                 fun.min = function(x) max(mean(x) - sd(x), 0), 
                 fun.max = function(x) mean(x) + sd(x), 
                 geom = "pointrange") +
    stat_smooth(method = "glm", formula = y ~ x) + 
    theme_minimal()
)

# Alpha analysis
summary(glm(
  (paternal_nonhomopolymer / (paternal_nonhomopolymer + maternal_nonhomopolymer)) ~ fatherAgeAtBirth,
  data = subset(str_data_by_indiv, 
                abs((str_data_by_indiv$motherAgeAtBirth - str_data_by_indiv$fatherAgeAtBirth) /
                      str_data_by_indiv$fatherAgeAtBirth) < 0.1),
  family = binomial(link = "logit")
))
# Binomial model is not a good fit here
summary(
  glm(
    (paternal / (paternal + maternal)) ~ fatherAgeAtBirth,
    data = subset(str_data_by_indiv, 
                  abs((str_data_by_indiv$motherAgeAtBirth - str_data_by_indiv$fatherAgeAtBirth) /
                        str_data_by_indiv$fatherAgeAtBirth) < 0.1),
    family = quasibinomial(link = "identity")
  )
)
# But quasibinomial with identity link is a good fit with positive correlation

ggsave(
  "../figures/alpha_by_paternal_age_quasibinomial_identity_20230729.pdf",
  ggplot(
    data = 
      subset(str_data_by_indiv, 
             abs((str_data_by_indiv$motherAgeAtBirth - str_data_by_indiv$fatherAgeAtBirth) /
                   str_data_by_indiv$fatherAgeAtBirth) < 0.1),
    aes(x = fatherAgeAtBirth, y = paternal_nonhomopolymer / (maternal_nonhomopolymer + paternal_nonhomopolymer))
  ) + 
    geom_jitter(height = 0.1, alpha = 0.2) +
    scale_y_continuous(limits = c(0, 1), oob = scales::squish) +
    theme_minimal() + 
    stat_smooth(method = "glm", formula = y ~ x, method.args = list(family = quasibinomial(link = "identity")))
)

# Maternal hotspot
maternal_hotspots_table <- read.delim("./nature24018-s3.txt", header = T)
maternal_hotspots_table <- subset(maternal_hotspots_table, maternal_hotspots_table$CG_enriched_region_indicator == 1)
maternal_hotspots_table$"start" <- maternal_hotspots_table$Pos_1Mb_window_number * 1000000
maternal_hotspots_table$"end" <- (maternal_hotspots_table$Pos_1Mb_window_number + 1) * 1000000 - 1

str_data$"is_in_maternal_hotspot" <-
  apply(
    str_data,
    1,
    function(x)
      nrow(
        subset(
          maternal_hotspots_table,
          maternal_hotspots_table$Chr == x[["chrom"]] &
            maternal_hotspots_table$start < as.numeric(x[["pos"]]) &
            maternal_hotspots_table$end >= as.numeric(x[["pos"]])
        )
      ) > 0
  )
str_data_by_indiv$"maternal_hotspot_maternal_muts" <-
  sapply(
    str_data_by_indiv$ssc_id,
    function(x)
      nrow(
        subset(
          str_data,
          str_data$is_in_maternal_hotspot == 1 &
            str_data$child == x &
            str_data$poocase == 3 &
            str_data$validation == "true_de_novo" &
            str_data$len_repeat_unit > 1
        )
      )
  )

ggsave(
  "../figures/maternal_hotspot_str_dnms_by_maternal_age_20230607.pdf",
  ggplot(
    data = str_data_by_indiv,
    aes(x = motherAgeAtBirth, y = maternal_hotspot_maternal_muts)
  ) +
    geom_jitter(width = 0.1, height = 0.1, alpha = 0.1) +
    theme_cowplot()
)

cor.test(
  str_data_by_indiv$maternal_hotspot_maternal_muts, str_data_by_indiv$motherAgeAtBirth,
)

summary(
  glm(
    maternal_hotspot_maternal_muts ~ motherAgeAtBirth,
    data = str_data_by_indiv,
    family = poisson(link = "identity")
  )
)
# No significant maternal age effect

# GC content
# Need to reweight
str_panel <- read.delim("./sites_annotated_reptiming_dnm_counts_ru.txt", header = F, na.strings = "NA")
names(str_panel) <-
  c("chrom", "pos", "repeat_unit", "len_ref_ru", "reptiming", "n_dnms")
str_panel$"any_gc" <-
  sapply(
    str_panel$repeat_unit,
    function(x)
      nchar(sub("[GC]", "", x)) < nchar(x)
  )
str_panel$"gc_content" <-
  sapply(str_panel$repeat_unit,
         function(x)
           nchar(gsub("[AT]", "", x)) /
           nchar(x))
ggsave(
  "../figures/nonhomopolymer_panel_gc_content_20230423.pdf",
  ggplot(
    data = str_panel,
    aes(x = gc_content)
  ) + geom_histogram(bins = 5, color = "white", fill = "gray") +
    theme_cowplot() + 
    geom_vline(xintercept = 0.125, linetype = "dashed")
)

n_any_gc <- nrow(subset(str_panel, nchar(str_panel$repeat_unit) > 1 & str_panel$any_gc))
n_at <- nrow(subset(str_panel, nchar(str_panel$repeat_unit) > 1 & !str_panel$any_gc))
rm(str_panel)

pat_str_dnm_by_nuc_content <-
  merge(
    with(
      subset(str_data, 
             str_data$validation == "true_de_novo" &
               !str_data$is_homopolymer &
               str_data$poocase == 2),
      as.data.frame.table(table(child, any_gc))
    ) %>%
      pivot_wider(., names_from = any_gc, values_from = Freq),
    str_data_by_indiv[, c("ssc_id", "fatherAgeAtBirth")],
    by.x = "child", by.y = "ssc_id", all.x = T, all.y = T
  ) %>%
  pivot_longer(., cols = c("FALSE", "TRUE"), names_to = "any_gc", values_to = "Freq") %>%
  mutate(., Freq = replace_na(Freq, 0))
pat_str_dnm_by_nuc_content$freq_reweight <-
  apply(
    pat_str_dnm_by_nuc_content,
    1,
    function(x)
      if(x[["any_gc"]]) as.numeric(x[["Freq"]]) * (n_at / n_any_gc)
    else as.numeric(x[["Freq"]])
  )

library(wesanderson)
ggsave(
  "../figures/nonhomopolymer_gc_content_reweight_by_paternal_age_20230607.pdf",
  ggplot(
    data = pat_str_dnm_by_nuc_content,
    aes(x = fatherAgeAtBirth, y = freq_reweight)
  ) + 
    geom_jitter(aes(color = any_gc), height = 0.2, alpha = 0.1) +
    stat_smooth(aes(color = any_gc), formula = y ~ x,
                method = "glm") +
    theme_cowplot() +
    scale_y_continuous(limits = c(0, 7), oob = scales::squish) +
    scale_color_manual(values = wes_palette("Royal1")),
  width = 7, height = 5
)

mat_str_dnm_by_nuc_content <-
  merge(
    with(
      subset(str_data, 
             str_data$validation == "true_de_novo" &
               !str_data$is_homopolymer & 
               str_data$poocase == 3),
      as.data.frame.table(table(child, any_gc))
    )  %>%
      pivot_wider(., names_from = any_gc, values_from = Freq),
    str_data_by_indiv[, c("ssc_id", "motherAgeAtBirth")],
    by.x = "child", by.y = "ssc_id", all.x = T, all.y = T
  ) %>%
  pivot_longer(., cols = c("FALSE", "TRUE"), names_to = "any_gc", values_to = "Freq") %>%
  mutate(., Freq = replace_na(Freq, 0))
mat_str_dnm_by_nuc_content$freq_reweight <-
  apply(
    mat_str_dnm_by_nuc_content,
    1,
    function(x)
      if(x[["any_gc"]]) as.numeric(x[["Freq"]]) * (n_at / n_any_gc)
    else as.numeric(x[["Freq"]])
  )

str_data_by_indiv <-
  merge(
    str_data_by_indiv,
    with(
      subset(str_data, 
             str_data$validation == "true_de_novo" &
               !str_data$is_homopolymer &
               str_data$poocase %in% c(2, 3)),
      as.data.frame.table(table(child, any_gc, poocase))
    )  %>%
      mutate(., any_gc_phase = paste(any_gc, poocase, sep = "_"), .keep = "unused") %>%
      pivot_wider(., names_from = any_gc_phase, values_from = Freq) %>%
      rename(., "AT_paternal" = "FALSE_2", "GC_paternal" = "TRUE_2", "AT_maternal" = "FALSE_3", "GC_maternal" = "TRUE_3"),
    by.x = "ssc_id", by.y = "child", all.x = T
  ) %>%
  mutate(., AT_paternal = replace_na(AT_paternal, 0),
         GC_paternal = replace_na(GC_paternal, 0),
         AT_maternal = replace_na(AT_maternal, 0),
         GC_maternal = replace_na(GC_maternal, 0))

summary(
  glm(
    n_muts ~ motherAgeAtBirth + gc_content + offset(log(denom)),
    data = 
      cbind(
        pivot_longer(str_data_by_indiv[, c("AT_maternal", "GC_maternal", "motherAgeAtBirth")], 
                     cols = c("AT_maternal", "GC_maternal"), 
                     names_to = "gc_content", 
                     values_to = "n_muts"),
        denom = rep(c(n_at, n_any_gc), times = nrow(str_data_by_indiv))
      ),
    family = poisson()
  )
)
summary(
  glm(
    n_muts ~ motherAgeAtBirth * gc_content + offset(log(denom)),
    data = 
      cbind(
        pivot_longer(str_data_by_indiv[, c("AT_maternal", "GC_maternal", "motherAgeAtBirth")], 
                     cols = c("AT_maternal", "GC_maternal"), 
                     names_to = "gc_content", 
                     values_to = "n_muts"),
        denom = rep(c(n_at, n_any_gc), times = nrow(str_data_by_indiv))
      ),
    family = poisson()
  )
)
summary(
  glm(
    n_muts ~ fatherAgeAtBirth * gc_content + offset(log(denom)),
    data = 
      cbind(
        pivot_longer(str_data_by_indiv[, c("AT_paternal", "GC_paternal", "fatherAgeAtBirth")], 
                     cols = c("AT_paternal", "GC_paternal"), 
                     names_to = "gc_content", 
                     values_to = "n_muts"),
        denom = rep(c(n_at, n_any_gc), times = nrow(str_data_by_indiv))
      ),
    family = poisson
  )
)
summary(
  glm(
    n_muts ~ fatherAgeAtBirth + gc_content + offset(log(denom)),
    data = 
      cbind(
        pivot_longer(str_data_by_indiv[, c("AT_paternal", "GC_paternal", "fatherAgeAtBirth")], 
                     cols = c("AT_paternal", "GC_paternal"), 
                     names_to = "gc_content", 
                     values_to = "n_muts"),
        denom = rep(c(n_at, n_any_gc), times = nrow(str_data_by_indiv))
      ),
    family = poisson
  )
)


ggsave(
  "../figures/nonhomopolymer_gc_content_reweight_by_maternal_age_read1_flank2_20230607.pdf",
  ggplot(
    data = mat_str_dnm_by_nuc_content,
    aes(x = motherAgeAtBirth, y = freq_reweight)
  ) + 
    geom_jitter(aes(color = any_gc), height = 0.2, alpha = 0.1) +
    stat_smooth(aes(color = any_gc), formula = y ~ x,
                method = "glm") +
    theme_cowplot() +
    scale_y_continuous(limits = c(0, 5), oob = scales::squish) +
    scale_color_manual(values = wes_palette("Royal1")),
  width = 7, height = 5
)

# Postzygotic analysis
# combi <- function(vec1)
# {
#   si <- length(vec1)
#   first <- rep(vec1, (si-1):0)
#   secR <- rev(vec1)
#   second <- secR[sequence(1:(si-1))]
#   second <- rev(second)
#   combi <- matrix(cbind(first, second), ncol = 2)
#   return(combi)
# }

paternal_diff_conditional_paternal_age <-
  bind_rows(
    lapply(
      unique(round(str_data_by_indiv$fatherAgeAtBirth)),
      function(pat_age){
        curr_rows = which(round(str_data_by_indiv$fatherAgeAtBirth) == pat_age)
        if(length(curr_rows) < 2) {return()}
        set.seed(0)
        curr_table =
          apply(
            matrix(
              c(
                sample(curr_rows, 
                       size = length(curr_rows)  - 
                         length(curr_rows) %% 2, replace = F)
              ), ncol = 2
            ),
            1,
            function(x){
              older = which(str_data_by_indiv$motherAgeAtBirth[unlist(x)] == 
                              max(str_data_by_indiv$motherAgeAtBirth[unlist(x)]))
              mat_age_diff = 
                str_data_by_indiv$motherAgeAtBirth[x[older]] - 
                str_data_by_indiv$motherAgeAtBirth[x[(older %% 2) + 1]]
              pat_mut_diff = 
                str_data_by_indiv$paternal_nonhomopolymer[x[older]] - 
                str_data_by_indiv$paternal_nonhomopolymer[x[(older %% 2) + 1]]
              mat_mut_diff = 
                str_data_by_indiv$maternal_nonhomopolymer[x[older]] - 
                str_data_by_indiv$maternal_nonhomopolymer[x[(older %% 2) + 1]]
              return(c(mat_age_diff, pat_mut_diff, mat_mut_diff))
            }
          )
        return(
          cbind.data.frame(
            pat_age,
            matrix(unlist(curr_table), byrow = T, ncol = 3)
          )
        )
      }
    )
  )
colnames(paternal_diff_conditional_paternal_age)[c(2:4)] <- c("mat_age_diff", "pat_mut_diff", "mat_mut_diff")

summary(
  glm(
    pat_mut_diff ~ mat_age_diff, 
    data = paternal_diff_conditional_paternal_age
  )
)
ggsave(
  "../figures/pat_mut_diff_mat_age_20230729.pdf",
  ggplot(
    data = 
      paternal_diff_conditional_paternal_age,
    aes(x = mat_age_diff, y = pat_mut_diff)
  ) +
    geom_jitter(alpha = 0.1, height = 0.5) +
    scale_x_continuous(limits = c(0, 17)) +
    theme_minimal() +
    stat_smooth(method = "glm", formula = y ~ x)
)
cor.test(
  paternal_diff_conditional_paternal_age$mat_age_diff,
  paternal_diff_conditional_paternal_age$pat_mut_diff,
  method = "spearman", alternative = "greater")
# No association

# Selection
sistr_data <-
  read.delim("./SSC_scores.txt", na.strings = c("N/A")) %>%
  mutate(., chrom = paste("chr", chrom, sep = ""))
sistr_data$"CI_lower" <- 
  sapply(sistr_data$ABC_s_95._CI, 
         function(x)
           as.numeric(strsplit(strsplit(x, split = " , ")[[1]][1], "\\(")[[1]][2])
  )
sistr_data$"CI_upper" <- 
  sapply(sistr_data$ABC_s_95._CI, 
         function(x)
           as.numeric(strsplit(strsplit(x, split = " , ")[[1]][2], "\\)")[[1]][1])
  )
str_data <-
  merge(
    str_data %>%
      mutate(., chrom_pos = paste(chrom, pos, sep = "_")),
    sistr_data[, c("chrom", "start", "ABC_s_median", "LRT_p_value", "CI_upper", "CI_lower")] %>%
      filter(., (CI_upper - CI_lower) < 0.3) %>%
      mutate(., chrom_pos = paste(chrom, start, sep = "_"), .keep = "unused"),
    by = "chrom_pos", all.x = T, all.y = F
  )
# str_data$"under_selection_called" <-
#   -log(str_data$ABC_s_median) > -log(median(str_data$ABC_s_median, na.rm = T))
# 
# str_data_by_indiv <-
#   merge(
#     str_data_by_indiv,
#     with(
#       subset(str_data, 
#              str_data$validation == "true_de_novo" &
#                !str_data$is_homopolymer &
#                str_data$poocase %in% c(2, 3)),
#       as.data.frame.table(table(child, poocase, under_selection_called))
#     ) %>% 
#       mutate(., selection_poocase = paste(under_selection_called, poocase, sep = "_"), .keep = "unused") %>%
#       pivot_wider(., names_from = selection_poocase, values_from = Freq) %>%
#       rename(., "paternal_constrained_called" = "TRUE_2", "paternal_neutral_called" = "FALSE_2",
#              "maternal_constrained_called" = "TRUE_3", "maternal_neutral_called" = "FALSE_3"),
#     by.x = "ssc_id", by.y = "child", all.x = T
#   ) %>%
#   mutate(., paternal_constrained_called = replace_na(paternal_constrained_called, 0),
#          paternal_neutral_called = replace_na(paternal_neutral_called, 0),
#          maternal_constrained_called = replace_na(maternal_constrained_called, 0),
#          maternal_neutral_called = replace_na(maternal_neutral_called, 0))
# 
# summary(
#   glm(
#     Freq ~ motherAgeAtBirth * selection,
#     data = 
#       str_data_by_indiv[, c("motherAgeAtBirth", "maternal_constrained_called", "maternal_neutral_called")] %>%
#       pivot_longer(., cols = c("maternal_neutral_called", "maternal_constrained_called"), names_to = "selection", values_to = "Freq"),
#     family = poisson(link = "identity")
#   )
# )
# summary(
#   glm(
#     Freq ~ fatherAgeAtBirth * selection,
#     data = 
#       str_data_by_indiv[, c("fatherAgeAtBirth", "paternal_constrained_called", "paternal_neutral_called")] %>%
#       pivot_longer(., cols = c("paternal_neutral_called", "paternal_constrained_called"), names_to = "selection", values_to = "Freq"),
#     family = poisson(link = "identity")
#   )
# )
# 
# ggplot(
#   data = 
#     str_data_by_indiv[, c("fatherAgeAtBirth", "paternal_constrained_called", "paternal_neutral_called")] %>%
#     pivot_longer(., cols = c("paternal_neutral_called", "paternal_constrained_called"), names_to = "selection", values_to = "Freq"),
#   aes(x = fatherAgeAtBirth, y = Freq)
# ) +
#   geom_jitter(aes(color = selection)) +
#   stat_smooth(aes(color = selection), formula = y ~ x, method = "glm")

ggsave(
  "../figures/constraint_by_maternal_age_decile_20230628.pdf",
  ggplot(
    data = 
      cbind(
        str_data,
        motherAgeDecile = 
          cut(str_data$motherAgeAtBirth, 
              breaks = quantile(str_data_by_indiv$motherAgeAtBirth, 
                                probs = seq(0, 1, by=0.1), na.rm = T))
      ) %>%
      filter(., poocase == 3 & validation == "true_de_novo" & !is_homopolymer & !is.na(motherAgeDecile)),
    aes(x = motherAgeDecile, y = -log(ABC_s_median))
  ) + 
    geom_boxplot()
)
ggsave(
  "../figures/constraint_by_paternal_age_decile_20230628.pdf",
  ggplot(
    data = 
      cbind(
        str_data,
        fatherAgeDecile = 
          cut(str_data$fatherAgeAtBirth, 
              breaks = quantile(str_data_by_indiv$fatherAgeAtBirth, 
                                probs = seq(0, 1, by=0.1), na.rm = T))
      ) %>%
      filter(., poocase == 2 & validation == "true_de_novo" & !is_homopolymer & !is.na(fatherAgeDecile)),
    aes(x = fatherAgeDecile, y = -log(ABC_s_median))
  ) + 
    geom_boxplot()
)
# maybe try adding some small amount to median selection coefficient to include the neutral ones?
summary(
  glm(
    -log(median_constraint) ~ motherAgeAtBirth,
    data = 
      str_data %>%
      filter(., poocase == 3 & validation == "true_de_novo" & !is_homopolymer & ABC_s_median > 0) %>%
      group_by(., child, motherAgeAtBirth) %>%
      summarise(., median_constraint = median(ABC_s_median))
  )
)
summary(
  glm(
    -log(median_constraint) ~ motherAgeAtBirth * phenotype,
    data = 
      str_data %>%
      filter(., poocase == 3 & validation == "true_de_novo" & !is_homopolymer & ABC_s_median > 0) %>%
      group_by(., child, motherAgeAtBirth, phenotype) %>%
      summarise(., median_constraint = median(ABC_s_median))
  )
)
summary(
  glm(
    -log(median_constraint) ~ fatherAgeAtBirth,
    data = 
      str_data %>%
      filter(., poocase == 2 & validation == "true_de_novo" & !is_homopolymer & ABC_s_median > 0) %>%
      group_by(., child, fatherAgeAtBirth) %>%
      summarise(., median_constraint = median(ABC_s_median))
  )
)
summary(
  glm(
    -log(median_constraint) ~ fatherAgeAtBirth + phenotype,
    data = 
      str_data %>%
      filter(., poocase == 2 & validation == "true_de_novo" & !is_homopolymer & ABC_s_median > 0) %>%
      group_by(., child, fatherAgeAtBirth, phenotype) %>%
      summarise(., median_constraint = median(ABC_s_median))
  )
)
summary(
  glm(
    -log(median_constraint) ~ fatherAgeAtBirth * phenotype,
    data = 
      str_data %>%
      filter(., poocase == 2 & validation == "true_de_novo" & !is_homopolymer & ABC_s_median > 0) %>%
      group_by(., child, fatherAgeAtBirth, phenotype) %>%
      summarise(., median_constraint = median(ABC_s_median))
  )
)

# can we find motifs enriched for maternal age signature?
complement <-
  function(x){
    if(x == "A"){return("T")}
    else if(x == "T"){return("A")}
    else if(x == "C"){return("G")}
    else if(x == "G"){return("C")}
    else {return(NA)}
  }

regularize_repeat_unit <-
  function(x){
    revcomp_x = paste0(sapply(rev(strsplit(x, "")[[1]]), complement), collapse = "")
    return(
      sort(
        c(sapply(
          c(1:nchar(x)),
          function(y)
            substr(paste0(c(x, x), collapse = ""), y, y + nchar(x) - 1)
        ),
        sapply(
          c(1:nchar(x)),
          function(y)
            substr(paste0(c(revcomp_x, revcomp_x), collapse = ""), y, y + nchar(x) - 1)
        ))
      )[1]
    )
  }
motif_has_string <-
  function(x, curr_string){
    revcomp_x = paste0(sapply(rev(strsplit(x, "")[[1]]), complement), collapse = "")
    return(
      max(
        unlist(
          lapply(
            strsplit(
              c(paste0(c(x, x), collapse = ""),
                paste0(c(revcomp_x, revcomp_x), collapse = "")
              ), split = curr_string
            ),
            length
          )
        )
      ) > 1
    )
  }
motif_h_dna <-
  function(x){
    return(min(
      nchar(c(gsub("[AG]+", "", x), 
              gsub("[CT]+", "", x)))
    ) == 0)
  }
motif_h_dna <-
  function(x){
    return(min(
      nchar(c(gsub("[AG]+", "", x), 
              gsub("[CT]+", "", x)))
    ) == 0)
  }

str_data$"regularized_motif" <-
  sapply(
    str_data$repeat_unit,
    regularize_repeat_unit
  )
str_data$"has_CAT" <-
  sapply(
    str_data$regularized_motif,
    function(x)
      motif_has_string(x = x, curr_string = "CAT")
  )
str_data$"h_dna" <-
  sapply(
    str_data$regularized_motif,
    motif_h_dna
  )
str_data$"z_dna" <-
  str_data$regularized_motif == "AC"

common_motif_df <-
  filter(
    as.data.frame.table(
      table(
        subset(str_data$regularized_motif, 
               str_data$validation == "true_de_novo" &
                 !str_data$is_homopolymer))), 
    Freq > 20) %>%
  rename(., "motif" = "Var1")

starting_vector <-
  rep(c(0.5, 0.51), times = nrow(common_motif_df) / 2)

# i need an offset
str_panel <- read.delim("./sites_annotated_reptiming_dnm_counts_ru.txt", header = F, na.strings = "NA")
names(str_panel) <-
  c("chrom", "pos", "repeat_unit", "len_ref_ru", "reptiming", "n_dnms")
str_panel$regularized_motif <-
  sapply(
    str_panel$repeat_unit,
    regularize_repeat_unit
  )
str_panel$"has_CAT" <-
  sapply(
    str_panel$regularized_motif,
    function(x)
      motif_has_string(x = x, curr_string = "CAT")
  )
str_panel$"h_dna" <-
  sapply(
    str_panel$regularized_motif,
    motif_h_dna
  )
str_panel$"z_dna" <-
  str_panel$regularized_motif == "AC"

maternal_age_interaction_by_motif_vector <-
  function(motif_vector) {
    in_denom = 
      str_panel %>%
      mutate(., motif_binary = 
               regularized_motif %in% common_motif_df$motif[
                 motif_vector > 0.5
               ]) %>%
      filter(., nchar(regularized_motif) > 1) %>%
      count(., motif_binary)
    return(
      (
        summary(
          glm(
            Freq ~ motherAgeAtBirth * motif_binary + offset(log(denom)),
            family = poisson(),
            data = 
              merge(
                str_data_by_indiv[, c("ssc_id", "motherAgeAtBirth")],
                str_data %>%
                  mutate(., motif_binary = 
                           regularized_motif %in% common_motif_df$motif[
                             motif_vector > 0.5
                           ]) %>%
                  filter(., poocase == 3 & !is_homopolymer & validation == "true_de_novo") %>%
                  count(., child, motif_binary) %>%
                  pivot_wider(., names_from = motif_binary, values_from = n) %>%
                  rename(., "in_motif" = "TRUE", "out_motif" = "FALSE") %>%
                  mutate(., in_motif = replace_na(in_motif, 0),
                         out_motif = replace_na(out_motif, 0)),
                by.x = "ssc_id", by.y = "child", all.x = T
              ) %>%
              pivot_longer(., cols = c("in_motif", "out_motif"), names_to = "motif_binary", values_to = "Freq") %>%
              mutate(., Freq = replace_na(Freq, 0)) %>%
              mutate(., denom = in_denom$n[as.numeric(motif_binary == "in_motif") + 1])
          )
        )$coefficients[4, 1])
    )
  }

run_motif_model <-
  optim(
    starting_vector,
    maternal_age_interaction_by_motif_vector,
    method = "L-BFGS-B",
    lower=rep(0, times = nrow(common_motif_df)),
    upper=rep(1, times = nrow(common_motif_df)),
    control = list(trace=4, maxit=100))

run_motif_model_10x <-
  lapply(
    c(1:10),
    function(x)
      optim(
        sample(starting_vector, size = length(starting_vector), replace = F),
        maternal_age_interaction_by_motif_vector,
        method = "L-BFGS-B",
        lower=rep(0, times = nrow(common_motif_df)),
        upper=rep(1, times = nrow(common_motif_df)),
        control = list(trace=4, maxit=100))
  )


in_denom_final = 
  str_panel %>%
  mutate(., motif_binary = 
           regularized_motif %in% common_motif_df$motif[
             run_motif_model$par > 0.5
           ]) %>%
  filter(., nchar(regularized_motif) > 1) %>%
  count(., motif_binary)

random_grid_search_output <-
  sapply(
    c(1:1000),
    function(x){
      v = sample(starting_vector, size = length(starting_vector), replace = F)
      return(
        list(
          output_val = maternal_age_interaction_by_motif_vector(v),
          curr_params = v
        )
      )
    }

  )

common_motif_df$"motif_slope_coef" <-
  sapply(
    common_motif_df$motif,
    function(x)
      summary(
        glm(
          n ~ motherAgeAtBirth,
          family = poisson(),
          data = 
            merge(
              str_data_by_indiv[, c("ssc_id", "motherAgeAtBirth")],
              str_data %>%
                filter(., 
                       regularized_motif == x & 
                         poocase == 3 & 
                         !is_homopolymer & 
                         validation == "true_de_novo") %>%
                count(., child),
              by.x = "ssc_id", by.y = "child", all.x = T
            ) %>%
            mutate(., n = replace_na(n, 0))
        )
      )$coefficients[2,1]
  )

common_motif_df$"motif_slope_pval" <-
  sapply(
    common_motif_df$motif,
    function(x)
      summary(
        glm(
          n ~ motherAgeAtBirth,
          family = poisson(),
          data = 
            merge(
              str_data_by_indiv[, c("ssc_id", "motherAgeAtBirth")],
              str_data %>%
                filter(., 
                       regularized_motif == x & 
                         poocase == 3 & 
                         !is_homopolymer & 
                         validation == "true_de_novo") %>%
                count(., child),
              by.x = "ssc_id", by.y = "child", all.x = T
            ) %>%
            mutate(., n = replace_na(n, 0))
        )
      )$coefficients[2,4]
  )

in_denom_h_dna = 
  str_panel %>%
  filter(., nchar(regularized_motif) > 1) %>%
  count(., h_dna)
summary(
  glm(
    Freq ~ motherAgeAtBirth * motif_binary + offset(log(denom)),
    family = poisson(),
    data = 
      merge(
        str_data_by_indiv[, c("ssc_id", "motherAgeAtBirth")],
        str_data %>%
          filter(., poocase == 3 & !is_homopolymer & validation == "true_de_novo") %>%
          count(., child, h_dna) %>%
          pivot_wider(., names_from = h_dna, values_from = n) %>%
          rename(., "in_motif" = "TRUE", "out_motif" = "FALSE") %>%
          mutate(., in_motif = replace_na(in_motif, 0),
                 out_motif = replace_na(out_motif, 0)),
        by.x = "ssc_id", by.y = "child", all.x = T
      ) %>%
      pivot_longer(., cols = c("in_motif", "out_motif"), names_to = "motif_binary", values_to = "Freq") %>%
      mutate(., Freq = replace_na(Freq, 0)) %>%
      mutate(., denom = in_denom_h_dna$n[as.numeric(motif_binary == "in_motif") + 1])
  )
)
summary(
  glm(
    maternal_nonhomopolymer ~ motherAgeAtBirth, data = str_data_by_indiv, family = poisson()
    )
  )
# h DNA model is a terrible fit

in_denom_z_dna = 
  str_panel %>%
  filter(., nchar(regularized_motif) > 1) %>%
  count(., z_dna)
summary(
  glm(
    Freq ~ motherAgeAtBirth * motif_binary + offset(log(denom)),
    family = poisson(),
    data = 
      merge(
        str_data_by_indiv[, c("ssc_id", "motherAgeAtBirth")],
        str_data %>%
          filter(., poocase == 3 & !is_homopolymer & validation == "true_de_novo") %>%
          count(., child, z_dna) %>%
          pivot_wider(., names_from = z_dna, values_from = n) %>%
          rename(., "in_motif" = "TRUE", "out_motif" = "FALSE") %>%
          mutate(., in_motif = replace_na(in_motif, 0),
                 out_motif = replace_na(out_motif, 0)),
        by.x = "ssc_id", by.y = "child", all.x = T
      ) %>%
      pivot_longer(., cols = c("in_motif", "out_motif"), names_to = "motif_binary", values_to = "Freq") %>%
      mutate(., Freq = replace_na(Freq, 0)) %>%
      mutate(., denom = in_denom_CAT$n[as.numeric(motif_binary == "in_motif") + 1])
  )
)


summary(
  glm(
    n ~ motherAgeAtBirth,
    family = poisson(),
    data = 
      merge(
        str_data_by_indiv[, c("ssc_id", "motherAgeAtBirth")],
        str_data %>%
          filter(., poocase == 3 & !is_homopolymer & validation == "true_de_novo" & !has_CAT) %>%
          count(., child),
        by.x = "ssc_id", by.y = "child", all.x = T
      ) %>%
      mutate(., n = replace_na(n, 0))
  )
)

summary(
  glm(
    Freq ~ motherAgeAtBirth * motif_binary + offset(log(denom)),
    family = poisson(),
    data = 
      merge(
        str_data_by_indiv[, c("ssc_id", "motherAgeAtBirth")],
        str_data %>%
          mutate(., motif_binary = 
                   regularized_motif %in% common_motif_df$motif[
                     run_motif_model$par > 0.5
                   ]) %>%
          filter(., poocase == 3 & !is_homopolymer & validation == "true_de_novo") %>%
          count(., child, motif_binary) %>%
          pivot_wider(., names_from = motif_binary, values_from = n) %>%
          rename(., "in_motif" = "TRUE", "out_motif" = "FALSE") %>%
          mutate(., in_motif = replace_na(in_motif, 0),
                 out_motif = replace_na(out_motif, 0)),
        by.x = "ssc_id", by.y = "child", all.x = T
      ) %>%
      pivot_longer(., cols = c("in_motif", "out_motif"), names_to = "motif_binary", values_to = "Freq") %>%
      mutate(., Freq = replace_na(Freq, 0)) %>%
      mutate(., denom = in_denom_final$n[as.numeric(motif_binary == "in_motif") + 1])
  )
)


ggplot(
  data = 
    merge(
      str_data_by_indiv[, c("ssc_id", "motherAgeAtBirth")],
      str_data %>%
        mutate(., motif_binary = 
                 regularized_motif %in% common_motif_df$motif[
                   run_motif_model$par > 0.5
                 ]) %>%
        filter(., poocase == 3 & !is_homopolymer & validation == "true_de_novo") %>%
        count(., child, motif_binary) %>%
        pivot_wider(., names_from = motif_binary, values_from = n) %>%
        rename(., "in_motif" = "TRUE", "out_motif" = "FALSE") %>%
        mutate(., in_motif = replace_na(in_motif, 0),
               out_motif = replace_na(out_motif, 0)),
      by.x = "ssc_id", by.y = "child", all.x = T
    ) %>%
    pivot_longer(., cols = c("in_motif", "out_motif"), names_to = "motif_binary", values_to = "Freq") %>%
    mutate(., Freq = replace_na(Freq, 0)),
  aes(x = motherAgeAtBirth, y = Freq)
) + geom_point(aes(color = motif_binary)) +
  stat_smooth(aes(color = motif_binary), formula = y ~ x, method = "glm")


    
  

