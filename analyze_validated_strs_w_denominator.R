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

ses_data <- read.delim("./SSC_SES_forMG.csv", sep = ",")
ses_data$"family" <- sub("\\..*", "", ses_data$individual)

str_data <- read.delim("./dnms_w_filter_info.csv", sep = ",")
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
str_data$"reptiming_early_late" <-
  sapply(
    str_data$reptiming_quartile,
    function(x) if(is.na(x)) NA else if(x > 2) "early" else "late"
  )
str_data <- merge(str_data, ancestry_data[, c("X2", "SuperPop", "AFR")], by.x = "child", by.y = "X2", all.x = T)
str_data$"len_repeat_unit_condensed" <- 
  sapply(str_data$len_repeat_unit,
         function(x) if(x >= 4) ">=4" else x)
str_data$len_repeat_unit_condensed <-
  factor(str_data$len_repeat_unit_condensed, levels = c("1", "2", "3", ">=4"))
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
str_data$"accurate_mutsize" <-
  str_data$mutate_from_homo_or_het == 1 | str_data$mat_gt == str_data$pat_gt
str_data$"accurate_expansion_or_deletion" <-
  apply(
    str_data,
    1,
    function(x)
      if(x[["accurate_mutsize"]]) return(TRUE)
    else{
      possible_alleles=c()
      if(x[["poocase"]] == 4) {
        possible_alleles = 
          as.integer(unlist(strsplit(c(x[["mat_gt"]], x[["pat_gt"]]), ",")))
      }
      else if(x[["poocase"]] == 2){
        possible_alleles = 
          as.integer(unlist(strsplit(x[["pat_gt"]], ",")))
      }
      else {
        possible_alleles = 
          as.integer(unlist(strsplit(x[["mat_gt"]], ",")))
      }
      return(x[["newallele"]] < min(possible_alleles) | 
               x[["newallele"]] > max(possible_alleles))
    }
  )
str_data$"quality_bin" <-
  cut(str_data$posterior, breaks = quantile(str_data$posterior, probs = seq(0, 1, by=0.5)))
str_data$"parents_share_allele" <-
  apply(
    str_data,
    1,
    function(x)
      any(unlist(strsplit(x[["mat_gt"]], ",")) %in% 
            unlist(strsplit(x[["pat_gt"]], ",")))
  )
str_data$"any_gc" <-
  sapply(
    str_data$repeat_unit,
    function(x)
      nchar(sub("[CG]", "", x)) < nchar(x)
  )
str_data$"abs_mut_size" <-
  abs(str_data$mutsize)

# table(
#   str_data$len_repeat_unit_condensed,
#   str_data$validation
# )


str_data_by_indiv <- 
  as.data.frame.table(
    table(
      subset(str_data[, c("child", "poocase")],
             str_data$validation == "true_de_novo")
    )
  ) %>% spread(poocase, Freq)
colnames(str_data_by_indiv) <- c("ssc_id", "paternal", "maternal", "unphased")
str_data_by_indiv$"total" <- rowSums(str_data_by_indiv[, c(2:4)])

str_data_by_indiv <-
  merge(str_data_by_indiv, 
        ped_file[, c("Sample", "father.id", "mother.id", "ssc_id", "sex")], 
        by="ssc_id")
# wilcox.test(subset(str_data_by_indiv$motherAgeAtBirth, str_data_by_indiv$sex == "M"), 
#             subset(str_data_by_indiv$motherAgeAtBirth, str_data_by_indiv$sex == "F"))
# no detected difference in parental age at birth between sexes of kids, not going to worry about excluding chrX from analyses
str_data_by_indiv <-
  merge(
    str_data_by_indiv,
    parental_ages,
    by.x="Sample", by.y="sample"
  )
str_data_by_indiv <-
  merge(
    str_data_by_indiv,
    as.data.frame.table(
      table(
        str_data[, c("child", "expansion_or_deletion")]
      )
    ) %>% spread(expansion_or_deletion, Freq),
    by.x="ssc_id", by.y="child"
  )

# phased deletion or expansions
str_data_by_indiv <-
  merge(
    str_data_by_indiv,
    as.data.frame.table(
      table(
        subset(str_data[, c("child", "expansion_or_deletion", "poocase")], str_data$poocase %in% c(2, 3) & str_data$validation == "true_de_novo")
      )
    ) %>% mutate(expansion_or_deletion_poocase=paste(expansion_or_deletion, poocase, sep = "_"), expansion_or_deletion=NULL, poocase=NULL) %>% 
      spread(expansion_or_deletion_poocase, Freq),
    by.x = "ssc_id", by.y = "child", all.x = T
  )
str_data_by_indiv$"motherAge_quantile" <-
  cut(str_data_by_indiv$motherAgeAtBirth, 
      breaks = quantile(str_data_by_indiv$motherAgeAtBirth, probs = seq(0, 1, by=0.2)))
str_data_by_indiv$"fatherAge_quantile" <-
  cut(str_data_by_indiv$fatherAgeAtBirth, 
      breaks = quantile(str_data_by_indiv$fatherAgeAtBirth, probs = seq(0, 1, by=0.2)))
str_data_by_indiv$"family" <-
  sub("\\..*", "", str_data_by_indiv$Sample)
str_data_by_indiv$"SuperPop" <- 
  sapply(str_data_by_indiv$Sample,
         function(x) ancestry_data$SuperPop[ancestry_data$eichler.id == x][1])
str_data_by_indiv$"SuperPop.mo" <- 
  sapply(str_data_by_indiv$family,
         function(x) ancestry_data$SuperPop[ancestry_data$family == x & ancestry_data$indiv == "mo"][1])
str_data_by_indiv$"SuperPop.fa" <- 
  sapply(str_data_by_indiv$family,
         function(x) ancestry_data$SuperPop[ancestry_data$family == x & ancestry_data$indiv == "fa"][1])
str_data_by_indiv$"AFR_fraction" <- 
  sapply(str_data_by_indiv$Sample,
         function(x) ancestry_data$AFR[ancestry_data$eichler.id == x][1])
str_data <-
  merge(str_data, str_data_by_indiv[, c("ssc_id", "SuperPop.mo", "SuperPop.fa")],
        by.x = "child", by.y = "ssc_id", all.x = T, all.y = F)
str_data_by_indiv <-
  merge(str_data_by_indiv, ses_data[, c("family", "bckgd_hx_annual_household")], by = "family", all.x = T)
str_data_by_indiv <- 
  merge(
    str_data_by_indiv,
    as.data.frame.table(
      table(
        subset(str_data[, c("child", "poocase")],
               str_data$validation == "true_de_novo" & 
                 str_data$len_repeat_unit > 1 & 
                 str_data$all_dnms_discoverable == 1)
      )
    ) %>% spread(poocase, Freq),
    by.x = "ssc_id", by.y = "child"
  ) %>% 
  rename(., c("paternal_nonhomopolymer_discoverable" = "2",
                        "maternal_nonhomopolymer_discoverable" = "3",
                        "unphased_nonhomopolymer_discoverable" = "4"))


as.data.frame.table(
  table(
    subset(str_data[, c("child", "poocase", "abs_mut_size")],
           str_data$validation == "true_de_novo" & 
             str_data$len_repeat_unit > 1 & 
             str_data$abs_mut_size < 4 &
             str_data$all_dnms_discoverable == 1)
  )
)


# denominators
str_denominators <- read.delim("./parents_callable_dnms_matrix.txt", header = F)
colnames(str_denominators) <-
  c("parent1", "parent2", "ru_len", "n_hets", "contains_gc", "dnm_size", "n")
str_denominators <- 
  merge(
    str_denominators,
    ped_file[, c("family", "ssc_id")],
    by.x = "parent1", by.y = "ssc_id"
  )

discoverable_str_data_by_indiv <-
  merge(
    str_denominators %>%
      filter(., ru_len > 1) %>%
      group_by(., family, dnm_size) %>%
      summarize(., n = sum(n)),
    merge(
      as.data.frame.table(
        table(
          subset(str_data[, c("child", "poocase", "abs_mut_size")],
                 str_data$validation == "true_de_novo" & 
                   str_data$len_repeat_unit > 1 & 
                   str_data$abs_mut_size < 4 &
                   str_data$all_dnms_discoverable == 1)
        )
      ),
      unique(str_data[, c("child", "family", "SuperPop")]), 
      by = "child", all.x = T, all.y = F
    ),
    by.x = c("family", "dnm_size"), by.y = c("family", "abs_mut_size")
  )

discoverable_mut_rate_by_indiv <-
  merge(
    discoverable_str_data_by_indiv %>%
      group_by(., child, family, dnm_size, SuperPop, poocase) %>%
      summarize(., dnm_rate = Freq / (2 * n)) %>%
      group_by(., child, family, SuperPop, poocase) %>%
      summarize(., mut_rate = sum(dnm_rate)),
    str_data_by_indiv[, c("ssc_id", "motherAgeAtBirth", "fatherAgeAtBirth")],
    by.x = "child", by.y = "ssc_id", all.y = F
  )

ggplot(
  data = 
    merge(
      discoverable_str_data_by_indiv %>%
        group_by(., child, family, dnm_size, SuperPop, poocase) %>%
        summarize(., dnm_rate = Freq / (2 * n)) %>%
        group_by(., child, family, SuperPop, poocase) %>%
        summarize(., mut_rate = sum(dnm_rate)),
      str_data_by_indiv[, c("ssc_id", "motherAgeAtBirth", "fatherAgeAtBirth")],
      by.x = "child", by.y = "ssc_id", all.y = F
    ) %>% 
    filter(., poocase == 3),
  aes(x = motherAgeAtBirth, y = mut_rate)
) + 
  geom_point() + 
  stat_smooth(formula = y ~ x, method = "glm")

summary(
  glm(
    mut_rate ~ motherAgeAtBirth,
    data = merge(
      discoverable_str_data_by_indiv %>%
        group_by(., child, family, dnm_size, SuperPop, poocase) %>%
        summarize(., dnm_rate = Freq / (2 * n)) %>%
        group_by(., child, family, SuperPop, poocase) %>%
        summarize(., mut_rate = sum(dnm_rate)),
      str_data_by_indiv[, c("ssc_id", "motherAgeAtBirth", "fatherAgeAtBirth")],
      by.x = "child", by.y = "ssc_id", all.y = F
    ) %>% 
      filter(., poocase == 3)
  )
)

summary(
  glm(
    mut_rate_across_phase ~ motherAgeAtBirth + fatherAgeAtBirth,
    data = discoverable_mut_rate_by_indiv %>%
        group_by(., child, family, SuperPop, motherAgeAtBirth, fatherAgeAtBirth) %>%
        summarize(., mut_rate_across_phase = sum(mut_rate)),
  )
)

summary(
  glm(
    mut_rate_across_phase ~  fatherAgeAtBirth,
    data = discoverable_mut_rate_by_indiv %>%
      group_by(., child, family, SuperPop, motherAgeAtBirth, fatherAgeAtBirth) %>%
      summarize(., mut_rate_across_phase = sum(mut_rate)),
  )
)

ggplot(
  data = discoverable_mut_rate_by_indiv %>%
    group_by(., child, family, SuperPop, motherAgeAtBirth, fatherAgeAtBirth) %>%
    summarize(., mut_rate_across_phase = sum(mut_rate)),
  aes(x = SuperPop, y = mut_rate_across_phase)
) + geom_boxplot()

summary(
  glm(
    mut_rate_across_phase ~ motherAgeAtBirth + fatherAgeAtBirth + SuperPop,
    data = discoverable_mut_rate_by_indiv %>%
      group_by(., child, family, SuperPop, motherAgeAtBirth, fatherAgeAtBirth) %>%
      summarize(., mut_rate_across_phase = sum(mut_rate)),
  )
)


summary(
  glm(
    mut_rate ~ fatherAgeAtBirth,
    data = merge(
      discoverable_str_data_by_indiv %>%
        group_by(., child, family, dnm_size, SuperPop, poocase) %>%
        summarize(., dnm_rate = Freq / (2 * n)) %>%
        group_by(., child, family, SuperPop, poocase) %>%
        summarize(., mut_rate = sum(dnm_rate)),
      str_data_by_indiv[, c("ssc_id", "motherAgeAtBirth", "fatherAgeAtBirth")],
      by.x = "child", by.y = "ssc_id", all.y = F
    ) %>% 
      filter(., poocase == 2)
  )
)






str_data$"abs_mutsize" <- abs(str_data$mutsize)
table(subset(str_data[, c("len_repeat_unit_condensed", "abs_mutsize")], str_data$child == "SSC00003" & 
         str_data$validation == "true_de_novo"))

str_denominators %>% 
  filter(., family == "11006") %>%
  group_by(ru_len, dnm_size) %>%
  summarize(n = sum(n))


# Reading in heterozygosity
couples_gts <-
  read.delim("./couples_gts.txt")
couples_gts$"total_0_hets" = rowSums(couples_gts[, c("X1_0_het", "X2_0_het", "X3_0_het", "X4_0_het")])
couples_gts$"total_1_hets" = rowSums(couples_gts[, c("X1_1_het", "X2_1_het", "X3_1_het", "X4_1_het")])
couples_gts$"total_2_hets" = rowSums(couples_gts[, c("X1_2_het", "X2_2_het", "X3_2_het", "X4_2_het")])
couples_gts <- 
  merge(couples_gts, ped_file[, c("family", "ssc_id")], by.x = "sample1", by.y = "ssc_id", all.x = T)

# chrX mutations do not have a phase assigned; this filtering is just a nice way to make sure we have consistency
str_data_by_indiv <-
  merge(
    str_data_by_indiv,
    unique(str_data[, c("child", "phase")]) %>%
      filter(., phase != ""),
    by.x = "ssc_id", by.y = "child", all.x = T
  )
# Currently (20220509) the best filtering strategy I have for dealing with dropout is insisting on mutsize == 1 and a stricter cutoff for quality
str_data_by_indiv <- 
  merge(
    str_data_by_indiv,
    as.data.frame.table(
      table(
        subset(str_data$child,
               str_data$mutsize == 1 &
                 str_data$posterior > 0.921 &
                 str_data$validation == "true_de_novo")
      )
    ) %>% rename(., c("ssc_id" = "Var1", "total_mutsize1_high_quality" = "Freq")),
    by = "ssc_id", all.x = T
  )
# It doesn't seem like the enclosing reads filter really helps me that much
# str_data_by_indiv <- 
#   merge(
#     str_data_by_indiv,
#     as.data.frame.table(
#       table(
#         subset(str_data$child,
#                str_data$mutsize == 1 &
#                  str_data$posterior > 0.921 &
#                  str_data$encl_child >= 10)
#       )
#     ) %>% rename(., c("ssc_id" = "Var1", "total_mutsize1_high_quality_encl10" = "Freq")),
#     by = "ssc_id"
#   )

str_data_by_indiv <- 
  merge(
    str_data_by_indiv,
    as.data.frame.table(
      table(
        subset(str_data[, c("child", "poocase")],
               str_data$mutsize == 1 &
                 str_data$posterior > 0.921 &
                 str_data$poocase != 4 &
                 str_data$validation == "true_de_novo")
      )
    ) %>% spread(poocase, Freq) %>% 
      rename(., c("ssc_id" = "child", "paternal_high_quality" = "2", "maternal_high_quality" = "3")),
    by = "ssc_id", all.x = T
  )

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

validation_data <- 
  read.delim("./mitra_et_al_2021_validations.txt")
validation_data <- 
  rename(validation_data, c("validation_family_number" = "Family",
                                      "chrom" = "Chrom",
                                      "pos" = "TR.start.position..hg38."))
validation_data$"family" <-
  sapply(validation_data$validation_family_number,
         function(x)
           if (x == 1) 14686
         else if (x == 2) 11963
         else if (x == 3) 12605
         else if (x == 4) 14329
         else if (x == 5) 11450)

str_data$"is_in_mitra_validation" <-
  sapply(
    c(1:nrow(str_data)),
    function(x)
      if(!str_data$family[x] %in% validation_data$family) FALSE
    else
      nrow(plyr::match_df(validation_data[, c("chrom", "pos", "family")], 
                    str_data[x, c("chrom", "pos", "family")])) > 0
  )
str_data$"pass_mitra_validation" <-
  sapply(
    c(1:nrow(str_data)),
    function(x)
      if(!str_data$family[x] %in% validation_data$family) NA
    else
      if(!str_data$pos[x] %in% subset(validation_data$pos, validation_data$family == str_data$family[x]))
        NA
    else
      str_data$pos[x] %in% subset(validation_data$pos, validation_data$Mutation.validated == "Y")
  )

ggsave(
  "../figures/mitra_vs_read1_flank2_validations_20220816.pdf",
  ggplot(
    data = subset(str_data, str_data$is_in_mitra_validation),
    aes(x = pass_mitra_validation)
  ) +
    geom_bar(aes(fill = validation), position = "fill")
)

summary(
  glm(
    (deletion_3_nonhomopolymer + expansion_3_nonhomopolymer) ~ motherAgeAtBirth, 
    data = str_data_by_indiv, 
    family = poisson(link = "identity")
  )
)
summary(
  glm(
    paternal ~ log(fatherAgeAtBirth), 
    data = str_data_by_indiv, 
    family = poisson(link = "identity")
  )
)


str_data_by_indiv <-
  merge(
    str_data_by_indiv,
    couples_gts[, c("family", "total_parental_heterozygosity")],
    by = "family"
  )

wilcox.test(
  subset(str_data_by_indiv$total_mutsize1_high_quality, str_data_by_indiv$SuperPop == "AFR"),
  subset(str_data_by_indiv$total_mutsize1_high_quality, str_data_by_indiv$SuperPop != "AFR")
)

ggplot(
  data = str_data_by_indiv,
  aes(x = total_parental_heterozygosity,
      y = total)
) + geom_point(alpha = 0.2)

summary(glm(
  total ~ SuperPop + total_parental_heterozygosity + motherAgeAtBirth + fatherAgeAtBirth,
  data = str_data_by_indiv,
  family = poisson()
)
)

ggsave(
  "../figures/mat_str_vs_mat_age_20220912.pdf",
  ggplot(
    data = 
      str_data_by_indiv,
    aes(x = motherAgeAtBirth, y = maternal)
  ) +
    geom_point(alpha = 0.5) +
    stat_smooth(method = "glm", formula = y ~ x, method.args = list(family = poisson(link = "identity"))),
  width = 7, height = 5
)
ggsave(
  "../figures/pat_str_vs_pat_age_20220912.pdf",
  ggplot(
    data = 
      str_data_by_indiv,
    aes(x = fatherAgeAtBirth, y = paternal)
  ) +
    geom_point(alpha = 0.5) +
    stat_smooth(method = "glm", formula = y ~ x, method.args = list(family = poisson(link = "identity"))),
  width = 7, height = 5
)
ggsave(
  "../figures/mat_str_vs_mat_age_20220912.pdf",
  ggplot(
    data = 
      str_data_by_indiv,
    aes(x = motherAgeAtBirth, y = maternal)
  ) +
    geom_point(alpha = 0.5) +
    stat_smooth(method = "glm", formula = y ~ x, method.args = list(family = poisson(link = "identity"))),
  width = 7, height = 5
)

ggsave(
  "../figures/mat_str_vs_mat_age_homopolymer_vs_nonhomopolymer_20220913.pdf",
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
        ),
        str_data_by_indiv[, c("ssc_id", "motherAgeAtBirth")],
        by.x = "child", by.y = "ssc_id", all.x = T), 
    aes(x = motherAgeAtBirth, y = Freq)
  ) +
    geom_point(aes(color = is_homopolymer), alpha = 0.5) +
    stat_smooth(aes(color = is_homopolymer), method = "glm", formula = y ~ x, method.args = list(family = poisson(link = "identity"))),
  width = 7, height = 5
)
ggsave(
  "../figures/pat_str_vs_pat_age_homopolymer_vs_nonhomopolymer_20220913.pdf",
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
        ),
        str_data_by_indiv[, c("ssc_id", "fatherAgeAtBirth")],
        by.x = "child", by.y = "ssc_id", all.x = T), 
    aes(x = fatherAgeAtBirth, y = Freq)
  ) +
    geom_point(aes(color = is_homopolymer), alpha = 0.5) +
    stat_smooth(aes(color = is_homopolymer), method = "glm", formula = y ~ x, method.args = list(family = poisson(link = "identity"))),
  width = 7, height = 5
)
ggsave(
  "../figures/mat_str_vs_mat_age_nonhomopolymer_20220916.pdf",
  ggplot(
    data = 
      str_data_by_indiv,
    aes(x = motherAgeAtBirth, y = expansion_3_nonhomopolymer + deletion_3_nonhomopolymer)
  ) +
    geom_point(alpha = 0.5) +
    theme_cowplot() +
    stat_smooth(method = "glm", formula = y ~ x, method.args = list(family = poisson(link = "identity"))),
  width = 7, height = 5
)
ggsave(
  "../figures/pat_str_vs_pat_age_nonhomopolymer_20220916.pdf",
  ggplot(
    data = 
      str_data_by_indiv,
    aes(x = fatherAgeAtBirth, y = expansion_2_nonhomopolymer + deletion_2_nonhomopolymer)
  ) +
    geom_point(alpha = 0.5) +
    theme_cowplot() +
    stat_smooth(method = "glm", formula = y ~ x, method.args = list(family = poisson(link = "identity"))),
  width = 7, height = 5
)

poisson_model_eqn <- function(m){
  eq <- substitute(lambda == a + b %.% italic(x)*","~~italic(P)~"="~p_val, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        p_val = format(summary(m)$coefficients[2,4], digits = 3)))
  as.character(as.expression(eq));
}
maternal_deletion_model <-
  glm(deletion_3 ~ motherAgeAtBirth, data = str_data_by_indiv,
      family=poisson(link = "identity"))
maternal_expansion_model <-
  glm(expansion_3 ~ motherAgeAtBirth, data = str_data_by_indiv,
      family=poisson(link = "identity"))
paternal_deletion_model <-
  glm(deletion_2 ~ fatherAgeAtBirth, data = str_data_by_indiv,
      family=poisson(link = "identity"))
paternal_expansion_model <-
  glm(expansion_2 ~ fatherAgeAtBirth, data = str_data_by_indiv,
      family=poisson(link = "identity"))

expansion_deletion_parental_age_matrix <-
  rbind(
    rbind(
      str_data_by_indiv[, c("ssc_id", "deletion_2_nonhomopolymer", "fatherAgeAtBirth")] %>%
        rename(., c("Freq" = "deletion_2_nonhomopolymer")),
      str_data_by_indiv[, c("ssc_id", "expansion_2_nonhomopolymer", "fatherAgeAtBirth")] %>%
        rename(., c("Freq" = "expansion_2_nonhomopolymer"))
    ) %>% rename(., c("parental_age" = "fatherAgeAtBirth")),
    rbind(
      str_data_by_indiv[, c("ssc_id", "deletion_3_nonhomopolymer", "motherAgeAtBirth")] %>%
        rename(., c("Freq" = "deletion_3_nonhomopolymer")),
      str_data_by_indiv[, c("ssc_id", "expansion_3_nonhomopolymer", "motherAgeAtBirth")] %>%
        rename(., c("Freq" = "expansion_3_nonhomopolymer"))
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

ggsave(
  "../figures/str_dnms_by_ancestry_filter_2bp_flank_1_read_20220809.pdf",
  ggplot(
    data = str_data_by_indiv,
    aes(x = SuperPop, y = total)
  ) +
    geom_boxplot()
)

ggsave(
  "../figures/frac_from_homo_filter_2bp_flank_1_read_20220809.pdf",
  ggplot(
    data = 
      merge(
        str_data_by_indiv[, c("ssc_id", "SuperPop")],
        as.data.frame.table(
          table(
            subset(str_data[, c("child", "n_parents_hets", "mutate_from_homo_or_het")],
                   str_data$n_parents_hets == 1 & 
                     str_data$mutate_from_homo_or_het != 0 &
                     str_data$validation == "true_de_novo")
          )
        ) %>% spread(mutate_from_homo_or_het, Freq) %>%
          rename(., c("homo" = "1", "het" = "2")) %>%
          mutate(., frac_homozygous = homo / (homo + het), .keep = "unused"), 
        by.x = "ssc_id", by.y = "child"
      ),
    aes(x = SuperPop, y = frac_homozygous)
  ) + geom_boxplot()
)

ggplot(
  data = 
    merge(
      with(
        subset(str_data, 
               str_data$poocase %in% c(3) & 
                 str_data$expansion_or_deletion == "deletion" &
                 str_data$validation == "true_de_novo"),
        as.data.frame.table(table(child, is_homopolymer))
      ),
      str_data_by_indiv[, c("ssc_id", "motherAgeAtBirth")],
      by.x = "child", by.y = "ssc_id", all.x = T
    ),
  aes(x = motherAgeAtBirth, y = Freq)
) + 
  geom_point(aes(color = is_homopolymer), alpha = 0.2) +
  stat_smooth(aes(color = is_homopolymer), formula = y ~ x,
              method = "glm", 
              method.args = list(family = poisson(link = "identity")))

# Is this still statistically significant?
summary(
  glm(
    Freq ~ motherAgeAtBirth * is_homopolymer,
    data = 
      merge(
        with(
          subset(str_data, 
                 str_data$poocase %in% c(3) & 
                   str_data$expansion_or_deletion == "deletion" &
                   str_data$validation == "true_de_novo"),
          as.data.frame.table(table(child, is_homopolymer))
        ),
        str_data_by_indiv[, c("ssc_id", "motherAgeAtBirth")],
        by.x = "child", by.y = "ssc_id", all.x = T
      ),
    family = poisson(link = "identity")
  )
)
summary(
  glm(
    Freq ~ motherAgeAtBirth * is_homopolymer * expansion_or_deletion + fatherAgeAtBirth + SuperPop,
    data = 
      merge(
        with(
          subset(str_data, 
                 str_data$poocase %in% c(3) & 
                   abs(str_data$mutsize) == 1 
                   ),
          as.data.frame.table(table(child, is_homopolymer, expansion_or_deletion))
        ),
        str_data_by_indiv[, c("ssc_id", "motherAgeAtBirth", "fatherAgeAtBirth", "SuperPop")],
        by.x = "child", by.y = "ssc_id", all.x = T
      ),
    family = poisson(link = "identity")
  )
)

summary(
  glm(
  total_mutsize1_high_quality ~ SuperPop + total_parental_heterozygosity + motherAgeAtBirth + fatherAgeAtBirth,
  data = str_data_by_indiv,
  family = poisson()
)
)
summary(
  glm(
    total ~ SuperPop + total_parental_heterozygosity + motherAgeAtBirth + fatherAgeAtBirth,
    data = str_data_by_indiv,
    family = poisson()
  )
)


ggplot(
  data = 
    merge(
      str_data_by_indiv[, c("ssc_id", "SuperPop")],
      as.data.frame.table(
        table(
          subset(
            str_data[, c("child", "n_parents_hets", "mutate_from_homo_or_het", "parents_share_allele")],
            str_data$n_parents_hets == 1 &
              str_data$mutate_from_homo_or_het != 0 &
              !str_data$parents_share_allele &
              abs(str_data$mutsize) == 1 &
              str_data$validation != "true_de_novo"
          )
        )
      ) %>% spread(mutate_from_homo_or_het, Freq) %>%
        rename(., c("homo" = "1", "het" = "2")) %>%
        mutate(., frac_homozygous = homo / (homo + het), .keep = "unused"), 
      by.x = "ssc_id", by.y = "child"
    ),
  aes(x = SuperPop, y = frac_homozygous)
) + geom_boxplot()


# does the fraction of mutations from homozygous change as a function of whether the mutation passed our read filter?

fisher.test(
  matrix(c(
    nrow(subset(str_data, str_data$SuperPop == "EUR" & str_data$n_parents_hets == 1 & !str_data$parents_share_allele & str_data$mutate_from_homo_or_het == 1 & str_data$validation == "true_de_novo")),
    nrow(subset(str_data, str_data$SuperPop == "EUR" & str_data$n_parents_hets == 1 & !str_data$parents_share_allele & str_data$mutate_from_homo_or_het == 1 & str_data$validation != "true_de_novo")),
    nrow(subset(str_data, str_data$SuperPop == "EUR" & str_data$n_parents_hets == 1 & !str_data$parents_share_allele & str_data$mutate_from_homo_or_het == 2 & str_data$validation == "true_de_novo")),
    nrow(subset(str_data, str_data$SuperPop == "EUR" & str_data$n_parents_hets == 1 & !str_data$parents_share_allele & str_data$mutate_from_homo_or_het == 2 & str_data$validation != "true_de_novo"))
  ), nrow = 2, ncol = 2)
)

fisher.test(
  matrix(c(
    nrow(subset(str_data, str_data$SuperPop == "EUR" & str_data$n_parents_hets == 1 & !str_data$parents_share_allele & str_data$mutate_from_homo_or_het == 1 & abs(str_data$mutsize) == 1)),
    nrow(subset(str_data, str_data$SuperPop == "EUR" & str_data$n_parents_hets == 1 & !str_data$parents_share_allele & str_data$mutate_from_homo_or_het == 1 & abs(str_data$mutsize) > 1)),
    nrow(subset(str_data, str_data$SuperPop == "EUR" & str_data$n_parents_hets == 1 & !str_data$parents_share_allele & str_data$mutate_from_homo_or_het == 2 & abs(str_data$mutsize) == 1)),
    nrow(subset(str_data, str_data$SuperPop == "EUR" & str_data$n_parents_hets == 1 & !str_data$parents_share_allele & str_data$mutate_from_homo_or_het == 2 & abs(str_data$mutsize) > 1))
  ), nrow = 2, ncol = 2)
)

fisher.test(
  matrix(c(
    nrow(subset(str_data, str_data$SuperPop == "EUR" & str_data$n_parents_hets == 1 & !str_data$parents_share_allele & str_data$mutate_from_homo_or_het == 1 & abs(str_data$mutsize) == 1 & str_data$validation == "true_de_novo")),
    nrow(subset(str_data, str_data$SuperPop == "EUR" & str_data$n_parents_hets == 1 & !str_data$parents_share_allele & str_data$mutate_from_homo_or_het == 1 & abs(str_data$mutsize) == 1 & str_data$validation != "true_de_novo")),
    nrow(subset(str_data, str_data$SuperPop == "EUR" & str_data$n_parents_hets == 1 & !str_data$parents_share_allele & str_data$mutate_from_homo_or_het == 2 & abs(str_data$mutsize) == 1 & str_data$validation == "true_de_novo")),
    nrow(subset(str_data, str_data$SuperPop == "EUR" & str_data$n_parents_hets == 1 & !str_data$parents_share_allele & str_data$mutate_from_homo_or_het == 2 & abs(str_data$mutsize) == 1 & str_data$validation != "true_de_novo"))
  ), nrow = 2, ncol = 2)
)

nrow(
  subset(
    str_data[, c("child", "n_parents_hets", "mutate_from_homo_or_het", "parents_share_allele")],
    str_data$n_parents_hets == 1 &
      str_data$mutate_from_homo_or_het != 0 &
      !str_data$parents_share_allele &
      abs(str_data$mutsize) == 1 &
      str_data$validation != "true_de_novo"
  )
)


merge(
  str_data_by_indiv[, c("ssc_id", "SuperPop")],
  as.data.frame.table(
    table(
      subset(
        str_data[, c("child", "n_parents_hets", "mutate_from_homo_or_het", "parents_share_allele")],
        str_data$n_parents_hets == 1 &
          str_data$mutate_from_homo_or_het != 0 &
          !str_data$parents_share_allele &
          abs(str_data$mutsize) == 1 &
          str_data$validation != "true_de_novo"
      )
    )
  ) %>% spread(mutate_from_homo_or_het, Freq) %>%
    rename(., c("homo" = "1", "het" = "2")) %>%
    mutate(., frac_homozygous = homo / (homo + het), .keep = "unused"), 
  by.x = "ssc_id", by.y = "child"
)

ggsave(
  "../figures/parental_heterozygosity_by_false_positive_20220803.pdf",
  ggplot(
    data = str_data,
    aes(x = n_parents_hets)
  ) + geom_bar(aes(fill=validation), position = "fill")
)

ggsave(
  "../figures/len_repeat_unit_by_false_positive_20220809.pdf",
  ggplot(
    data = str_data,
    aes(x = len_repeat_unit_condensed)
  ) + geom_bar(aes(fill=validation), position = "fill")
)
table(str_data$len_repeat_unit_condensed)


str_data_by_indiv <-
  merge(
    str_data_by_indiv,
    as.data.frame.table(
      table(subset(
        str_data[, c("child", "poocase", "any_gc")],
        str_data$validation == "true_de_novo"
      ))
    ) %>% mutate(any_gc_poocase=paste(any_gc, poocase, sep = "_"), poocase = NULL, any_gc = NULL) %>%
      spread(any_gc_poocase, Freq),
    by.x = "ssc_id", by.y = "child", all.x = T
  )

summary(
  glm(
    TRUE_2 ~ fatherAgeAtBirth, 
    data = str_data_by_indiv,
    family = poisson(link = "identity")
  )
)

summary(
  glm(
    FALSE_3 ~ motherAgeAtBirth, 
    data = str_data_by_indiv,
    family = poisson(link = "identity")
  )
)

ggplot(
  data = 
    merge(
      with(
        subset(str_data, 
               str_data$poocase %in% c(2) & 
                 str_data$validation == "true_de_novo" &
                 !str_data$is_homopolymer),
        as.data.frame.table(table(child, any_gc))
      ),
      str_data_by_indiv[, c("ssc_id", "fatherAgeAtBirth")],
      by.x = "child", by.y = "ssc_id", all.x = T
    ),
  aes(x = fatherAgeAtBirth, y = Freq)
) + 
  geom_point(aes(color = any_gc), alpha = 0.2) +
  stat_smooth(aes(color = any_gc), formula = y ~ x,
              method = "glm", 
              method.args = list(family = poisson(link = "identity")))
ggsave(
  "../figures/nonhomopolymer_gc_content_by_maternal_age_read1_flank2_20220901.pdf",
  ggplot(
    data = 
      merge(
        with(
          subset(str_data, 
                 str_data$poocase %in% c(3) & 
                   str_data$validation == "true_de_novo" &
                   !str_data$is_homopolymer),
          as.data.frame.table(table(child, any_gc))
        ),
        str_data_by_indiv[, c("ssc_id", "motherAgeAtBirth")],
        by.x = "child", by.y = "ssc_id", all.x = T
      ),
    aes(x = motherAgeAtBirth, y = Freq)
  ) + 
    geom_point(aes(color = any_gc), alpha = 0.2) +
    stat_smooth(aes(color = any_gc), formula = y ~ x,
                method = "glm", 
                method.args = list(family = poisson(link = "identity")))
  
)

ggsave(
  "../figures/nonhomopolymer_gc_content_by_paternal_age_read1_flank2_20220901.pdf",
  ggplot(
    data = 
      merge(
        with(
          subset(str_data, 
                 str_data$validation == "true_de_novo" &
                   !str_data$is_homopolymer),
          as.data.frame.table(table(child, any_gc))
        ),
        str_data_by_indiv[, c("ssc_id", "fatherAgeAtBirth")],
        by.x = "child", by.y = "ssc_id", all.x = T
      ),
    aes(x = fatherAgeAtBirth, y = Freq)
  ) + 
    geom_point(aes(color = any_gc), alpha = 0.2) +
    stat_smooth(aes(color = any_gc), formula = y ~ x,
                method = "glm", 
                method.args = list(family = poisson(link = "identity")))
  
)
str_data$"short" <- ((str_data$newallele + str_data$mutsize) * str_data$period) <= 15
pat_str_dnm_by_nuc_content <-
  merge(
    with(
      subset(str_data, 
             str_data$validation == "true_de_novo" &
               !str_data$is_homopolymer &
               str_data$poocase == 2),
      as.data.frame.table(table(child, any_gc))
    ),
    str_data_by_indiv[, c("ssc_id", "fatherAgeAtBirth")],
    by.x = "child", by.y = "ssc_id", all.x = T
  )
pat_str_dnm_by_nuc_content$freq_reweight <-
  apply(
    pat_str_dnm_by_nuc_content,
    1,
    function(x)
      if(x[["any_gc"]]) as.numeric(x[["Freq"]]) * (205303 / 477663)
    else as.numeric(x[["Freq"]])
  )

library(wesanderson)
ggsave(
  "../figures/nonhomopolymer_gc_content_reweight_by_paternal_age_read1_flank2_20220907.pdf",
  ggplot(
    data = pat_str_dnm_by_nuc_content,
    aes(x = fatherAgeAtBirth, y = freq_reweight)
  ) + 
    geom_point(aes(color = any_gc), alpha = 0.2) +
    stat_smooth(aes(color = any_gc), formula = y ~ x,
                method = "glm") +
    theme_cowplot() +
    scale_color_manual(values = wes_palette("Royal1")),
  width = 7, height = 5
)

str_data$"gc_content" <-
  sapply(str_data$repeat_unit,
         function(x)
           nchar(gsub("[AT]", "", x)) /
           nchar(x))

panel_statistics <-
  read.delim("./sites_annotated_reptiming_dnm_counts_ru.txt", header = F)
names(panel_statistics) <-
  c(
    "chrom", "pos", "repeat_unit", "ref_length", "reptiming", "dnm_counts"
  )
panel_statistics$"gc_content" <-
  sapply(panel_statistics$repeat_unit,
         function(x)
           nchar(gsub("[AT]", "", x)) /
           nchar(x))
ggsave(
  "../figures/nonhomopolymer_panel_gc_content_20220916.pdf",
  ggplot(
    data = panel_statistics,
    aes(x = gc_content)
  ) + geom_histogram(bins = 5) +
    theme_cowplot() + 
    geom_vline(xintercept = 0.125, linetype = "dashed")
)

summary(
  glm(
    freq_reweight ~ fatherAgeAtBirth * any_gc,
    data = pat_str_dnm_by_nuc_content
  )
)

mat_str_dnm_by_nuc_content <-
  merge(
    with(
      subset(str_data, 
             str_data$validation == "true_de_novo" &
               !str_data$is_homopolymer & 
               str_data$poocase == 3),
      as.data.frame.table(table(child, any_gc))
    ),
    str_data_by_indiv[, c("ssc_id", "motherAgeAtBirth")],
    by.x = "child", by.y = "ssc_id", all.x = T
  )
mat_str_dnm_by_nuc_content$freq_reweight <-
  apply(
    mat_str_dnm_by_nuc_content,
    1,
    function(x)
      if(x[["any_gc"]]) as.numeric(x[["Freq"]]) * (205303 / 477663)
    else as.numeric(x[["Freq"]])
  )
ggsave(
  "../figures/nonhomopolymer_gc_content_reweight_by_maternal_age_read1_flank2_20220907.pdf",
  ggplot(
    data = mat_str_dnm_by_nuc_content,
    aes(x = motherAgeAtBirth, y = freq_reweight)
  ) + 
    geom_point(aes(color = any_gc), alpha = 0.2) +
    stat_smooth(aes(color = any_gc), formula = y ~ x,
                method = "glm") +
    theme_cowplot() +
    scale_color_manual(values = wes_palette("Royal1")),
  width = 7, height = 5
)

summary(
  glm(
    freq_reweight ~ motherAgeAtBirth * any_gc,
    data = mat_str_dnm_by_nuc_content
  )
)



summary(
  glm(
    Freq ~ fatherAgeAtBirth * any_gc,
    data = 
      merge(
        with(
          subset(str_data, 
                 str_data$poocase %in% c(2) & 
                   str_data$validation == "true_de_novo" &
                   !str_data$is_homopolymer),
          as.data.frame.table(table(child, any_gc))
        ),
        str_data_by_indiv[, c("ssc_id", "fatherAgeAtBirth")],
        by.x = "child", by.y = "ssc_id", all.x = T
      ),
    family = poisson(link = "identity")
  )
)

summary(
  glm(
    Freq ~ motherAgeAtBirth * any_gc,
    data = 
      merge(
        with(
          subset(str_data, 
                 str_data$poocase %in% c(3) & 
                   str_data$validation == "true_de_novo" &
                   !str_data$is_homopolymer),
          as.data.frame.table(table(child, any_gc))
        ),
        str_data_by_indiv[, c("ssc_id", "motherAgeAtBirth")],
        by.x = "child", by.y = "ssc_id", all.x = T
      ),
    family = poisson(link = "identity")
  )
)

repeat_units_sequenced <-
  read.delim("./repeat_units_sequenced.txt", header = F)
colnames(repeat_units_sequenced) <- 
  c("repeat_unit", "N")
repeat_units_sequenced$"any_gc" <-
  sapply(
    repeat_units_sequenced$repeat_unit,
    function(x)
      nchar(sub("[CG]", "", x)) < nchar(x)
  )
repeat_units_sequenced$"is_homopolymer" <-
  nchar(repeat_units_sequenced$repeat_unit) == 1
sum(subset(repeat_units_sequenced$N, 
           repeat_units_sequenced$any_gc & 
             !repeat_units_sequenced$is_homopolymer))
sum(subset(repeat_units_sequenced$N, 
           !repeat_units_sequenced$any_gc & 
             !repeat_units_sequenced$is_homopolymer))





