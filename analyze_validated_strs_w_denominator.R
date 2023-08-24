# checking out the validated str dnms

options(stringsAsFactors = F)
setwd("~/Documents/Harris/Simons/for_michael/")

library(ggplot2)
library(tidyr)
library(dplyr)
library(ggridges)
library(cowplot)

parental_ages <- read.table("../SSC_parental_age_at_birth_in_years.txt", header = T)
ped_file <- read.delim("./ssc_copy.ped", header = F, sep = " ")
colnames(ped_file) <- c("family", "Sample", "father.id", "mother.id", "asd_state", "sex", "ssc_id")
parental_ages <-
  merge(
    parental_ages,
    ped_file[, c("Sample", "ssc_id")],
    by.x = "sample", by.y = "Sample",
    all.x = F, all.y = F
  )

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
str_data$"abs_mut_size" <-
  abs(str_data$mutsize)


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

# this needs work!

discoverable_str_data_by_indiv <-
  merge(
    str_denominators %>%
      filter(., ru_len > 1) %>%
      group_by(., family, dnm_size) %>%
      summarize(., n = sum(n)),
    merge(
      as.data.frame.table(
        table(
          str_data[, c("child", 
                       "poocase", 
                       "abs_mut_size", 
                       "validation", 
                       "len_repeat_unit", 
                       "all_dnms_discoverable")]
        )) %>%
          filter(., validation == "true_de_novo" &
                   len_repeat_unit %in% c(2, 3, 4) &
                   abs_mut_size %in% c(1, 2, 3) &
                   all_dnms_discoverable == 1) %>%
        mutate(., len_repeat_unit = as.numeric(len_repeat_unit),
               abs_mut_size = as.numeric(abs_mut_size)) %>%
        group_by(., child, poocase, abs_mut_size) %>%
        summarize(., Freq = sum(Freq)),
      unique(str_data[, c("child", "family")]), 
      by = "child", all.x = T, all.y = T
    ),
    by.x = c("family", "dnm_size"), by.y = c("family", "abs_mut_size"),
    all.x = F, all.y = F
  )

# only 833 maternally phased mutations that are fully discoverable and not at homopolymers
# nrow(subset(str_data, str_data$poocase == 3 & 
#               str_data$all_dnms_discoverable == 1 & 
#               str_data$validation == "true_de_novo" & 
#               str_data$len_repeat_unit > 1))

discoverable_mut_rate_by_indiv <-
  merge(
    discoverable_str_data_by_indiv %>%
      group_by(., child, family, dnm_size, poocase) %>%
      summarize(., dnm_rate = Freq / (2 * n)) %>%
      group_by(., child, family, poocase) %>%
      summarize(., mut_rate = sum(dnm_rate)),
    parental_ages[, c("ssc_id", "motherAgeAtBirth", "fatherAgeAtBirth")],
    by.x = "child", by.y = "ssc_id", all.y = F, all.x = F
  )

ggplot(
  data = discoverable_mut_rate_by_indiv %>% filter(., poocase == 3),
  aes(x = motherAgeAtBirth, y = mut_rate)
) + 
  geom_point() + 
  stat_smooth(formula = y ~ x, method = "glm")

summary(
  glm(
    mut_rate ~ motherAgeAtBirth,
    data = discoverable_mut_rate_by_indiv %>% 
      filter(., poocase == 3)
  )
)
summary(
  glm(
    mut_rate ~ fatherAgeAtBirth,
    data = discoverable_mut_rate_by_indiv %>% 
      filter(., poocase == 2)
  )
)
summary(
  glm(
    mut_rate_across_phase ~ motherAgeAtBirth + fatherAgeAtBirth,
    data = discoverable_mut_rate_by_indiv %>%
        group_by(., child, family, motherAgeAtBirth, fatherAgeAtBirth) %>%
        summarize(., mut_rate_across_phase = sum(mut_rate)),
  )
)
summary(
  glm(
    mut_rate_across_phase ~  fatherAgeAtBirth,
    data = discoverable_mut_rate_by_indiv %>%
      group_by(., child, family, motherAgeAtBirth, fatherAgeAtBirth) %>%
      summarize(., mut_rate_across_phase = sum(mut_rate)),
  )
)

summary(
  glm(
    Freq ~ motherAgeAtBirth, offset(log(n)), family = poisson(link = "log"),
    data = 
      merge(
        discoverable_str_data_by_indiv %>%
          group_by(., child, family, poocase) %>%
          summarize(., n = sum(n) * 2, Freq = sum(Freq)),
        parental_ages[, c("ssc_id", "motherAgeAtBirth")],
        by.x = "child", by.y = "ssc_id", all.y = F
      ) %>%
      filter(., poocase == 3)
  )
)
summary(
  glm(
    Freq ~ fatherAgeAtBirth, offset(log(n)), family = poisson(link = "log"),
    data = 
      merge(
        discoverable_str_data_by_indiv %>%
          group_by(., child, family, poocase) %>%
          summarize(., n = sum(n) * 2, Freq = sum(Freq)),
        parental_ages[, c("ssc_id", "fatherAgeAtBirth")],
        by.x = "child", by.y = "ssc_id", all.y = F
      ) %>%
      filter(., poocase == 2)
  )
)


