## Mitochondrial genes and lung cancer: manuscript figures

library(fs)
library(tidyverse)
library(ggpubr)
library(rstatix)

#### expression subplot function ####

exp_subplot <- function(fpkm, genes) {
  subplot <- fpkm %>%
    dplyr::select(genes$mt_gene, shortLetterCode) %>%
    filter(shortLetterCode == "NT" | shortLetterCode == "TP") %>%
    pivot_longer(cols = -shortLetterCode,
                 names_to = "gene", values_to = "expression") %>%
    ggplot() +
      geom_boxplot(aes(shortLetterCode, expression, color = shortLetterCode)) +
      facet_wrap(~ gene, ncol = 6, scales = "free") +
      labs(y = "log2 expression") +
      scale_x_discrete(labels = c("NT"="norm", "TP"="tum")) +
      scale_color_brewer(type = "qual") +
      geom_text(data = genes, aes(x = -Inf, y = Inf, 
                                  label = rounded),
                hjust = -0.1, vjust = 1.5, size = 1.5) +
      theme_bw() +
      theme(axis.title.y = element_blank(),
            axis.text.y = element_text(size = 6),
            axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            legend.position = "none",
            strip.text.x = element_text(size = 6, margin = margin(1, 1, 1, 1)))
}
luad_ttest_a <- exp_subplot(fpkm_luad, luad_not_sig) +
  theme(strip.background = element_rect(fill="gray"))
luad_ttest_a 

#### collapse stages function ####

collapse_stages <- function(dat, stage_column) {
  stages <- dat %>%
    filter(shortLetterCode == "TP") %>%
    # remove subcategories
    mutate(collapsed_stages = str_remove({{ stage_column }}, "[ABCD]")) %>%
    # combine III and IV (to Stage III)
    mutate(collapsed_stages = str_replace_all(collapsed_stages, "Stage IV", "Stage III")) %>%
    filter(!is.na(collapsed_stages))
}
luad_stages <- collapse_stages(fpkm_luad, ajcc_pathologic_stage)

#### stages stats function ####

stages_stats <- function(fpkm, genes) {
  fpkm_stats <- fpkm %>%
    select(genes$mt_gene, collapsed_stages) %>%
    pivot_longer(cols = -collapsed_stages,
                 names_to = "gene", values_to = "expression") %>%
    group_by(gene) %>%
    t_test(expression ~ collapsed_stages) %>%
    adjust_pvalue(method = "bonferroni") %>%
    add_significance() %>%
    add_xy_position(x = "collapsed_stages")
}
luad_stages_stats_a <- stages_stats(luad_stages, luad_sig_down)

#### stages subplot function ####

stages_subplot <- function(fpkm, genes, stats) {
  subplot <- fpkm %>%
    select(genes$mt_gene, collapsed_stages) %>%
    pivot_longer(cols = genes$mt_gene,
                 names_to = "gene", values_to = "expression") %>%
    ggplot() +
      geom_boxplot(aes(collapsed_stages, expression, color = collapsed_stages)) +
      facet_wrap("gene", ncol = 6, scales = "free") +
      labs(y = "log2 expression") +
      scale_x_discrete(
           labels = c("Stage I"="I", "Stage II"="II", "Stage III"="III/IV")) +
      scale_color_brewer(type = "qual") +
      stat_pvalue_manual(stats, size = 2, hide.ns = TRUE) +
      theme_bw() +
      theme(axis.title.y = element_blank(),
          axis.text.y = element_text(size = 6),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = "none",
          strip.text.x = element_text(size = 6, margin = margin(1, 1, 1, 1)))
}
lusc_stages_a <- stages_subplot(lusc_stages, luad_sig_down, luad_stages_stat_a) +
  theme(strip.background = element_rect(fill="#ffff99"))
lusc_stages_a 

#### smoking stats function ####

smoking_stats <- function(fpkm, genes) {
  fpkm_stats <- fpkm %>%
    filter(shortLetterCode == "TP") %>%
    filter(paper_Smoking.Status != "[Not Available]") %>% 
    filter(paper_Smoking.Status != "N/A") %>%
    select(genes$mt_gene, paper_Smoking.Status) %>%
    pivot_longer(cols = -paper_Smoking.Status,
                 names_to = "gene", values_to = "expression") %>%
    group_by(gene) %>%
    t_test(expression ~ paper_Smoking.Status) %>%
    adjust_pvalue(method = "bonferroni") %>%
    add_significance() %>%
    add_xy_position(x = "paper_Smoking.Status")
}
luad_smoking_stats_a <- smoking_stats(fpkm_luad, luad_sig_down)

#### smoking subplot function ####

smoking_subplot <- function(fpkm, genes, stats) {
  subplot <- fpkm %>%
    filter(shortLetterCode == "TP") %>%
    filter(paper_Smoking.Status != "[Not Available]") %>% 
    filter(paper_Smoking.Status != "N/A") %>%
    arrange(paper_Smoking.Status) %>%
    mutate(smoking = factor(paper_Smoking.Status, levels=c("Lifelong Non-smoker", "Current reformed smoker for > 15 years", "Current reformed smoker for < or = 15 years", "Current smoker"))) %>%
    select(genes$mt_gene, smoking) %>%
    pivot_longer(cols = genes$mt_gene,
                 names_to = "gene", values_to = "expression") %>%
    ggplot() +
    geom_boxplot(aes(smoking, expression, color = smoking)) +
    facet_wrap("gene", ncol = 6, scales = "free") +
    labs(y = "log2 expression") +
    scale_x_discrete(
      labels = c("Lifelong Non-smoker"="Non", 
                 "Current reformed smoker for > 15 years"=">15", 
                 "Current reformed smoker for < or = 15 years"="<=15", 
                 "Current smoker"="current")) +
    scale_color_brewer(type = "qual") +
    stat_pvalue_manual(stats, size = 2, hide.ns = TRUE) +
    theme_bw() +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_text(size = 6),
          axis.title.x = element_blank(),
          #axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = "none",
          strip.text.x = element_text(size = 6, margin = margin(1, 1, 1, 1)))
}
lusc_smoke_a <- smoking_subplot(fpkm_luad, luad_sig_down, luad_smoking_stats_a) +
  theme(strip.background = element_rect(fill="#ffff99"))
lusc_smoke_a 

#### survival plot function ####

survival_plot <- function(fpkm, gene, stats) {
  fpkmTum <- fpkm %>%
    select(gene, shortLetterCode, vital_status, days_to_death) %>%
    filter(shortLetterCode == "NT" | shortLetterCode == "TP") %>%
    # code vital status for survival analysis 
    mutate(vital = as.numeric(as_factor(vital_status)) - 1)
  # rename target gene column
  names(fpkmTum)[1] <- "test_gene"
  # adding column for gene expression categories
  fpkmTum <- fpkmTum %>%
    mutate(gene_hl = test_gene)
  # determine cutpoint
  cut <- surv_cutpoint(fpkmTum, time = "days_to_death", 
                       event = "vital", variables = "gene_hl")
  # categorize based on cutpoint
  fpkmSurv <- surv_categorize(cut)
  # fit model
  death_fit <- survfit(Surv(days_to_death, vital) ~ gene_hl, 
                       data = fpkmSurv)
  # plot
  survPlot <- ggsurvplot(death_fit, 
                         data = fpkmSurv, 
                         conf.int = TRUE,
                         break.time.by = 500,
                         xlab = "days to death",
                         ggtheme = theme_bw(),
                         legend = "none",
                         # green=hi, purple=lo
                         palette = c("#7fc97f", "#beaed4"))
  survPlot$plot <- survPlot$plot +
    ggplot2::annotate(geom = "text", 
                      label = paste0(gene, ", ", 
                                     stats$rounded[stats$mt_gene == gene]), 
                      x = Inf, y = Inf, 
                      hjust = 1.1, vjust = 1.5)
}
surv_plot <- survival_plot(fpkm_luad, "ACSF2", surv_luad)
surv_plot

#### subtypes stats function ####

subtype_stats <- function(fpkm, genes) {
  fpkm_stats <- fpkm %>%
    filter(!is.na(subtypes)) %>%
    select(genes$mt_gene, subtypes) %>%
    pivot_longer(cols = -subtypes,
                 names_to = "gene", values_to = "expression") %>%
    group_by(gene) %>%
    t_test(expression ~ subtypes) %>%
    adjust_pvalue(method = "bonferroni") %>%
    add_significance() %>%
    add_xy_position(x = "subtypes")
}
luad_subtype_stats_a <- subtype_stats(fpkm_luad, luad_sig_down)

#### subtypes plot function ####

subtype_subplot <- function(fpkm, genes, stats) {
  subplot <- fpkm %>%
    filter(!is.na(subtypes)) %>%
    select(genes$mt_gene, subtypes) %>%
    pivot_longer(cols = genes$mt_gene,
                 names_to = "gene", values_to = "expression") %>%
    ggplot() +
      geom_boxplot(aes(subtypes, expression, color = subtypes)) +
      facet_wrap("gene", ncol = 6, scales = "free") +
      labs(y = "log2 expression") +
      #scale_x_discrete(
      #  labels = c("Stage I"="I", "Stage II"="II", "Stage III"="III/IV")) +
      scale_color_brewer(type = "qual") +
      stat_pvalue_manual(stats, size = 2, hide.ns = TRUE) +
      theme_bw() +
      theme(axis.title.y = element_blank(),
            axis.text.y = element_text(size = 6),
            axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            legend.position = "none",
            strip.text.x = element_text(size = 6, margin = margin(1, 1, 1, 1)))
}
luad_subtype_a <- subtype_subplot(fpkm_luad, luad_sig_down, luad_subtype_stats_a) +
  theme(strip.background = element_rect(fill="#ffff99"))
luad_subtype_a 

#### collapse race function ####

collapse_race <- function(fpkm) {
  stages <- fpkm %>%
    filter(shortLetterCode == "TP") %>%
    mutate(race_cat = ifelse(race == "black or african american" & ethnicity == "not hispanic or latino", "AA",
                      ifelse(race == "white" & ethnicity == "hispanic or latino", "HIS",
                      ifelse(race == "white" & ethnicity == "not hispanic or latino", "CA", NA))))
}
luad_race <- collapse_race(fpkm_luad)

#### race stats function ####

race_stats <- function(fpkm, genes) {
  fpkm_race <- fpkm %>%
    filter(!is.na(race_cat), race_cat != "HIS") %>%
    select(genes$mt_gene, race_cat) %>%
    pivot_longer(cols = -race_cat,
                 names_to = "gene", values_to = "expression") %>%
    group_by(gene) %>%
    t_test(expression ~ race_cat) %>%
    adjust_pvalue(method = "bonferroni") %>%
    add_significance() %>%
    add_xy_position(x = "race_cat")
}
luad_race_stats_a <- race_stats(luad_race, luad_sig_down)

#### Data prep: LUAD ####

# read in data
fpkm_luad <- read_csv("TCGA-LUAD.csv")
ttest_luad <- read_csv("mt-LUAD-ttest-sig.csv")

# not significant
luad_not_sig <- ttest_luad %>%
  dplyr::filter(ttest_sig == FALSE) %>%
  dplyr::mutate(rounded = paste("p = ", signif(pval_ttest_adj, 3)),
                gene = mt_gene)

# significant, downregulated in cancer
luad_sig_down <- ttest_luad %>%
  dplyr::filter(ttest_sig == TRUE, overexpressed == FALSE) %>%
  dplyr::mutate(rounded = paste("p = ", signif(pval_ttest_adj, 3)),
                gene = mt_gene)

# significant, upregulated in cancer
luad_sig_up <- ttest_luad %>%
  dplyr::filter(ttest_sig == TRUE, overexpressed == TRUE) %>%
  dplyr::mutate(rounded = paste("p = ", signif(pval_ttest_adj, 3)),
                gene = mt_gene)

#### Data prep: LUSC ####

# read in data
fpkm_lusc <- read_csv("TCGA-lusc.csv")
ttest_lusc <- read_csv("mt-lusc-ttest-sig.csv")

# not significant
lusc_not_sig <- ttest_lusc %>%
  dplyr::filter(ttest_sig == FALSE) %>%
  dplyr::mutate(rounded = paste("p = ", signif(pval_ttest_adj, 3)),
                gene = mt_gene)

# significant, downregulated in cancer
lusc_sig_down <- ttest_lusc %>%
  dplyr::filter(ttest_sig == TRUE, overexpressed == FALSE) %>%
  dplyr::mutate(rounded = paste("p = ", signif(pval_ttest_adj, 3)),
                gene = mt_gene)

# significant, upregulated in cancer
lusc_sig_up <- ttest_lusc %>%
  dplyr::filter(ttest_sig == TRUE, overexpressed == TRUE) %>%
  dplyr::mutate(rounded = paste("p = ", signif(pval_ttest_adj, 3)),
                gene = mt_gene)

#### Fig 1 Expression graphs for all genes in LUAD ####

nrow(filter(fpkm_luad, shortLetterCode == "NT"))
nrow(filter(fpkm_luad, shortLetterCode == "TP"))

# A: not significant
luad_ttest_a <- exp_subplot(fpkm_luad, luad_not_sig) +
  theme(strip.background = element_rect(fill="gray"))
ggsave(plot = luad_ttest_a, 
       filename = "luad_expression_a.jpg", 
       height = 1, width = 5)
# B: significant, downregulated in cancer
luad_ttest_b <- exp_subplot(fpkm_luad, luad_sig_down) +
  theme(strip.background = element_rect(fill="#ffff99"))
ggsave(plot = luad_ttest_b, 
       filename = "luad_expression_b.jpg", 
       height = 1, width = 6)
# C: significant, upregulated in cancer
luad_ttest_c <- exp_subplot(fpkm_luad, luad_sig_up) +
  theme(strip.background = element_rect(fill="#fdc086"))
ggsave(plot = luad_ttest_c, 
       filename = "luad_expression_c.jpg", 
       height = 6, width = 6)

#### Fig 2 Expression graphs for all genes in LUSC ####

nrow(filter(fpkm_lusc, shortLetterCode == "NT"))
nrow(filter(fpkm_lusc, shortLetterCode == "TP"))

# A: not significant
lusc_ttest_a <- exp_subplot(fpkm_lusc, lusc_not_sig) +
  theme(strip.background = element_rect(fill="gray"))
ggsave(plot = lusc_ttest_a, 
       filename = "lusc_expression_a.jpg", 
       height = 1, width = 3)

# B: significant, downregulated in cancer
lusc_ttest_b <- exp_subplot(fpkm_lusc, lusc_sig_down) +
  theme(strip.background = element_rect(fill="#ffff99"))
ggsave(plot = lusc_ttest_b, 
       filename = "lusc_expression_b.jpg", 
       height = 2, width = 6)

# C: significant, upregulated in cancer
lusc_ttest_c <- exp_subplot(fpkm_lusc, lusc_sig_up) +
  theme(strip.background = element_rect(fill="#fdc086"))
ggsave(plot = lusc_ttest_c, 
       filename = "lusc_expression_c.jpg", 
       height = 5, width = 6)

#### Fig 3 Tumor stage and gene expression in LUAD ####

# prep stages data
luad_stages <- collapse_stages(fpkm_luad, ajcc_pathologic_stage)

table(luad_stages$collapsed_stages)

# A: expression significant, downregulated in cancer
luad_stages_stat_a <- stages_stats(luad_stages, luad_sig_down)
luad_stages_a <- stages_subplot(luad_stages, luad_sig_down, luad_stages_stat_a) +
  theme(strip.background = element_rect(fill="#ffff99"))
ggsave(plot = luad_stages_a, 
       filename = "luad_stages_a.jpg", 
       height = 1.5, width = 6)

# B: expression significant, upregulated in cancer
luad_stages_stat_b <- stages_stats(luad_stages, luad_sig_up)
luad_stages_b <- stages_subplot(luad_stages, luad_sig_up, luad_stages_stat_b) +
  theme(strip.background = element_rect(fill="#fdc086"))
ggsave(plot = luad_stages_b, 
       filename = "luad_stages_b.jpg", 
       height = 6, width = 6)

#### Fig 4 Tumor stage and gene expression in LUSC ####

# prep stages data
lusc_stages <- collapse_stages(fpkm_lusc, ajcc_pathologic_stage)

table(lusc_stages$collapsed_stages)

# A: expression significant, downregulated in cancer
lusc_stages_stat_a <- stages_stats(lusc_stages, lusc_sig_down)
lusc_stages_a <- stages_subplot(lusc_stages, lusc_sig_down, lusc_stages_stat_a) +
  theme(strip.background = element_rect(fill="#ffff99"))
ggsave(plot = lusc_stages_a, 
       filename = "lusc_stages_a.jpg", 
       height = 2, width = 6)

# B: expression significant, upregulated in cancer
lusc_stages_stat_b <- stages_stats(lusc_stages, lusc_sig_up)
lusc_stages_b <- stages_subplot(lusc_stages, lusc_sig_up, lusc_stages_stat_b) +
  theme(strip.background = element_rect(fill="#fdc086"))
ggsave(plot = lusc_stages_b, 
       filename = "lusc_stages_b.jpg", 
       height = 5, width = 6)

#### Fig 5 Survival for significant (expression) genes in LUAD ####

# import pvalues
surv_luad <- read_csv("mt-LUAD-surv-sig.csv")
# extract only significant genes
surv_luad <- surv_luad %>%
  filter(surv_sig == TRUE) %>%
  mutate(rounded = paste("p = ", signif(pval_surv_adj, 3)),
         gene = mt_gene)
surv_luad

surv_ACSF2 <- survival_plot(fpkm_luad, "ACSF2", surv_luad)
surv_ACSF2
ggsave(plot = surv_ACSF2, 
       filename = "luad_surv_ACSF2.jpg",
       height = 2, width = 4)

surv_GCAT <- survival_plot(fpkm_luad, "GCAT", surv_luad)
surv_GCAT
ggsave(plot = surv_GCAT, 
       filename = "luad_surv_GCAT.jpg",
       height = 2, width = 4)

surv_PDK2 <- survival_plot(fpkm_luad, "PDK2", surv_luad)
surv_PDK2
ggsave(plot = surv_PDK2, 
       filename = "luad_surv_PDK2.jpg",
       height = 2, width = 4)

surv_SLC25A4 <- survival_plot(fpkm_luad, "SLC25A4", surv_luad)
surv_SLC25A4
ggsave(plot = surv_SLC25A4, 
       filename = "luad_surv_SLC25A4.jpg",
       height = 2, width = 4)

#### Fig 6 Survival for significant (expression) genes in LUSC ####

# import pvalues
surv_lusc <- read_csv("mt-LUSC-surv-sig.csv")
# extract only significant genes
surv_lusc <- surv_lusc %>%
  filter(surv_sig == TRUE) %>%
  mutate(rounded = paste("p = ", signif(pval_surv_adj, 3)),
         gene = mt_gene)
sig_surv_lusc

surv_MACROD1 <- survival_plot(fpkm_lusc, "MACROD1", surv_lusc)
surv_MACROD1
ggsave(plot = surv_MACROD1, 
       filename = "lusc_surv_MACROD1.jpg",
       height = 2, width = 4)

surv_SLC25A4 <- survival_plot(fpkm_lusc, "SLC25A4", surv_lusc)
surv_SLC25A4
ggsave(plot = surv_SLC25A4, 
       filename = "lusc_surv_SLC25A4.jpg",
       height = 2, width = 4)

#### Fig 7 Smoking and expression in LUAD ####

fpkm_luad %>%
  filter(shortLetterCode == "TP") %>%
  count(paper_Smoking.Status)
  
# A: expression significant, downregulated in cancer
luad_smoking_stats_a <- smoking_stats(fpkm_luad, luad_sig_down)
luad_smoke_a <- smoking_subplot(fpkm_luad, luad_sig_down, luad_smoking_stats_a) +
  theme(strip.background = element_rect(fill="#ffff99"),
        axis.text.x = element_blank())
ggsave(plot = luad_smoke_a, 
       filename = "luad_smoking_a.jpg", 
       height = 1.25, width = 6) 

# B: expression significant, upregulated in cancer
luad_smoking_stats_b <- smoking_stats(fpkm_luad, luad_sig_up)
luad_smoke_b <- smoking_subplot(fpkm_luad, luad_sig_up, luad_smoking_stats_b) +
  theme(strip.background = element_rect(fill="#fdc086"),
        axis.text.x = element_blank())
ggsave(plot = luad_smoke_b, 
       filename = "luad_smoking_b.jpg", 
       height = 6, width = 6)

#### Fig 8 Smoking and expression in LUSC ####

fpkm_lusc %>%
  filter(shortLetterCode == "TP") %>%
  count(paper_Smoking.Status)

# A: expression significant, downregulated in cancer
lusc_smoking_stat_a <- smoking_stats(fpkm_lusc, lusc_sig_down)
lusc_smoke_a <- smoking_subplot(fpkm_lusc, lusc_sig_down, lusc_smoking_stat_a) +
  theme(strip.background = element_rect(fill="#ffff99"),
        axis.text.x = element_blank())
ggsave(plot = lusc_smoke_a, 
       filename = "lusc_smoking_a.jpg", 
       height = 2, width = 6)

# B: expression significant, upregulated in cancer
lusc_smoking_stat_b <- smoking_stats(fpkm_lusc, lusc_sig_up)
lusc_smoke_b <- smoking_subplot(fpkm_lusc, lusc_sig_up, lusc_smoking_stat_b) +
  theme(strip.background = element_rect(fill="#fdc086"),
        axis.text.x = element_blank())
ggsave(plot = lusc_smoke_b, 
       filename = "lusc_smoking_b.jpg", 
       height = 5, width = 6)


#### Fig 7 Pathways for significant genes (expression and survival) ####


#### Fig 8 IHC analysis of SLC25A4 in both LUAD and LUSC and correlation analysis with stage and survival ####
#### Fig 9 expression subtype LUAD ####

fpkm_luad <- fpkm_luad %>%
  rename(subtypes = paper_expression_subtype)
table(fpkm_luad$subtypes)

# A: expression significant, downregulated in cancer
luad_subtype_stats_a <- subtype_stats(fpkm_luad, luad_sig_down)
luad_subtype_a <- subtype_subplot(fpkm_luad, luad_sig_down, luad_subtype_stats_a) +
  theme(strip.background = element_rect(fill="#ffff99"))
luad_subtype_a 
ggsave(plot = luad_subtype_a, 
       filename = "luad_subtype_a.jpg", 
       height = 1.5, width = 6)

# B: expression significant, upregulated in cancer
luad_subtype_stats_b <- subtype_stats(fpkm_luad, luad_sig_up)
luad_subtype_b <- subtype_subplot(fpkm_luad, luad_sig_up, luad_subtype_stats_b) +
  theme(strip.background = element_rect(fill="#ffff99"))
luad_subtype_b 
ggsave(plot = luad_subtype_b, 
       filename = "luad_subtype_b.jpg", 
       height = 6, width = 6)

#### Fig 10 expression subtype LUSC ####

fpkm_lusc <- fpkm_lusc %>%
  rename(subtypes = paper_Expression.Subtype)
table(fpkm_lusc$subtypes)

# A: expression significant, downregulated in cancer
lusc_subtype_stats_a <- subtype_stats(fpkm_lusc, lusc_sig_down)
lusc_subtype_a <- subtype_subplot(fpkm_lusc, lusc_sig_down, lusc_subtype_stats_a) +
  theme(strip.background = element_rect(fill="#ffff99"))
lusc_subtype_a 
ggsave(plot = lusc_subtype_a, 
       filename = "lusc_subtype_a.jpg", 
       height = 2, width = 6)

# B: expression significant, upregulated in cancer
lusc_subtype_stats_b <- subtype_stats(fpkm_lusc, lusc_sig_up)
lusc_subtype_b <- subtype_subplot(fpkm_lusc, lusc_sig_up, lusc_subtype_stats_b) +
  theme(strip.background = element_rect(fill="#ffff99"))
lusc_subtype_b
ggsave(plot = lusc_subtype_b, 
       filename = "lusc_subtype_b.jpg", 
       height = 5, width = 6)

#### Fig 11 race LUAD ####

# overlap between race and ethnicity
fpkm_luad %>%
  filter(shortLetterCode == "TP") %>%
  group_by(race, ethnicity) %>%
  tally()
#only hispanic/latino are white

# A: expression significant, downregulated in cancer
luad_race <- collapse_race(fpkm_luad)
count(luad_race, race_cat)

luad_race_stats_a <- race_stats(luad_race, luad_sig_down)
luad_race_a <- race_subplot(luad_race, luad_sig_down, luad_race_stats_a) +
  theme(strip.background = element_rect(fill="#ffff99"))
luad_race_a 
ggsave(plot = luad_race_a, 
       filename = "luad_race_a.jpg", 
       height = 1.5, width = 6)

# B: expression significant, upregulated in cancer
luad_race <- collapse_race(fpkm_luad)
luad_race_stats_b <- race_stats(luad_race, luad_sig_up)
luad_race_b <- race_subplot(luad_race, luad_sig_up, luad_race_stats_b) +
  theme(strip.background = element_rect(fill="#ffff99"))
luad_race_b
ggsave(plot = luad_race_b, 
       filename = "luad_race_b.jpg", 
       height = 6, width = 6)

#### Fig 12 race LUSC ####

# overlap between race and ethnicity
fpkm_luad %>%
  group_by(race, ethnicity) %>%
  tally()
#only hispanic/latino are white or not reported

# A: expression significant, downregulated in cancer
lusc_race <- collapse_race(fpkm_lusc, lusc_sig_down)


lusc_race_stats_a <- subtype_stats(fpkm_lusc, lusc_sig_down)
lusc_subtype_a <- subtype_subplot(fpkm_lusc, lusc_sig_down, lusc_subtype_stats_a) +
  theme(strip.background = element_rect(fill="#ffff99"))
lusc_subtype_a 
ggsave(plot = lusc_subtype_a, 
       filename = "lusc_subtype_a.jpg", 
       height = 2, width = 6)

# B: expression significant, upregulated in cancer
lusc_subtype_stats_b <- subtype_stats(fpkm_lusc, lusc_sig_up)
lusc_subtype_b <- subtype_subplot(fpkm_lusc, lusc_sig_up, lusc_subtype_stats_b) +
  theme(strip.background = element_rect(fill="#ffff99"))
lusc_subtype_b
ggsave(plot = lusc_subtype_b, 
       filename = "lusc_subtype_b.jpg", 
       height = 5, width = 6)
