## Mitochondrial genes and lung cancer driver scripts

library(rmarkdown)
library(fs)

#### Data preparation ####
rmarkdown::render("lung-cancer-mt-TCGAdata.Rmd", output_file = "LUAD_genes/mt-LUAD.pdf",
                    params = list(cancer = "LUAD"))
rmarkdown::render("lung-cancer-mt-TCGAdata.Rmd", output_file = "LUSC_genes/mt-LUSC.pdf",
                  params = list(cancer = "LUSC"))

# create list of genes
gene_list <- read_csv("mt-genes.csv")

#### create report for each gene: LUAD ####

# create list of reports
reports_LUAD <- tibble(
  each_gene = gene_list$gene,
  filename = stringr::str_c("LUAD-", each_gene, ".pdf"),
  params = purrr::map(each_gene, ~ list(gene = ., cancer = "LUAD"))
)

# create p-value tables
write_csv(tibble(mt_gene = "", pvalue_ttest = "", overexpressed = "", pvalue_surv = ""), "mt-LUAD-pvalue.csv")

# run reports
reports_LUAD %>% 
  select(output_file = filename, params) %>% 
  purrr::pwalk(rmarkdown::render, input = "lung-cancer-mt-gene.Rmd", output_dir = "LUAD_genes")

#### create report for each gene: LUSC ####

# create list of reports
reports_LUSC <- tibble(
  each_gene = gene_list$gene,
  filename = stringr::str_c("LUSC-", each_gene, ".pdf"),
  params = purrr::map(each_gene, ~ list(gene = ., cancer = "LUSC"))
)

# create header for p-value table
write_csv(tibble(mt_gene = "", pvalue_ttest = "", overexpressed = "", pvalue_surv = ""), "mt-LUSC-pvalue.csv")

# run reports
reports_LUSC %>% 
  select(output_file = filename, params) %>% 
  purrr::pwalk(rmarkdown::render, input = "lung-cancer-mt-gene.Rmd", output_dir = "LUSC_genes")

#### summary reports ####

rmarkdown::render("lung-cancer-mt-summary.Rmd", output_file = "summary-mt-LUAD.pdf",
                  params = list(cancer = "LUAD"))
rmarkdown::render("lung-cancer-mt-summary.Rmd", output_file = "summary-mt-LUSC.pdf",
                  params = list(cancer = "LUSC"))

rmarkdown::render("LUAD-mt-sig.Rmd")
rmarkdown::render("LUSC-mt-sig.Rmd")

