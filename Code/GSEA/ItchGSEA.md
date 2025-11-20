Itch GSEA
==============
Erik Larsen
2025-02-10

# Overview

This markdown was developed to formally modernize the gene set enrichment analysis of bulk RNA-seq data for **Dr. Martha Bhattacharya's lab** at the **University of Arizona** for the [Itch paper](https://journals.lww.com/pain/abstract/2022/05000/transmembrane_protein_tmem184b_is_necessary_for.18.aspx)

# Set up Environment

## Load packages

+ (install and) attach necessary packages


```r
packages <- c('BiocManager', 'tidyverse', 'stringr', 'GO.db', 'org.Mm.eg.db', 'forcats', 'ggplot2', 'ggrepel')
for ( package in packages ) {
  if (!require(package, character.only = T)) {
    install.packages(package)
    library(package, character.only = T)
  }
}
```

## Write Functions for Importing and Wrangling Data

### aDRG import function


```r
read_in_aDRG_GO_and_PANTHER_analysis_files <- function(category){
  if ( !grepl(category, pattern = 'PANTHER') |> any() ) {
    
    assign(
    paste0('aDRG_GO_',
           category,
           '_results'),
    value = read.table(
    paste0('path/to/file',
           'aDRG_GO_analysis_',
           stringr::str_to_upper(category),
           '.txt'
           ),
    sep = '\t', skip = 11, header = T, quote = ''
    ),
    envir = .GlobalEnv)
    
  } else if (grepl(category, pattern = 'PANTHER') |> any() ) {
    
    assign(
      paste0('aDRG_PANTHER_Pathways'),
    value = read.table(
      'path/to/file/aDRG_PANTHER_Pathways.txt',
      sep = '\t',
      skip = 11,
      header = T,
      quote = '') |>
      dplyr::rename_with(
        ~gsub(
          .x,
          pattern = '\\.\\.',
          replacement = '')) %>% # remove the .'s  from the GO BP column
      dplyr::rename_with(
        ~gsub(
          .x,
          pattern = 'Mus.musculus.|upload_1|\\.$',
          replacement = '')) %>% # rm appended/excess column name strings
      dplyr::rename(Pathway.Size = REFLIST21836) |>
      dplyr::rename(Overlap = `364`) %>%
      dplyr::mutate(
        dplyr::across(c(2:4,6:8),
                      ~ as.numeric(.x))) |> # co-erce to numeric for operations
      tidyr::separate(
        PANTHER.Pathways,
        into = c("PANTHER.Pathway", "Panther.Pathway.ID"),
        sep = '\\s\\(P') |> # split to preserve names & IDs
      dplyr::mutate(
        Panther.Pathway.ID = paste0(
          'P',
          gsub(Panther.Pathway.ID,
               pattern = '\\)',
               replacement = ''))),
    envir = .GlobalEnv)
  }
}
```

### aDRG wrangling function


```r
aDRG_GO_and_PANTHER_wrangling <- function(category){

  if ( grepl(category, pattern = 'BP|Biological Process', ignore.case = T) |> any() ) {
    go_term = "GO.biological.process.complete"
  } else if ( grepl(category, pattern = 'CC|Cellular Component', ignore.case = T) |> any() ) {
    go_term = "GO.cellular.component.complete"
  } else if ( grepl(category, pattern = 'MF|Molecular Function', ignore.case = T) |> any() ) {
    go_term = "GO.molecular.function.complete"
  }
  
  if ( !grepl(category, pattern = 'PANTHER') |> any() ) {
    
    assign(
    x = paste0('aDRG_GO_', category, '_results'),
    
    value = get(x = paste0('aDRG_GO_', category, '_results'))%>%
      dplyr::rename_with(
        ~gsub(
          .x,
          pattern = '\\.\\.',
          replacement = '')) %>% # remove the .'s  from the GO BP column
      dplyr::rename_with(
        ~gsub(
          .x,
          pattern = 'Mus.musculus.|upload_1|\\.$',
          replacement = '')) |> # remove excess strings
      dplyr::rename(
        GO_Term = go_term) |> # rename the GO category column
      dplyr::rename(
        GO_Term_Size = REFLIST21836) |> # rename the GO Term Size column
      dplyr::rename(
        Overlap = 3) %>% # rename the overlap column
      dplyr::mutate(
        dplyr::across(c(2:4,7,8), ~as.numeric(.x))) |> # convert to numeric
      tidyr::separate(
        GO_Term,
        into = c("GO_Term", "GO_Term_ID"),
        sep = '\\s\\(GO:') |> # split the Term ID from the Term
      dplyr::mutate(
        GO_Term_ID = paste0('GO:',
                            gsub(GO_Term_ID,
                                 pattern = '\\)',
                                 replacement = ''))) |>
      dplyr::mutate(
        GO_Type = category, .after = GO_Term_ID) |> # create the category column
      dplyr::mutate(
        fold.Enrichment = stringr::str_extract(
          fold.Enrichment,
          pattern = '[:alnum:]+') |>
          as.numeric(fold.Enrichment)) |> # convert the FE column to numeric
      dplyr::arrange(
        FDR, desc(fold.Enrichment)
        ), # sort by the highest FE and lowest adjusted p-value
    
    envir = .GlobalEnv
    )
    
  }
}
```

### eDRG import function


```r
read_in_eDRG_GO_and_PANTHER_analysis_files <- function(age, category){
  if ( !grepl(category, pattern = 'PANTHER') |> any() ) {
    
    assign(
    paste0(age,
           '_DRG_GO_',
           category,
           '_results'),
    value = read.table(
    paste0('path/to/file',
           age,
           '_DRG_GO_analysis_',
           stringr::str_to_upper(category),
           '.txt'
           ),
    sep = '\t', skip = 11, header = T, quote = ''
    ),
    envir = .GlobalEnv)
    
  } else if (grepl(category, pattern = 'PANTHER') |> any() ) {
    
    assign(
    paste0(age,
           '_PANTHER_Pathways'),
    value = read.table(
    paste0('path/to/file',
           age,
           '_PANTHER_Pathways.txt'
           ),
    sep = '\t', skip = 11, header = T, quote = ''
    ),
    envir = .GlobalEnv)
    
  }
  
}
```

### eDRG wrangling function


```r
eDRG_GO_and_PANTHER_wrangling <- function(age, category){

  if ( grepl(category, pattern = 'BP|Biological Process', ignore.case = T) |> any() ) {
    go_term = "GO.biological.process.complete"
  } else if ( grepl(category, pattern = 'CC|Cellular Component', ignore.case = T) |> any() ) {
    go_term = "GO.cellular.component.complete"
  } else if ( grepl(category, pattern = 'MF|Molecular Function', ignore.case = T) |> any() ) {
    go_term = "GO.molecular.function.complete"
  }
  
  if ( !grepl(category, pattern = 'PANTHER') |> any() ) {
    
    assign(
    x = paste0(age, '_DRG_GO_', category, '_results'),
    
    value = get(x = paste0(age, '_DRG_GO_', category, '_results'))|>
      dplyr::rename_with(
        ~gsub(
          .x,
          pattern = '\\.\\.',
          replacement = '')) |> # remove the .'s  from the GO BP column
      dplyr::rename_with(
        ~gsub(
          .x,
          pattern = 'Mus.musculus.|upload_1|\\.$',
          replacement = '')) |> # remove excess strings
      dplyr::rename(
        GO_Term = go_term) |> # rename the GO category column
      dplyr::rename(
        GO_Term_Size = REFLIST21836) |> # rename the GO Term Size column
      dplyr::rename(
        Overlap = 3) |> # rename the overlap column
      dplyr::mutate(
        dplyr::across(c(2:4,7,8), ~as.numeric(.x))) |> # convert to numeric
      tidyr::separate(
        GO_Term,
        into = c("GO_Term", "GO_Term_ID"),
        sep = '\\s\\(GO:') |> # split the Term ID from the Term
      dplyr::mutate(
        GO_Term_ID = paste0('GO:',
                            gsub(GO_Term_ID,
                                 pattern = '\\)',
                                 replacement = ''))) |>
      dplyr::mutate(
        GO_Type = category, .after = GO_Term_ID) |> # create the category column
      dplyr::mutate(
        fold.Enrichment = stringr::str_extract(
          fold.Enrichment,
          pattern = '[:alnum:]+') |>
          as.numeric(fold.Enrichment)) |> # convert the FE column to numeric
      dplyr::arrange(
        FDR, desc(fold.Enrichment)
        ), # sort by the highest FE and lowest adjusted p-value
    
    envir = .GlobalEnv
    )
    
  } else {
    
    assign(
      x = paste0(age, '_PANTHER_Pathways'),
      value = get(x = paste0(age, '_PANTHER_Pathways')) |> 
      dplyr::mutate(devStage = age, .before = 1) |>
      dplyr::rename_with(
        ~gsub(
          .x,
          pattern = '\\.\\.',
          replacement = '')) %>% # remove the .'s  from the GO BP column
      dplyr::rename_with(
        ~gsub(
          .x,
          pattern = 'Mus.musculus.|upload_1|\\.$',
          replacement = '')) |> # rm appended/excess column name strings
      dplyr::rename(Pathway.Size = 4) |>
      dplyr::rename(Overlap = 3) %>%
      dplyr::mutate(
        dplyr::across(c(3:5,7:9), ~
                        as.numeric(.x))) |> # co-erce to numeric for operations
      tidyr::separate(PANTHER.Pathways,
                      into = c("PANTHER.Pathway", "Panther.Pathway.ID"),
                      sep = '\\s\\(P') |> # split to preserve names & IDs
      dplyr::mutate(Panther.Pathway.ID = paste0('P',
                                                gsub(Panther.Pathway.ID,
                                                     pattern = '\\)',
                                                     replacement = ''))), 
      envir = .GlobalEnv
    )
    
  }
}
```

# Import and Harmonize **aDRG** Data



+ Import the `aDRG GSEA results`

## Read in **aDRG** GSEA Data

+ `GO_{category}_results` = gene.ontology GSEA results


```r
aDRG_categories <- c('BP', 'CC', 'MF', 'PANTHER')

for (i in 1:length(aDRG_categories)) {
  read_in_aDRG_GO_and_PANTHER_analysis_files(category = aDRG_categories[i])
}
```

## Harmonize and Join **aDRG** GSEA Data

+ Harmonize and join the `aDRG GSEA results`


```r
for (i in 1:length(aDRG_categories)) {
  aDRG_GO_and_PANTHER_wrangling(category = aDRG_categories[i])
}

aDRG_GO_results <- aDRG_GO_BP_results |> 
  dplyr::full_join(aDRG_GO_CC_results) |> 
  dplyr::full_join(aDRG_GO_MF_results)

rm(list = ls()[which(grepl(ls(), pattern = 'aDRG_GO_[A-Z]') == T)])
```

# aDRG GSEA Results

+ Plot the `aDRG GSEA results`

## **aDRG** GO GSEA


```r
aDRG_GO_results |>
  dplyr::filter(!is.na(GO_Term)) |>
  dplyr::mutate(fold.Enrichment = as.numeric(fold.Enrichment)) |>
  dplyr::arrange(desc(fold.Enrichment), FDR) |>
  dplyr::filter(!grepl(GO_Term, pattern = "Unclassified", ignore.case = T)) |>
  dplyr::mutate(sig = dplyr::case_when(FDR <= 0.001 ~ "***",
                                dplyr::between(FDR, 0.001, 0.0499) ~ "**",
                                dplyr::between(FDR, 0.01, 0.05) ~ "*",
                                TRUE ~ "NS"), .after = FDR
  ) |>
  dplyr::filter(sig != 'NS') |>
  dplyr::slice(c(1:30)) |> 
  ggplot(aes(x = fold.Enrichment,
             y = forcats::fct_rev(forcats::fct_infreq(GO_Term, fold.Enrichment)),
             color = GO_Type,
             fill = GO_Type,
             size = Overlap,
             alpha = -log10(FDR))) +
  geom_point() +
  labs(title = expression(bold(paste("GO Analysis (GSEA) of ",
                                     italic("Tmem184b")^bold("GT/GT"),
                                     " DRG mRNA Exp. Changes")
                               )
                          ),
       # x = bquote(~ -log[10] ~ '(Adj.P-value)'),
       x = 'fold enrichment',
       y = 'GO Term',
       color = 'GO\nCategory',
       fill = 'GO\nCategory',
       size = '# Genes in List',
       # alpha = 'fold\nEnrichment'
       alpha = bquote(~ -log[10] ~ '(Adj.P-value)')
       ) +
  # geom_vline(xintercept = -log10(0.05), linetype = 'dashed', color = 'firebrick') +
  theme_bw() +
  theme(plot.title = element_text(face = 'bold', size = 13),
        axis.title = element_text(face = 'bold', size = 12),
        axis.text = element_text(face = 'bold', size = 12, vjust = 0.5),
        axis.ticks = element_line(linetype = 'solid'),
        legend.text = element_text(face = 'bold', size = 12),
        legend.title = element_text(face = 'bold', size = 12)) +
  guides(fill = guide_legend(override.aes = list(size = 5,
                                                 shape = 19,
                                                 fill = c('navy',
                                                          'darkgoldenrod3',
                                                          'firebrick'))),
         alpha = guide_legend(override.aes = list(size = 5, shape = 19)),
         size = guide_legend(override.aes = list(shape = 19))) +
  scale_color_manual(values = c('navy', 'darkgoldenrod3', 'firebrick')) +
  scale_fill_manual(values = c('navy', 'darkgoldenrod3', 'firebrick')) +
  geom_text(aes(x = fold.Enrichment+ 0.2,
                y = GO_Term,
                color = GO_Type,
                label = sig,
                fontface = 'bold'),
            size = 5.5,
            angle = 90,
            vjust = 1,
            alpha = 1)
```

![](ItchGSEA_files/figure-html/Plot aDRG GO GSEA results-1.png)<!-- -->

## **aDRG** PANTHER Pathways GSEA


```r
aDRG_PANTHER_Pathways |>
  dplyr::mutate(sig = dplyr::case_when(FDR <= 0.001 ~ "***",
                                dplyr::between(FDR, 0.001, 0.0499) ~ "**",
                                dplyr::between(FDR, 0.01, 0.05) ~ "*",
                                TRUE ~ "NS"), .after = FDR
  ) |>
  dplyr::filter(sig != 'NS') |>
  dplyr::arrange(desc(fold.Enrichment), FDR) |> 
  dplyr::slice(c(1:10)) |> 
  ggplot(aes(x = -log10(FDR),
             y = forcats::fct_rev(forcats::fct_infreq(PANTHER.Pathway, -log10(FDR))),
             color = -log10(FDR),
             fill = -log10(FDR))) +
  geom_bar(stat = 'identity') +
  labs(title = expression(bold(paste("PANTHER Pathway Analysis (GSEA) of ",
                                     italic("Tmem184b")^bold("GT/GT"),
                                     " aDRG mRNA Exp. Changes")
                               )
                          ),
       x = bquote(~ -log[10] ~ '(Adj.P-value)'),
       y = 'PANTHER Pathway',
       alpha = bquote(~ -log[10] ~ '(Adj.P-value)')
       ) +
  geom_vline(xintercept = -log10(0.05), linetype = 'dashed', color = 'firebrick') +
  theme_bw() +
  theme(plot.title = element_text(face = 'bold', size = 13),
        axis.title = element_text(face = 'bold', size = 12),
        axis.text = element_text(face = 'bold', size = 12, vjust = 0.5),
        axis.ticks = element_line(linetype = 'solid'),
        legend.position = 'none'
        ) +
  # guides(alpha = guide_legend(override.aes = list(size = 5, shape = 19)),
  #        size = guide_legend(override.aes = list(shape = 19))
  #        ) +
  geom_text(aes(label = 'FDR\n\u2264 0.05'),
            size = 5,
            color = 'firebrick',
            x = 3,
            y = 2,
            show.legend = FALSE) +
  geom_text(aes(x = -log10(FDR)+0.1,
                label = sig,
                fontface = 'bold'),
            size = 6.5,
            color = 'black',
            angle = 90,
            vjust = 1,
            alpha = 1)
```

![](ItchGSEA_files/figure-html/Plot aDRG PANTHER Pathways GSEA-1.png)<!-- -->


# Import and Harmonize **eDRG** Data




+ Import the `eDRG GSEA results`

## Read in **eDRG** GSEA Data

+ `({dev_stage_group})_GO_{category}_results` = gene.ontology results
  
+ Call the function to read in all `eDRG GSEA` datasets


```r
ages <- c('e13', 'p0', 'p10')
categories <- c('BP', 'CC', 'MF', 'PANTHER')

for (i in 1:length(ages)) {
  for (j in 1:length(categories)) {
   read_in_eDRG_GO_and_PANTHER_analysis_files(age = ages[i],
                                              category = categories[j]) 
  }
}
```

## Harmonize and Join **eDRG** GSEA Data

+ Harmonize the `eDRG GSEA results`


```r
for (i in 1:length(ages)) {
  for (j in 1:length(categories)) {
   eDRG_GO_and_PANTHER_wrangling(age = ages[i],
                                 category = categories[j]) 
  }
}
```

+ Join all of the dataframes by their developmental stage


```r
e13_DRG_GO_results <- e13_DRG_GO_BP_results |>
  dplyr::full_join(e13_DRG_GO_CC_results) |>
  dplyr::full_join(e13_DRG_GO_MF_results) |>
  dplyr::mutate(devStage = 'e13', .before = 1)

p0_DRG_GO_results <- p0_DRG_GO_BP_results |>
  dplyr::full_join(p0_DRG_GO_CC_results) |>
  dplyr::full_join(p0_DRG_GO_MF_results) |>
  dplyr::mutate(devStage = 'p0', .before = 1)

p10_DRG_GO_results <- p10_DRG_GO_BP_results |>
  dplyr::full_join(p10_DRG_GO_CC_results) |>
  dplyr::full_join(p10_DRG_GO_MF_results) |>
  dplyr::mutate(devStage = 'p10', .before = 1)

rm(list = ls()[which(grepl(ls(), pattern = '_DRG_GO_[A-Z]|aDRG_GO_[A-Z]') == T)])
```

# eDRG GSEA Results

+ Plot the `eDRG GSEA results`

## **eDRG** GO GSEA Results


```r
eDRG_GO_results <- e13_DRG_GO_results |>
  dplyr::full_join(p0_DRG_GO_results) |>
  dplyr::full_join(p10_DRG_GO_results)

eDRG_GO_results |>
  dplyr::filter(!is.na(GO_Term)) |>
  dplyr::mutate(fold.Enrichment = as.numeric(fold.Enrichment)) |>
  dplyr::arrange(desc(fold.Enrichment), FDR) |>
  dplyr::filter(!grepl(GO_Term, pattern = "Unclassified", ignore.case = T)) |>
  dplyr::mutate(sig = dplyr::case_when(FDR <= 0.001 ~ "***",
                                dplyr::between(FDR, 0.001, 0.0499) ~ "**",
                                dplyr::between(FDR, 0.01, 0.05) ~ "*",
                                TRUE ~ "NS"), .after = FDR
  ) |>
  dplyr::filter(sig != 'NS') |>
  dplyr::slice(c(1:30)) |> 
  dplyr::mutate(
    GO_Term = dplyr::case_when(
      grepl(GO_Term,
            pattern = 'adaptive immune response based on somatic') ~ 
        'adaptive immune response based on somatic recombination of\nimmune receptors built from immunoglobulin superfamily domains',
      TRUE ~ GO_Term)) |>
  ggplot(aes(x = fold.Enrichment,
             y = forcats::fct_rev(forcats::fct_infreq(GO_Term, fold.Enrichment)),
             color = GO_Type,
             fill = GO_Type,
             shape = devStage,
             size = Overlap,
             alpha = -log10(FDR))) +
  geom_point(aes(shape = devStage)) +
  labs(title = expression(
    bold(
      paste("GO Analysis (GSEA) of ",
            italic("Tmem184b")^bold("GT/GT"),
            " eDRG mRNA Exp. Changes"))),
    
       x = 'fold Enrichment',
       # x = bquote(~ -log[10] ~ '(Adj.P-value)'),
       y = 'GO Term',
       color = 'GO\nCategory',
       fill = 'GO\nCategory',
       shape = 'Dev Stage',
       size = 'Overlap',
       # alpha = 'fold\nEnrichment'
       alpha = '-log10(FDR)'
       ) +
  # geom_vline(xintercept = -log10(0.05), linetype = 'dashed', color = 'firebrick') +
  theme_bw() +
  theme(plot.title = element_text(face = 'bold', size = 13),
        axis.title = element_text(face = 'bold', size = 12),
        axis.text = element_text(face = 'bold', size = 12, vjust = 0.5),
        axis.ticks = element_line(linetype = 'solid'),
        legend.text = element_text(face = 'bold', size = 12),
        legend.title = element_text(face = 'bold', size = 12),
        strip.background = element_rect(color = 'black', fill = 'white'),
        strip.text = element_text(face = 'bold', size = 13)) +
  guides(shape = guide_legend(override.aes = list(size = 4, shape = c(19, 17))),
         fill = guide_legend(override.aes = list(size = 4, shape = 19, color = c('navy', 'darkgoldenrod3', 'firebrick'))),
         alpha = guide_legend(override.aes = list(size = 5, shape = 19)),
         size = guide_legend(override.aes = list(shape = 19))
         # size = 'none'
         ) +
  scale_color_manual(values = c('navy', 'darkgoldenrod3', 'firebrick')) +
  scale_fill_manual(values = c('navy', 'darkgoldenrod3', 'firebrick')) +
  # geom_text(aes(label = 'FDR\n\u2264 0.05'),
  #           size = 5,
  #           color = 'firebrick',
  #           x = 2,
  #           y = 4,
  #           show.legend = FALSE) +
  geom_text(aes(x = fold.Enrichment+1,
                label = sig,
                fontface = 'bold',
                color = GO_Type),
            size = 6.5,
            # color = 'black',
            angle = 90,
            vjust = 1,
            alpha = 1) +
  facet_grid(~devStage)
```

![](ItchGSEA_files/figure-html/Plot eDRG Results-1.png)<!-- -->

## **eDRG** PANTHER Pathway GSEA Results


```r
eDRG_PANTHER_Pathways <- e13_PANTHER_Pathways  |>
  dplyr::full_join(p0_PANTHER_Pathways) |>
  dplyr::full_join(p10_PANTHER_Pathways) |>
  dplyr::filter(!grepl(PANTHER.Pathway, pattern = 'Unclassified')) |>
  dplyr::arrange(
    desc(fold.Enrichment),
    FDR
    ) |> # sort by the highest FE and lowest adjusted p-value
  
  dplyr::mutate(sig = dplyr::case_when(FDR <= 0.001 ~ "***",
                                dplyr::between(FDR, 0.001, 0.0499) ~ "**",
                                dplyr::between(FDR, 0.01, 0.05) ~ "*",
                                TRUE ~ "NS"), .after = FDR
  )

eDRG_PANTHER_Pathways |>
  dplyr::filter(sig != 'NS') |>
  dplyr::arrange(desc(fold.Enrichment), FDR) |> 
  dplyr::slice(c(1:30)) |> 
  # dplyr::group_by(devStage, PANTHER.Pathway) |>
  ggplot(aes(x = fold.Enrichment,
             y = forcats::fct_rev(forcats::fct_infreq(PANTHER.Pathway, fold.Enrichment)),
             color = devStage,
             fill = devStage,
             size = Overlap,
             alpha = -log10(FDR))) +
  geom_point() +
  # geom_bar(stat = 'identity') +
  labs(title = expression(bold(paste("Pathway Analysis (GSEA) of ",
                                     italic("Tmem184b")^bold("GT/GT"),
                                     " eDRG mRNA Exp. Changes"))),
       x = 'fold Enrichment',
       # x = bquote(~ -log[10] ~ '(Adj.P-value)'),
       y = 'PANTHER Pathway',
       color = 'devStage',
       fill = 'devStage',
       size = 'Overlap',
       # alpha = 'fold\nEnrichment'
       alpha = '-log10(FDR)'
       ) +
  # geom_vline(xintercept = -log10(0.05), linetype = 'dashed', color = 'firebrick') +
  theme_bw() +
  theme(plot.title = element_text(face = 'bold', size = 13),
        axis.title = element_text(face = 'bold', size = 12),
        axis.text = element_text(face = 'bold', size = 12, vjust = 0.5),
        axis.ticks = element_line(linetype = 'solid'),
        legend.text = element_text(face = 'bold', size = 12),
        legend.title = element_text(face = 'bold', size = 12),
        strip.background = element_rect(color = 'black', fill = 'white'),
        strip.text = element_text(face = 'bold', size = 12)) +
  guides(fill = guide_legend(override.aes = list(size = 4,
                                                 shape = 19,
                                                 color = c('navy',
                                                           'darkgoldenrod3'))),
         alpha = guide_legend(override.aes = list(size = 5, shape = 19)),
         size = guide_legend(override.aes = list(shape = 19))
         ) +
  scale_color_manual(values = c('navy', 'darkgoldenrod3')) +
  scale_fill_manual(values = c('navy', 'darkgoldenrod3')) +
  # geom_text(aes(label = 'FDR\n\u2264 0.05'),
  #           size = 5,
  #           color = 'firebrick',
  #           x = 2,
  #           y = 4,
  #           show.legend = FALSE) +
  geom_text(aes(x = fold.Enrichment+0.1,
                label = sig,
                fontface = 'bold',
                color = devStage),
            size = 6.5,
            angle = 90,
            vjust = 1,
            alpha = 1) +
  facet_grid(~devStage)
```

![](ItchGSEA_files/figure-html/Plot eDRG PANTHER Results-1.png)<!-- -->

# Paper Fig 5B

The `GeneOntology.org` db has changed since 2018 when we did our original analyses, thus some terms are different in terms of their significance.

+ The terms in the paper have each been filtered directly

## **Old** Fig 5B


```r
eDRG_GO_results |>
  dplyr::filter(devStage == 'e13', GO_Type == 'BP') |>
  dplyr::arrange(FDR, desc(fold.Enrichment)) |> 
  dplyr::mutate(sig = dplyr::case_when(FDR <= 0.001 ~ "***",
                                dplyr::between(FDR, 0.001, 0.0499) ~ "**",
                                dplyr::between(FDR, 0.01, 0.05) ~ "*",
                                TRUE ~ "NS"), .after = FDR) |>
  dplyr::filter(
    
    GO_Term %in% c('axonogenesis',
                   'nervous system development',
                   'regulation of calcium-dependent activation of synaptic vesicle fusion', # re-categorized from extant 'reg. of calcium-dependent exocytosis'
                   'axon development', 'chromosome localization',
                   'neuron differentiation', 'ganglion development',
                   'G1/S transition of mitotic cell cycle', # re-categorized from G1/S transition
                   'epithelial cell-cell adhesion', # re-categorized from epithelial cell adhesion
                   'generation of neurons' # re-categorized from neuron generation
                   )
    ) |>
  
  dplyr::mutate(GO_Term = dplyr::case_when(

    grepl(GO_Term, pattern = 'calcium-dependent') ~ 'regulation of calcium-dependent\nactivation of synaptic vesicle fusion',
                                    TRUE ~ GO_Term)
    ) |> 
  ggplot2::ggplot(aes(x = forcats::fct_rev(forcats::fct_infreq(GO_Term, FDR)),
                      y = -log10(FDR),
                      color = -log10(FDR),
                      fill = -log10(FDR))) +
  ggplot2::geom_bar(stat = 'identity') +
  ggplot2::theme_bw()+
  ggplot2::theme(plot.title = element_text(face = 'bold', size = 13),
                 axis.title = element_text(face = 'bold', size = 12),
                 axis.text.y = element_text(face = 'bold', size = 12, vjust = 0.5),
                 axis.text.x = element_text(face = 'bold', size = 12, vjust = 0.5,
                                            angle = 90, hjust = 1),
                 axis.ticks = element_line(linetype = 'solid'),
                 legend.position = 'none',
                 # legend.text = element_text(face = 'bold', size = 12),
                 # legend.title = element_text(face = 'bold', size = 12),
                 panel.grid = element_blank()) +
  ggplot2::labs(title = expression(
    bold(
      paste("GO Analysis (GSEA) of ",
            italic("Tmem184b")^bold("GT/GT"), " e13 DRG mRNA Exp. Changes")
      )
    ),
       x = 'GO Biological Process',
       y = bquote(~ -log[10] ~ '(Adj.P-value)')
       ) +
  ggplot2::geom_hline(yintercept = -log10(0.05), color = 'firebrick', linetype = 'dashed') +
  ggplot2::geom_text(aes(x = GO_Term,
                         y = -log10(FDR)+0.75,
                         label = sig,
                fontface = 'bold'),
                size = 6.5,
                color = 'black',
                vjust = 0.5,
                hjust = 0.5,
                alpha = 1) +
  ggplot2::geom_text(aes(x = 9.5,
                         y = 5,
                         label = 'FDR \u2264 0.05'),
                     fontface = 'bold',
                     size = 4.5,
                     color = 'firebrick')
```

![](ItchGSEA_files/figure-html/Retro Paper Pathway fig-1.png)<!-- -->

## **New** Fig 5B


```r
eDRG_GO_results |>
  dplyr::filter(devStage == 'e13', GO_Type == 'BP') |>
  dplyr::arrange(FDR, desc(fold.Enrichment)) |> 
  dplyr::mutate(sig = dplyr::case_when(FDR <= 0.001 ~ "***",
                                dplyr::between(FDR, 0.001, 0.0499) ~ "**",
                                dplyr::between(FDR, 0.01, 0.05) ~ "*",
                                TRUE ~ "NS"), .after = FDR) |>
  dplyr::slice(c(1:15)) |> 
  ggplot2::ggplot(aes(x = forcats::fct_rev(forcats::fct_infreq(GO_Term, FDR)), y = -log10(FDR), color = -log10(FDR), fill = -log10(FDR))) +
  ggplot2::geom_bar(stat = 'identity') +
  ggplot2::theme_bw()+
  ggplot2::theme(plot.title = element_text(face = 'bold', size = 13),
                 axis.title = element_text(face = 'bold', size = 12),
                 axis.text.y = element_text(face = 'bold', size = 12, vjust = 0.5),
                 axis.text.x = element_text(face = 'bold', size = 12, vjust = 0.5, angle = 90, hjust = 1),
                 axis.ticks = element_line(linetype = 'solid'),
                 legend.text = element_text(face = 'bold', size = 12),
                 legend.title = element_text(face = 'bold', size = 12),
                 panel.grid = element_blank()) +
  ggplot2::labs(title = expression(bold(paste("GO Analysis (GSEA) of ", italic("Tmem184b")^bold("GT/GT"), " e13 DRG mRNA Exp. Changes"))),
                x = 'GO Biological Process',
                y = bquote(~ -log[10] ~ '(Adj.P-value)')
                ) +
  ggplot2::geom_hline(yintercept = -log10(0.05), color = 'firebrick', linetype = 'dashed') +
  ggplot2::geom_text(aes(x = GO_Term,
                         y = -log10(FDR)+0.75,
                         label = sig,
                         fontface = 'bold'),
                     size = 6.5,
                     color = 'black',
                     vjust = 0.5,
                     hjust = 0.5,
                     alpha = 1)
```

![](ItchGSEA_files/figure-html/Current Paper Pathway fig-1.png)<!-- -->


