Itch DEA
================
Erik Larsen
2025-02-10

# Overview

This markdown was developed to formally modernize the differential
expression data analysis of bulk RNA-seq data for **Dr. Martha
Bhattacharya’s lab** at the **University of Arizona** for the [Itch
paper](https://journals.lww.com/pain/abstract/2022/05000/transmembrane_protein_tmem184b_is_necessary_for.18.aspx)

## Set up Environment

#### **Load packages**

``` r
packages <- c('BiocManager',
              'tidyverse',
              'stringr',
              'GO.db',
              'org.Mm.eg.db',
              'forcats',
              'ggplot2',
              'ggrepel',
              'ggbiplot',
              'pheatmap')
for ( package in packages ) {
  if (!require(package, character.only = T)) {
    install.packages(package)
    library(package, character.only = T)
  }
}
```

# Read in and Prep aDRG Data

#### **Read in aDRG Data**

Not shown; variable names:
+ `aDRG` = adult mouse DRG DEA results
  + `GO_{category}_results` = gene.ontology gene set enrichment analysis
results

# Plot aDRG Results

## **DEA**

### **Volcano**

``` r
  ## create the handy notin pipe operator
"%notin%" = Negate("%in%")

  ## list the differentially expressed genes known to be involved in itch
itch_DEGs = c("Il31ra", "Cysltr2", "Npy2r", "Sst",
              "Htr1a", "P2rx3", "Lpar3", "Lpar5",
              "Scn11a", "Scn10a", "Mrgprd", "Trpc6",
              "Trpc3", "F2rl2", "Htr1f", "Osmr",
              "Fxyd2", "Htr4", "Mrgprx1", "Ptgdr",
              "Trpa1", "Trpm6", "Hrh1", "Mrgpra3", "Tmem184b")

  ## Create grouping and labeling variables
aDRG <- aDRG |> 
  dplyr::mutate(goi = case_when(-log10(AdjP) >= -log10(0.05) &
                                  GeneID %in% itch_DEGs ~ 'Itch-related DEG',
                                -log10(AdjP) >= -log10(0.05) &
                                  GeneID %notin% itch_DEGs ~ 'DEG',
                                TRUE ~ 'Non-DEG'),
                labs = case_when(GeneID %in% c('Tmem184b',
                                               'Il31ra',
                                               'Cysltr2',
                                               'Npy2r',
                                               'Sst',
                                               'Htr1a',
                                               'Htr1f',
                                               'Lpar5',
                                               'Lpar3',
                                               'Mrgprd',
                                               'Mgrpra3',
                                               'Trpa1') ~ GeneID,
                                 TRUE ~ ''))

  ## volcano plot
ggplot(data = aDRG) +
  geom_point(data = subset(aDRG, goi == 'Non-DEG'),
             aes(x = log2.FC.,
                 y = -log10(AdjP),
                 color = goi),
             size = 3,
             alpha = 0.1) +
  geom_point(data = subset(aDRG, goi == 'DEG'),
             aes(x = log2.FC.,
                 y = -log10(AdjP),
                 color = goi),
             size = 4,
             alpha = 0.5) +
  geom_point(data = subset(aDRG, goi == 'Itch-related DEG'),
             aes(x = log2.FC.,
                 y = -log10(AdjP),
                 color = goi),
             size = 5,
             alpha = 1) +
  coord_cartesian(ylim = c(0,53), xlim = c(-3.2,2)) +
  theme_classic() +
  theme(plot.title = element_text(face = 'bold', size = 15),
        axis.title = element_text(face = 'bold', size = 12),
        axis.text = element_text(face = 'bold', size = 12, vjust = 0.5),
        axis.ticks = element_line(linetype = 'solid'),
        legend.position = 'none'
        ) +
  labs(title = expression(bold(paste("DEA of ",
                                     italic("Tmem184b")^bold("GT/GT"),
                                     " DRG mRNA Exp. Changes"))),
       x = bquote(~ log[2] ~ '(fold change mRNA Relative to WT)'),
       y = bquote(~ -log[10] ~ '(Adj.P-value)'),
       col = 'Gene\nData\nType'
       ) +
  geom_hline(yintercept = -log10(0.05),
             linetype = 'dashed',
             color = 'firebrick') +
  scale_color_manual(values = c('darkgoldenrod4', 'navy', 'gray48'))+
  geom_text(aes(x = -1.9, 
                y = 2.5,
                label = 'FDR \U2264 0.05'),
            color = 'firebrick',
            size = 5) + 
  geom_text_repel(data = aDRG,
                  aes(label = aDRG$labs),
                  x = aDRG$log2.FC.,
                  y = -log10(aDRG$AdjP),
                  color = 'black',
                  max.overlaps = Inf,
                  size = 5)
```

![](ItchDEAgit_files/figure-gfm/Plot%20aDRG%20DEA%20results-1.png)<!-- -->

## **Trx Profiling**

### **Z-score the Normalized Counts**

#### **Z-score Function**

``` r
  ## new dataframe should not contain any NAs in p-value columns
aDRG_filtered <- aDRG |>
  dplyr::filter(!is.na(AdjP)) |> # remove any NAs in p-value columns
  dplyr::filter(!grepl(GeneID, # remove rRNAs, mitochondrial tRNAs, pseudogenes 
                       pattern = 'Rps.+.?$|RP.+.?$|Rpl.+.?$|MRPL.+.?$|Mrpl.+.?$|MRPS.+.?$|Mrps.+.?$|.*Rik.+$|.*Rik$|Gm.+.?$|^[A-Z]+[A-Z].+.?$|^[0-9]+.+.?$|mt.+.?$'))

  ## Rename columns
colnames(aTPM) = c("GeneID",
                   "WT1", "WT2", "WT3", "WT4",
                   "Mut1", "Mut2", "Mut3", "Mut4")

  ## Pass profile through Zscore function
Find_Row_Z(Expression_Profile = aTPM[,c(1:9)])
```

### **PCA**

``` r
  ## Create a covariance matrix of the centered and normalized data
Z_cov_mat = cov(Z[,c(2:9)])
  ## Pass the covariance matrix into R's PCA function to obtain the PCA object

PCA = prcomp(Z_cov_mat, center = FALSE, scale. = FALSE)
summary(PCA)
```

    ## Importance of components:
    ##                           PC1    PC2    PC3    PC4    PC5    PC6    PC7
    ## Standard deviation     0.4671 0.4284 0.3736 0.3606 0.3491 0.3423 0.3222
    ## Proportion of Variance 0.2151 0.1810 0.1376 0.1282 0.1202 0.1155 0.1024
    ## Cumulative Proportion  0.2151 0.3961 0.5337 0.6619 0.7821 0.8976 1.0000
    ##                              PC8
    ## Standard deviation     2.442e-17
    ## Proportion of Variance 0.000e+00
    ## Cumulative Proportion  1.000e+00

``` r
ggbiplot(PCA, 
         ellipse = FALSE,
         obs.scale = 1,
         var.scale = 1,
         labels = colnames(Z[c(2:9)]),
         groups = c(colnames(Z[c(2:9)])[c(1:4)],
                    colnames(Z[c(2:9)])[c(5:8)]),
         labels.size = 9,
         # point.size = 8,
         # varname.size = 20,
         var.axes = FALSE) +
  labs(title = expression(paste(italic("Tmem184b")^"GT/GT", " aDRG mRNA Exp. Changes PCA"))) +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5, size = 20),
        legend.position = 'none',
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 14),
        legend.title = element_blank()) +
  scale_color_manual(values = c(rep("darkgoldenrod4", 4), rep("gray25", 4)))
```

![](ItchDEAgit_files/figure-gfm/aDRG%20PCA-1.png)<!-- -->

### **Hierarchical Clustering the Z-scores**

``` r
  ## Don't remove non-DEGs
Z <- Z |>
  dplyr::filter(GeneID %in% c(aDRG |> 
                                dplyr::select(GeneID) |>
                                unlist() |>
                                as.character()
                              )
                )

  ## Create a list of gene names ordered by euclidean distance by which to re-order the dataframe
Euclid_dist_order = hclust(dist(Z[,c(2:9)], method = "euclidean"))$order
  ## The names (not the numbers)
Euclid_dist_ord_Genes = c(Z$GeneID[Euclid_dist_order])

  ## Transform again to order by clusters (Euclidean distances)
Za = Z |>
  dplyr::mutate(GeneID = factor(GeneID, levels = Euclid_dist_ord_Genes)) |>
  dplyr::arrange(GeneID)

  ## Find where the "Itch-related DEGs" are in the clustered matrix subsetted to DEGs for the heatmap
itch_DEGs = c("Il31ra", "Cysltr2", "Npy2r", "Sst",
              "Htr1a", "P2rx3", "Lpar3", "Lpar5",
              "Scn11a", "Scn10a", "Mrgprd", "Trpc6",
              "Trpc3", "F2rl2", "Htr1f", "Osmr",
              "Fxyd2", "Htr4", "Mrgprx1", "Ptgdr",
              "Trpa1", "Trpm6", "Hrh1", "Mrgpra3", "Tmem184b")

temp = vector()
for( i in 1:length(itch_DEGs)){
  temp[i] = which(Euclid_dist_ord_Genes == itch_DEGs[i])
}
itch_index = temp
#itch_index

  ## Confirm by indexing; copy and paste to customize the clustered heatmap
#Euclid_dist_ord_Genes[c(itch_index)]
  ## Confirm the location/order of the transformed matrix/dataframe
  ## Should be "Il31ra"
Za[itch_index[1],]
```

    ##      GeneID       WT1       WT2      WT3      WT4       Mut1       Mut2
    ## 7285 Il31ra 0.8348239 0.4522143 1.156546 1.199576 -0.9437066 -0.8612782
    ##            Mut3      Mut4
    ## 7285 -0.9282107 -0.909964

### **Plot Heatmap**

- View full heatmap

``` r
 ## Visualize the Z-scored, Euclidean-clustered and ordered gene TPMs across replicates with "pheatmap"
pheatmap(mat = Za[,2:9],
         color = colorRampPalette(c("navy", "white", "darkgoldenrod4"))((50)),
         clustsering_distance_rows = "euclidean",
         angle_col = 0,
         cutree_rows = 2,
         treeheight_row = 35,
         treeheight_col = 7,
         show_rownames = F)
```

![](ItchDEAgit_files/figure-gfm/Plot%20aDRG%20heatmap-1.png)<!-- -->

- Zoom to cluster near Tmem184b

``` r
  ## Visualize around the Tmem184b cluster
pheatmap(mat = Za[1135:1175,2:9],
         color = colorRampPalette(c("navy", "white", "darkgoldenrod4"))((50)),
         clustsering_distance_rows = "euclidean",
         angle_col = 0,
         treeheight_row = 35,
         treeheight_col = 9,
         labels_row = Euclid_dist_ord_Genes[1135:1175])
```

![](ItchDEAgit_files/figure-gfm/Plot%20aDRG%20Tmem%20Zoom%20heatmap-1.png)<!-- -->

- Zoom to itch-related genes

``` r
  ## Visualize the "Itch-related DEGs"
pheatmap(mat = Za[c(itch_index),2:9],
         color = colorRampPalette(c("navy", "white", "darkgoldenrod4"))((50)),
         clustsering_distance_rows = "euclidean",
         angle_col = 0,
         treeheight_row = 35,
         treeheight_col = 9,
         cutree_rows = 4,
         labels_row = Euclid_dist_ord_Genes[c(itch_index)])
```

![](ItchDEAgit_files/figure-gfm/Plot%20itch%20DEG%20heatmap-1.png)<!-- -->

# Read in and Prep eDRG Data

#### **Read in eDRG Data**

Not shown; variable names: + `({dev_stage_group})_DRG` = DRG DEA
results + `({dev_stage_group})_GO_{category}_results` = gene.ontology
results

``` r
neuron_differentiation_DEGs <- GO_INFO_by_TERM_df |>
  dplyr::filter(GO_Term == 'neuron differentiation') |>
  dplyr::select(Gene_IDs_from_List) |>
  unlist() |>
  as.character() %>%
  str_split(., pattern = ';') |>
  unlist()

eDRG <- e13_DRG |>
  dplyr::mutate(devStage = 'e13', .before = 1) |> 
  dplyr::full_join(p0_DRG |>
                     dplyr::mutate(devStage = 'p0', .before = 1) |>
                     dplyr::rename(Base.Mean = Base.mean) |>
                     dplyr::rename(Pval = PVal)) |>
  dplyr::full_join(p10_DRG |>
                     dplyr::mutate(devStage = 'p10', .before = 1) |>
                     dplyr::rename(Base.Mean = Base.mean) |>
                     dplyr::rename(Pval = PVal)) |>
  dplyr::filter(!is.na(AdjP)) |> # remove any NAs in p-value columns
  dplyr::filter(!grepl(GeneID, # remove rRNAs, mitochondrial tRNAs, pseudogenes
                       pattern = 'Rps.+.?$|RP.+.?$|Rpl.+.?$|MRPL.+.?$|Mrpl.+.?$|MRPS.+.?$|Mrps.+.?$|.*Rik.+$|.*Rik$|Gm.+.?$|^[A-Z]+[A-Z].+.?$|^[0-9]+.+.?$|mt.+.?$')
                ) |>
  dplyr::mutate(goi = case_when(AdjP <= 0.05 &
                                  grepl(GeneID,
                                        pattern = 'Wnt.+$|Pou.+$|Fzd.+$|Neuro.+$|Mapk.+$')
                                ~ "Neuron Differentiation DEG",
                                AdjP <= 0.05 &
                                  !grepl(GeneID,
                                         pattern = 'Wnt.+$|Pou.+$|Fzd.+$|Neuro.+$|Mapk.+$')
                                ~ "DEG",
                                TRUE ~ "Non-DEG"))
```

# Plot Embryonic DRG DEA Results

#### **Volcano**

``` r
ggplot(data = eDRG) +
  geom_point(data = subset(eDRG, goi == 'Non-DEG'),
             aes(x = log2.FC., y = -log10(AdjP), color = goi), alpha = 0.3) +
  geom_point(data = subset(eDRG, goi == 'DEG'),
             aes(x = log2.FC., y = -log10(AdjP), color = goi), alpha = 0.5) +
  geom_point(data = subset(eDRG, goi == 'Neuron Differentiation DEG'),
             aes(x = log2.FC., y = -log10(AdjP), color = goi), alpha = 0.5) +
  labs(title = expression(
    bold(paste("Differential Expression Analysis (DEA) of ",
               italic("Tmem184b")^bold("GT/GT"),
               " eDRG mRNA Exp. Changes"))),
       x = bquote(~ log[2] ~ 'fold change (relative to WT mRNA expression)'),
       y = bquote(~ -log[10] ~ '(Adj.P-value)')
       ) +
  theme_bw() +
  theme(plot.title = element_text(face = 'bold', size = 13),
        axis.title = element_text(face = 'bold', size = 12),
        axis.text = element_text(face = 'bold', size = 12, vjust = 0.5),
        axis.ticks = element_line(linetype = 'solid'),
        legend.text = element_text(face = 'bold', size = 12),
        legend.title = element_text(face = 'bold', size = 12),
        strip.background = element_rect(color = 'black', fill = 'white'),
        strip.text = element_text(face = 'bold', size = 13)) +
  geom_hline(yintercept = -log10(0.05), color = 'firebrick', linetype = 'dashed') +
  scale_color_manual(values = c('darkgoldenrod4', 'navy', 'gray48')) +
  geom_text(aes(label = 'FDR\n\u2264 0.05'),
            x = -5,
            y = -log10(0.05),
            size = 5,
            color = 'firebrick',
            show.legend = FALSE) +
  facet_grid(~devStage)
```

![](ItchDEAgit_files/figure-gfm/Plot%20eDRG%20Results-1.png)<!-- -->

#### **MA**

``` r
ggplot(data = eDRG) +
  geom_point(data = subset(eDRG, goi == 'Non-DEG'),
             aes(x = log(Base.Mean), y = log2.FC., color = goi), alpha = 0.3) +
  geom_point(data = subset(eDRG, goi == 'DEG'),
             aes(x = log(Base.Mean), y = log2.FC., color = goi), alpha = 0.5) +
  geom_point(data = subset(eDRG, goi == 'Neuron Differentiation DEG'),
             aes(x = log(Base.Mean), y = log2.FC., color = goi), alpha = 1) +
  labs(title = expression(
    bold(paste("Differential Expression Analysis (DEA) of ",
               italic("Tmem184b")^bold("GT/GT"),
               " eDRG mRNA Exp. Changes"))),
       x = bquote(~ ln ~ '(Mean of TPM)'),
       y = bquote(~ -log[2] ~ 'fold change (relative to WT mRNA expression)')
       ) +
  theme_bw() +
  theme(plot.title = element_text(face = 'bold', size = 13),
        axis.title = element_text(face = 'bold', size = 12),
        axis.text = element_text(face = 'bold', size = 12, vjust = 0.5),
        axis.ticks = element_line(linetype = 'solid'),
        legend.text = element_text(face = 'bold', size = 12),
        legend.title = element_text(face = 'bold', size = 12),
        strip.background = element_rect(color = 'black', fill = 'white'),
        strip.text = element_text(face = 'bold', size = 13)) +
  geom_hline(yintercept = 0, color = 'black', linetype = 'solid') +
  scale_color_manual(values = c('darkgoldenrod4', 'navy', 'gray48')) +
  scale_fill_manual(values = c('darkgoldenrod4', 'navy', 'gray48')) +
  facet_grid(~devStage)
```

![](ItchDEAgit_files/figure-gfm/Plot%20eDRG%20MA%20plots-1.png)<!-- -->

#### **Plot Embryonic DRG Heatmap (Paper Fig 5C)**

Generate the Z-scored normalized counts matrix

``` r
Find_Row_Z(Expression_Profile = TPM |> dplyr::select(c(GeneID, dplyr::starts_with('e13'))))
```

``` r
  ## Create a list of gene names ordered by euclidean distance by which to re-order the dataframe
Euclid_dist_order = hclust(dist(Z[,c(2:7)], method = "euclidean"))$order
  ## The names (not the numbers)
Euclid_dist_ord_Genes = c(Z$GeneID[Euclid_dist_order])

  ## Transform again to order by clusters (Euclidean distances)
Za = Z |>
  dplyr::mutate(GeneID = factor(GeneID, levels = Euclid_dist_ord_Genes)) |>
  dplyr::arrange(GeneID)


temp = vector()
for( i in 1:length(neuron_differentiation_DEGs)){
  temp[i] = which(Euclid_dist_ord_Genes == neuron_differentiation_DEGs[i])
}
diff_index = temp
#diff_index

pheatmap(mat = Za[c(diff_index),2:7],
         color = colorRampPalette(c("navy", "white", "darkgoldenrod4"))((50)),
         clustsering_distance_rows = "euclidean",
         angle_col = 0,
         treeheight_row = 35,
         treeheight_col = 9,
         cutree_rows = 4,
         labels_row = Euclid_dist_ord_Genes[c(diff_index)])
```

![](ItchDEAgit_files/figure-gfm/Plot%20Paper%20eDRG%20Heatmap-1.png)<!-- -->

#### **Paper Fig 5B**

The GeneOntology.org db has changed since 2018 when we did our original
analyses, thus some terms are different in terms of their significance.

- The terms in the paper have each been filtered directly

![](ItchDEAgit_files/figure-gfm/Retro%20Paper%20Pathway%20fig-1.png)<!-- -->

Here is the current analysis of the same data

``` r
eDRG_GO_results |>
  dplyr::filter(devStage == 'e13', GO_Type == 'BP') |>
  dplyr::arrange(FDR, desc(fold.Enrichment)) |> 
  dplyr::mutate(sig = case_when(FDR <= 0.001 ~ "***",
                                between(FDR, 0.001, 0.0499) ~ "**",
                                between(FDR, 0.01, 0.05) ~ "*",
                                TRUE ~ "NS"), .after = FDR) |>
  dplyr::slice(c(1:15)) |> 
  ggplot(aes(x = forcats::fct_rev(forcats::fct_infreq(GO_Term, FDR)),
             y = -log10(FDR),
             color = -log10(FDR),
             fill = -log10(FDR))) +
  geom_bar(stat = 'identity') +
  theme_bw()+
  theme(plot.title = element_text(face = 'bold', size = 13),
        axis.title = element_text(face = 'bold', size = 12),
        axis.text.y = element_text(face = 'bold', size = 12, vjust = 0.5),
        axis.text.x = element_text(face = 'bold',
                                   size = 12,
                                   vjust = 0.5,
                                   angle = 90,
                                   hjust = 1),
        axis.ticks = element_line(linetype = 'solid'),
        legend.text = element_text(face = 'bold', size = 12),
        legend.title = element_text(face = 'bold', size = 12),
        panel.grid = element_blank()) +
  labs(title = expression(bold(paste("GO Analysis (GSEA) of ",
                                     italic("Tmem184b")^bold("GT/GT"),
                                     " e13 DRG mRNA Exp. Changes"))),
       x = 'GO Biological Process',
       y = bquote(~ -log[10] ~ '(Adj.P-value)')
       ) +
  geom_text(aes(x = GO_Term,
                y = -log10(FDR)+0.75,
                label = sig,
                fontface = 'bold'),
            size = 6.5,
            color = 'black',
            vjust = 0.5,
            hjust = 0.5,
            alpha = 1)
```

![](ItchDEAgit_files/figure-gfm/Current%20Paper%20Pathway%20fig-1.png)<!-- -->
