Gasperini (2019) data documentation
================
Tim Barry
2022-07-04

``` r
library(tidyverse)
library(ondisc)
```

# Overview

This repository contains code to import and process the [Gasperini
2019](https://pubmed.ncbi.nlm.nih.gov/30612741/) data. Gasperini et
al. developed a high MOI, single-cell CRISPR screen assay to map
putative enhancers at genome-wide scale in a population of \>200,000
K562 cells.

The `gasperini-2019-v2` directory structure is as follows:

    └── at-scale
        ├── intermediate
        ├── processed
        │   ├── gene
        │   │   ├── matrix.odm
        │   │   └── metadata.rds
        │   ├── grna_assignment
        │   │   ├── matrix.odm
        │   │   └── metadata.rds
        │   └── grna_expression
        │       ├── matrix.odm
        │       └── metadata.rds
        └── raw

The contents of the `raw` and `intermediate` directories are suppressed,
as they are unimportant. We set file paths to the `gene`,
`grna_assignment`, and `grna_expression` modalities of the `processed`
directory below.

``` r
gasp_offsite <- paste0(.get_config_path("LOCAL_GASPERINI_2019_V2_DATA_DIR"), "at-scale/")
processed_data_dir <- paste0(gasp_offsite, "processed/")
processed_gene_dir <- paste0(processed_data_dir, "gene/")
processed_gRNA_expression_dir <- paste0(processed_data_dir, "grna_expression/")
processed_gRNA_assignment_dir <- paste0(processed_data_dir, "grna_assignment/")
```

# Gene

``` r
gene_odm <- read_odm(paste0(processed_gene_dir, "matrix.odm"),
                     paste0(processed_gene_dir, "metadata.rds"))
gene_odm
```

    ## A covariate_ondisc_matrix with the following components:
    ##  An integer-valued ondisc_matrix with 13135 features and 207324 cells.
    ##  A cell covariate matrix with columns n_nonzero, n_umis, p_mito, batch.
    ##  A feature covariate matrix with columns mean_expression, coef_of_variation, n_nonzero.

The gene data consist of 13,135 genes measured across 207,324 cells. The
cell covariates include `n_nonzero` (number of genes expressed in cell),
`n_umis` (cell sequencing depth or library size), `p_mito`, and `batch`.
There are two batches: `prep_batch_1` and `prep_batch_2`.

``` r
gene_odm |> get_cell_covariates() |> head()
```

    ##                                  n_nonzero n_umis      p_mito        batch
    ## AAACCTGAGAGGTACC-1_1A_1_SI-GA-E2      3549  17566 0.058786706 prep_batch_1
    ## AAACCTGAGTCAATAG-1_1A_1_SI-GA-E2      2543   8917 0.036086518 prep_batch_1
    ## AAACCTGCAAACAACA-1_1A_1_SI-GA-E2      3191  14626 0.069823051 prep_batch_1
    ## AAACCTGCACTTCTGC-1_1A_1_SI-GA-E2      4539  22783 0.026186508 prep_batch_1
    ## AAACCTGCATGTAGTC-1_1A_1_SI-GA-E2      2605  10124 0.007991318 prep_batch_1
    ## AAACCTGGTAGCGCAA-1_1A_1_SI-GA-E2      2187   9743 0.022356681 prep_batch_1

# gRNA (expression)

The gRNA (expression) modality contains data on 13,189 gRNAs and 207,324
cells.

``` r
gRNA_expression_odm <- read_odm(paste0(processed_gRNA_expression_dir, "matrix.odm"),
                                paste0(processed_gRNA_expression_dir, "metadata.rds"))
gRNA_expression_odm
```

    ## A covariate_ondisc_matrix with the following components:
    ##  An integer-valued ondisc_matrix with 13189 features and 207324 cells.
    ##  A cell covariate matrix with columns n_nonzero, n_umis.
    ##  A feature covariate matrix with columns mean_expression, n_nonzero, gRNA_group, target_type, target, target_gene.

The cell covariate matrix contains columns `n_nonzero` and `n_umis`. The
feature covariate matrix, meanwhile, is a bit more complicated. Let us
take a look.

``` r
gRNA_expression_odm |> get_feature_covariates() |> head()
```

    ##                      mean_expression n_nonzero   gRNA_group target_type
    ## AAACCGCTCCCGAGCACGGG      0.08524339      1450 SH3BGRL3_TSS    gene_tss
    ## AAATAGTGGGAAGATTCGTG      0.02863151       556 MTRNR2L8_TSS    gene_tss
    ## AACACACCACGGAGGAGTGG      0.06421350       987   FAM83A_TSS    gene_tss
    ## AACAGCCCGGCCGGCCAAGG      0.07196465      1192   ZNF593_TSS    gene_tss
    ## AACGAGAGACTGCTTGCTGG      0.03283749       693   ATPIF1_TSS    gene_tss
    ## AACGGCTCGGAAGCCTAGGG      0.07587158      1251    TIPRL_TSS    gene_tss
    ##                                        target     target_gene
    ## AAACCGCTCCCGAGCACGGG   chr1:26605667-26605668 ENSG00000142669
    ## AAATAGTGGGAAGATTCGTG  chr11:10530735-10530736 ENSG00000255823
    ## AACACACCACGGAGGAGTGG chr8:124191287-124191288 ENSG00000147689
    ## AACAGCCCGGCCGGCCAAGG   chr1:26496362-26496363 ENSG00000142684
    ## AACGAGAGACTGCTTGCTGG   chr1:28562620-28562621 ENSG00000130770
    ## AACGGCTCGGAAGCCTAGGG chr1:168148171-168148172 ENSG00000143155

The column `target_type` indicates the target type of a given gRNA, one
of `candidate_enhancer`, `non-targeting`, `known_enhancer`, and
`gene_tss`.

``` r
gRNA_expression_odm |>
  get_feature_covariates() |>
  pull(target_type) |> table()
```

    ## 
    ## candidate_enhancer      non-targeting     known_enhancer           gene_tss 
    ##              12200                101                 14                762

`candidate_enhancer` indicates a region of the genome that is a putative
enhancer; `non-targeting` indicates a non-targeting gRNA;
`known_enhancer` indicates a region of the genome that is known to be an
enhancer (based on previous research); and `gene_tss` is the
transcription start site of a gene.

Next, the `target` column indicates the chromosomal region that a given
gRNA targets:

``` r
target_str <- gRNA_expression_odm |>
  get_feature_covariates() |>
  pull(target)
```

Most chromosomal locations are targeted by two or four gRNAs.
Non-targeting gRNAs have the string “non-targeting” in the `target`
column:

``` r
gRNA_expression_odm |>
  get_feature_covariates() |>
  filter(target_type == "non-targeting") |> head()
```

    ##                      mean_expression n_nonzero   gRNA_group   target_type
    ## AACACAACACACCAAAACTG      0.06001717      1064    random_24 non-targeting
    ## AAGTTGACTCTACATAGCAG      0.03666725       728     random_7 non-targeting
    ## AATATTCTCCCTCATTCTGG      0.04364666       818    random_13 non-targeting
    ## AATCATGGTGGAAGGTGAAG      0.03595821       524    random_19 non-targeting
    ## AATCCTCTAATGGACGAAGA      0.03273138       702    random_10 non-targeting
    ## AATGAGGAGCAAACGAAAAT      0.03384558       665 scrambled_20 non-targeting
    ##                             target target_gene
    ## AACACAACACACCAAAACTG non-targeting        <NA>
    ## AAGTTGACTCTACATAGCAG non-targeting        <NA>
    ## AATATTCTCCCTCATTCTGG non-targeting        <NA>
    ## AATCATGGTGGAAGGTGAAG non-targeting        <NA>
    ## AATCCTCTAATGGACGAAGA non-targeting        <NA>
    ## AATGAGGAGCAAACGAAAAT non-targeting        <NA>

Next, `target_gene` indicates the gene targeted by TSS-targeting (i.e.,
type `gene_tss`) gRNAs. gRNAs that do not target gene TSSs have an
`<NA>` in this column.

``` r
# tss-targeting gRNAs
gRNA_expression_odm |>
  get_feature_covariates() |>
  filter(target_type == "gene_tss") |> head()
```

    ##                      mean_expression n_nonzero   gRNA_group target_type
    ## AAACCGCTCCCGAGCACGGG      0.08524339      1450 SH3BGRL3_TSS    gene_tss
    ## AAATAGTGGGAAGATTCGTG      0.02863151       556 MTRNR2L8_TSS    gene_tss
    ## AACACACCACGGAGGAGTGG      0.06421350       987   FAM83A_TSS    gene_tss
    ## AACAGCCCGGCCGGCCAAGG      0.07196465      1192   ZNF593_TSS    gene_tss
    ## AACGAGAGACTGCTTGCTGG      0.03283749       693   ATPIF1_TSS    gene_tss
    ## AACGGCTCGGAAGCCTAGGG      0.07587158      1251    TIPRL_TSS    gene_tss
    ##                                        target     target_gene
    ## AAACCGCTCCCGAGCACGGG   chr1:26605667-26605668 ENSG00000142669
    ## AAATAGTGGGAAGATTCGTG  chr11:10530735-10530736 ENSG00000255823
    ## AACACACCACGGAGGAGTGG chr8:124191287-124191288 ENSG00000147689
    ## AACAGCCCGGCCGGCCAAGG   chr1:26496362-26496363 ENSG00000142684
    ## AACGAGAGACTGCTTGCTGG   chr1:28562620-28562621 ENSG00000130770
    ## AACGGCTCGGAAGCCTAGGG chr1:168148171-168148172 ENSG00000143155

``` r
# all other gRNAs
gRNA_expression_odm |>
  get_feature_covariates() |>
  filter(target_type != "gene_tss") |> head()
```

    ##                      mean_expression n_nonzero         gRNA_group
    ## TAGAGGGGTGAAACAAAAAG      0.01610040       291 chr10.1007_top_two
    ## CCATCCAGCAAACTTAGCAG      0.04669020       864 chr10.1007_top_two
    ## GTCTTGCTCTGCTTAAAGTA      0.04682526       990 chr10.1018_top_two
    ## ATAGATAAGGAGGAAACAGG      0.02062472       384 chr10.1018_top_two
    ## GCGAAGATGATCAAACCCCC      0.06437750      1100 chr10.1019_top_two
    ## AATGACTGCAAAGGTCTGGG      0.03107214       593 chr10.1019_top_two
    ##                             target_type                  target target_gene
    ## TAGAGGGGTGAAACAAAAAG candidate_enhancer chr10:29023920-29024458        <NA>
    ## CCATCCAGCAAACTTAGCAG candidate_enhancer chr10:29023920-29024458        <NA>
    ## GTCTTGCTCTGCTTAAAGTA candidate_enhancer chr10:29111303-29112331        <NA>
    ## ATAGATAAGGAGGAAACAGG candidate_enhancer chr10:29111303-29112331        <NA>
    ## GCGAAGATGATCAAACCCCC candidate_enhancer chr10:29115307-29115488        <NA>
    ## AATGACTGCAAAGGTCTGGG candidate_enhancer chr10:29115307-29115488        <NA>

Finally, `gRNA_group` contains the group to which each gRNA belongs, *as
defined by the original authors*. Negative control gRNAs are grouped
into pairs (with the exception of group “bassik_mch”, which contains
only a single gRNA).

``` r
gRNA_expression_odm |>
  get_feature_covariates() |>
  filter(target_type == "non-targeting") |> 
  pull(gRNA_group) |>
  table()
```

    ## 
    ##   bassik_mch     random_1    random_10    random_11    random_12    random_13 
    ##            1            2            2            2            2            2 
    ##    random_14    random_15    random_16    random_17    random_18    random_19 
    ##            2            2            2            2            2            2 
    ##     random_2    random_20    random_21    random_22    random_23    random_24 
    ##            2            2            2            2            2            2 
    ##    random_25     random_3     random_4     random_5     random_6     random_7 
    ##            2            2            2            2            2            2 
    ##     random_8     random_9  scrambled_1 scrambled_10 scrambled_11 scrambled_12 
    ##            2            2            2            2            2            2 
    ## scrambled_13 scrambled_14 scrambled_15 scrambled_16 scrambled_17 scrambled_18 
    ##            2            2            2            2            2            2 
    ## scrambled_19  scrambled_2 scrambled_20 scrambled_21 scrambled_22 scrambled_23 
    ##            2            2            2            2            2            2 
    ## scrambled_24 scrambled_25  scrambled_3  scrambled_4  scrambled_5  scrambled_6 
    ##            2            2            2            2            2            2 
    ##  scrambled_7  scrambled_8  scrambled_9 
    ##            2            2            2

# gRNA assignment

Finally, the gRNA assignment modality is a thresholded version of the
gRNA expression modality. The threshold used was 5 gRNA UMIs/cell, same
as in the original study.

``` r
gRNA_assignment_odm <- read_odm(paste0(processed_gRNA_assignment_dir, "matrix.odm"),
                                paste0(processed_gRNA_assignment_dir, "metadata.rds"))
gRNA_assignment_odm
```

    ## A covariate_ondisc_matrix with the following components:
    ##  A logical ondisc_matrix with 13189 features and 207324 cells.
    ##  A cell covariate matrix with columns n_nonzero.
    ##  A feature covariate matrix with columns n_nonzero, gRNA_group, target_type, target, target_gene.

The feature covariates of the gRNA assignment ODM coincide with those of
the gRNA expression ODM. The cell covariates, by contrast, are
`n_nonzero`.
