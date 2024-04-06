Gasperini (2019) data documentation
================
Gene Katsevich
2024-04-06

``` r
library(sceptre)
library(dplyr)
library(conflicted)
conflict_prefer("filter", "dplyr")
```

# Overview

This repository contains code to import and process the [Gasperini
2019](https://pubmed.ncbi.nlm.nih.gov/30612741/) data. Gasperini et
al. developed a high MOI, single-cell CRISPR screen assay to map
putative enhancers at genome-wide scale in a population of \>200,000
K562 cells.

The `gasperini-2019-v3` directory structure is as follows:

    └── at-scale
        ├── intermediate
        ├── processed
        │   ├── grna.odm
        │   ├── response.odm
        │   ├── sceptre_object.rds
        │   └── sceptre_object_in_memory.rds
        └── raw

The contents of the `raw` and `intermediate` directories are suppressed,
as they are unimportant. The `grna.odm`, `response.odm`, and
`sceptre_object.rds` files are the three pieces of an ondisc-backed
`sceptre_object`. The `sceptre_object_mem.rds` object is an in-memory
`sceptre_object` containing the same information. It requires about 2GB
on disk and 8GB in memory. Let us load the disk-backed `sceptre_object`:

``` r
gasp_offsite <- paste0(.get_config_path("LOCAL_GASPERINI_2019_V3_DATA_DIR"), "at-scale/")
processed_data_dir <- paste0(gasp_offsite, "processed/")
sceptre_object <- read_ondisc_backed_sceptre_object(
  sceptre_object_fp = paste0(processed_data_dir, "sceptre_object.rds"),
  response_odm_file_fp = paste0(processed_data_dir, "response.odm"),
  grna_odm_file_fp = paste0(processed_data_dir, "grna.odm")
)
```

We can print the `sceptre_object` for more information:

``` r
sceptre_object
```

    ## An object of class sceptre_object.
    ## 
    ## Attributes of the data:
    ##  • 207324 cells
    ##  • 13135 responses
    ##  • High multiplicity-of-infection 
    ##  • 12976 targeting gRNAs (distributed across 6105 targets) 
    ##  • 101 non-targeting gRNAs 
    ##  • 6 covariates (grna_n_nonzero, grna_n_umis, prep_batch, response_n_nonzero, response_n_umis, response_p_mito)

Let us take a look at the cell covariates:

``` r
sceptre_object |> get_cell_covariates() |> head()
```

    ##                                  response_n_nonzero response_n_umis
    ## AAACCTGAGAGGTACC-1_1A_1_SI-GA-E2               3549           17566
    ## AAACCTGAGTCAATAG-1_1A_1_SI-GA-E2               2543            8917
    ## AAACCTGCAAACAACA-1_1A_1_SI-GA-E2               3191           14626
    ## AAACCTGCACTTCTGC-1_1A_1_SI-GA-E2               4539           22783
    ## AAACCTGCATGTAGTC-1_1A_1_SI-GA-E2               2605           10124
    ## AAACCTGGTAGCGCAA-1_1A_1_SI-GA-E2               2187            9743
    ##                                  response_p_mito grna_n_nonzero grna_n_umis
    ## AAACCTGAGAGGTACC-1_1A_1_SI-GA-E2      0.05880679             83         980
    ## AAACCTGAGTCAATAG-1_1A_1_SI-GA-E2      0.03611080             44         347
    ## AAACCTGCAAACAACA-1_1A_1_SI-GA-E2      0.06987556             81         919
    ## AAACCTGCACTTCTGC-1_1A_1_SI-GA-E2      0.02620375             55         568
    ## AAACCTGCATGTAGTC-1_1A_1_SI-GA-E2      0.00800079             69        1018
    ## AAACCTGGTAGCGCAA-1_1A_1_SI-GA-E2      0.02237504             78        1275
    ##                                    prep_batch
    ## AAACCTGAGAGGTACC-1_1A_1_SI-GA-E2 prep_batch_1
    ## AAACCTGAGTCAATAG-1_1A_1_SI-GA-E2 prep_batch_1
    ## AAACCTGCAAACAACA-1_1A_1_SI-GA-E2 prep_batch_1
    ## AAACCTGCACTTCTGC-1_1A_1_SI-GA-E2 prep_batch_1
    ## AAACCTGCATGTAGTC-1_1A_1_SI-GA-E2 prep_batch_1
    ## AAACCTGGTAGCGCAA-1_1A_1_SI-GA-E2 prep_batch_1

We can see that there are two batches in these data:

``` r
sceptre_object |> get_cell_covariates() |> count(prep_batch)
```

    ##     prep_batch      n
    ## 1 prep_batch_1 107947
    ## 2 prep_batch_2  99377

Let us take a look at the gRNA-target data frame:

``` r
sceptre_object@grna_target_data_frame
```

    ## # A tibble: 13,077 × 5
    ##    grna_id              grna_target     chr       start       end
    ##    <chr>                <chr>           <chr>     <int>     <int>
    ##  1 AAGCCCGGCCAGCCGACCGG ENSG00000075624 chr7    5570339   5570340
    ##  2 ACGGGCGGCGGATCGGCAAG ENSG00000075624 chr7    5570339   5570340
    ##  3 AGAGTGCGAGAGCTGGCAGG ENSG00000184009 chr17  79479815  79479816
    ##  4 CACCGGCAGAGAAACGCGAG ENSG00000184009 chr17  79479815  79479816
    ##  5 CTGGACCGAGCTGCAGAGTG ENSG00000119640 chr14  75530722  75530723
    ##  6 TCCGGGACCACCACGCCAAG ENSG00000119640 chr14  75530722  75530723
    ##  7 CTACATCCCGCGGCGCTGGG ENSG00000159346 chr1  202927391 202927392
    ##  8 GGGAGGGCGCTGAAGATCGG ENSG00000159346 chr1  202927391 202927392
    ##  9 CCAAGGCGTCCTCAGACCAG ENSG00000128918 chr15  58306298  58306299
    ## 10 GCCGGGTGTCCCTAGCCCGG ENSG00000128918 chr15  58306298  58306299
    ## # ℹ 13,067 more rows

The gRNA id is the gRNA barcode. Different gRNAs have different target
types.

Some target genes:

``` r
sceptre_object@grna_target_data_frame |>
  filter(grepl("ENSG", grna_target))
```

    ## # A tibble: 754 × 5
    ##    grna_id              grna_target     chr       start       end
    ##    <chr>                <chr>           <chr>     <int>     <int>
    ##  1 AAGCCCGGCCAGCCGACCGG ENSG00000075624 chr7    5570339   5570340
    ##  2 ACGGGCGGCGGATCGGCAAG ENSG00000075624 chr7    5570339   5570340
    ##  3 AGAGTGCGAGAGCTGGCAGG ENSG00000184009 chr17  79479815  79479816
    ##  4 CACCGGCAGAGAAACGCGAG ENSG00000184009 chr17  79479815  79479816
    ##  5 CTGGACCGAGCTGCAGAGTG ENSG00000119640 chr14  75530722  75530723
    ##  6 TCCGGGACCACCACGCCAAG ENSG00000119640 chr14  75530722  75530723
    ##  7 CTACATCCCGCGGCGCTGGG ENSG00000159346 chr1  202927391 202927392
    ##  8 GGGAGGGCGCTGAAGATCGG ENSG00000159346 chr1  202927391 202927392
    ##  9 CCAAGGCGTCCTCAGACCAG ENSG00000128918 chr15  58306298  58306299
    ## 10 GCCGGGTGTCCCTAGCCCGG ENSG00000128918 chr15  58306298  58306299
    ## # ℹ 744 more rows

Some target enhancers:

``` r
sceptre_object@grna_target_data_frame |>
  filter(grepl("chr", grna_target))
```

    ## # A tibble: 12,222 × 5
    ##    grna_id              grna_target             chr      start      end
    ##    <chr>                <chr>                   <chr>    <int>    <int>
    ##  1 CGGCCGGGTGAGCTCCTCAG chr16:1470800-1470801   chr16  1470800  1470801
    ##  2 GGCGTCAGTCGAGGAGTCAG chr16:1470800-1470801   chr16  1470800  1470801
    ##  3 CTCGCAGCCGTCATGGCAGG chr1:85725315-85725316  chr1  85725315 85725316
    ##  4 TTGTCGCTCGCAGCCGTCAG chr1:85725315-85725316  chr1  85725315 85725316
    ##  5 GGTGGCCCGGGTCTATAATG chr21:33984895-33984896 chr21 33984895 33984896
    ##  6 TGGGATTCCGGACTGTGAGG chr21:33984895-33984896 chr21 33984895 33984896
    ##  7 GCCCGACTGTTGAGTGAGGG chr6:31802385-31802386  chr6  31802385 31802386
    ##  8 GCGAAGTGGGCTCGTGGGTG chr6:31802385-31802386  chr6  31802385 31802386
    ##  9 TAGAGGGGTGAAACAAAAAG chr10:29023920-29024458 chr10 29023920 29024458
    ## 10 CCATCCAGCAAACTTAGCAG chr10:29023920-29024458 chr10 29023920 29024458
    ## # ℹ 12,212 more rows

Some are non-targeting:

``` r
sceptre_object@grna_target_data_frame |>
  filter(grna_target == "non-targeting")
```

    ## # A tibble: 101 × 5
    ##    grna_id              grna_target   chr   start   end
    ##    <chr>                <chr>         <chr> <int> <int>
    ##  1 CTCACCTGTGGGAGTAACGG non-targeting <NA>     NA    NA
    ##  2 AATCCTCTAATGGACGAAGA non-targeting <NA>     NA    NA
    ##  3 AGATACCTATGGCCATATAG non-targeting <NA>     NA    NA
    ##  4 AGAGGTAACCAAAATAGCAA non-targeting <NA>     NA    NA
    ##  5 ATATGTAACCTCCAGAATGA non-targeting <NA>     NA    NA
    ##  6 TACAACTGCATTACATGCCA non-targeting <NA>     NA    NA
    ##  7 TCCGCAGTCAAAAGACCGAG non-targeting <NA>     NA    NA
    ##  8 AATATTCTCCCTCATTCTGG non-targeting <NA>     NA    NA
    ##  9 TTAAAATTGATTCTGCCACT non-targeting <NA>     NA    NA
    ## 10 ATTGGCAGTCTCTAAGAAGT non-targeting <NA>     NA    NA
    ## # ℹ 91 more rows
