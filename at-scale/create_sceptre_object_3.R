library(sceptre)
library(monocle)
library(dplyr)

# set paths
gasp_offsite <- paste0(.get_config_path("LOCAL_GASPERINI_2019_V3_DATA_DIR"), "at-scale/")
raw_data_dir <- paste0(gasp_offsite, "raw/")
processed_data_dir <- paste0(gasp_offsite, "processed/")
intermediate_data_dir <- paste0(.get_config_path("LOCAL_GASPERINI_2019_V3_DATA_DIR"), "at-scale/intermediate/")
if (!dir.exists(processed_data_dir)) dir.create(path = processed_data_dir, recursive = TRUE)

# load data
monocle_obj <- readRDS(paste0(raw_data_dir, "/GSE120861_at_scale_screen.cds.rds"))
cell_metadata <- pData(monocle_obj)
gene_expression_matrix <- exprs(monocle_obj)
grna_expression_matrix <- readRDS(paste0(intermediate_data_dir, "grna_count_matrix.rds"))
grna_feature_covariates <- readRDS(paste0(intermediate_data_dir, "grna_feature_covariates.rds"))

rm(monocle_obj)

# subset grna expression matrix to grnas contained in grna_feature_covariates
grna_expression_matrix <- grna_expression_matrix[grna_feature_covariate_df$barcode,]

# construct grna_target_data_frame
grna_target_data_frame <- grna_feature_covariates |>
  dplyr::mutate(grna_id = barcode,
                grna_target = if_else(!is.na(target_gene), target_gene, target)) |>
  dplyr::select(grna_id, grna_target)

# extract extra covariates
extra_covariates <- cell_metadata[,c("prep_batch", "percent.mito")]

# import into sceptre_object
sceptre_object <- import_data(response_matrix = gene_expression_matrix,
                              grna_matrix = grna_expression_matrix,
                              grna_target_data_frame = grna_target_data_frame,
                              moi = "high",
                              extra_covariates = extra_covariates,
                              use_ondisc = TRUE, 
                              directory = processed_data_dir)