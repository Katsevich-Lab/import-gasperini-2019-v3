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
grna_expression_matrix <- readRDS(paste0(intermediate_data_dir, "gRNA_count_matrix.rds"))
grna_feature_covariate_df <- readRDS(paste0(intermediate_data_dir, "gRNA_feature_covariates.rds"))
all_results <- readr::read_tsv(file = paste0(raw_data_dir, "GSE120861_all_deg_results.at_scale.txt"),
                               col_names = TRUE, col_types = "c")
gene_names <- monocle_obj@featureData@data$gene_short_name

# remove monocle object to save memory
rm(monocle_obj)

# subset grna expression matrix to grnas contained in grna_feature_covariates
grna_expression_matrix <- grna_expression_matrix[grna_feature_covariate_df$barcode,]

# extract gRNA target location information
grna_info <- all_results |> 
  dplyr::select(gRNA_group, target_site.chr, target_site.start, target_site.stop) |> 
  dplyr::rename(grna_group = gRNA_group) |>
  unique()

# construct grna_target_data_frame
grna_target_data_frame <- grna_feature_covariate_df |>
  dplyr::left_join(grna_info, by = "grna_group") |>
  dplyr::mutate(grna_id = barcode,
                grna_target = if_else(!is.na(target_gene), target_gene, target),
                chr = target_site.chr,
                start = target_site.start,
                end = target_site.stop) |>
  dplyr::select(grna_id, grna_target, chr, start, end)

# extract extra covariates
extra_covariates <- cell_metadata[, "prep_batch", drop = FALSE]

# import into ondisc-backed sceptre_object
sceptre_object_disk <- import_data(
  response_matrix = gene_expression_matrix,
  grna_matrix = grna_expression_matrix,
  grna_target_data_frame = grna_target_data_frame,
  moi = "high",
  extra_covariates = extra_covariates,
  response_names = gene_names,
  use_ondisc = TRUE,
  directory = processed_data_dir
)

# write sceptre_object to disk
write_ondisc_backed_sceptre_object(
  sceptre_object = sceptre_object_disk,
  directory_to_write = processed_data_dir
)

# import into in-memory sceptre_object
sceptre_object_mem <- import_data(
  response_matrix = gene_expression_matrix,
  grna_matrix = grna_expression_matrix,
  grna_target_data_frame = grna_target_data_frame,
  moi = "high",
  extra_covariates = extra_covariates,
  response_names = gene_names,
  use_ondisc = FALSE
)

# save to disk
saveRDS(sceptre_object_mem, 
        file = paste0(processed_data_dir, "sceptre_object_mem.rds"))