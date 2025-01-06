# github: https://github.com/Joker-Jerome/UTMOST

# step 1 Prepare gwas sumstats data
# Column names that need to be included： SNP effect_allele other_allele se beta pval

# step 2 input file
# 用于单组织分析的基因权重文件夹路径
model_db_path <- "./sample_data/weight_db_GTEx"

# 用于单组织分析的协方差矩阵covariance matrices文件夹路径
covariance <- "./sample_data/covariance_tissue"

# gwas 文件夹
gwas_folder <- "./sample_data/gwas"

# 文件名
gwas_file_pattern <- "gwas_file_pattern"
covariance_joint <- "./sample_data/covariance_joint/"


# step 3 single tissue
# 输出结果文件
results <- "./results"
if (!dir.exists(results)) {
  dir.create(results, recursive = T)
}

# 获取可分析的组织名称，也可直接写对应组织名，例如：  c("Adipose_Subcutaneous","Adipose_Visceral_Omentum")
single_tissue <- tools::file_path_sans_ext(dir(model_db_path, pattern = ".db"))
for (tissue in single_tissue) {
  out_file <- paste0(results, "/", tissue, ".csv")
  args <- c(
    "single_tissue_association_test.py",
    "--model_db_path ", shQuote(paste0(model_db_path, "/", tissue, ".db")),
    "--covariance ", shQuote(paste0(covariance, "/", tissue, ".txt.gz")),
    "--gwas_folder ", shQuote(gwas_folder),
    "--gwas_file_pattern ", gwas_file_pattern,
    "--snp_column ", "SNP",
    "--effect_allele_column ", "effect_allele",
    "--non_effect_allele_column ", "other_allele",
    "--beta_column ", "b",
    "--se_column ", "se",
    "--pvalue_column ", "p",
    "--output_file ", out_file
  )
  system2("python.exe", args = args)
}


# step 4 joint GBJ
results_gbj <- "./results_gbj"
if (!dir.exists(results_gbj)) {
  dir.create(results_gbj, recursive = T)
}

# 生成基因文件
results_file <- dir(results, ".csv$", full.names = T)
genes <- c()
for (f in results_file) {
  if (!file.exists(f)) {
    next
  }
  dat <- data.table::fread(f, data.table = F, fill = T)
  genes <- append(genes, dat$gene)
}
genes <- unique(genes)

gene_path <- paste0(results_gbj, "/gene_info.txt")
data.table::fwrite(data.frame(gene_ensg = genes), gene_path, sep = ",", quote = F, row.names = F, na = NA)

args <- c(
  "joint_GBJ_test.py",
  "--weight_db", shQuote(model_db_path),
  "--output_dir", shQuote(results_gbj),
  "--cov_dir", shQuote(covariance_joint),
  "--input_folder", shQuote(results),
  "--gene_info", shQuote(gene_path),
  "--output_name", "joint_GBJ",
  "--start_gene_index", 1,
  "--end_gene_index", length(genes),
)

system2("python.exe", args = args)
