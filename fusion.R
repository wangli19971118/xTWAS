# Official website: http://gusevlab.org/projects/fusion/
# github: https://github.com/gusevlab/fusion_twas

# step 1 Prepare gwas sumstats data
# Column names that need to be included：  SNP  A1  A2  Z

# step 2 input file
# gwas文件
sumstats <- "PGC2.SCZ.sumstats"
# 权重文件pos
weights <- "./WEIGHTS/GTEx.Whole_Blood.pos"
# 权重文件dir
weights_dir <- "./WEIGHTS/"
# 参考文件
ref_ld_chr <- "./LDREF/1000G.EUR."

# step 3 Performing the expression imputation

# 输出结果文件夹
assoc_test_result <- "./assoc_test_result"
if (!dir.exists(assoc_test_result)) {
  dir.create(assoc_test_result, recursive = T)
}

for (i in seq_len(22)) {
  # 输出文件
  out_file <- paste0(assoc_test_result, "/fusion_assoc_", i, ".dat")
  args <- c(
    "FUSION.assoc_test.R",
    "--sumstats", shQuote(sumstats),
    "--weights", shQuote(weights),
    "--weights_dir", shQuote(weights_dir),
    "--ref_ld_chr", shQuote(ref_ld_chr),
    "--chr", i,
    "--out", out_file
  )
  system2("Rscript.exe", args = args)
}

# 汇总结果数据
fusion_sum <- data.frame()
for (i in dir(assoc_test_result, ".dat$", full.names = T)) {
  if (!file.exists(i)) {
    next
  }
  temp_dat <- data.table::fread(i, data.table = F, fill = T)
  fusion_sum <- rbind(fusion_sum, temp_dat)
}
fn_sumstats <- paste0(assoc_test_result, "/fusion_assoc_summary.dat")
data.table::fwrite(fusion_sum, file = fn_sumstats, row.names = F, sep = "\t", quote = F)


# step 4 Joint/conditional tests and plots

assoc_join_result <- "./assoc_join_result"
if (!dir.exists(assoc_join_result)) {
  dir.create(assoc_join_result, recursive = T)
}

for (i in seq_len(22)) {
  out_file <- paste0(result, "/fusion_assoc_", i, ".dat")
  args <- c(
    "FUSION.post_process.R",
    "--sumstats", shQuote(sumstats),
    "--input", shQuote(fn_sumstats),
    "--ref_ld_chr", shQuote(ref_ld_chr),
    "--chr", i,
    "--out", out_file,
    "--plot",
    "--locus_win", 100000
  )
  system2("Rscript.exe", args = args)
}