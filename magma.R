#  Official website: https://cncr.nl/research/magma/

# step 1 Prepare gwas sumstats data
# Column names that need to be included：SNP  chr pos samplesize

# step 2 input file
# magma可执行程序
magma_exe <- "e:/magma_v1.10_win/magma.exe"
# gwas文件
gwas_path <- "e:/data_prepare/ieu-a-30.txt"
# 参考文件
bfile <- "y:/data_ref/1kg.v3/EUR"
# 含基因位置信息的文件路径
geneloc_file <- "e:/magma_v1.10_win/MAGMA_gene_boundary_files/ENSGv110.coding.genes.txt"
# 含基因集的文件路径
geneset_file <- "e:/magma_v1.10_win/Gene_set_files/MSigDB_20231Hs_MAGMA.txt"
# 含基因表达水平基因集的文件路径
genecovar_file <- "e:/magma_v1.10_win/Gene_expression_files/gtex_v8_ts_avg_log2TPM_54tissues.txt"
# 是否进行基因集分析
geneset_anno <- T
# 输出文件前缀
out_prefix <- "magma"


# step 3 Prepare input file
gwas_dat <- data.table::fread(gwas_path, data.table = F, fill = T)
bim <- gwas_dat[, c("SNP", "chr", "pos")]
fn_snp_loc <- paste0(out_prefix, "_snp_loc.txt")
write.table(bim, file = fn_snp_loc, quote = F, row.names = F, sep = "\t")
snp_p <- data.frame(SNP = gwas_dat$SNP, P = gwas_dat$pval, samplesize = gwas_dat$samplesize)
fn_snp_p <- paste0(out_prefix, "_snp_pval.txt")
write.table(snp_p, file = fn_snp_p, quote = F, row.names = F, sep = "\t")

# step 4

fn_step1 <- paste0(out_prefix, "_step1", "_snp_loc")
system2(
  command = magma_exe,
  args = c(
    "--annotate",
    "--snp-loc",
    shQuote(fn_snp_loc),
    "--gene-loc",
    shQuote(geneloc_file),
    "--out",
    shQuote(fn_step1)
  )
)

# step 5

fn_step1_file <- paste0(fn_step1, ".genes.annot")
if (!file.exists(fn_step1_file)) {
  stop("not exist file")
}

fn_step2 <- paste0(out_prefix, "_step2", "_snp_loc")
system2(
  command = magma_exe,
  args = c(
    "--bfile",
    shQuote(bfile),
    "--pval",
    fn_snp_p,
    "ncol=",
    col_samplesize,
    "--gene-annot",
    shQuote(fn_step1_file),
    ifelse(geneset_anno == T, "", "--genes-only"),
    "--out",
    shQuote(fn_step2)
  )
)

# step 6
fn_step2_file <- paste0(fn_step2, ".genes.out")

gene <- data.table::fread(fn_step2_file, data.table = F)
geneloc <- data.table::fread(geneloc_file, data.table = F)
gene_mapped <- merge(geneloc[, c("V1", "V6")], gene, by.x = "V1", by.y = "GENE")
colnames(gene_mapped)[colnames(gene_mapped) == "V6"] <- "GENE_SYMBOL"
colnames(gene_mapped)[colnames(gene_mapped) == "V1"] <- "GENE_Ensembl_ID"

fn_step2_file_mapped <- paste0(out_prefix, "_step2_snp_loc.genes.out.mapped.csv")
write.csv(gene_mapped, file = fn_step2_file_mapped, row.names = F)

# step 7
if (geneset_anno == T) {
  fn_step3 <- paste0(out_prefix, "_step3_snp_loc_geneSET")
  system2(
    command = magma_exe,
    args = c(
      "--gene-results",
      shQuote(paste0(fn_step2, ".genes.raw")),
      "--set-annot",
      shQuote(geneset_file),
      "--out",
      shQuote(fn_step3)
    )
  )
}

# step 8
if (geneset_anno == T) {
  fn_step4 <- paste0(out_prefix, "_step4_snp_loc_geneCOVAR")

  system2(
    command = magma_exe,
    args = c(
      "--gene-results",
      shQuote(paste0(fn_step2, ".genes.raw")),
      "--gene-covar",
      shQuote(genecovar_file),
      "--model direction-covar=greater condition-hide=Average",
      "--out",
      shQuote(fn_step4)
    )
  )
}