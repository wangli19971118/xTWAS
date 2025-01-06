

library(QTLMR)

setwd("F:/QTLMR")
# 配置FOCUS软件环境，如何已经本地配置过，勿需再在线配置python环境

TWAS_FOCUS_install(force = TRUE)


#########  芬兰12 版本
a <- data_info_finngen()


# 1 数据转换 --------------------------------------------------------------------
format_data_FinnGen(
  GWASfile = "./data_raw/finngen_R11_N14_PROSTHYPERPLA.gz",
  GWAS_name = "BPH_finn_37",
  save_path = "./1.data/R11/",
  type = "outcome",
  min_pval = 1e-200,
  build_to_hg19 = T,
  Twosample_dat = T,
  SMR_dat = TRUE,
  MTAG_dat = T,
  METAL_dat = T)


# 3.2 单组织分析：进行单个组织的TWAS-UTMOST分析 ------------------------------------------
data_utmost <- TWAS_UTMOST_single_tissue_test(
  test_help = F,  # 是否显示帮助文档信息
  GTEX_tissue = "all",  # 使用所有组织数据
  model_db_path = "./2.UTMOST/database_normalized_pruned/",  # 归一化后的模型数据库路径
  covariance = "./2.UTMOST/covariance_GTEx8_normalized_pruned/",  # 归一化后的协方差矩阵路径
  gwas_folder = "./1.data/R12/",  # GWAS数据的文件夹路径
  gwas_file_pattern = "BPH_finn_METAL.txt",  # GWAS的文件名
  snp_column = "SNP",  # SNP列的名称
  effect_allele_column = "effect_allele",  # 效应等位基因列的名称
  non_effect_allele_column = "other_allele",  # 非效应等位基因列的名称
  beta_column = "beta",  # beta系数列的名称
  pvalue_column = "pval",  # p值列的名称
  opt_arguments = NULL,  # 额外的参数
  save_path = "./2.UTMOST/UTMOST_result/",  # 保存结果的路径
  cores = 14  # 运行使用电脑线程数
)









#5.1 使用MAGMA软件进行基因层面和基因集合层面的关联分析  -----------------------------------------------------------------

gene_based_dat <- MAGMA_gene_based(
  GWAS_file = "./1.data/R11/BPH_finn_37_METAL.txt",  # 输入的GWAS文件路径
  bfile_1000G = "./6.MAGMA/g1000_eur/g1000_eur",  # 1000G参考数据路径
  gene_loc = "./6.MAGMA/ENSGv110.coding.genes.txt",  # 基因位置文件路径
  set_annot = "./6.MAGMA/MSigDB_20231Hs_MAGMA.txt",  # 基因集合注释文件路径
  SNP_P_col = c(3, 10),  # SNP p值列的位置
  samplesize_col = "N",  # 样本量列的名称
  save_name = "BPH",  # 保存的名称
  save_path = "./6.MAGMA"  # 保存路径
)




# 5.2 使用MAGMA进行组织特异性分析 ----------------------------------------------------

Tissue_specific_dat <- MAGMA_Tissue_specific(
  genes_raw = "./6.MAGMA/magma_step2_snp_loc.genes.raw",  # 原始基因数据路径
  gene_covar = "./6.MAGMA/gtex_v8_ts_avg_log2TPM.txt",  # 组织特异性表达数据路径
  save_name = "BPH",  # 保存的名称
  save_path = "./6.MAGMA/MAGMA_result/")  # 保存路径


# 5.3 对MAGMA基因层面分析结果进行基因注释及曼哈顿图绘制 -----------------------------------------

MAGMA_genes <- MAGMA_genes_Manhattanplot(
  genes_out = "./6.MAGMA/first_result/magma_step2_snp_loc.genes.out.txt",  # 基因分析结果文件路径
  gene_loc = "./6.MAGMA/ENSGv110.coding.genes.txt",  # 基因位置文件路径
  Manhtn_plot = TRUE,  # 是否绘制曼哈顿图
  threshold_sig = "bonferroni",  # 显著性阈值调整方法
  Manhtn_gene_sig = 10,  # 曼哈顿图的基因标记标签的数目
  signal_cex = 1,  # 显著性基因图形圆点的大小
  width = 9,  # 图形宽度
  height = 7,  # 图形高度
  save_name = "MAGMA_gsa"  # 保存的名称
)

# 筛选FDR值小于0.05的基因
MAGMA_genes <- subset(MAGMA_genes, MAGMA_genes$P.fdr < 0.05)

write.csv(MAGMA_genes, "./6.MAGMA/MAGMA_result/MAGMA_genes_fdr.csv", row.names = FALSE)


# #对MAGMA通路富集分析结果图形可视化 ----------------------------------------------------
plot_dat <- MAGMA_gsa_barplot(
  gsa_out = "./6.MAGMA/first_result/magma_step3_snp_loc_geneSET.gsa.out.txt",  # GSA结果文件路径
  showNum = 50,  # 显示的前50个通路
  X_text_size = 5,  # X轴文本大小
  X_text_angle = 80,  # X轴文本角度
  Y_text_size = 10,  # Y轴文本大小
  save_plot = TRUE,  # 是否保存图形
  pdf_name = "gsa_plot",  # PDF文件名
  width = 9,  # 图形宽度
  height = 7,  # 图形高度
  save_path = "./6.MAGMA/"  # 保存路径
)


# # 对MAGMA组织特异性分析结果图形可视化 --------------------------------------------------

Tissue_specific_dat <- MAGMA_gsa_barplot(
  gsa_out = "./6.MAGMA/MAGMA_result/BPH_gai.Tissue_specific.gsa.out.txt",  # 组织特异性GSA结果文件路径
  showNum = 54,  # 显示的54个组织柱状图
  X_text_size = 7,  # X轴文本大小
  X_text_angle = 70,  # X轴文本角度
  Y_text_size = 10,  # Y轴文本大小
  save_plot = TRUE,  # 是否保存图形
  pdf_name = "Tissue_specific_gsa_plot",  # PDF文件名
  width = 9,  # 图形宽度
  height = 7,  # 图形高度
  save_path = "./6.MAGMA/"  # 保存路径
)
















# FOCUS数据格式转换

TWAS_FOCUS_format_data(GWASfile = "./data_prepare/BPH_finn_TwosampleMR.txt",
                       N = NULL,
                       FOCUS_munge = TRUE,
                       save_name = "BPH_finn",
                       save_path = "./FOCUS/")


# 权重数据转换

# FOCUS_WEIGHTS文件夹是已经帮大家转换好了的GTExv8权重数据，大家可以不用再转换!!!

TWAS_FOCUS_import(import_help = F,
                  path = "./WEIGHTS/GTExv8.EUR.Whole_Blood.nofilter.pos",
                  resource = "fusion",
                  tissue = "Whole_Blood",
                  name = "Whole_Blood",
                  assay = "rnaseq",
                  opt_arguments = c("--use-ens-id","--from-gencode"),
                  save_name = "GTExv8.EUR.Whole_Blood.nofilter.pos",
                  save_path = "./FOCUS_WEIGHTS")

#    FOCUS分析循环

file_list <- list.files(path = "./FOCUS_WEIGHTS/",pattern = ".pos.db$")

FOCUS_dat_res <- c()
for (i in file_list) {

  FOCUS_dat <- TWAS_FOCUS_test(finemap_help = F,
                               Sumstatsfile = "./BPH_finn.sumstats.gz",
                               ref_ld_dir = "./EUR",
                               weights_file_db = paste0("./FOCUS_WEIGHTS/",i),
                               start_chr = 1,
                               end_chr = 22,
                               plot = TRUE,
                               locations = "38:EUR",
                               prior_prob = "gencode38",
                               p_threshold = 5e-08,
                               opt_arguments = NULL,
                               pip_sig = 0.8,
                               save_name = "BPH",
                               save_path = "./BPH_FOCUS",
                               cores = 1)

  FOCUS_dat_res <- rbind(FOCUS_dat_res,FOCUS_dat)
}



############ 不循环
TWAS_FOCUS_test(finemap_help = F,
                Sumstatsfile = "./BPH_finn.sumstats.gz",
                ref_ld_dir = "./EUR",
                weights_file_db = "./FOCUS_WEIGHTS/GTExv8.EUR.Esophagus_Muscularis.nofilter.pos.db",
                start_chr = 1,
                end_chr = 22,
                plot = TRUE,
                locations = "38:EUR",
                prior_prob = "gencode38",
                p_threshold = 5e-08,
                opt_arguments = NULL,
                pip_sig = 0.8,
                save_name = "BPH",
                save_path = "./BPH_FOCUS/Esophagus_Muscularis",
                cores = 2)


##########  读入数据后矫正
library(data.table)


# 获取所有csv文件名
files <- list.files(path = "./result/", pattern = "*.csv", full.names = TRUE)

# 读取并合并
data_list <- lapply(files, read.csv)
combined_data <- do.call(rbind, data_list)

# 写出合并后的文件
write.csv(combined_data, "./result/combined_output.csv", row.names = FALSE)

############# 进行fdr矫正
colnames(combined_data)
###筛选位于90%置信区间且pip>0.8的基因
pip_th=0.8  #阈值

focusall_sig = combined_data[combined_data$in_cred_set_pop1 == 1 & combined_data$pips_pop1> pip_th &combined_data$ens_gene_id!="NULL.MODEL"]

write.table(focusall_sig, "./result/focusall_sig20.txt", row.names = FALSE)




focusall_sig$FDR =p.adjust(Adipose_Subcutaneous$cv.R2.pval_pop1,method = "fdr")




?p.adjust

write.table(Insomnia, "./LDSC_SEG/Insomnia.cell_type_results_fdr.txt", row.names = FALSE)
write.table(BPH, "./LDSC_SEG/BPH.cell_type_results_fdr.txt", row.names = FALSE)







