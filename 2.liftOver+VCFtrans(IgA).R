





library(data.table)

# 1. 读取转换后的 hg38 坐标文件 (4列：chr, start, end, ID)
# V1 = CHROM_hg38, V3 = POS_hg38, V4 = ID
new_coords <- fread("/users/zhangjh/202601_GWAS/tools/QC/IgA_raw_hg38.bed", 
                    select = c(1, 3, 4), 
                    col.names = c("CHROM_hg38", "POS_hg38", "ID"))

# 2. 去掉染色体编号前的 "chr" 前缀 (例如把 chr1 变成 1)，使其与原始数据格式统一
new_coords[, CHROM_hg38 := gsub("chr", "", CHROM_hg38)]

# 3. 读取原始的 VCF 数据
raw_vcf <- fread("/media/DS1/zhangjh/202509_GWAS/input/kidney_disease/IgA_nephropathy/ebi-a-GCST90018866.vcf.gz", 
                 sep="\t", header=TRUE, skip="#CHROM")

# 4. 按 ID 进行合并 (将 hg38 坐标信息加入到原始数据中)
# 我们使用 merge 仅保留转换成功的位点
data <- merge(new_coords, raw_vcf, by.x = "ID", by.y = "ID")

#5.VCF转txt
# 提取并重命名列
data[, CHR_hg38 := CHROM_hg38]
data[, BP_hg38 := POS_hg38]
data[, SNP := ID]
data[, effect_allele := ALT]
data[, other_allele := REF]

# 分割 FORMAT 列和实际数据列(拆分耗时比较久，因为数据太大啦!)
format_cols <- strsplit(data$FORMAT[1], ":")[[1]]
data_values <- data[, tstrsplit(`ebi-a-GCST90018866`, ":", fixed = TRUE)]

# 获取 AF 对应的列索引
af_idx <- which(format_cols == "AF")

# 提取 AF 作为 eaf
data[, eaf := as.numeric(data_values[[af_idx]])]

# 获取其他对应列名的位置
es_idx <- which(format_cols == "ES")
se_idx <- which(format_cols == "SE")
lp_idx <- which(format_cols == "LP")
#ss_idx <- which(format_cols == "SS")

# 创建新列
data[, beta := as.numeric(data_values[[es_idx]])]
data[, se := as.numeric(data_values[[se_idx]])]
data[, pval := 10^(-as.numeric(data_values[[lp_idx]]))]  # -log10 to p-value
#data[, n := as.numeric(data_values[[ss_idx]])]

# 选择并重命名所需的列
standard_data <- data[, .(SNP, CHR_hg38, BP_hg38, effect_allele, other_allele, eaf, beta, se, pval)]

# 保存或进一步处理
fwrite(standard_data, "/media/DS1/zhangjh/202509_GWAS/input/kidney_disease/IgA_nephropathy/IgA_nephropathy_standard_gwas_data.txt", sep = "\t", quote = FALSE)
