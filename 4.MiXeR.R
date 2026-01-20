
###1.根据1000g生成22个染色体的.bim / .bed 文件[.sh]
# 假设原始文件是 g1000_eur.bed/bim/fam
for i in {1..22}
do
/users/zhangjh/202509_gwas/tools/MR/plink/plink \
--bfile /media/DS1/zhangjh/202509_GWAS/Ref/MR/LD/1kg.v3/EUR \
--chr $i \
--make-bed \
--out /media/DS1/zhangjh/202601_GWAS/ref/MiXeR/1000G_EUR_QC/1000G_EUR_QC_$i
done

###2.用 MiXeR生成每个染色体的.ld 文件,约2h（Docker容器中）[.sh]
for i in {1..22}
do
python /tools/mixer/precimed/mixer.py ld \
--bfile /home/ref/MiXeR/1000G_EUR_QC/1000G_EUR_QC_$i \ 
--out /home/ref/MiXeR/1000G_EUR_QC_ld/1000G_EUR_QC_$i.run4.ld \ 
--r2min 0.05 \
--ldscore-r2min 0.05 \
--ld-window-kb 30000
done


###3.生成.snp文件（Docker容器中）[.sh]
#全基因组生成一个
python /tools/mixer/precimed/mixer.py snps \
--lib /tools/mixer/src/build/lib/libbgmg.so \
--bim-file /home/ref/MiXeR/1000G_EUR_QC/1000G_EUR_QC_@.bim \
--ld-file  /home/ref/MiXeR/1000G_EUR_QC_ld/1000G_EUR_QC_@.run4.ld \
--out /home/ref/MiXeR/1000G_EUR_QC_snp/1000G_EUR_QC_all.snps \
--maf 0.05 \
--r2 0.8 \
--subset 10000000 \
--seed 1


###4.准备gwas数据（hg19）
##hg38转为hg19
#1.提取.bed文件【.R】
library(data.table)
# 读取你的 summary stats
gwas <- fread("/media/DS1/zhangjh/202601_GWAS/input/kidney_cardiovascular/CKDIorII_qc.txt", header = TRUE, sep = "\t")
head(gwas)
# 生成 BED 文件数据框( BED 格: chr, start(0-based), end(1-based), name(rsID))
bed_df <- gwas[, .(
  chrom = paste0("chr", chromosome),  # 加 chr 前缀，hg38ToHg19.over.chain.gz 链文件一般使用 chr10 而不是 10
  start = as.integer(base_pair_location - 1),  # 强制转换为整数
  end   = as.integer(base_pair_location),
  name  = rsid
)]
head(bed_df)
# 写出 BED 文件
fwrite(bed_df, "/users/zhangjh/202601_GWAS/tools/QC/bed/CKDIorII_qc_hg38.bed",
       sep = "\t", col.names = FALSE)
#2.运行坐标转换.【.sh】
/users/zhangjh/202601_GWAS/tools/QC/liftOver \
/users/zhangjh/202601_GWAS/tools/QC/bed/CKDIorII_qc_hg38.bed \
/users/zhangjh/202601_GWAS/tools/QC/hg38ToHg19.over.chain.gz \
/users/zhangjh/202601_GWAS/tools/QC/bed/CKDIorII_qc_hg19.bed \
/users/zhangjh/202601_GWAS/tools/QC/bed/CKDIorII_unmapped.txt
#内存不够时可以分染色体进行再合并，意义同上
#分开进行
for chr in {1..22}; do
echo "Processing chr$chr"
awk -v c="$chr" '$1=="chr"c {print $0}' /users/zhangjh/202601_GWAS/tools/QC/bed/CKDIorII_qc_hg38.bed \
> /users/zhangjh/202601_GWAS/tools/QC/bed/CKDIorII_qc_hg38_chr${chr}.bed

/users/zhangjh/202601_GWAS/tools/QC/liftOver \
/users/zhangjh/202601_GWAS/tools/QC/bed/CKDIorII_qc_hg38_chr${chr}.bed \
/users/zhangjh/202601_GWAS/tools/QC/hg38ToHg19.over.chain.gz \
/users/zhangjh/202601_GWAS/tools/QC/bed/CKDIorII_qc_hg19_chr${chr}.bed \
/users/zhangjh/202601_GWAS/tools/QC/bed/CKDIorII_unmapped_chr${chr}.txt
done
#合并结果（如果需要全基因组 BED）：
cat /users/zhangjh/202601_GWAS/tools/QC/bed/CKDIorII_qc_hg19_chr*.bed > /users/zhangjh/202601_GWAS/tools/QC/bed/CKDIorII_qc_hg19.bed
cat /users/zhangjh/202601_GWAS/tools/QC/bed/CKDIorII_unmapped_chr*.txt > /users/zhangjh/202601_GWAS/tools/QC/bed/CKDIorII_unmapped.txt
#查看并去除异常数据（有的话）
awk 'NF!=4' /users/zhangjh/202601_GWAS/tools/QC/bed/CKDIorII_qc_hg19.bed | head
awk 'NF==4' /users/zhangjh/202601_GWAS/tools/QC/bed/CKDIorII_qc_hg19.bed > /users/zhangjh/202601_GWAS/tools/QC/bed/CKDIorII_qc_hg19.clean.bed
#删除
rm -f /users/zhangjh/202601_GWAS/tools/QC/bed/CKDIorII_qc_hg38_chr*.bed
rm -f /users/zhangjh/202601_GWAS/tools/QC/bed/CKDIorII_qc_hg19_chr*.bed
rm -f /users/zhangjh/202601_GWAS/tools/QC/bed/CKDIorII_unmapped_chr*.txt

#3.合并hg19信息并提取所需列
# 读取转换后的 hg38 坐标文件 (4列：chr, start, end, ID)
# V1 = CHROM_hg19, V3 = POS_hg19, V4 = ID
new_coords <- fread("/users/zhangjh/202601_GWAS/tools/QC/bed/CKDIorII_qc_hg19.clean.bed", 
                    select = c(1, 3, 4), 
                    col.names = c("CHROM_hg19", "POS_hg19", "ID"))

# 2. 去掉染色体编号前的 "chr" 前缀 (例如把 chr1 变成 1)，使其与原始数据格式统一
new_coords[, CHROM_hg19 := gsub("chr", "", CHROM_hg19)]

# 3. 读取原始gwas数据
raw_vcf <- fread("/media/DS1/zhangjh/202601_GWAS/input/kidney_cardiovascular/CKDIorII_qc.txt", 
                 sep="\t", header=TRUE)
# 4. 按 ID 进行合并 (将 hg19 坐标信息加入到原始数据中)
# 我们使用 merge 仅保留转换成功的位点
head(raw_vcf)
trait_file <- merge(new_coords, raw_vcf, by.x = "ID", by.y = "rsid")
#① 读入数据（你已经完成）
#library(data.table)
#trait_file <- fread(
#  "/media/DS1/zhangjh/202601_GWAS/input/kidney_cardiovascular/IgA_qc.txt",
#  sep = "\t",
#  header = TRUE
#)
head(trait_file)
#② 计算 Z 值（MiXeR 强烈推荐）
#MiXeR 优先用 Z，而不是 beta + se。
trait_file[, Z := beta / standard_error]

#trait_file$Z <- qnorm(trait_file$p_value / 2, lower.tail = FALSE) * sign(trait_file$beta)# 缺少se时使用 p_value 计算 Z-score 的绝对值，并结合 beta 的符号确定方向
#③ 去掉 Z 无法计算或非法的 SNP（非常关键）
#MiXeR 不能接受 NA / Inf
trait_file <- trait_file[
  is.finite(Z) &
    !is.na(ID) &
    !is.na(effect_allele) &
    !is.na(other_allele)
]

#④ 构造 MiXeR 所需的 7 列
# 保留标准染色体上的snp
trait_file <- trait_file[
  CHROM_hg19 %in% as.character(1:22) | CHROM_hg19 %in% c("X","Y","MT")
]
#⚠️ A1 = effect_allele（Z 的正方向）
mixer_dt <- trait_file[, .(
  SNP = ID,
  CHR = as.integer(CHROM_hg19),
  BP  = as.integer(POS_hg19),
  A1  = toupper(effect_allele),
  A2  = toupper(other_allele),
  Z   = Z
)]
#⑤ 加上样本量 N（必须）
#如果 AF GWAS 样本量是固定的（例如 N = 1,030,836）：
N_AF <- 431894  # ← 请替换为你的真实 N
mixer_dt[, N := N_AF]
#⑥ 额外质量控制（强烈建议）
#1️⃣ 去除非 A/C/G/T 等位基因
mixer_dt <- mixer_dt[
  A1 %in% c("A","C","G","T") &
    A2 %in% c("A","C","G","T")
]
#2️⃣ 去除 A1 == A2 的异常 SNP
mixer_dt <- mixer_dt[A1 != A2]
#3️⃣ 去重（保留第一条）
setkey(mixer_dt, SNP)
mixer_dt <- unique(mixer_dt, by = "SNP")
head(mixer_dt)
#⑦ 写出 MiXeR 输入文件（TSV）
fwrite(
  mixer_dt,
  file = "/media/DS1/zhangjh/202601_GWAS/input/MiXeR/CKDIorII_mixer.tsv",
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

###5.MiXeR
###1️⃣ 配置变量（Docker容器中）[.sh]
#!/bin/bash
# 目录配置
WKDIR=/home/output/MiXeR
PYTHON=/usr/local/bin/python3
MIXER=/tools/mixer/precimed/mixer.py
LIBBGMG=/tools/mixer/src/build/lib/libbgmg.so
LDSR=/home/ref/MiXeR/1000G_EUR_QC
SNPS=/home/ref/MiXeR/1000G_EUR_QC_snp/1000G_EUR_QC_all.snps

# GWAS 文件
TRAIT1=/home/input/MiXeR/AF_mixer.tsv
TRAIT2=/home/input/MiXeR/CKDIorII_mixer.tsv
# 输出文件夹
mkdir -p ${WKDIR}/fit1
mkdir -p ${WKDIR}/fit2
mkdir -p ${WKDIR}/test2
mkdir -p ${WKDIR}/results

###2️⃣ Fit1：每染色体单性状建模【.sh】
#第一步：split GWAS summary statistics（必须）
${PYTHON} ${MIXER} split_sumstats \
--trait1-file /home/input/MiXeR/CKDIorII_mixer.tsv \
--out /home/input/MiXeR/CKDIorII_mixer.chr@.sumstats.gz
#第二步：用「每条染色体自己的 GWAS 文件」跑 fit1
for i in {1..22}; do
python /tools/mixer/precimed/mixer.py fit1 \
--trait1-file /home/input/MiXeR/CKDIorII_mixer.chr${i}.sumstats.gz \
--bim-file /home/ref/MiXeR/1000G_EUR_QC/1000G_EUR_QC_${i}.bim \
--ld-file  /home/ref/MiXeR/1000G_EUR_QC_ld/1000G_EUR_QC_${i}.run4.ld \
--extract  /home/ref/MiXeR/1000G_EUR_QC_snp/1000G_EUR_QC_all.snps \
--chr      ${i} \
--out      /home/output/MiXeR/fit1/CKDIorII.fit.rep${i} \
--lib      /tools/mixer/src/build/lib/libbgmg.so
done
###3️⃣ Fit2：两性状联合建模
for i in {1..22}; do
echo "Fit2 chr$i"
${PYTHON} ${MIXER} fit2 \
--trait1-file ${TRAIT1} \
--trait2-file ${TRAIT2} \
--trait1-params-file ${WKDIR}/fit1/AF.fit.rep${i}.json.json \
--trait2-params-file ${WKDIR}/fit1/CKDIorII.fit.rep${i}.json \
--out ${WKDIR}/fit2/AF_vs_CKDIorII.fit.rep${i} \
--extract ${SNPS} \
--bim-file ${LDSR}/1000G_EUR_QC_${i}.bim \
--ld-file ${LDSR}_ld/1000G_EUR_QC_${i}.run4.ld \
--chr ${i} \
--lib ${LIBBGMG}
done

###4️⃣ Test2：预测性能评估
for i in {1..22}; do
echo "Test2 chr$i"
${PYTHON} ${MIXER} test2 \
--trait1-file ${TRAIT1} \
--trait2-file ${TRAIT2} \
--load-params-file ${WKDIR}/fit2/AF_vs_CKDIorII.fit.rep${i}.json \
--out ${WKDIR}/test2/AF_vs_CKDIorII.test.rep${i}\
--bim-file ${LDSR}/1000G_EUR_QC_${i}.bim \
--ld-file ${LDSR}_ld/1000G_EUR_QC_${i}.run4.ld \
--chr ${i} \
--lib ${LIBBGMG}
done

###5️⃣ Combine：合并所有染色体结果
MFIGURE=/tools/mixer/precimed/mixer_figures.py
ID=AF_vs_CKDIorII
# 合并 fit
${PYTHON} ${MFIGURE} combine \
--json ${WKDIR}/fit2/${ID}.fit.rep@.json \
--out ${WKDIR}/results/${ID}.fit.json
# 合并 test
${PYTHON} ${MFIGURE} combine \
--json ${WKDIR}/test2/${ID}.test.rep@.json \
--out ${WKDIR}/results/${ID}.test

###6️⃣ Two：绘制联合图 / 计算统计量
${PYTHON} ${MFIGURE} two \
--json-fit ${WKDIR}/results/${ID}.fit.json \
--json-test ${WKDIR}/results/${ID}.test.json \
--trait1 AF \
--trait2 CKDIorII \
--out ${WKDIR}/results/${ID} \
--ext png svg \
--statistic mean std












