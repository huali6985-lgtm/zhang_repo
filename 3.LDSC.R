##基础知识：https://mp.weixin.qq.com/s/gDY5SqCjrGiOvEdgcSZnWg
##实战教程：https://github.com/mglev1n/ldscr
#devtools::install_github("mglev1n/ldscr")
setwd("/media/DS1/zhangjh/202601_GWAS/input/kidney_cardiovascular_LDSC/")
#### LDSC方法二 ####
#devtools::install_github("GenomicSEM/GenomicSEM")
library(GenomicSEM)
library(data.table)
library(dplyr)
library(stringr)
###1.整理含有必须列的标准数据
# convert_or_to_beta函数：将 odds_ratio 转换为 beta
convert_or_to_beta <- function(trait_files){
  for (file_path in trait_files){
    # 读取文件
    dt <- fread(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    # 检查是否存在 odds_ratio 列
    if(!"odds_ratio" %in% names(dt)){
      warning(paste("File", file_path, "does not contain 'odds_ratio'. Skipped."))
      next
    }
    # 将 odds_ratio 转为 beta
    dt[, beta := log(odds_ratio)]
    # 保存回原路径
    fwrite(dt, file_path, sep = "\t", quote = FALSE)
    message(paste("Processed and saved:", file_path))
  }
}
#输入
traits <- c("CKDIorII", "CKDIII", "CKDIV", "ESRD", "IHD")
trait_files <- paste0("/media/DS1/zhangjh/202601_GWAS/input/kidney_cardiovascular/", traits, "_qc.txt")
#运行convert_or_to_beta函数
convert_or_to_beta(trait_files)

###2.加载munge_data函数
munge_data<-function(file,hm3,N,trait.name){
  ref <- fread(hm3, header = T, data.table = F)
  if ("MAF" %in% colnames(file)) {
    file$MAF <- ifelse(file$MAF <= 0.5, file$MAF, (1 - file$MAF))
  }
  file$A1 <- factor(toupper(file$A1), c("A", "C", "G", "T"))
  file$A2 <- factor(toupper(file$A2), c("A", "C", "G", "T"))
  file <- merge(ref, file, by = "SNP", all.x = F, all.y = F)
  file <- subset(file, !(is.na(file$P)))
  file <- subset(file, !(is.na(file$effect)))
  a1 <- file$effect[[1]]
  file$effect <- ifelse(rep(round(median(file$effect, na.rm = T)) == 
                              1, nrow(file)), log(file$effect), file$effect)
  a2 <- file$effect[[1]]
  file$effect <- ifelse(file$A1.x != (file$A1.y) & file$A1.x == 
                          (file$A2.y), file$effect * -1, file$effect)
  file <- subset(file, !(file$A1.x != (file$A1.y) & file$A1.x != 
                           (file$A2.y)))
  file <- subset(file, !(file$A2.x != (file$A2.y) & file$A2.x != 
                           (file$A1.y)))
  file$Z <- sign(file$effect) * sqrt(qchisq(file$P, 1, lower = F))
  if ("N" %in% colnames(file)) {
    output <- cbind.data.frame(file$SNP, file$N, file$Z, 
                               file$A1.x, file$A2.x)
  }else {
    output <- cbind.data.frame(file$SNP, N, file$Z, file$A1.x, 
                               file$A2.x)
  }
  output$N<-N
  colnames(output) <- c("SNP", "N", "Z", "A1", "A2")
  trait.name <- str_replace_all(trait.name, fixed(" "), "")
  fwrite(x = output, file = paste0(trait.name, ".sumstats"), 
         sep = "\t", quote = FALSE, row.names = F)
  R.utils::gzip(paste0(trait.name, ".sumstats"), overwrite = overwrite)
  return()
}

###3.生成sumstats.gz文件
process_trait <- function(trait_name, file_path, hm3_snplist_path, N1) {
  # 读取数据
  trait_file <- fread(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  # 重命名列，适配 munge_data
  trait_file <- trait_file %>% 
    dplyr::rename(
      SNP = rsid,
      A1 = effect_allele,
      A2 = other_allele,
      effect = beta,
      P = p_value
    )
  # 调用 munge_data
  munge_data(trait_file, hm3 = hm3_snplist_path, N1, trait_name)
}
# 定义参数
trait_info <- list(
  "CRF"      = 458440,
  "CYSTKD"   = 458440,
  "NS"       = 458440,
  "HF"       = 458440,
  "SBP"      = 431014,
  "DBP"      = 431018,
  "VTE"      = 458440,
  "AMI"      = 458440,
  "PVD"      = 458440,
  "HTN"      = 458440
)
base_path <- "/media/DS1/zhangjh/202601_GWAS/input/kidney_cardiovascular/"
hm3_snplist <- "/media/DS1/zhangjh/202509_GWAS/Ref/LDSC/eur_w_ld_chr/w_hm3.snplist"
# 循环处理每个 trait
for (trait_name in names(trait_info)) {
  N1 <- trait_info[[trait_name]]
  file_path <- paste0(base_path, trait_name, "_qc.txt")
  
  process_trait(trait_name, file_path, hm3_snplist, N1)
}

###4.LDSC分析
run_pairwise_ldsc <- function(
    trait1_list,
    trait2_list,
    sumstats_dir,
    output_dir,
    ld_dir = "/media/DS1/zhangjh/202509_GWAS/Ref/LDSC/eur_w_ld_chr/"
) {
  library(GenomicSEM)
  # 确保输出目录存在
  if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  # 遍历 trait1
  for (trait1 in trait1_list) {
    trait1_file <- paste0(sumstats_dir, "/", trait1, ".sumstats.gz")
    # 遍历 trait2
    for (trait2 in trait2_list) {
      trait2_file <- paste0(sumstats_dir, "/", trait2, ".sumstats.gz")
      message("====== 正在运行 LDSC：", trait1, " vs ", trait2, " ======")
      # traits 文件名
      traits <- c(trait1_file, trait2_file)
      sample.prev <- c(NA, NA)
      population.prev <- c(NA, NA)
      trait.names <- c(trait1, trait2)
      # 运行 LDSC
      LDSCoutput <- GenomicSEM::ldsc(
        traits = traits,
        sample.prev = sample.prev,
        population.prev = population.prev,
        ld = ld_dir,
        wld = ld_dir,
        trait.names = trait.names,
        stand = TRUE
      )
      # 整理结果
      h2_trait1     <- as.numeric(LDSCoutput$S[1,1])
      h2_trait1_SE  <- sqrt(LDSCoutput$V[1,1])
      int_trait1    <- LDSCoutput$I[1,1]
      int_trait2    <- LDSCoutput$I[2,2]
      h2_trait2     <- as.numeric(LDSCoutput$S[2,2])
      h2_trait2_SE  <- sqrt(LDSCoutput$V[3,3])
      rgcov         <- as.numeric(LDSCoutput$S[1,2])
      k             <- nrow(LDSCoutput$S)
      SE_matrix     <- matrix(0, k, k)
      SE_matrix[lower.tri(SE_matrix, diag = TRUE)] <- sqrt(diag(LDSCoutput$V_Stand))
      rg_SE         <- SE_matrix[2,1]
      rg            <- rgcov / sqrt(h2_trait1 * h2_trait2)
      rg_P          <- 2 * pnorm(abs(rg / rg_SE), lower.tail = FALSE)
      rgcov_SE      <- sqrt(LDSCoutput$V[2,2])
      # 构建结果数据框
      ldsc_result <- data.frame(
        trait1,
        trait2,
        h2_trait1,
        h2_trait1_SE,
        int_trait1,
        int_trait2,
        h2_trait2,
        h2_trait2_SE,
        rgcov,
        rgcov_SE,
        rg,
        rg_SE,
        rg_P,
        stringsAsFactors = FALSE
      )
      # 赋列名
      colnames(ldsc_result) <- c(
        "trait1",
        "trait2",
        paste0("h2_", trait1),
        paste0("h2_", trait1, "_se"),
        paste0("intercept_", trait1),
        paste0("intercept_", trait2),
        paste0("h2_", trait2),
        paste0("h2_", trait2, "_se"),
        "rgcov",
        "rgcov_se",
        "rg",
        "rg_se",
        "rg_p"
      )
      # 输出文件
      out_file <- paste0(output_dir, "/", trait1, "&", trait2, "_ldsc_result.txt")
      write.table(
        ldsc_result,
        file = out_file,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE,
        col.names = TRUE
      )
      message("✔ 输出完成：", out_file)
    }
  }
}
#输入参数
trait1_list <- c("CKDIorII", "CKDIII", "CKDIV", "ESRD", "CRF", "CYSTKD", "NS", 
                 "eGFRcr_cys", "IgA", "HF", "AF", "IS", "IHD", "SBP", "DBP", 
                 "VTE", "AMI", "PVD", "HTN")
trait2_list <- trait1_list
#输入路径
sumstats_dir <- "/media/DS1/zhangjh/202601_GWAS/input/kidney_cardiovascular_LDSC"
#输出路径
output_dir   <- "/users/zhangjh/202601_GWAS/output/LDSC"
#运行
run_pairwise_ldsc(trait1_list, trait2_list, sumstats_dir, output_dir)

###5.合并所有LDSC结果
merge_ldsc_results <- function(
    input_dir,
    output_file = "AA_ldsc_result.csv",
    file_pattern = "\\.txt$"
) {
  # 进入结果目录
  setwd(input_dir)
  # 列出所有 LDSC 结果文件
  files <- list.files(pattern = file_pattern)
  if (length(files) == 0) {
    stop("❌ 未在目录中找到匹配的 LDSC 结果文件")
  }
  ldsc_result <- NULL
  for (i in seq_along(files)) {
    # 1️⃣ 读取单个 LDSC 结果
    ldsc_result0 <- read.table(
      files[i],
      header = TRUE,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    # 2️⃣ 统一第 3–8 列名称，确保 rbind 不错位
    colnames(ldsc_result0)[3:8] <- c(
      "h2_trait1", 
      "h2_trait1_se", 
      "intercept_trait1", 
      "intercept_trait2", 
      "h2_trait2", 
      "h2_trait2_se"
    )
    # 3️⃣ 从文件名中提取 trait1（保持你原来的逻辑）
    ldsc_result0$source_file <- strsplit(files[i], "_")[[1]][1]
    # 4️⃣ 合并
    ldsc_result <- rbind(ldsc_result, ldsc_result0)
  }
  # 5️⃣ 输出合并后的结果
  write.table(
    ldsc_result,
    file = output_file,
    sep = ",",          # 真正 CSV
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
  )
}
#填写输入输出并运行
merge_ldsc_results(
  input_dir  = "/users/zhangjh/202601_GWAS/output/LDSC/LDSC_result",
  output_file = "AA_ldsc_result.csv"
)

###6.所有疾病的遗传力（h2）柱状图（Bar Plot）
plot_ldsc_heritability <- function(
    input_file,
    output_pdf,
    kidney_traits,
    cardio_traits
) {
  # ===============================
  # 2️⃣ 读取数据并整理
  # ===============================
  library(data.table)
  library(ggplot2)
  library(ggbreak)
  library(dplyr)
  # 读取文件
  df <- fread(input_file)
  # 取前 19 行 & 第 2、7 列
  plot_df <- df[1:19, .(trait2, h2_trait2, h2_trait2_se)]
  # 添加疾病类别
  plot_df[, category := fifelse(
    trait2 %in% kidney_traits, "Kidney",
    fifelse(trait2 %in% cardio_traits, "Cardiovascular", NA)
  )]
  # 去掉不在两类中的（保险起见）
  plot_df <- plot_df[!is.na(category)]
  # ===============================
  # 排序
  # Cardiovascular：从小到大
  # Kidney：从大到小
  # ===============================
  plot_df <- plot_df %>%
    group_by(category) %>%
    arrange(
      category,
      if_else(category == "Cardiovascular",
              h2_trait2,        # 心血管：从小到大
              -h2_trait2        # 肾脏：从大到小
      )
    ) %>%
    ungroup() %>%
    mutate(
      trait2 = factor(trait2, levels = trait2)
    )
  # ===============================
  # 3️⃣ Nature 风格柱形图（核心）
  # ===============================
  p <- ggplot(plot_df, aes(x = trait2, y = h2_trait2, fill = category)) +
    # ① 绘制柱形图（geom_col）
    geom_col(
      width = 0.7,        # 柱子宽度（0~1），0.7 适中
      color = "black",    # 柱子边框颜色
      linewidth = 0.3     # 边框线粗细，ggplot2 新版用 linewidth 替代 size 避免警告
    ) +
    # ② 添加标准误（误差条）
    geom_errorbar(
      aes(
        ymin = h2_trait2 - h2_trait2_se,  # 下界 = 估计值 - 标准误
        ymax = h2_trait2 + h2_trait2_se   # 上界 = 估计值 + 标准误
      ),
      width = 0.35,       # 横线宽度（误差条横线）
      linewidth = 0.4     # 误差条线宽
    ) +
    # ③ 断轴（核心部分）
    scale_y_break(
      c(0.035, 0.1),     # 断轴区间：0.035 到 0.1 被跳过
      scales = 0.3,      # 断轴后上下部分比例（可调整柱子视觉比例）
      space = 0.39       # 调节中间空白大小（断轴间距）
    ) +
    # ④ 自定义填充颜色
    scale_fill_manual(
      values = c(
        "Kidney" = "#FC757B",           # 肾脏类为红色
        "Cardiovascular" = "#3C9BC9"    # 心血管类为蓝色
      )
    ) +
    # ⑤ 标签设置
    labs(
      x = NULL,  # x 轴标题为空
      y = expression("Heritability (h"^2*")"),  # y 轴用数学公式显示 h²
      fill = NULL  # 图例标题为空
    ) +
    
    # ⑥ Nature 风格主题
    theme_classic(base_size = 15) +  # 基础字体大小 15，背景白色、无网格线
    
    # ⑦ 进一步自定义主题元素
    theme(
      axis.text.x = element_text(
        angle = 45,   # x 轴标签旋转 45°，避免重叠
        hjust = 1,    # 水平对齐
        vjust = 1     # 垂直对齐
      ),
      axis.line = element_line(linewidth = 0.6),   # 坐标轴线粗细
      axis.ticks = element_line(linewidth = 0.6),  # 坐标轴刻度线粗细
      legend.position = "top",                     # 图例放在顶部
      legend.text = element_text(size = 15),       # 图例字体大小
      plot.margin = margin(10, 10, 10, 10)        # 图边距，上右下左各 10
    )
  # 打印绘图对象
  print(p)
  # ===============================
  # 保存为 PDF 文件
  # ===============================
  ggsave(
    filename = output_pdf,  # 输出文件名
    plot = p,               # ggplot 对象
    width = 8,              # PDF 宽度（inch）
    height = 6.5,           # PDF 高度（inch）
    units = "in",           # 单位为英寸
    dpi = 300               # 分辨率，PDF 矢量图不受影响，但用于兼容性
  )
  # 返回 ggplot 对象，方便后续继续修改
  invisible(p)
}
#设置可调参数
kidney_traits <- c("CKDIorII","CKDIII","CKDIV","ESRD","CRF",
                   "CYSTKD","NS","eGFRcr_cys","IgA")
cardio_traits <- c("HF","AF","IS","IHD","SBP",
                   "DBP","VTE","AMI","PVD","HTN")
#输入输出并运行
plot_ldsc_heritability(
  input_file = "/users/zhangjh/202601_GWAS/output/LDSC/LDSC_result/AA_ldsc_result.csv",
  output_pdf = "/users/zhangjh/202601_GWAS/output/LDSC/LDSC_plot/heritability_plot.pdf",
  kidney_traits = kidney_traits,
  cardio_traits = cardio_traits
)

###7.遗传相关性弦图(Genetic Correlation Circos Plot)
#' 绘制肾脏 vs 心血管疾病的遗传相关性弦图
#' @param ldsc_file  LDSC 结果 CSV 文件路径
#' @param output_pdf 输出 PDF 文件路径
#' @param kidney_traits 肾脏疾病 trait 列表（可自定义）
#' @param cardio_traits 心血管疾病 trait 列表（可自定义）
draw_genetic_correlation_circos <- function(
    ldsc_file,
    output_pdf,
    kidney_traits, 
    cardio_traits 
) {
  # 加载包
  library(data.table)
  library(dplyr)
  library(circlize)
  library(ComplexHeatmap)
  library(grid)
  
  #1️⃣ 读取 LDSC 结果并筛选
  df <- fread(ldsc_file)
  
  rg_sig <- df %>%
    filter(
      trait1 %in% cardio_traits,
      trait2 %in% kidney_traits,
      rg_p < 0.05
    ) %>%
    transmute(
      from  = trait2,   # Kidney
      to    = trait1,   # Cardiovascular
      value = rg        # 遗传相关性
    )
  
  #2️⃣ 彩带颜色 & 宽度映射
  col_fun <- colorRamp2(c(-1, 0, 1), c("#3182BD", "white", "#DE2D26"))  # 蓝 → 白 → 红
  link_col <- col_fun(rg_sig$value)   # 彩带颜色固定 -1~1
  link_lwd <- 0.6 + abs(rg_sig$value) * 4  # 彩带宽度 ∝ |rg|
  
  #3️⃣ 扇区顺序与颜色
  sector_order <- c(kidney_traits, cardio_traits)
  grid_col <- c(
    setNames(rep("#FC757B", length(kidney_traits)), kidney_traits),
    setNames(rep("#3C9BC9", length(cardio_traits)), cardio_traits)
  )
  
  #4️⃣ 输出 PDF
  pdf(output_pdf, width = 12, height = 10)
  
  circos.clear()
  circos.par(start.degree = 90)
  
  chordDiagram(
    x = rg_sig,
    order = sector_order,
    grid.col = grid_col,
    col = link_col,        # 彩带颜色
    link.lwd = link_lwd,   # 彩带宽度
    transparency = 0.25,
    annotationTrack = "grid",
    preAllocateTracks = list(track.height = 0.15)
  )
  
  #5️⃣ 添加扇区名称）
  circos.trackPlotRegion(
    track.index = 1,
    panel.fun = function(x, y) {
      sector = get.cell.meta.data("sector.index")
      xlim   = get.cell.meta.data("xlim")
      ylim   = get.cell.meta.data("ylim")
      
      circos.text(
        mean(xlim),
        ylim[1] + 0.1,
        sector,
        facing = "clockwise",
        niceFacing = TRUE,
        adj = c(0, 0.5),
        cex = 1.2   # 字体相对大小
      )
    },
    bg.border = NA
  )
  
  #6️⃣ 添加连续色标（heatmap-style bar）
  lgd <- Legend(col_fun = col_fun,
                title = expression(r[g]),
                direction = "vertical",
                at = c(-1, 0, 1),
                labels = c("-1","0","1"))
  
  # 放置在弦图右上角
  pushViewport(viewport(x=0.85, y=0.95, width=0.03, height=0.3, just=c("right","top")))
  grid.draw(lgd)
  popViewport()
  
  dev.off()
}

# 输入参数并运行
draw_genetic_correlation_circos(
  ldsc_file = "/users/zhangjh/202601_GWAS/output/LDSC/LDSC_result/AA_ldsc_result.csv",
  output_pdf = "/users/zhangjh/202601_GWAS/output/LDSC/LDSC_plot/Genetic_Correlation_Circos_plot.pdf",
  kidney_traits = c("CKDIorII","CKDIII","CKDIV","ESRD","CRF",
                    "CYSTKD","NS","eGFRcr_cys","IgA"),
  cardio_traits = c("HF","AF","IS","IHD","SBP",
                    "DBP","VTE","AMI","PVD","HTN")
)




























