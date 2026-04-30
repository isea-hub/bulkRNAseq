library(DESeq2)
library(dplyr)
library(readr)
library(readxl) 
library(tidyverse) 
# 读取数据
data <- read_excel("GSE276500_gene_expression.xlsx")
data <- as.data.frame(data)
head(data)
str(data)
summary(data)
#数据清洗只保留需要的数据,并设置行名
rownames(data) <- make.unique(data$gene_symbol)
data <- data[, -c(1:2,4:15)]
data <- data[, -c(7:8)]
#删除gene_symbol一列，便于后续分析
data$gene_symbol <- NULL  
head(data)

# 计算每行缺失值数量
row_zero_percentage <- rowMeans(data == 0)  
summary(row_zero_percentage)  # 查看 0 值比例的分布
# 计算每列缺失值数量
col_zero_percentage <- colMeans(data == 0)  
summary(col_zero_percentage)  # 查看 0 值比例的分布
total_zero_percentage <- mean(data == 0)  
paste("数据整体 0 值比例：", round(total_zero_percentage * 100, 2), "%")
library(tidyverse)
install.packages('DMwR2')
library(DMwR2)  # 用于 KNN 插补
# 对于超过20%缺失值的行删除
threshold <- 0.2  # 设定阈值：20%
# 仅保留 0 值比例 ≤ 20% 的基因
filtered_data <- data[row_zero_percentage <= threshold, ] 
paste("删除了", nrow(data) - nrow(filtered_data))
data <- filtered_data  # 更新数据
# 使用KNN进行插补
data[data == 0] <- NA  
# 使用 KNN 进行缺失值插补（k=5 表示使用 5个最近邻）
data_imputed <- knnImputation(data, k = 5)
# 检查是否仍然存在缺失值（NA）
sum(is.na(data_imputed))
# 更新数据
data <- data_imputed
#定义group
# 如果你的列名就是 xxx 等
control_samples <- c("read_count_NC1", "read_count_NC2", "read_count_NC3")
siSPIN1_samples <- c("read_count_SI1", "read_count_SI2", "read_count_SI3")

# 创建分组
group <- ifelse(colnames(data) %in% control_samples, 
                "control", "siSPIN1")
group <- factor(group, levels = c("control", "siSPIN1"))
head(group) 

library(ggplot2)
library(ggrepel)
library(ggsci)      # 期刊配色，如未安装则运行 install.packages("ggsci")
library(scales)     # 用于格式调整
# 1. 确认数据结构
cat("原始数据维度（行=基因，列=样本）:", dim(data), "\n")

# 2. 找出方差为零的基因（行）
row_var <- apply(data, 1, var, na.rm = TRUE)
cat("方差为零的基因数:", sum(row_var == 0, na.rm = TRUE), "\n")

# 3. 过滤
data_filt <- data[row_var > 0, ]

# 4. PCA
pca_result <- prcomp(t(data_filt), scale. = TRUE)

# 5. 查看结果
summary(pca_result)
# 计算方差解释比例
var_explained <- summary(pca_result)$importance[2, 1:2] * 100

# 提取 PC 得分
pca_df <- as.data.frame(pca_result$x[, 1:2])
pca_df$Sample <- rownames(pca_df)
pca_df$Group <- group  # 需提前定义，例如 c("Control","Control","siSPIN1","siSPIN1")


# ==================== 2. 提取 PCA 结果 ====================
# 假设 pca_result 已存在，group 已定义
#________________________________定义grop

head(group) 
#——————————————————————————————————————————————————————————————
pca_df <- as.data.frame(pca_result$x[, 1:2])   # 取前两个主成分
colnames(pca_df) <- c("PC1", "PC2")
pca_df$Sample <- rownames(pca_df)
pca_df$Group <- group

# 计算方差解释百分比
var_exp <- round(summary(pca_result)$importance[2, 1:2] * 100, 1)
pc1_lab <- paste0("PC1 (", var_exp[1], "%)")
pc2_lab <- paste0("PC2 (", var_exp[2], "%)")

# ==================== 3. 主题设置 ====================
# 定义统一的主题元素
theme_publication <- function(base_size = 12, base_family = "sans") {
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      # 面板边框（黑色实线，这是你要求的框线）
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
      # 网格线（轻微灰色，保留主网格线）
      panel.grid.major = element_line(color = "gray85", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      # 坐标轴线
      axis.line = element_line(color = "black", linewidth = 0.3),
      # 坐标轴刻度
      axis.ticks = element_line(color = "black", linewidth = 0.3),
      axis.ticks.length = unit(0.15, "cm"),
      # 坐标轴标题
      axis.title = element_text(size = base_size + 1, face = "bold"),
      axis.text = element_text(size = base_size, color = "black"),
      # 图例
      legend.position = "bottom",
      legend.title = element_text(size = base_size, face = "bold"),
      legend.text = element_text(size = base_size - 1),
      legend.key.size = unit(0.6, "cm"),
      legend.background = element_rect(fill = "white", color = "black", linewidth = 0.2),
      # 标题
      plot.title = element_text(size = base_size + 2, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = base_size, hjust = 0.5),
      plot.caption = element_text(size = base_size - 2, face = "italic", hjust = 0)
    )
}

# ==================== 4. 绘图 ====================
p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group, fill = Group)) +
  geom_point(size = 3.5, alpha = 0.9) +
  # 95% 置信椭圆（基于正态分布）
  stat_ellipse(aes(group = Group), type = "norm", level = 0.95, linewidth = 0.6) +
  # 样本标签（可选，若样本多可去掉 geom_text_repel）
  geom_text_repel(aes(label = Sample), 
                  size = 3, 
                  show.legend = FALSE,
                  box.padding = 0.4, 
                  point.padding = 0.2,
                  segment.color = "gray50",
                  segment.size = 0.2) +
  # 坐标轴标签（带上百分比）
  labs(x = pc1_lab, y = pc2_lab, 
       title = "PCA of expression profiles",
       subtitle = paste0("Groups: ", paste(levels(group), collapse = " vs ")),
       color = "Group", fill = "Group") +
  # 应用自定义主题
  theme_publication(base_size = 12) +
  # （可根据喜好更换）
  scale_color_npg() +
  scale_fill_npg() +
  # 确保坐标轴从原点开始（可选，有时不强制）
  coord_equal(ratio = 1)   # 使 PC1 和 PC2 比例一致，避免扭曲

# 显示图形
print(p)

# ==================== 5. 保存文件 ====================
# PDF（矢量格式，推荐用于期刊，可无限放大）
ggsave("hcc827sispin1PCA.pdf", plot = p, width = 6, height = 5, dpi = 300, device = "pdf")
# TIFF（位图，某些期刊要求，300 dpi 以上）
ggsave("hcc827sispin1PCA.tiff", plot = p, width = 6, height = 5, dpi = 300, compression = "lzw")
# PNG（备选，用于快速预览）
ggsave("hcc827sispin1PCA.png", plot = p, width = 6, height = 5, dpi = 300)


#edgeR分析edgeR 要求输入 原始 read counts（整数矩阵）
#且每个基因每个样本对应的值必须是未经过任何标准化或对数变换的原始计数
library(edgeR)
# 1) 创建 DGEList
dge <- DGEList(counts = data, group = group)
# 2) 计算归一化因子 (TMM)
dge <- calcNormFactors(dge)
# 3) 估计离散度
dge <- estimateDisp(dge)
# 4) 精确检验 (两组比较)
et <- exactTest(dge, pair = c("control", "siSPIN1"))
# 6) 提取结果 (FDR < 0.05 为显著)
results<- topTags(et, n = Inf, adjust.method = "BH")
results<- as.data.frame(results)
# 7) 查看显著基因数量
table(results$FDR < 0.05)
sig_genes<-results$FDR < 0.05
# 保存所有差异基因结果到 CSV
write.csv(results, file = "hcc827sispin1edgeR_DEG_results_all.csv")
# 保存显著基因列表
write.csv(sig_genes, file = "hcc827sispin1edgeR_DEG_significant.csv")

#绘制火山图
# 添加显著性标记
results$Significance <- "Not Significant"
results$Significance[results$logFC > 1 & results$FDR < 0.05] <- "Up"
results$Significance[results$logFC < -1 & results$FDR < 0.05] <- "Down"
# 指定文件路径
output_path <- "hcc827sispin1differential_expression_results.csv"
write.csv(results, file = output_path, row.names = TRUE)
library(ggrepel)
# 选择最显著上调和下调的基因
top_up <- results %>% 
  filter(Significance == "Up") %>% 
  top_n(10, logFC)

top_down <- results %>% 
  filter(Significance == "Down") %>% 
  top_n(-10, logFC)
# 合并需要标注的基因
top_genes <- bind_rows(top_up, top_down)


# 加载必要的包
library(ggplot2)
library(ggrepel)

# 假设 results 数据框包含 logFC, FDR, Significance 列
# 假设 top_genes 数据框是需要标注的显著差异基因（行名为基因名）
max_abs_logFC <- max(abs(results$logFC), na.rm = TRUE)
# 设定阈值
p_cutoff <- 0.05
lfc_cutoff <- 1

# 绘制火山图（带边框和优化标题）
p <-ggplot(results, aes(x = logFC, y = -log10(FDR), color = Significance)) +
  geom_point(alpha = 0.7, size = 1.5) +               # 调整点的大小和透明度
  scale_color_manual(
    values = c("Up" = "#CC3333", "Down" = "#0066CC", "Not Significant" = "gray80"),
    name = "Regulation"                               # 图例标题更清晰
  ) +
  labs(
    title = "Volcano Plot: HCC827 Response to SPIN1 Knockdown",
    subtitle = "Differential gene expression (GSE163110 dataset)",
    x = expression(log[2] ~ "Fold Change"),
    y = expression(-log[10] ~ "FDR"),
    caption = "Red: up-regulated, Blue: down-regulated, Gray: not significant"
  ) +
  theme_minimal(base_size = 12) +                     # 基础字号
  xlim(-max_abs_logFC, max_abs_logFC) +
  # 参考线
  geom_hline(yintercept = -log10(p_cutoff), linetype = "dashed", color = "gray30", linewidth = 0.4) +
  geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed", color = "gray30", linewidth = 0.4) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.3, linetype = "solid") +
  theme(
    # 添加黑色边框（框线）
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    # 调整网格线（可选，让边框更突出）
    panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    # 图例位置和背景
    legend.position = "bottom",
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.2),
    # 标题居中（可选）
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    plot.caption = element_text(hjust = 0, face = "italic", size = 9)
  ) +
  # 用带框标签标注 top_genes（默认黑色字体、白色填充）
  geom_label_repel(
    data = top_genes,
    aes(label = rownames(top_genes)),
    size = 3.5,
    color = "black",           # 文字颜色
    fill = "white",            # 背景填充色（形成“框框”）
    box.padding = 0.35,
    point.padding = 0.2,
    segment.color = "black",   # 连接线颜色
    segment.size = 0.3,
    label.r = 0.15,            # 圆角程度
    label.size = 0.25          # 标签边框粗细
  )
p
ggsave("HCC827sispin1volcano_plot.png", plot = p, width = 8, height = 6, dpi = 300, bg = "white")
ggsave("HCC827sispin1volcano_plot.pdf", plot = p, width = 8, height = 8, dpi = 300)

