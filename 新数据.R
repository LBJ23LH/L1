# 设置CRAN镜像以加速下载
options(repos = c(CRAN = "https://cloud.r-project.org/"))

# 检查并安装基础工具包
if (!require("devtools")) install.packages("devtools")
if (!require("remotes")) install.packages("remotes")

# 1. 安装CRAN包 -------------------------------------------------------------
cran_packages <- c(
  "Seurat", "data.table", "stringr", "tibble", "ggplot2", "patchwork",
  "dplyr", "tidyverse", "ggrepel", "pryr", "clustree", "jsonlite"
)

# 检查并安装缺失的CRAN包
new_cran <- cran_packages[!cran_packages %in% installed.packages()[,"Package"]]
if(length(new_cran)) install.packages(new_cran)

# 2. 安装Bioconductor包 -----------------------------------------------------
if (!require("BiocManager")) install.packages("BiocManager")
bioc_packages <- c("SingleR", "ComplexHeatmap", "clusterProfiler", "org.Hs.eg.db", "enrichplot")

# 检查并安装缺失的Bioconductor包
new_bioc <- bioc_packages[!bioc_packages %in% installed.packages()[,"Package"]]
if(length(new_bioc)) BiocManager::install(new_bioc)

# 3. 安装GitHub包 ----------------------------------------------------------
# 检查并安装DoubletFinder
if (!require("DoubletFinder")) {
  remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
}

# 检查并安装ROGUE
if (!require("ROGUE")) {
  devtools::install_github("PaulingLiu/ROGUE")
}

# 检查并安装CellChat
if (!require("CellChat")) {
  devtools::install_github("sqjin/CellChat")
}

# 检查并安装presto
if (!require("presto")) {
  devtools::install_github("immunogenomics/presto")
}

# 4. 加载所有包 ------------------------------------------------------------
suppressPackageStartupMessages({
  library(Seurat)        # 单细胞分析核心工具
  library(data.table)    # 高效数据处理
  library(stringr)       # 字符串处理
  library(tibble)        # 数据框操作
  library(ggplot2)       # 可视化
  library(patchwork)     # 拼图工具
  library(dplyr)         # 数据操作
  library(tidyverse)     # 数据科学生态
  library(ggrepel)       # 标签防重叠
  library(clustree)      # 聚类树可视化
  library(harmony)       # 批次校正
  library(SingleR)       # 细胞类型注释
  library(CellChat)      # 细胞通讯分析
  library(ROGUE)         # 细胞纯度评估
  library(presto)        # 快速差异分析
  library(org.Hs.eg.db)  # 人类基因注释
  library(clusterProfiler) # 富集分析
  library(enrichplot)    # 富集结果可视化
})

# 验证Harmony加载
if (!require("harmony")) install.packages("harmony"); library(harmony)

# 内存管理
cat("Current memory usage:", pryr::mem_used()/1024^2, "MB\n")





# 1. 数据导入与整合

# 设置工作目录
setwd("D:/lhRes/code/newCode/")

# 定义数据路径结构
data_dir <- "data/"
groups <- list(
  Control = "Control",  # 对照组
  Disease = c("Subacute", "Moderate", "Severe")  # 患病组及其亚型
)

# 创建空列表存储所有样本的Seurat对象
all_seurat_list <- list()

# 遍历每个分组（Control + Disease各亚型）
for (group in names(groups)) {
  if (group == "Control") {
    # 处理对照组
    control_path <- file.path(data_dir, groups[[group]])
    samples <- list.dirs(control_path, full.names = TRUE, recursive = FALSE)
    
    for (sample_path in samples) {
      sample_name <- basename(sample_path)
      
      # 读取10x数据并创建Seurat对象
      seurat_data <- Read10X(data.dir = sample_path)
      seurat_obj <- CreateSeuratObject(
        counts = seurat_data,
        project = sample_name,
        min.features = 200,
        min.cells = 3
      )
      
      # 添加metadata：分组标签 + 疾病亚型
      seurat_obj$group <- "Control"        # 组别标签（Disease/Control）
      seurat_obj$disease_stage <- "None"   # 疾病亚型（仅患病组有）
      
      # 重命名细胞ID（添加样本名前缀）
      seurat_obj <- RenameCells(seurat_obj, new.names = paste0(sample_name, "_", colnames(seurat_obj)))
      
      all_seurat_list[[length(all_seurat_list) + 1]] <- seurat_obj
    }
    
  } else {
    # 处理患病组各亚型（Subacute/Moderate/Severe）
    for (subtype in groups[[group]]) {
      subtype_path <- file.path(data_dir, subtype)
      samples <- list.dirs(subtype_path, full.names = TRUE, recursive = FALSE)
      
      for (sample_path in samples) {
        sample_name <- basename(sample_path)
        
        # 读取10x数据并创建Seurat对象
        seurat_data <- Read10X(data.dir = sample_path)
        seurat_obj <- CreateSeuratObject(
          counts = seurat_data,
          project = sample_name,
          min.features = 200,
          min.cells = 3
        )
        
        # 添加metadata：分组标签 + 疾病亚型
        seurat_obj$group <- "Disease"          # 组别标签
        seurat_obj$disease_stage <- subtype    # 疾病亚型
        
        # 重命名细胞ID（添加样本名前缀）
        seurat_obj <- RenameCells(seurat_obj, new.names = paste0(sample_name, "_", colnames(seurat_obj)))
        
        all_seurat_list[[length(all_seurat_list) + 1]] <- seurat_obj
      }
    }
  }
}

# 合并所有样本为一个Seurat对象
combined_seurat <- merge(all_seurat_list[[1]], y = all_seurat_list[-1])
combined_seurat <- JoinLayers(combined_seurat)

# 检查metadata是否正确
table(combined_seurat$group)      # 统计Disease vs Control细胞数
table(combined_seurat$disease_stage)  # 统计各疾病亚型细胞数

# 保存
setwd("D:/lhRes/code/newCode/results/1. 数据导入与整合/")
saveRDS(combined_seurat, file = "combined_seurat.rds")

# 加载
combined_seurat <- readRDS("combined_seurat.rds")





# 2. 质量控制与过滤

# 设置结果目录
setwd("D:/lhRes/code/newCode/results/2. 质量控制与过滤/")

# 1. 查看样本分布 -----------------------------------------------------------------
# 获取样本分布情况
sample_distribution <- table(combined_seurat@meta.data$orig.ident)

# 将样本分布转换为数据框
sample_distribution_df <- as.data.frame(sample_distribution)

# 设置列名
colnames(sample_distribution_df) <- c("Sample", "Cell_Count")

# 设置保存路径
output_file <- "sample_distribution.csv"

# 将数据框存储为CSV文件
write.csv(sample_distribution_df, file = output_file, row.names = FALSE)

# 输出质控前细胞数
message("质控前细胞数: ", ncol(combined_seurat))

# 打印样本分布表格
print(sample_distribution_df)



# 2. 计算质控指标 -----------------------------------------------------------------
# 线粒体基因比例
combined_seurat[["percent.mt"]] <- PercentageFeatureSet(combined_seurat, pattern = "^MT-")

# 血红蛋白基因比例（检查基因是否存在于数据中）
HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB.genes <- CaseMatch(HB.genes, rownames(combined_seurat))  # 确保大小写匹配
if(length(HB.genes) > 0) {
  combined_seurat[["percent.HB"]] <- PercentageFeatureSet(combined_seurat, features = HB.genes)
} else {
  warning("No hemoglobin genes found in the dataset.")
}

# 3. 可视化质控指标相关性 ----------------------------------------------------------
# 生成散点图并直接保存到当前目录
scatter_plots <- list()

# 散点图1：线粒体 vs UMI总数
scatter_plots[[1]] <- FeatureScatter(combined_seurat, "nCount_RNA", "percent.mt", 
                                     group.by = "orig.ident") +
  ggtitle("Mitochondrial vs UMI Count")
cor(combined_seurat$percent.mt, combined_seurat$nCount_RNA, method = "pearson")


# 散点图2：基因数 vs UMI总数
scatter_plots[[2]] <- FeatureScatter(combined_seurat, "nCount_RNA", "nFeature_RNA", 
                                     group.by = "orig.ident") +
  ggtitle("Gene Count vs UMI Count")
cor(combined_seurat$nFeature_RNA, combined_seurat$nCount_RNA, method = "pearson")

# 可视化
combined_scatter <- scatter_plots[[1]] | scatter_plots[[2]]  # 使用patchwork横向排列
ggsave(
  filename = "scatter_QC.pdf",
  plot = combined_scatter,
  width = 24,  # 宽度调大以适应两图横向排列
  height = 12
)

# 4. 质控前指标分布（按样本分组） -------------------------------------------------
theme.set2 <- theme(axis.title.x = element_blank())
plot.features <- c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.HB")

# 定义每个指标的最大值（根据质控前数据动态计算）
y_limits <- list(
  nFeature_RNA = c(0, max(combined_seurat$nFeature_RNA)),
  nCount_RNA = c(0, max(combined_seurat$nCount_RNA)),
  percent.mt = c(0, max(combined_seurat$percent.mt)),
  percent.HB = c(0, ifelse("percent.HB" %in% colnames(combined_seurat@meta.data), 
                           max(combined_seurat$percent.HB), 1))
)

# 绘图
plot_violin <- function(seurat_obj, y_limits) {
  wrap_plots(
    lapply(plot.features, function(feature) {
      VlnPlot(seurat_obj, group.by = "orig.ident", features = feature, pt.size = 0.001) + 
        theme.set2 + 
        NoLegend() +
        ylim(y_limits[[feature]])  # 强制统一纵坐标范围
    }),
    nrow = 4
  )
}

# 生成质控前小提琴图
violin_before <- plot_violin(combined_seurat, y_limits)

# 保存到当前目录
ggsave("vlnplot_before_qc.pdf", plot = violin_before, width = 24, height = 24)

# 5. 动态设定质控阈值 -------------------------------------------------------------
# 检查分位数（根据合并数据调整阈值）
quantile(combined_seurat$nFeature_RNA, c(0.01, 0.99))  # 基因数范围
quantile(combined_seurat$nCount_RNA, c(0.01, 0.99))     # UMI总数范围
quantile(combined_seurat$percent.mt, c(0.95, 1))        # 线粒体基因上限

# 设定阈值（根据分位数结果调整）
qc_thresholds <- list(
  minGene = 250,
  maxGene = 4500,
  minUMI = 650,
  maxUMI = 22000,
  pctMT = 25,
  pctHB = 1
)

# 执行质控过滤
combined_seurat <- subset(
  combined_seurat,
  subset = nFeature_RNA > qc_thresholds$minGene &
    nFeature_RNA < qc_thresholds$maxGene &
    nCount_RNA > qc_thresholds$minUMI &
    nCount_RNA < qc_thresholds$maxUMI &
    percent.mt < qc_thresholds$pctMT
)
if("percent.HB" %in% colnames(combined_seurat@meta.data)) {
  combined_seurat <- subset(combined_seurat, percent.HB < qc_thresholds$pctHB)
}

# 6. 质控后验证 -------------------------------------------------------------------
# 生成质控后小提琴图
violin_after <- plot_violin(combined_seurat, y_limits)

# 保存到当前目录
ggsave("vlnplot_after_qc.pdf", plot = violin_after, width = 24, height = 24)
message("质控后细胞数: ", ncol(combined_seurat))

# 7. 保存质控后的对象 -------------------------------------------------------------
saveRDS(combined_seurat, "combined_seurat_qc.rds")

# 加载
combined_seurat <- readRDS("combined_seurat_qc.rds")





# 3. 数据标准化 -----------------------------------------------------------

# 设置结果目录
setwd("D:/lhRes/code/newCode/results/3. 数据标准化/")

# 1. 数据标准化 -----------------------------------------------------------------
# 使用LogNormalize方法标准化（结果存储在RNA assay的data层）
combined_seurat <- NormalizeData(
  combined_seurat,
  normalization.method = "LogNormalize",
  scale.factor = 10000  # 根据文库大小标准化
)

# 2. 检测高变基因 ---------------------------------------------------------------
# 识别2000个高变基因用于下游分析
combined_seurat <- FindVariableFeatures(
  combined_seurat,
  selection.method = "vst",  # 方差稳定转换方法
  nfeatures = 2000,          # 选择top 2000高变基因
  verbose = FALSE            # 关闭冗长输出
)

# 可视化高变基因
top20 <- head(VariableFeatures(combined_seurat), 20)
var_plot <- LabelPoints(
  plot = VariableFeaturePlot(combined_seurat), 
  points = top20, 
  repel = TRUE
)
ggsave("variable_features.pdf", plot = var_plot, width = 12, height = 8)

# 3. 细胞周期评分 --------------------------------------------------------------
# 加载细胞周期基因集
cc_genes <- list(
  S.genes = cc.genes.updated.2019$s.genes,
  G2M.genes = cc.genes.updated.2019$g2m.genes
)

# 计算细胞周期评分
combined_seurat <- CellCycleScoring(
  combined_seurat,
  s.features = cc_genes$S.genes,
  g2m.features = cc_genes$G2M.genes,
  set.ident = FALSE
)

# 4. 数据缩放与回归 -------------------------------------------------------------
# 缩放数据并回归技术性因素（结果存储在RNA assay的scale.data层）
combined_seurat <- ScaleData(
  combined_seurat,
  vars.to.regress = c("S.Score", "G2M.Score", "percent.mt", "nCount_RNA"), # 回归细胞周期和质控指标
  verbose = FALSE
)

# 5. 保存处理后的对象 -----------------------------------------------------------
saveRDS(combined_seurat, "combined_seurat_normalized.rds")





# 4. 降维聚类与批次校正 ---------------------------------------------------------

# 设置结果目录
setwd("D:/lhRes/code/newCode/results/4. 降维聚类与批次校正/")

# 1. 加载标准化后的数据 --------------------------------------------------------
combined_seurat <- readRDS("../3. 数据标准化/combined_seurat_normalized.rds")

# 2. PCA降维 -------------------------------------------------------------------
combined_seurat <- RunPCA(
  combined_seurat,
  features = VariableFeatures(combined_seurat),
  verbose = FALSE,
  npcs = 50  # 计算前50个主成分
)

# 可视化PCA相关图形
pca_plots <- list()

# PCA散点图
pdf("PCA_scatter.pdf", width=8, height=6)
print(DimPlot(combined_seurat, reduction = "pca", group.by = "orig.ident") +
        ggtitle("PCA - Before Harmony"))
dev.off()

# 主成分负荷图
pdf("PC_loadings.pdf", width=10, height=8)
print(VizDimLoadings(combined_seurat, dims = 1:4, reduction = "pca") +
        theme(axis.text = element_text(size=6)))
dev.off()

# 主成分热图
pdf("PC_heatmap.pdf", width=12, height=9) 
DimHeatmap(
  combined_seurat,
  dims = 1:4,
  cells = 500,
  nfeatures = 30,
  balanced = TRUE
)
dev.off()

# 3. 确定主成分数量 ------------------------------------------------------------
# 计算主成分贡献率
pct <- combined_seurat[["pca"]]@stdev^2 / sum(combined_seurat[["pca"]]@stdev^2) * 100
cum_pct <- cumsum(pct)

# 自动选择主成分数（累计贡献>85% 且 相邻主成分差异<2%）
selected_pcs <- which(cum_pct > 85 & abs(diff(c(cum_pct,100))) < 2)[1]
message("Selected PCs: ", selected_pcs)

# 肘部法则图
elbow_plot <- ElbowPlot(combined_seurat, ndims = 50) +
  geom_vline(xintercept = selected_pcs, linetype = "dashed", color = "red") +
  ggtitle("Elbow Plot for PCA")
ggsave("elbow_plot.pdf", plot = elbow_plot, width = 8, height = 6)

# 4. 批次校正（Harmony） -------------------------------------------------------
# 运行Harmony
combined_seurat <- RunHarmony(
  object = combined_seurat,
  group.by.vars = "orig.ident",
  reduction.use = "pca",
  dims.use = 1:selected_pcs,
  theta = 2,
  lambda = 0.5,
  max_iter = 30,
  verbose = FALSE
)

# 5. 非线性降维 ---------------------------------------------------------------
# 运行UMAP和tSNE
combined_seurat <- RunUMAP(
  combined_seurat,
  reduction = "harmony",
  dims = 1:selected_pcs,
  umap.method = "uwot",
  min.dist = 0.3
) %>% RunTSNE(
  reduction = "harmony",
  dims = 1:selected_pcs,
  perplexity = 30
)

# 6. 聚类分析 ------------------------------------------------------------------
# 测试不同分辨率
resolutions <- seq(0.2, 2.0, by = 0.2)
combined_seurat <- FindNeighbors(
  combined_seurat,
  reduction = "harmony",
  dims = 1:selected_pcs
)

for (res in resolutions) {
  combined_seurat <- FindClusters(combined_seurat, resolution = res)
}

# 聚类树可视化
clustree_plot <- clustree(combined_seurat, prefix = "RNA_snn_res.") +
  scale_color_brewer(palette = "Set2") +
  ggtitle("Cluster Resolution Selection")
ggsave("clustree_plot.pdf", plot = clustree_plot, width = 14, height = 10)

# 确定最终分辨率
combined_seurat <- FindClusters(combined_seurat, resolution = 1.0)

# 7. 可视化分析 ----------------------------------------------------------------
dim_plots <- list()

# UMAP可视化
dim_plots[[1]] <- DimPlot(combined_seurat, reduction = "umap", label = TRUE) +
  ggtitle("UMAP - Cell Clusters")

dim_plots[[2]] <- DimPlot(combined_seurat, reduction = "umap", group.by = "orig.ident") +
  ggtitle("UMAP - Batch Distribution")

# tSNE可视化
dim_plots[[3]] <- DimPlot(combined_seurat, reduction = "tsne", label = TRUE) +
  ggtitle("tSNE - Cell Clusters")

dim_plots[[4]] <- DimPlot(combined_seurat, reduction = "tsne", group.by = "orig.ident") +
  ggtitle("tSNE - Batch Distribution")

# 保存可视化结果
ggsave("dim_reduction_plots.pdf", 
       plot = wrap_plots(dim_plots, ncol = 2), 
       width = 18, height = 12)

# 8. 保存最终结果 -------------------------------------------------------------
# 添加重要元数据注释
combined_seurat$cell_id <- colnames(combined_seurat)
combined_seurat$combined_cluster <- paste0("C", combined_seurat$seurat_clusters)

# 输出关键元数据
write.csv(combined_seurat@meta.data, "cell_metadata.csv", row.names = FALSE)
saveRDS(combined_seurat, "combined_seurat_final.rds")





# 5. 细胞注释 -------------------------------------------------------------------

# 设置结果目录
setwd("D:/lhRes/code/newCode/results/5. 细胞注释/")

# 1. 加载参考数据库（确保路径正确）
ref_path <- "D:/lhRes/code/newCode/data/ref_Human_all.RData"
if (file.exists(ref_path)) {
  load(ref_path)
} else {
  stop("参考数据库文件不存在，请检查路径: ", ref_path)
}

# 2. 准备输入数据 -------------------------------------------------------------
# 获取标准化后的表达矩阵（logcounts）
testdata <- GetAssayData(
  object = combined_seurat,
  assay = "RNA",
  layer = "data"  # 使用标准化后的数据
)

# 获取聚类信息（使用resolution = 1.0的聚类结果）
clusters <- combined_seurat$seurat_clusters

# 3. 运行SingleR注释 ---------------------------------------------------------
# 使用主细胞类型标签（label.main）
cellpred <- SingleR(
  test = testdata,
  ref = ref_Human_all,
  clusters = clusters,
  labels = ref_Human_all$label.main,  # 使用主分类标签
  assay.type.test = "logcounts",
  assay.type.ref = "logcounts"
)

# 4. 整合注释结果到metadata ---------------------------------------------------
# 生成注释数据框
celltype <- data.frame(
  ClusterID = rownames(cellpred),
  celltype = cellpred$labels,
  stringsAsFactors = FALSE
)

# 添加注释到Seurat对象
combined_seurat$SingleR_main <- "Unassigned"
for (i in seq_len(nrow(celltype))) {
  cluster_cells <- which(combined_seurat$seurat_clusters == celltype$ClusterID[i])
  combined_seurat$SingleR_main[cluster_cells] <- celltype$celltype[i]
}

# 5. 可视化分析 ---------------------------------------------------------------
# 注释热图
p_heatmap <- plotScoreHeatmap(cellpred)
ggsave("SingleR_score_heatmap.pdf", p_heatmap, width = 14, height = 8)

# UMAP可视化注释结果
p_umap <- DimPlot(
  combined_seurat,
  group.by = "SingleR_main",
  label = TRUE,
  repel = TRUE,
  reduction = "umap"
) + ggtitle("Cell Type Annotation (Main)")
ggsave("SingleR_umap.pdf", p_umap, width = 12, height = 8)

# 7. 高级注释 ----------------------------------------------------------
# 使用细粒度注释（label.fine）
cellpred_fine <- SingleR(
  test = testdata,
  ref = ref_Human_all,
  clusters = clusters,
  labels = ref_Human_all$label.fine,  # 使用精细分类标签
  assay.type.test = "logcounts",
  assay.type.ref = "logcounts"
)

# 整合精细注释
celltype_fine <- data.frame(
  ClusterID = rownames(cellpred_fine),
  celltype = cellpred_fine$labels,
  stringsAsFactors = FALSE
)

combined_seurat$SingleR_fine <- "Unassigned"
for (i in seq_len(nrow(celltype_fine))) {
  cluster_cells <- which(combined_seurat$seurat_clusters == celltype_fine$ClusterID[i])
  combined_seurat$SingleR_fine[cluster_cells] <- celltype_fine$celltype[i]
}

# UMAP可视化
p_umap_fine <- DimPlot(
  combined_seurat,
  group.by = "SingleR_fine",
  label = TRUE,
  repel = TRUE,              # 避免标签重叠
  label.size = 3,            # 调小标签字号
  reduction = "umap",
  shuffle = TRUE             # 随机打点提升可视化效果
) + 
  ggtitle("Cell Type Annotation (Fine)") +
  theme(legend.position = "right", legend.text = element_text(size=7))
ggsave("SingleR_umap_fine.pdf", p_umap_fine, width = 14, height = 10)

# 保存最终对象
saveRDS(combined_seurat, "combined_seurat_annotated.rds")






# 6. 组间差异分析 ----------------------------------------------------------------

# 设置结果目录
setwd("D:/lhRes/code/newCode/results/6. 组间差异分析/")

# 1. 准备分析数据 -------------------------------------------------------------
# 加载注释后的数据
combined_seurat <- readRDS("../5. 细胞注释/combined_seurat_annotated.rds")

# 确认分组信息
table(combined_seurat$group)  # 应显示Control和Disease的细胞数
table(combined_seurat$disease_stage) # 查看疾病亚型分布

# 2. 整体差异分析（所有细胞）-------------------------------------------------

# Wilcoxon秩和检验
markers_all <- FindMarkers(
  combined_seurat,
  ident.1 = "Disease",
  ident.2 = "Control",
  group.by = "group",
  test.use = "wilcox",  # 默认方法
  logfc.threshold = 0.25,
  min.pct = 0.1,
  only.pos = FALSE
)

# 保存结果
write.csv(markers_all, "DEGs_all_cells_wilcox.csv")

# 添加基因名称列
deg_all_filtered$gene <- rownames(deg_all_filtered)

# 火山图 - 高亮基因
highlight_genes <- c("ISG15", "IFIT1", "OAS1", "MX1", "IL1B", "CXCL10", "TNF", "IRF7")  # 需要高亮的基因

# 绘制火山图
volcano_plot <- ggplot(deg_all_filtered, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(aes(color = ifelse(avg_log2FC > 0, "Up", "Down")), alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Up" = "red", "Down" = "blue")) +
  # 高亮显示特定基因
  geom_text_repel(
    data = filter(deg_all_filtered, gene %in% highlight_genes),
    aes(label = gene),
    size = 3, max.overlaps = 10, box.padding = 0.5
  ) +
  labs(title = "Volcano Plot for All Cells (Highlighted Genes)",
       x = "Log2 Fold Change",
       y = "-log10 Adjusted P-value") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none")

# 保存火山图为PDF
ggsave("Volcano_Plot_all_cells.pdf", plot = volcano_plot, width = 8, height = 6)

colnames(deg_all_filtered)


# 3. 细胞类型特异性差异分析 --------------------------------------------------
# 获取所有细胞类型
celltypes <- unique(combined_seurat$SingleR_fine)
message("需要分析的细胞类型：", paste(celltypes, collapse = ", "))

# 并行计算加速
library(future)
plan("multisession", workers = 4)  # 根据CPU核心数调整

# 循环分析每个细胞类型
deg_results <- list()
for (ct in celltypes) {
  message("\n正在分析：", ct)
  subset_obj <- subset(combined_seurat, SingleR_fine == ct)
  
  # 执行差异分析
  markers <- FindMarkers(
    subset_obj,
    ident.1 = "Disease",
    ident.2 = "Control",
    group.by = "group",
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1
  )
  
  # 添加基因符号和注释
  markers$gene <- rownames(markers)
  markers$celltype <- ct
  markers$direction <- ifelse(markers$avg_log2FC > 0, "Up", "Down")
  
  deg_results[[ct]] <- markers
  write.csv(markers, file = paste0("DEGs_", gsub("/", "_", ct), ".csv"))
}

# 合并所有结果
combined_degs <- bind_rows(deg_results)
write.csv(combined_degs, "Combined_DEGs_all_celltypes.csv")

# 4. 可视化分析 --------------------------------------------------------------
# 高亮基因
highlight_genes <- c("ISG15", "IFIT1", "OAS1", "MX1", "IL1B", "CXCL10", "TNF", "IRF7")

# 创建列表存放图形
volcano_list <- list()

# 限制每个图显示的最多基因数
min_genes <- 20

# 遍历每个细胞类型，生成火山图并存入列表
for (ct in unique(deg_filtered$celltype)) {
  ct_deg <- filter(deg_filtered, celltype == ct)
  if (nrow(ct_deg) < min_genes) next
  
  p <- ggplot(ct_deg, aes(avg_log2FC, -log10(p_val_adj))) +
    geom_point(aes(color = direction), alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c(Up = "red", Down = "blue")) +
    geom_text_repel(
      data = filter(ct_deg, gene %in% highlight_genes),
      aes(label = gene),
      size = 2.5, max.overlaps = 10
    ) +
    labs(
      title = ct,
      x = "Log2 Fold Change",
      y = "-log10 Adjusted P-value"
    ) +
    theme_minimal(base_size = 10) +
    theme(legend.position = "none")
  
  volcano_list[[ct]] <- p
}

# 使用 patchwork 组合所有图
combined_plot <- wrap_plots(volcano_list, ncol = 3) +
  plot_annotation(title = "Volcano Plots of DEGs across Cell Types")

# 保存为PDF
ggsave("All_Volcano_Plots_combined.pdf", plot = combined_plot, width = 18, height = 12)

# 5. 保存最终结果 ------------------------------------------------------------
saveRDS(combined_seurat, "combined_seurat_with_DEGs.rds")





# 7. 功能富集分析 -------------------------------------------------------------

# 设置结果目录
setwd("D:/lhRes/code/newCode/results/7. 功能富集分析/")

# 加载差异分析结果
combined_degs <- read.csv("../6. 组间差异分析/Combined_DEGs_all_celltypes.csv")

# 创建分析目录结构
dir.create("GO", showWarnings = FALSE)
dir.create("KEGG", showWarnings = FALSE)

# 定义富集分析函数
run_enrichment <- function(genes, celltype, direction) {
  # 1. ID转换（SYMBOL -> ENTREZID）-----------------------------------------
  ids <- clusterProfiler::bitr(
    genes,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )
  
  # 若基因映射失败，返回空结果
  if (nrow(ids) == 0) {
    warning("No genes mapped to ENTREZID.")
    return(list(go = NULL, kegg = NULL))
  }
  
  # 2. GO富集分析 ---------------------------------------------------------
  go_res <- enrichGO(
    gene = ids$ENTREZID,
    OrgDb = org.Hs.eg.db,
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    readable = TRUE  # 将ENTREZID转换为基因名显示
  )
  
  # 3. KEGG富集分析（修复ENTREZID映射）-------------------------------------
  kegg_res <- enrichKEGG(
    gene = ids$ENTREZID,
    organism = 'hsa',
    keyType = 'ncbi-geneid',  # 明确指定输入ID类型为ENTREZID
    pvalueCutoff = 0.05
  )
  
  # 4. 计算通路相似性矩阵（修复emapplot报错）-----------------------------
  if (!is.null(kegg_res) && nrow(kegg_res) > 0) {
    kegg_res <- pairwise_termsim(kegg_res)  # 生成termsim矩阵
  }
  
  # 5. 保存结果 ----------------------------------------------------------
  output_dir_go <- "GO"
  output_dir_kegg <- "KEGG"
  
  safe_ct <- gsub("[:/ ]", "_", ct)  # 清理路径非法字符
  
  if (!is.null(go_res) && nrow(go_res@result) > 0) {
    write.csv(go_res@result, file.path(output_dir_go, paste0(safe_ct, "_", direction, "_GO.csv")))
  }
  if (!is.null(kegg_res) && nrow(kegg_res@result) > 0) {
    write.csv(kegg_res@result, file.path(output_dir_kegg, paste0(safe_ct, "_", direction, "_KEGG.csv")))
  }
  
  
  # 6. 返回结果 ----------------------------------------------------------
  return(list(go = go_res, kegg = kegg_res))
}

# 执行富集分析
enrichment_results <- list()

# 按细胞类型和方向分组处理
deg_groups <- combined_degs %>%
  filter(p_val_adj < 0.05) %>%
  group_by(celltype, direction) %>%
  group_split()

for (group in deg_groups) {
  ct <- unique(group$celltype)
  dir <- unique(group$direction)
  message("Processing: ", ct, " - ", dir)
  
  tryCatch({
    res <- run_enrichment(group$gene, ct, dir)
    enrichment_results[[paste(ct, dir, sep = "_")]] <- res
  }, error = function(e) {
    message("Error in ", ct, " ", dir, ": ", e$message)
  })
}

# 可视化分析，绘制所有细胞类型方向的富集结果
pdf("GO/All_GO_Plots.pdf", width = 10, height = 8)

for (key in names(enrichment_results)) {
  res <- enrichment_results[[key]]
  dir <- sub(".*_([A-Za-z]+)$", "\\1", key)
  ct <- sub("_Up$|_Down$", "", key)
  
  if (!is.null(res$go) && nrow(res$go@result) > 0) {
    p_go <- dotplot(res$go, showCategory = min(15, nrow(res$go@result))) +
      ggtitle(paste("GO Enrichment:", ct, "-", dir))
    print(p_go)
  }
}

dev.off()

pdf("KEGG/All_KEGG_Plots.pdf", width = 10, height = 8)

for (key in names(enrichment_results)) {
  res <- enrichment_results[[key]]
  dir <- sub(".*_([A-Za-z]+)$", "\\1", key)
  ct <- sub("_Up$|_Down$", "", key)
  
  if (!is.null(res$kegg) && nrow(res$kegg@result) > 1) {
    tryCatch({
      p_kegg <- emapplot(res$kegg, showCategory = min(15, nrow(res$kegg@result))) +
        ggtitle(paste("KEGG Enrichment:", ct, "-", dir))
      print(p_kegg)
    }, error = function(e) {
      message("跳过 ", ct, " ", dir, " 的 KEGG 图：", e$message)
    })
  }
}

dev.off()





# 8. 细胞间通讯分析 -----------------------------------------------------------
# 【设置工作目录和加载数据】
setwd("D:/lhRes/code/newCode/results/")

combined_seurat <- readRDS("5. 细胞注释/combined_seurat_annotated.rds")

# 创建目录
com_dir <- "8. 细胞间通讯分析"
dir.create(file.path(com_dir, "疾病组"),    recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(com_dir, "对照组"),    recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(com_dir, "比较分析"), recursive = TRUE, showWarnings = FALSE)

# 加载 CellChat 数据库
data("CellChatDB.human")
data("PPI.human")
CellChatDB <- CellChatDB.human
ppi        <- PPI.human
options(future.globals.maxSize = 8 * 1024^3)

# 【分组】
pbmc_disease <- subset(combined_seurat, subset = group == "Disease")
pbmc_control <- subset(combined_seurat, subset = group == "Control")
data.input_disease <- GetAssayData(pbmc_disease, assay = "RNA", layer = "data")
data.input_control <- GetAssayData(pbmc_control, assay = "RNA", layer = "data")

# 【构建CellChat对象 - 疾病组】
cellchat_disease <- createCellChat(object = data.input_disease, meta = pbmc_disease@meta.data, group.by = "SingleR_fine")
cellchat_disease@DB <- CellChatDB
saveRDS(cellchat_disease, file.path(com_dir, "疾病组", "cellchat_disease_step1.rds"))
cellchat_disease <- readRDS(file.path(com_dir, "疾病组", "cellchat_disease_step1.rds"))

cellchat_disease <- subsetData(cellchat_disease)
cellchat_disease <- identifyOverExpressedGenes(cellchat_disease)
cellchat_disease <- identifyOverExpressedInteractions(cellchat_disease)
cellchat_disease <- projectData(cellchat_disease, ppi)
cellchat_disease <- computeCommunProb(cellchat_disease, type = "truncatedMean", trim = 0.1, raw.use = FALSE)
cellchat_disease <- filterCommunication(cellchat_disease, min.cells = 5)
cellchat_disease <- computeCommunProbPathway(cellchat_disease)
cellchat_disease <- netAnalysis_computeCentrality(cellchat_disease)
cellchat_disease <- aggregateNet(cellchat_disease)
saveRDS(cellchat_disease, file.path(com_dir, "疾病组", "cellchat_disease_final.rds"))
cellchat_disease <- readRDS(file.path(com_dir, "疾病组", "cellchat_disease_final.rds"))

# 【构建CellChat对象 - 对照组】
cellchat_control <- createCellChat(object = data.input_control, meta = pbmc_control@meta.data, group.by = "SingleR_fine")
cellchat_control@DB <- CellChatDB
saveRDS(cellchat_control, file.path(com_dir, "对照组", "cellchat_control_step1.rds"))
cellchat_control <- readRDS(file.path(com_dir, "对照组", "cellchat_control_step1.rds"))

cellchat_control <- subsetData(cellchat_control)
cellchat_control <- identifyOverExpressedGenes(cellchat_control)
cellchat_control <- identifyOverExpressedInteractions(cellchat_control)
cellchat_control <- projectData(cellchat_control, ppi)
cellchat_control <- computeCommunProb(cellchat_control, type = "truncatedMean", trim = 0.1, raw.use = FALSE)
cellchat_control <- filterCommunication(cellchat_control, min.cells = 5)
cellchat_control <- computeCommunProbPathway(cellchat_control)
cellchat_control <- netAnalysis_computeCentrality(cellchat_control)
cellchat_control <- aggregateNet(cellchat_control)
saveRDS(cellchat_control, file.path(com_dir, "对照组", "cellchat_control_final.rds"))
cellchat_control <- readRDS(file.path(com_dir, "对照组", "cellchat_control_final.rds"))



# 4. 绘制总体网络圈图 - 疾病组与对照组 ------------------------------------------

## 4.1 疾病组
groupSize_disease <- as.numeric(table(cellchat_disease@idents))
vertex.weight_disease <- pmin(pmax(groupSize_disease / max(groupSize_disease) * 5, 1), 5)

# 绘制疾病组总体网络圈图
pdf(file.path(com_dir, "疾病组", "Disease_Overall_Network_Strength.pdf"), width = 8, height = 8)
netVisual_circle(cellchat_disease@net$weight, vertex.weight = vertex.weight_disease, title.name = "Disease - Strength")
dev.off()

## 4.2 对照组
groupSize_control <- as.numeric(table(cellchat_control@idents))
vertex.weight_control <- pmin(pmax(groupSize_control / max(groupSize_control) * 5, 1), 5)

# 绘制对照组总体网络圈图
pdf(file.path(com_dir, "对照组", "Control_Overall_Network_Strength.pdf"), width = 8, height = 8)
netVisual_circle(cellchat_control@net$weight, vertex.weight = vertex.weight_control, title.name = "Control - Strength")
dev.off()

# 5. 修改提取函数支持默认10个通路 -----------------------------------------------

# 提取前N个活跃通路的函数
extract_top_pathways <- function(cellchat_obj, top_n = 10) {
  pathways <- cellchat_obj@netP$pathways
  pathway_strength <- sapply(pathways, function(path) {
    sum(cellchat_obj@netP$prob[, , path], na.rm = TRUE)
  })
  top_pathways <- names(head(sort(pathway_strength, decreasing = TRUE), top_n))
  return(top_pathways)
}

# 6. 疾病组分析 ---------------------------------------------------------------

# 提取疾病组前10个活跃通路
top10_disease <- extract_top_pathways(cellchat_disease, 10)
cat("疾病组前10个活跃通路:\n")
print(top10_disease)

# 绘制气泡图
pdf(file.path(com_dir, "疾病组", "Disease_Top10_Pathways_Bubble.pdf"), width = 10, height = 6)
netAnalysis_signalingRole_scatter(cellchat_disease, signaling = top10_disease)
dev.off()

# 绘制角色分析网络图
pdf(file.path(com_dir, "疾病组", "Disease_Role_Analysis.pdf"), width = 10, height = 8)
netAnalysis_signalingRole_network(cellchat_disease, signaling = top10_disease)
dev.off()

# 7. 对照组分析 ---------------------------------------------------------------

# 提取对照组前10个活跃通路
top10_control <- extract_top_pathways(cellchat_control, 10)
cat("对照组前10个活跃通路:\n")
print(top10_control)

# 绘制气泡图
pdf(file.path(com_dir, "对照组", "Control_Top10_Pathways_Bubble.pdf"), width = 10, height = 6)
netAnalysis_signalingRole_scatter(cellchat_control, signaling = top10_control)
dev.off()

# 绘制角色分析网络图
pdf(file.path(com_dir, "对照组", "Control_Role_Analysis.pdf"), width = 10, height = 8)
netAnalysis_signalingRole_network(cellchat_control, signaling = top10_control)
dev.off()

# 8. 组间比较分析 ---------------------------------------------------------------

## 8.1 合并CellChat对象
cellchat_merged <- mergeCellChat(
  list(Control = cellchat_control, Disease = cellchat_disease),
  add.names = c("Control", "Disease")
)
saveRDS(cellchat_merged, file.path(com_dir, "比较分析", "cellchat_merged.rds"))

## 8.2 加载合并对象
cellchat_merged <- readRDS(file.path(com_dir, "比较分析", "cellchat_merged.rds"))

## 8.3 绘制比较Strength差异图
pdf(file.path(com_dir, "比较分析", "Differential_Strength_Comparison.pdf"), width = 8, height = 8)
netVisual_diffInteraction(
  cellchat_merged, 
  measure = "weight", 
  comparison = c("Control", "Disease")
)
dev.off()

## 8.4 绘制比较通路活性排行图
pdf(file.path(com_dir, "比较分析", "Pathway_Activity_Ranking.pdf"), width = 10, height = 20)
rankNet(
  cellchat_merged,
  mode = "comparison",
  stacked = TRUE
) + theme(axis.text.x = element_text(size = 6, angle = 45, hjust = 1))
dev.off()





