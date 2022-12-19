############################################
### Hierarchical clustering of apoptosis expression levels
### If you use this code or input files, please cite:
### Fitzgerald et al., 2023 XXXXXX
### Contact niamhmconnolly@rcsi.com for more information
############################################

library(ComplexHeatmap)

# Set working directory (UPDATE!!!)
setwd("XXX/YYY")

# Read in protein levels (normalised to HeLa)
proteins <- read.csv("Protein expression_individual.csv") 
rownames(proteins) <- proteins[,1]
proteins <- proteins[,-1]

log2proteins <- log2(proteins)

#zscore across proteins
z_log2proteins <- as.data.frame(scale(log2proteins, center = TRUE, scale = TRUE))

# Add subgroups
subgroup <- c("SHH","SHH", "SHH","Group3","Group3","Group4","Group4","Group4")

# Set cell line annotation
cell_line_annotation <- 
  HeatmapAnnotation(subgroup=subgroup, 
                    which="col", 
                    col=list(subgroup=c("Group3"="navy", "Group4"="blue", "SHH"="cyan")))

# Cluster heatmap
Heatmap(as.matrix(t(z_log2proteins)), 
        top_annotation = cell_line_annotation,
        heatmap_legend_param=list(title="zscore\nexpr."))
