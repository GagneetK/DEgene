# load required packages

req_packages = c("Biobase", "cluster", "clusterProfiler", "cowplot",
                 "data.table", "edgeR", "ggpubr", 
                 "ggrepel", 
                 "GO.db", "goseq", "grid", "gridExtra", "lattice", 
                 "pheatmap", "qvalue", "RColorBrewer", 
                 "Rmisc", "RUVSeq","splitstackshape", "statmod", "stringr", "tidyverse",
                 "VennDiagram", "viridis")
## load them, quietly
invisible(suppressWarnings(suppressMessages(
    lapply(req_packages, require, character.only = TRUE)
)))

## The Cowplot package changes the default themes of ggplot2. Set the bw theme with larger font sizes like so:
theme_set(theme_bw(base_size = 16))
## ... or set the default theme
# theme_set(theme_gray())

## suppress excessive VennDiagram log files
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
     
