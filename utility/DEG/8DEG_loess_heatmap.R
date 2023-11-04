suppressMessages({
  library("plyr")
  library("dplyr")
  library("tidyverse")
  library("ggrepel")
  library("ggthemes")
  library("ggpubr")
  library("Seurat")
  library("factoextra")
  library("scales")
  library("patchwork")
  library("clusterProfiler")
  library("pheatmap")
  library("colorspace")
  library("GenomicRanges")
  library("edgeR")
})

set.seed(123)

#non_gene_cols=c("PC1","PC2","PC3","age", "race", "gender")
non_gene_cols=c("age")
run_loess <- function(gene, data) {
  # Extract gene specific data
  data <- data[, c(gene, non_gene_cols)]
  
  # Generate LOESS model and save to lo_data
  lo <- loess(get(gene) ~ age, data, span = 0.75)
  
  # compute p-value
 # p=permutate_p( get(gene), data$age)
 
 # colnames(p) = gene

  # Generate predicted values
  lo_predict <- predict(lo, data.frame(age = seq(age_min, age_max, 1))) %>%
    # lo_predict <- predict(lo, data.frame(age = seq(age_min, age_max, 0.0833))) %>%
    as.data.frame()
  colnames(lo_predict) <- gene

 # return(lo_predict,p)
  
 return(lo_predict)
}


heatmap_of_loess <- function(lo_predict, output_dir, cell_type_label,cluster=NULL,modal=NULL) {
  
  
  # Convert wide to long
  if(modal == "rna"){
  lo_predict_long <- gather(lo_predict, gene, expression, 2:ncol(lo_predict))
  }else if(modal == "atac"){
  lo_predict_long <- gather(lo_predict, gene, chromatin, 2:ncol(lo_predict))

  }
  if(is.null(cluster)){ 
  # Order genes by hclustering
  clust_order <- hclust(dist(t(lo_predict[,2:ncol(lo_predict)])), method = "ward.D2" )$order
 }else{
  clust_order=cluster
  }
  lo_predict_long$gene <- factor(lo_predict_long$gene, levels = unique(lo_predict_long$gene)[clust_order])
  gene_order=clust_order #data.frame(gene=unique(lo_predict_long$gene), order=clust_order)
 
  # Add limits for heatmap plotting
  zscore_limit <- 0.5
#  zscore_limit <- 1.5

  lo_predict_long_plotting <- lo_predict_long
  if(modal == "rna"){
  lo_predict_long_plotting[ which(lo_predict_long_plotting$expression < -zscore_limit), "expression"] <- -zscore_limit
  lo_predict_long_plotting[ which(lo_predict_long_plotting$expression > zscore_limit), "expression"] <- zscore_limit
  }else if(modal == "atac"){
 lo_predict_long_plotting[ which(lo_predict_long_plotting$expression < -zscore_limit), "chromatin"] <- -zscore_limit
  lo_predict_long_plotting[ which(lo_predict_long_plotting$expression > zscore_limit), "chromatin"] <- zscore_limit


  }
  genenum=length(unique(lo_predict_long_plotting$gene))
  # Generate heatmap
  if(modal == "rna"){
  p <- ggplot(lo_predict_long_plotting , aes(x = age, y = gene, fill = expression)) +
    theme_Publication_blank() +
    geom_raster(interpolate=TRUE) +
    theme(axis.text.y=element_blank()) +
    scale_fill_gradient2(low = "cyan", mid = "black", high = "yellow",
                         breaks = c(-zscore_limit, zscore_limit),
                         limits = c(-zscore_limit, zscore_limit))+ggtitle(paste0(genenum," genes changed over aging"))


#                         limits = c(-zscore_limit, zscore_limit))+ggtitle(paste0(genenum," genes correlated with aging"))
  }else if(modal == "atac"){
  p <- ggplot(lo_predict_long_plotting , aes(x = age, y = gene, fill = chromatin)) +
    theme_Publication_blank() +
    geom_raster(interpolate=TRUE) +
    theme(axis.text.y=element_blank()) +
    scale_fill_gradient2(low = "cyan", mid = "black", high = "yellow",
                         breaks = c(-zscore_limit, zscore_limit),
                         limits = c(-zscore_limit, zscore_limit))+ggtitle(paste0(genenum," genes changed over aging"))

#                         limits = c(-zscore_limit, zscore_limit))+ggtitle(paste0(genenum," genes correlated with aging"))


}
  # Export heatmap
  ggplot2::ggsave(paste0(output_dir, "/", cell_type_label, "_",modal,"_loess_predicted_heatmap.pdf"), p, width = 5, height = 4)
  return(gene_order)
}


line_plot_of_loess <- function(data, color, ylim = NULL, data_compare = NULL,
                               alpha = 0.05, size = 0.3, plot_ind_genes = TRUE,cluster=NULL) {
  # Find average expression over age of all genes in cluster
  avg <- data %>%
    group_by(age) %>%
   summarize_each(funs(mean, sd, se=sd(.)/sqrt(n())), expression)

#    summarize_each(funs(mean, sd, se=sd(.)/sqrt(n())), expression)
  
  # Generate LOESS model and save to lo_data
  lo <- loess(mean ~ age, avg, span = 0.75)
  lo_predict <- predict(lo, data.frame(age = avg$age))
  avg <- add_column(avg, lo = lo_predict)
  
  # Define gene number annotation df
  annotations <- data.frame(
    xpos = c(-Inf),
    ypos =  c(Inf),
    annotateText = c(paste0("c", cluster ," n=", length(unique(data$gene)))),
    hjustvar = c(-0.2) ,
    vjustvar = c(1))
  
  if (is.null(data_compare)) {
    # Generate plot
    plot <- ggplot(data, aes(x = age, y = expression)) +
      theme_Publication_blank() + 
      geom_line(stat="smooth", method = "loess", span = 0.75, se = FALSE,
                aes(group = gene),
                alpha = .05, color = color, size = size) +
      scale_y_continuous(expand = c(0,0)) +
      scale_x_continuous(expand = c(0,0)) +
      geom_line(data = avg,
                aes(x = age, y = mean), #avg_exp),
                stat="smooth", method = "loess", span = 0.75, se = FALSE,
                color = darken(color, amount = 0.3), size = size * 3) +
      theme(legend.position = "none",
            aspect.ratio = 1,
            text = element_text(size = 15),
            panel.margin = unit(c(0, 0, 0, 0), "null")) +
      geom_text(data=annotations,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText)) +
      {if(!is.null(ylim))ylim(ylim)}
  } else {
    # Find average expression over age of all genes in cluster
    avg_compare <- data_compare %>%
      group_by(age) %>%
      summarize_each(funs(mean, sd, se=sd(.)/sqrt(n())), expression)
    
    # Generate LOESS model and save to lo_data
    lo <- loess(mean ~ age, avg_compare, span = 0.75)
    lo_predict <- predict(lo, data.frame(age = avg_compare$age))
    avg_compare <- add_column(avg_compare, lo = lo_predict)
    
    # Find shared age range
    min_age_compare <- max(min(data$age), min(data_compare$age))
    max_age_compare <- min(max(data$age), max(data_compare$age))
    
    # Extract regression values
    p_tmp <- ggplot(data, aes(x = age, y = expression)) +
      geom_line(data = avg, aes(x = age, y = mean),
                stat="smooth", method = "loess", span = 0.75)
    stats <- ggplot_build(p_tmp)$data[[1]][seq(1,80,by=8),]
    
    p_tmp <- ggplot(data_compare, aes(x = age, y = expression)) +
      geom_line(data = avg_compare, aes(x = age, y = mean),
                stat="smooth", method = "loess", span = 0.75)
    compare_stats <- ggplot_build(p_tmp)$data[[1]][seq(1,80,by=8),]
    
    # Generate plot
    plot <- ggplot(data, aes(x = age, y = expression)) +
      theme_Publication_blank() +
      {if(plot_ind_genes)
        geom_line(stat="smooth", method = "loess", span = 0.75, se = FALSE,
                  aes(group = gene),
                  alpha = alpha, color = "gray", size = size)} +
      {if(plot_ind_genes)
        geom_line(data = data_compare,
                  aes(x = age, y = expression, group = gene),
                  stat="smooth", method = "loess", span = 0.75, se = FALSE,
                  alpha = alpha, color = color, size = size)} +
      geom_line(data = avg,
                aes(x = age, y = mean),
                stat="smooth", method = "loess", span = 0.75, se = FALSE,
                color = "black", size = size * 3) + #gray40
      geom_line(data = avg_compare,
                aes(x = age, y = mean),
                stat="smooth", method = "loess", span = 0.75, se = FALSE,
                color = color, size = size * 3) +
      geom_errorbar(inherit.aes = FALSE, data = stats, color = "black",
                    mapping = aes(x = x, ymin = y-se, ymax=y+se), width = 0.5, size = 0.3) +
      geom_errorbar(inherit.aes = FALSE, data = compare_stats, color = color,
                    mapping = aes(x = x, ymin = y-se, ymax=y+se), width = 0.5, size = 0.3) +
      theme(legend.position = "none",
            aspect.ratio = 1,
            text = element_text(size = 8),
            panel.margin = unit(c(0, 0, 0, 0), "null")) +
      coord_cartesian(xlim=c(min_age_compare, max_age_compare)) +
      geom_text(data=annotations,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText)) +
      {if(!is.null(ylim))ylim(ylim)}
  }
  return(plot)
}


cluster_loess <- function( data_scaled, output_dir, cell_type_label,
                          cluster_nums = 1:5, lo_predict = NULL, hclust_cut_merged = NULL) {
  # If not performing comparison
  if (!is.null(lo_predict)) {
    # Perform distance calculation on all trajectories (removing age column)
    dist_mat <- dist(t(lo_predict[,2:ncol(lo_predict)]))
    
    # Plot hierarchal clusters
    pdf(paste0(output_dir, "/", cell_type_label, "_dendogram_loess_predicted.pdf"))
    clust <- hclust(dist_mat, method = "complete")
    #clust <- hclust(dist_mat, method = "average")

    plot(clust, labels = FALSE)
    dev.off()
    
    # Initialize hclust_cut_list
    hclust_cut_list <- list()
  }
  
  # Generate clustering with k = cluster_nums
  for (k in cluster_nums) {
    # Check if we have enough genes
    if (!is.null(lo_predict)) {
      if (k > (ncol(lo_predict) - 1)) {
        print(paste0("Insufficient gene number for ", k, " clusters"))
        next
      }
      
      # Generate hierarchical clusters
      hclust_cut <- data.frame(cutree(clust, k = k))
    } else {
      # Isolate hierarchical clusters
      hclust_cut <- hclust_cut_merged[,match(k, cluster_nums),drop = FALSE]
    }
    
    # Generate color vector
    color_vector <- hue_pal()(k)
    
    # Initialize plot list
    hclust_plot_list <- list()
    
    for (i in seq(k)) {
      # Subset genes of interest
      cols <- c(rownames(hclust_cut[ which(hclust_cut[,1] == i), ,drop = FALSE]), "age")
      data_scaled_subset <- data_scaled[ ,cols ]
      
      if (!is.null(lo_predict)) {
        data_long <- gather(data_scaled_subset, gene, expression, -age)
        
        # Run plotting
        plot <- line_plot_of_loess(data_long, color = color_vector[i], cluster=i)
      } else {
        # Extract comparison parameters of interest
        group <- c("clonal", "Diagnosis", "sex")[grep("compare", c(clonal, Diagnosis, sex))]
        if (group == "clonal") {
          group1 <- "NC"
          group2 <- "C"
        } else if (group == "Diagnosis") {
          group1 <- "HC"
          group2 <- "MCI/AD"
        } else if (group == "sex") {
          group1 <- "m"
          group2 <- "f"
        }
        
        # Subset two groups
        data_tmp <- data_scaled_subset[which(data_scaled[,group] == group1),]
        data_tmp_compare <- data_scaled_subset[which(data_scaled[,group] == group2),]
        data_long <- gather(data_tmp, gene, expression, -age)
        data_long_compare <- gather(data_tmp_compare, gene, expression, -age)
        
        # Run plotting
        plot <- line_plot_of_loess(data_long, color = color_vector[i],
                                   data_compare = data_long_compare, cluster=i)
      }
      
      # Append plot to list
      hclust_plot_list[[i]] <- plot
    }
    
    # Generate figure of clustered trajectories
    if (k == 1) {
      pdf(paste0(output_dir, "/", cell_type_label, "_hclust_k", k, ".pdf"))
      print(hclust_plot_list[[1]])
      dev.off()
    } else {
      patchwork <- wrap_plots(hclust_plot_list, guides = "collect")
      if (is.null(lo_predict)) {
        patchwork <- patchwork + plot_annotation(subtitle = paste0(cell_type_label, ": gray = ", group1, ", color = ", group2))
      } else {
        patchwork <- patchwork + plot_annotation(subtitle = cell_type_label)
      }
      pdf(paste0(output_dir, "/", cell_type_label, "_hclust_k", k, ".pdf"))
      print(patchwork)
      dev.off()
    }
    
    # Append cut to list
    if (!is.null(lo_predict)) {hclust_cut_list[[k]] <- hclust_cut}
  }
  # Write out cluster cut
  if (!is.null(lo_predict)) {
    hclust_cut_merged <- as.data.frame(do.call(cbind, hclust_cut_list[cluster_nums]))
    write.csv(hclust_cut_merged, file = paste0(output_dir, "/", cell_type_label, "_hclust_cut.csv"))
  }
}

theme_Publication_blank <- function(base_size=12, base_family="") {
  (theme_foundation(base_size=base_size, base_family=base_family)
   + theme(plot.title = element_text(size = rel(1.2), hjust = 0.5),
           text = element_text(),
           panel.background = element_rect(fill = "transparent",colour = NA),
           plot.background = element_rect(fill = "transparent",colour = NA),
           panel.border = element_rect(colour = NA, fill = "transparent"),
           axis.title = element_text(size = rel(1)),
           axis.title.y = element_text(angle=90,margin=margin(0,10,0,0)),
           axis.title.x = element_text(margin=margin(10,0,0,0)),
           axis.text = element_text(), 
           axis.line = element_line(colour="black"),
           axis.ticks = element_line(size = 0.3),
           axis.line.x = element_line(size = 0.3, linetype = "solid", colour = "black"),
           axis.line.y = element_line(size = 0.3, linetype = "solid", colour = "black"),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           legend.key = element_rect(colour = NA, fill="transparent"),
           legend.position = "bottom",
           legend.margin = margin(t = 10, unit='pt'),
           plot.margin=unit(c(10,5,5,5),"mm"),
           strip.background=element_rect(colour="#d8d8d8",fill="#d8d8d8")
   ))
} 

args <- commandArgs(trailingOnly = TRUE)

cell_type_label=args[1]
out_dir=args[2]
in_dir=args[3]
exp_dir=args[4]
atac_rds=args[5]
atac_gtf=args[6]
atac_meta=args[7]
#clu_num=args[8]

atac_cn=args[8]
atac_count_ct=args[9]
qval=args[10]
seq=args[11]
seq1=args[12]
#out_dir="DEG_loess_ge_ATAC_clu_genebodyPro_donor"
#in_dir="DEG_loess_ge_ATAC_genebodyPro_donor"
#exp_dir="genexp_donor_cell_raw_new"
#atac_rds="geneBodyProFragCount_Group_celltype_donor"
#atac_gtf="cellranger_gex-GRCh38_tss_str_genebody_pro.gtf"
#atac_meta="human_meta/data/scATAC_sample_list_all_new2_demograph_donor"
#clu_num="6:12"


output_dir=paste0("/storage/chenlab/Users/junwang/human_meta/data/",out_dir,"/",cell_type_label)
dir.create(output_dir,recursive=T)
setwd(output_dir)

#####gene exp##########
data=read.table(paste0("/storage/chenlab/Users/junwang/human_meta/data/region_DESeq2_batch_ageNum_dream/",cell_type_label,"_interval_cpm01_snRNA_clean_young_cpm07_",seq1,"/",cell_type_label,"_DEG_res_cpm_age"),header=T)

gene_list=rownames(data[data$qval< qval ,])

data_scaled=read.table(paste0("/storage/chenlab/Users/junwang/human_meta/data/",exp_dir,"/exp_",cell_type_label,"_logNorm_",seq,"_DEswan"),header=T)

age_min=min(data_scaled$age)
age_min1=53
age_max=max(data_scaled$age)



data_scaled1=scale(data_scaled[, -c(1:6) ],center=TRUE,scale=TRUE)
data_scaled2=cbind(data_scaled[,c(1:6)],data_scaled1)
genes <- colnames(data_scaled2[, -c(1:6) ])

lo_predict_list <- lapply(genes, run_loess, data = data_scaled2)

lo_predict <- as.data.frame(lo_predict_list)


saveRDS(lo_predict,paste0(output_dir,"/exp_",cell_type_label,"_logNorm_",seq,"_lo_DEswan_predict.rds"))

lo_predict <- as.data.frame(lo_predict_list)


colnames(lo_predict)=make.unique(colnames(lo_predict))

m=match(gene_list,colnames(lo_predict),nomatch=0)

lo_predict=lo_predict[,m]



lo_predict$age <- seq(age_min, age_max, 1)

lo_predict <- dplyr::relocate(lo_predict, age)


cluster=NULL
cluster=heatmap_of_loess(lo_predict, output_dir, cell_type_label,cluster,"rna")



write.csv(lo_predict, paste0(output_dir, "/", cell_type_label, "_lo_predict.csv"))





