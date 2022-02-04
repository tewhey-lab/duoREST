head(standard_res$GM12878)

all_data_single_df <- data.frame()
for(celltype in names(standard_res)){
  message(celltype)
  
  id_split <- colsplit(standard_res[[celltype]]$ID, "\\^",c("enhancer","silencer"))
  
  for(enh in unique(id_split$enhancer)){
    message(enh)
    
    temp <- id_split$silencer[which(id_split$enhancer==enh)]
    temp2 <- standard_res[[celltype]][which(id_split$enhancer==enh), c("project","log2FoldChange","lfcSE")]
    
    temp <- cbind(temp, temp2)
    
    colnames(temp) <- c("silencer","project",paste0(celltype,".",enh,".l2FC"), paste0(celltype,".",enh,".l2FC_SE"))
    
    if(celltype == "GM12878" & enh == "En02"){
      all_data_single_df <- temp
    }
    else{
      all_data_single_df <- merge(all_data_single_df, temp, by=c("silencer","project"), all=T)
    }
  }
}

head(all_data_single_df)

write.table(all_data_single_df,"results/silencer_project_l2FC_l2FC.SE_all_cell_enh.txt", quote=F, sep = "\t", row.names = F)


A_above <- cell_enh_sub$log2FoldChange[which(cell_enh_sub$REST=="bound" & cell_enh_sub$silencer %in% temp_sil_list & cell_enh_sub$silencer %in% oligo_split_list[[celltype]][[enh]]$above)]
A_above_med <- median(A_above, na.rm = T)
A_below <- cell_enh_sub$log2FoldChange[which(cell_enh_sub$REST=="bound" & cell_enh_sub$silencer %in% temp_sil_list & cell_enh_sub$silencer %in% oligo_split_list[[celltype]][[enh]]$below)]
A_below_med <- median(A_below, na.rm = T)
C_above <- cell_enh_sub$log2FoldChange[which(cell_enh_sub$REST=="bound" & cell_enh_sub$silencer %notin% temp_sil_list & cell_enh_sub$silencer %in% oligo_split_list[[celltype]][[enh]]$above)]
C_above_med <- median(C_above, na.rm = T)
C_below <- cell_enh_sub$log2FoldChange[which(cell_enh_sub$REST=="bound" & cell_enh_sub$silencer %notin% temp_sil_list & cell_enh_sub$silencer %in% oligo_split_list[[celltype]][[enh]]$below)]
C_below_med <- median(C_below, na.rm = T)


### REDO with color relating to significance and positivity of odds ratio as opposed to the actual odds ratio value
or_tf_df <- data.frame()
for(celltype in unique(med_comp_plot_df$celltype)){
  message(celltype)
  for(enh in enh_all){
    message(enh)
    cell_enh_sub <- med_comp_plot_df[which(med_comp_plot_df$celltype==celltype & med_comp_plot_df$enh==enh),]
    
    cell_enh_sub2 <- bound_factor[which(bound_factor$celltype==celltype & bound_factor$enhancer==enh),1:8]
    
    cell_enh_sub$OR <- NA
    cell_enh_sub$OR_sig <- NA
    for(tf in names(tf_sil_list[[celltype]])){
      temp_sil_list <- tf_sil_list[[celltype]][[tf]]
      A <- cell_enh_sub2$log2FoldChange[which(cell_enh_sub2$REST=="bound" & cell_enh_sub2$silencer %in% temp_sil_list & cell_enh_sub2$silencer %in% oligo_split_list[[celltype]][[enh]]$above)]
      
      B <- cell_enh_sub2$log2FoldChange[which(cell_enh_sub2$REST=="unbound" & cell_enh_sub2$silencer %in% temp_sil_list & cell_enh_sub2$silencer %notin% oligo_split_list[[celltype]][[enh]]$above)]
      
      C <- cell_enh_sub2$log2FoldChange[which(cell_enh_sub2$REST=="bound" & cell_enh_sub2$silencer %notin% temp_sil_list & cell_enh_sub2$silencer %in% oligo_split_list[[celltype]][[enh]]$above)]
      
      D <- cell_enh_sub2$log2FoldChange[which(cell_enh_sub2$REST=="unbound" & cell_enh_sub2$silencer %notin% temp_sil_list & cell_enh_sub2$silencer %notin% oligo_split_list[[celltype]][[enh]]$above)]
      
      OR_above_mat <- matrix(c(length(A),length(B),length(C),length(D)), nrow = 2)
      ft_above <- fisher.test(OR_above_mat)
      OR_above_est <- as.numeric(ft_above$estimate)
      OR_above_pval <- as.numeric(ft_above$p.value)
      cell_enh_sub$OR[which(cell_enh_sub$tf==tf)] <- OR_above_est
      cell_enh_sub$OR_sig[which(cell_enh_sub$tf==tf)] <- OR_above_pval
    }
    cell_enh_sub$log2OR <- log2(cell_enh_sub$OR)
    or_tf_df <- rbind(or_tf_df, cell_enh_sub)
    
    max_ac <- max(cell_enh_sub$delta_AC, na.rm=T)
    min_ac <- min(cell_enh_sub$delta_AC, na.rm=T)
    max_bd <- max(cell_enh_sub$delta_BD, na.rm = T)
    min_bd <- min(cell_enh_sub$delta_BD, na.rm = T)
    x_max <- max(abs(max_ac), abs(min_ac))
    y_max <- max(abs(max_bd), abs(min_bd))
    pdf(paste0("plots/tf_analysis/scatter_plots/odds_ratio/",celltype,"_",enh,"_delta_bound_v_unbound_centered.pdf"), height=55/25.4, width = 55/25.4)
    print(ggplot(cell_enh_sub, aes(x=delta_AC, y=delta_BD, col=log2OR, shape=significance, key=tf)) + geom_point(aes(alpha=0.3), size=1) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
      xlim(-x_max,x_max) + ylim(-y_max,y_max) + xlab("") + ylab("") +
      theme_light(base_size = 6) + scale_shape_manual(values = c("both"=15,"REST+ only"=18,"REST- only"=17, "neither"=16),name="Significance") + scale_color_gradient(low="blue",high="orange") +
      theme(legend.position = "bottom", legend.margin=margin(t = -0.5, unit='cm'),legend.justification = "left", legend.key.width = unit(0.7, "mm"), legend.box = "vertical") +
      guides(shape=guide_legend(override.aes = list(size = 0.5)),color = guide_legend(override.aes = list(size = 0.5)),
             alpha=F, size=F))
    dev.off()
    
    pdf(paste0("plots/tf_analysis/scatter_plots/odds_ratio/",celltype,"_",enh,"_delta_bound_v_unbound.pdf"), height=55/25.4, width = 55/25.4)
    print(ggplot(cell_enh_sub, aes(x=delta_AC, y=delta_BD, col=log2OR, shape=significance, key=tf)) + geom_point(aes(alpha=0.3), size=1) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
            xlim(-2.5,7.25) + ylim(-2.5,7.25) + xlab("") + ylab("") +
            theme_light(base_size = 6) + scale_shape_manual(values = c("both"=15,"REST+ only"=18,"REST- only"=17, "neither"=16),name="Significance") + scale_color_gradient(low="blue",high="orange") +
            theme(legend.position = "bottom", legend.margin=margin(t = -0.5, unit='cm'),legend.justification = "left", legend.key.width = unit(0.7, "mm"), legend.box = "vertical") +
            guides(shape=guide_legend(override.aes = list(size = 0.5)),color = guide_legend(override.aes = list(size = 0.5)),
                   alpha=F, size=F))
    dev.off()
  }
}


or_tf_df$quadrant <- 0
or_tf_df$quadrant[which(or_tf_df$delta_AC > 0 & or_tf_df$delta_BD > 0)] <- 1
or_tf_df$quadrant[which(or_tf_df$delta_AC < 0 & or_tf_df$delta_BD > 0)] <- 2
or_tf_df$quadrant[which(or_tf_df$delta_AC < 0 & or_tf_df$delta_BD < 0)] <- 3
or_tf_df$quadrant[which(or_tf_df$delta_AC > 0 & or_tf_df$delta_BD < 0)] <- 4

or_tf_df$sig_or <- ">=0.05"
or_tf_df$sig_or[which(or_tf_df$OR_sig < 0.05)] <- "<0.05"

for(celltype in unique(or_tf_df$celltype)){
  message(celltype)
  for(enh in enh_all){
    message(enh)
    to_plot <- or_tf_df[which(or_tf_df$celltype==celltype & or_tf_df$enh==enh & or_tf_df$significance=="both"),]
    to_plot$quadrant <- factor(to_plot$quadrant, levels=c("1","2","3","4"))
    # pdf(paste0("plots/tf_analysis/jitter/",celltype,"_",enh,"_OR_sig_quadrant_both_only.pdf"), height=55/25.4, width = 55/25.4)
    
    tmp <- ggplot(to_plot, aes(x=quadrant, y=log2OR, col=sig_or, key=tf)) + geom_jitter(alpha=0.4, size=0.6) + # geom_point(alpha=0.4, size=0.6) +
            scale_color_manual(values=c(">=0.05"="black","<0.05"="red"), name="OR Significance") + theme_light(base_size = 6) + theme(legend.position = "bottom")
    tmply <- ggplotly(tmp)
    
    htmlwidgets::saveWidget(tmply,paste0("plots/tf_analysis/jitter/plotly/",celltype,"_",enh,"_OR_sig_quadrant_both_only.html"))
    # dev.off()
  }
}
