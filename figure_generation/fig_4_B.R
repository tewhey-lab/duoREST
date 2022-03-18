tf_plots <- data.frame()
for(celltype in names(standard_res)){
  temp <- standard_res[[celltype]][standard_res[[celltype]]$project=="motif",c("ID","log2FoldChange")]
  temp$celltype <- celltype
  tf_plots <- rbind(tf_plots,temp)
}

tf_id_split <- colsplit(tf_plots$ID,"\\^",c("enhancer","silencer"))
tf_plots$enhancer <- tf_id_split$enhancer
tf_plots$silencer <- tf_id_split$silencer

tf_plots$tf_stat <- 0
for(celltype in unique(tf_plots$celltype)){
  message(celltype)
  for(sil in unique(tf_plots$silencer)){
    if(any(chip_tf_celltypes[[celltype]][which(chip_tf_celltypes[[celltype]]$name==sil),5:ncol(chip_tf_celltypes[[celltype]])]==1)){
      tf_plots$tf_stat[which(tf_plots$silencer==sil & tf_plots$celltype==celltype)] <- 1
    }
  }
}

tf_plots$tf_stat <- as.factor(tf_plots$tf_stat)
tf_plots$celltype <- as.factor(tf_plots$celltype)

accession_tf_name_map <- read.delim("ENCODE_ChIP/complete_histone_tf_list.txt", stringsAsFactors = F)
colnames(accession_tf_name_map) <- c("Accession","Target","Celltype")

accession_celltype_list <- list()
accession_tf_name_map$Celltype[which(accession_tf_name_map$Celltype=="SK-N-SH")] <- "SKNSH"
accession_tf_name_map$Celltype[which(accession_tf_name_map$Celltype=="HEPG2")] <- "HepG2"
for(celltype in unique(accession_tf_name_map$Celltype)){
  message(celltype)
  accession_celltype_list[[celltype]] <- list()
  for(tf in unique(accession_tf_name_map$Target[which(accession_tf_name_map$Celltype==celltype)])){
    # if(tf=="REST") next
    accession_celltype_list[[celltype]][[tf]] <- accession_tf_name_map$Accession[which(accession_tf_name_map$Target==tf & accession_tf_name_map$Celltype==celltype)]
  }
}

bound_factor <- tf_plots[,1:5]
bound_factor$REST <- "unbound"
for(celltype in unique(bound_factor$celltype)){
  message(celltype)
  for(sil in unique(bound_factor$silencer[which(bound_factor$celltype==celltype)])){
    if(isTRUE(chip_tf_celltypes[[celltype]][which(chip_tf_celltypes[[celltype]]$name==sil),accession_celltype_list[[celltype]]$REST[1]]==1) | isTRUE(chip_tf_celltypes[[celltype]][which(chip_tf_celltypes[[celltype]]$name==sil),accession_celltype_list[[celltype]]$REST[2]]==1)){
      bound_factor$REST[which(bound_factor$silencer==sil & bound_factor$celltype==celltype)] <- "bound"
    }
  }
}

scramble_skew <- list()
for(celltype in names(emVAR_all)){
  message(celltype)
  intersect_list <- bound_factor$ID[bound_factor$ID[bound_factor$celltype==celltype] %in% emVAR_all[[celltype]]$ID]
  scramble_skew[[celltype]] <- emVAR_all[[celltype]][emVAR_all[[celltype]]$ID %in% intersect_list,c("ID","B_log2FC","LogSkew")]
}

bound_factor$scr_l2fc <- NA
bound_factor$LogSkew <- NA
for(celltype in names(scramble_skew)){
  message(celltype)
  # bound_factor[which(bound_factor$celltype==celltype),] <- merge(bound_factor[which(bound_factor$celltype==celltype),], scramble_skew[[celltype]], by="ID", all=T)
  for(id in scramble_skew[[celltype]]$ID){
    bound_factor$scr_l2fc[which(bound_factor$celltype==celltype & bound_factor$ID==id)] <- scramble_skew[[celltype]]$B_log2FC[which(scramble_skew[[celltype]]$ID==id)]
    bound_factor$LogSkew[which(bound_factor$celltype==celltype & bound_factor$ID==id)] <- scramble_skew[[celltype]]$LogSkew[which(scramble_skew[[celltype]]$ID==id)]
  }
}

bound_factor$skew_na <- is.na(bound_factor$LogSkew)

bound_factor$num_bound <- NA
bound_tf_names <- list()
for(celltype in names(chip_tf_celltypes)){
  message(celltype)
  bound_tf_names[[celltype]] <- list()
  for(sil in unique(bound_factor$silencer[which(bound_factor$celltype==celltype)])){
    acc_bound <- colnames(chip_tf_celltypes[[celltype]])[chip_tf_celltypes[[celltype]][which(chip_tf_celltypes[[celltype]]$name==sil),]==1]
    tf_bound <- unique(accession_tf_name_map$Target[which(accession_tf_name_map$Accession %in% acc_bound)])
    bound_tf_names[[celltype]][[sil]] <- tf_bound
    bound_factor$num_bound[which(bound_factor$silencer==sil & bound_factor$celltype==celltype)] <- length(tf_bound)
  }
}

# accession_tf_name_map <- read.delim("ENCODE_ChIP/complete_histone_tf_list.txt", stringsAsFactors = F)
# colnames(accession_tf_name_map) <- c("Accession","Target","Celltype")
# 
# accession_celltype_list <- list()
# accession_tf_name_map$Celltype[which(accession_tf_name_map$Celltype=="SK-N-SH")] <- "SK.N.SH"
# accession_tf_name_map$Celltype[which(accession_tf_name_map$Celltype=="HEPG2")] <- "HepG2"
# for(celltype in unique(accession_tf_name_map$Celltype)){
#   message(celltype)
#   accession_celltype_list[[celltype]] <- list()
#   for(tf in unique(accession_tf_name_map$Target[which(accession_tf_name_map$Celltype==celltype)])){
#     # if(tf=="REST") next
#     accession_celltype_list[[celltype]][[tf]] <- accession_tf_name_map$Accession[which(accession_tf_name_map$Target==tf & accession_tf_name_map$Celltype==celltype)]
#   }
# }

tf_sil_list <- list()
for(celltype in names(accession_celltype_list)){
  message(celltype)
  tf_sil_list[[celltype]] <- list()
  temp <- chip_tf_celltypes[[celltype]]
  for(tf in names(accession_celltype_list[[celltype]])){
    temp_sil_list <- c()
    for(acc in accession_celltype_list[[celltype]][[tf]]){
      temp_sil_list <- c(temp_sil_list, chip_tf_celltypes[[celltype]]$name[which(chip_tf_celltypes[[celltype]][,acc]==1)])
    }
    tf_sil_list[[celltype]][[tf]] <- temp_sil_list
  }
}


med_comp_list <- list()
med_comp_list$A_B <- list()
med_comp_list$A_D <- list()
med_comp_list$C_B <- list()
med_comp_list$C_D <- list()
med_comp_list$other <- list()
for(celltype in names(tf_sil_list)){
  message(celltype)
  med_comp_list$A_B[[celltype]] <- list()
  med_comp_list$A_D[[celltype]] <- list()
  med_comp_list$C_B[[celltype]] <- list()
  med_comp_list$C_D[[celltype]] <- list()
  med_comp_list$other[[celltype]] <- list()
  cell_sub <- bound_factor[which(bound_factor$celltype==celltype),1:8]
  for(enh in unique(cell_sub$enhancer)){
    message(enh)
    med_comp_list$A_B[[celltype]][[enh]] <- list()
    med_comp_list$A_D[[celltype]][[enh]] <- list()
    med_comp_list$C_B[[celltype]][[enh]] <- list()
    med_comp_list$C_D[[celltype]][[enh]] <- list()
    med_comp_list$other[[celltype]][[enh]] <- list()
    med_comp_list$A_B[[celltype]][[enh]]$TFs <- c()
    med_comp_list$A_D[[celltype]][[enh]]$TFs <- c()
    med_comp_list$C_B[[celltype]][[enh]]$TFs <- c()
    med_comp_list$C_D[[celltype]][[enh]]$TFs <- c()
    med_comp_list$other[[celltype]][[enh]]$TFs <- c()
    cell_enh_sub <- cell_sub[which(cell_sub$enhancer==enh),]
    for(tf in names(tf_sil_list[[celltype]])){
      temp_sil_list <- tf_sil_list[[celltype]][[tf]]
      RB_TFB <- median(cell_enh_sub$log2FoldChange[which(cell_enh_sub$REST=="bound" & cell_enh_sub$silencer %in% temp_sil_list)], na.rm = T)
      RU_TFB <- median(cell_enh_sub$log2FoldChange[which(cell_enh_sub$REST=="unbound" & cell_enh_sub$silencer %in% temp_sil_list)], na.rm = T)
      RB_TFU <- median(cell_enh_sub$log2FoldChange[which(cell_enh_sub$REST=="bound" & cell_enh_sub$silencer %notin% temp_sil_list)], na.rm = T)
      RU_TFU <- median(cell_enh_sub$log2FoldChange[which(cell_enh_sub$REST=="unbound" & cell_enh_sub$silencer %notin% temp_sil_list)], na.rm = T)
      if(any(is.na(c(RB_TFB, RU_TFB)))){
        med_comp_list$other[[celltype]][[enh]]$TFs <- c(med_comp_list$other[[celltype]][[enh]]$TFs,tf)
        next
      }
      if(RB_TFB > RB_TFU & RU_TFB > RU_TFU){
        med_comp_list$A_B[[celltype]][[enh]]$TFs <- c(med_comp_list$A_B[[celltype]][[enh]]$TFs,tf)
        med_comp_list$A_B[[celltype]][[enh]]$med_skew <- median(cell_enh_sub$LogSkew[which(cell_enh_sub$silencer %in% temp_sil_list)], na.rm = T)
      }
      else if(RB_TFB > RB_TFU & RU_TFB < RU_TFU){
        med_comp_list$A_D[[celltype]][[enh]]$TFs <- c(med_comp_list$A_D[[celltype]][[enh]]$TFs,tf)
        med_comp_list$A_D[[celltype]][[enh]]$med_skew <- median(cell_enh_sub$LogSkew[which((cell_enh_sub$REST=="bound" & cell_enh_sub$silencer %in% temp_sil_list) | 
                                                                                             (cell_enh_sub$REST=="unbound" & cell_enh_sub$silencer %notin% temp_sil_list))], na.rm = T)
      }
      else if(RB_TFB < RB_TFU & RU_TFB > RU_TFU){
        med_comp_list$C_B[[celltype]][[enh]]$TFs <- c(med_comp_list$C_B[[celltype]][[enh]]$TFs,tf)
        med_comp_list$C_B[[celltype]][[enh]]$med_skew <- median(cell_enh_sub$LogSkew[which((cell_enh_sub$REST=="bound" & cell_enh_sub$silencer %notin% temp_sil_list) | 
                                                                                             (cell_enh_sub$REST=="unbound" & cell_enh_sub$silencer %in% temp_sil_list))], na.rm = T)
      }
      else if(RB_TFB < RB_TFU & RU_TFB < RU_TFU){
        med_comp_list$C_D[[celltype]][[enh]]$TFs <- c(med_comp_list$C_D[[celltype]][[enh]]$TFs,tf)
        med_comp_list$C_D[[celltype]][[enh]]$med_skew <- median(cell_enh_sub$LogSkew[which(cell_enh_sub$silencer %notin% temp_sil_list)], na.rm = T)
      }
    }
  }
}

med_comp_plot_df <- data.frame()
for(celltype in names(tf_sil_list)){
  message(celltype)
  cell_sub <- bound_factor[which(bound_factor$celltype==celltype),1:8]
  for(enh in unique(cell_sub$enhancer)){
    message(enh)
    cell_enh_sub <- cell_sub[which(cell_sub$enhancer==enh),]
    for(tf in names(tf_sil_list[[celltype]])){
      temp_sil_list <- tf_sil_list[[celltype]][[tf]]
      A <- cell_enh_sub$log2FoldChange[which(cell_enh_sub$REST=="bound" & cell_enh_sub$silencer %in% temp_sil_list)]
      A_med <- median(A, na.rm = T)
      B <- cell_enh_sub$log2FoldChange[which(cell_enh_sub$REST=="unbound" & cell_enh_sub$silencer %in% temp_sil_list)]
      B_med <- median(B, na.rm = T)
      C <- cell_enh_sub$log2FoldChange[which(cell_enh_sub$REST=="bound" & cell_enh_sub$silencer %notin% temp_sil_list)]
      C_med <- median(C, na.rm = T)
      D <- cell_enh_sub$log2FoldChange[which(cell_enh_sub$REST=="unbound" & cell_enh_sub$silencer %notin% temp_sil_list)]
      D_med <- median(D, na.rm = T)
      if(length(A[is.na(A)])==length(A)){
        A <- c()
      }
      if(length(B[is.na(B)])==length(B)){
        B <- c()
      }
      if(length(C[is.na(C)])==length(C)){
        C <- c()
      }
      if(length(D[is.na(D)])==length(D)){
        D <- c()
      }
      if(length(A) > 0 & length(C) > 0){
        A_C_sig <- wilcox.test(A,C)$p.value
      }
      if(length(A) == 0 | length(C) == 0){
        A_C_sig <- NA
      }
      if(length(B) > 0 & length(D) > 0){
        B_D_sig <- wilcox.test(B,D)$p.value
      }
      if(length(B) == 0 | length(D) == 0){
        B_D_sig <- NA
      }
      temp <- data.frame(celltype, enh, tf, A_med, B_med, C_med, D_med, A_C_sig, B_D_sig)
      med_comp_plot_df <- rbind(med_comp_plot_df, temp)
    }
  }
}

med_comp_plot_df$delta_AC <- med_comp_plot_df$A_med - med_comp_plot_df$C_med
med_comp_plot_df$delta_BD <- med_comp_plot_df$B_med - med_comp_plot_df$D_med

med_comp_plot_df$A_C_adj <- p.adjust(med_comp_plot_df$A_C_sig, method = "hochberg")
med_comp_plot_df$B_D_adj <- p.adjust(med_comp_plot_df$B_D_sig, method = "hochberg")

med_comp_plot_df$significance <- NA
med_comp_plot_df$significance[which(med_comp_plot_df$A_C_adj < 0.05 & med_comp_plot_df$B_D_adj < 0.05)] <- "both"
med_comp_plot_df$significance[which(med_comp_plot_df$A_C_adj < 0.05 & med_comp_plot_df$B_D_adj > 0.05)] <- "REST+ only"
med_comp_plot_df$significance[which(med_comp_plot_df$A_C_adj > 0.05 & med_comp_plot_df$B_D_adj < 0.05)] <- "REST- only"
med_comp_plot_df$significance[which(med_comp_plot_df$A_C_adj > 0.05 & med_comp_plot_df$B_D_adj > 0.05)] <- "neither"
med_comp_plot_df$significance <- factor(med_comp_plot_df$significance, levels = c("both","REST+ only","REST- only","neither"))

for(celltype in unique(med_comp_plot_df$celltype)){
  message(celltype)
  for(enh in enh_all){
    message(enh)
    cell_enh_sub <- med_comp_plot_df[which(med_comp_plot_df$celltype==celltype & med_comp_plot_df$enh==enh),]
    max_ac <- max(cell_enh_sub$delta_AC, na.rm=T)
    min_ac <- min(cell_enh_sub$delta_AC, na.rm=T)
    max_bd <- max(cell_enh_sub$delta_BD, na.rm = T)
    min_bd <- min(cell_enh_sub$delta_BD, na.rm = T)
    x_max <- max(abs(max_ac), abs(min_ac))
    y_max <- max(abs(max_bd), abs(min_bd))
    # pdf(paste0("plots/tf_analysis/scatter_plots/",celltype,"_",enh,"_delta_bound_v_unbound_centered_bhadj.pdf"), height=55/25.4, width = 55/25.4)
    a <- ggplot(cell_enh_sub, aes(x=delta_AC, y=delta_BD, col=significance, key=tf)) + geom_point(aes(alpha=0.3), size=2) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
      xlim(-x_max,x_max) + ylim(-y_max,y_max) + xlab("") + ylab(NULL) +
      theme_light(base_size = 6) + scale_color_manual(values = c("red","dodgerblue","green3","black"),name="Significance") +
      theme(legend.position = "bottom", legend.margin=margin(t = -0.6, unit='cm'),legend.justification = "left", legend.key.width = unit(0.6, "mm")) + 
      guides(fill=guide_legend(nrow = 1),shape=guide_legend(override.aes = list(size = 0.5)),color = guide_legend(override.aes = list(size = 0.5)),
             alpha=F, size=F)
    # print(a)
    # dev.off()
    
    a_ly <- ggplotly(a)
    htmlwidgets::saveWidget(as_widget(a_ly), paste0("plots/tf_analysis/scatter_plots/",celltype,"_",enh,"_delta_bound_v_unbound_centered_bhadj.html"))
    
    # pdf(paste0("plots/tf_analysis/scatter_plots/",celltype,"_",enh,"_delta_bound_v_unbound.pdf"), height=55/25.4, width = 55/25.4)
    # print(ggplot(cell_enh_sub, aes(x=delta_AC, y=delta_BD, col=significance)) + geom_point(aes(alpha=0.3), size=2) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
    #         xlim(-2.5,7.25) + ylim(-2.5,7.25) + xlab("") + ylab("") +
    #         theme_light(base_size = 6) + scale_color_manual(values = c("red","dodgerblue","green3","black"),name="Significance") +
    #         theme(legend.position = "bottom", legend.margin=margin(t = -0.5, unit='cm'),legend.justification = "left", legend.key.width = unit(0.7, "mm")) + 
    #         guides(fill=guide_legend(nrow = 1),shape=guide_legend(override.aes = list(size = 0.5)),color = guide_legend(override.aes = list(size = 0.5)),
    #                alpha=F, size=F)) 
    # dev.off()
    
  }
}
