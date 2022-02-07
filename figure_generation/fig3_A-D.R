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
for(celltype in names(chip_tf_celltypes)){
  message(celltype)
  for(sil in unique(bound_factor$silencer[which(bound_factor$celltype==celltype)])){
    acc_bound <- colnames(chip_tf_celltypes[[celltype]])[chip_tf_celltypes[[celltype]][which(chip_tf_celltypes[[celltype]]$name==sil),]==1]
    tf_bound <- unique(accession_tf_name_map$Target[which(accession_tf_name_map$Accession %in% acc_bound)])
    bound_tf_names[[celltype]][[sil]] <- tf_bound
    bound_factor$num_bound[which(bound_factor$silencer==sil & bound_factor$celltype==celltype)] <- length(tf_bound)
  }
}

neg_ctrl_all <- data.frame()
for(celltype in names(standard_res)){
  message(celltype)
  neg_ctrl_temp <- standard_res[[celltype]][which(standard_res[[celltype]]$project=="NegCtrl"),c("ID","log2FoldChange","project")]
  neg_ctrl_temp$celltype <- celltype
  neg_ctrl_temp$scr_l2fc <- NA
  neg_ctrl_temp$LogSkew <- NA
  neg_ctrl_temp$REST <- "NegCtrl"
  neg_ctrl_all <- rbind(neg_ctrl_all, neg_ctrl_temp)
}

neg_id_split <- colsplit(neg_ctrl_all$ID,"\\^",c("enhancer","silencer"))
neg_ctrl_all$enhancer <- neg_id_split$enhancer
neg_ctrl_all$silencer <- neg_id_split$silencer
neg_ctrl_all$REST <- 0

bound_factor <- merge(bound_factor, RESTscreen_attr_comb[,c("ID","project")], by="ID", all.x=T)
fig3 <- bound_factor[,c("ID","enhancer","silencer","celltype","log2FoldChange","scr_l2fc","LogSkew","project","REST")]
fig3 <- rbind(neg_ctrl_all, fig3)

fig3$RCOR1 <- "-"
fig3$REST <- "-"
fig3$HDAC2 <- "-"
fig3$MIER <- "-"
fig3$EHMT2 <- "-"
for(celltype in unique(fig3$celltype)){
  message(celltype)
  fig3$REST[which(fig3$celltype==celltype & fig3$silencer %in% tf_sil_list[[celltype]]$REST)] <- "+"
  fig3$RCOR1[which(fig3$celltype==celltype & fig3$silencer %in% tf_sil_list[[celltype]]$RCOR1)] <- "+"
  fig3$HDAC2[which(fig3$celltype==celltype & fig3$silencer %in% tf_sil_list[[celltype]]$HDAC2)] <- "+"
  fig3$MIER[which(fig3$celltype==celltype & fig3$silencer %in% tf_sil_list[[celltype]]$MIER1)] <- "+"
  fig3$EHMT2[which(fig3$celltype==celltype & fig3$silencer %in% tf_sil_list[[celltype]]$EHMT2)] <- "+"
}

fig3$all <- interaction(fig3$REST,fig3$MIER,fig3$EHMT2,fig3$RCOR1,fig3$HDAC2, sep="/")

fig3$REST_RCOR_HDAC <- interaction(fig3$REST,fig3$RCOR1,fig3$HDAC2, sep = "/")
fig3$REST_MIER_EHMT2 <- interaction(fig3$REST,fig3$MIER,fig3$EHMT2, sep = "/")

fig3 <- remove.factors(fig3)
fig3$all[which(fig3$project=="NegCtrl")] <- "NegCtrl"
fig3$REST_MIER_EHMT2[which(fig3$project=="NegCtrl")] <- "NegCtrl"
fig3$REST_RCOR_HDAC[which(fig3$project=="NegCtrl")] <- "NegCtrl"

# inter_stats <- list()
for(celltype in c("K562")){
  message(celltype)
  # inter_stats[[celltype]] <- list()
  for(enh in c("En19")){
    message(enh)
    # inter_stats[[celltype]][[enh]] <- list()
    plot_temp <- fig3[which(fig3$enhancer==enh & fig3$celltype==celltype),]
    if(nrow(plot_temp[which(plot_temp$project=="motif"),])==0) next
    
    ## Cycle through the sets that I need to plot along the x-axis and subset each plot_temp and skew_temp to include only that
    ## Rename the columns to have the relevant column include only the immediately relevant columns: celltype, enhancer, silencer, ID, log2FoldChange, LogSkew and relevant binding
    ## re-name columns to normalize the name of the binding column across all different sets.
    ## Make sure to get the names to match up properly for xlab and legend title
    
    # for(inter in c("REST_RCOR_HDAC", "REST_MIER_EHMT2")){
    # message(inter)
    # skew_levels <- c("-/-/-","+/-/-","+/+/-","+/-/+","+/+/+","-/+/-","-/-/+","-/+/+")
    # l2fc_levels <- c("NegCtrl","+/-/-","+/+/-","+/-/+","+/+/+","-/-/-","-/+/-","-/-/+","-/+/+")
    
    # if(inter=="REST_RCOR_HDAC"){
    #   lab_name <- "REST/RCOR/HDAC binding"
    # }
    # if(inter=="REST_MIER_EHMT2"){
    #   lab_name <- "REST/MIER1/EHMT2"
    # }
    
    
    l2fc_levels <- c("NegCtrl","+/-/-/-/-","+/+/-/-/-","+/-/+/-/-","+/-/-/+/-","+/-/-/-/+","+/+/+/-/-","+/+/-/+/-","+/+/-/-/+","+/-/+/+/-","+/-/+/-/+","+/-/-/+/+","+/-/+/+/+","+/+/-/+/+","+/+/+/-/+","+/+/+/+/-","+/+/+/+/+",
                     "-/-/-/-/-","-/+/-/-/-","-/-/+/-/-","-/-/-/+/-","-/-/-/-/+","-/+/+/-/-","-/+/-/+/-","-/+/-/-/+","-/-/+/+/-","-/-/+/-/+","-/-/-/+/+","-/-/+/+/+","-/+/+/-/+","-/+/-/+/+","-/+/+/+/-","-/+/+/+/+")
    lab_name <- "REST/MIER/EHMT2/RCOR/HDAC2"
    
    red_plot <- plot_temp[,c("ID", "enhancer","silencer", "celltype", "log2FoldChange", "scr_l2fc", "LogSkew", "project","all")] #inter instead of all if redoing things
    colnames(red_plot) <- c(colnames(red_plot)[1:8],"binding")
    # skew_temp <- red_plot[complete.cases(red_plot$LogSkew),]
    
    red_plot$binding <- factor(red_plot$binding, levels = l2fc_levels)
    
    # skew_temp$binding <- factor(skew_temp$binding, levels = skew_levels)
    
    l2fc_min <- -2.5
    l2fc_max <- 10.5
    # skew_min <- -3
    # skew_max <- 5.5
    
    neg_ctrl_med <- median(red_plot$log2FoldChange[which(red_plot$project=="NegCtrl")], na.rm = T)
    
    pdf(paste0("plots/tf_analysis/box_plots/neg_ctrl_motif_comp/",enh,"_",celltype,"_all_yrange.pdf"), width=10, height = 5) #width=40/25.4, height=45/25.4)
    print(ggplot(red_plot, aes(x=binding,y=log2FoldChange, fill=binding)) + geom_hline(yintercept = neg_ctrl_med,lwd=0.2) + geom_boxplot(outlier.shape = NA,lwd=0.2) + theme_light(base_size = 6) +
            coord_cartesian(ylim = c(l2fc_min,l2fc_max)) + xlab(lab_name) + # scale_fill_manual(values=c("grey",brewer_paired[c(1,3,5,7,2,4,6,8)])) + xlab(lab_name) +
            theme(legend.position = "none",axis.text.x = element_text(angle = 45, hjust = 1)))
    dev.off()
    # pdf(paste0("plots/tf_analysis/box_plots/skew/",enh,"_",celltype,"_",inter,"_yrange.pdf"), width=40/25.4, height=45/25.4)
    # print(ggplot(skew_temp, aes(x=binding,y=LogSkew, fill=binding)) + geom_hline(yintercept = 0,lwd=0.2) + geom_boxplot(outlier.shape = NA, lwd=0.2) + theme_light(base_size = 6) +
    #         coord_cartesian(ylim = c(skew_min,skew_max)) + scale_fill_manual(values=c(brewer_paired[c(2,1,3,5,7,4,6,8)])) + xlab(lab_name) +
    #         theme(legend.position = "none",axis.text.x = element_text(angle = 45, hjust = 1)))
    # dev.off()
    
    # }
    
    # message("l2fc sig")
    # inter_stats[[celltype]][[enh]][[inter]]$l2fc <- red_plot %>%
    #   wilcox_test(formula = log2FoldChange ~ binding, ref.group = "NegCtrl") %>%
    #   adjust_pvalue(method = "BH") %>%
    #   add_significance()
    # 
    # message("skew sig")
    # if(inter=="CoREST"){
    #   inter_stats[[celltype]][[enh]][[inter]]$skew <- skew_temp %>%
    #     wilcox_test(formula = LogSkew ~ binding, ref.group = "-/-") %>%
    #     adjust_pvalue(method = "BH") %>%
    #     add_significance()
    # }
    # if(inter!="CoREST"){
    #   inter_stats[[celltype]][[enh]][[inter]]$skew <- skew_temp %>%
    #     wilcox_test(formula = LogSkew ~ binding, ref.group = "-/-/-") %>%
    #     adjust_pvalue(method = "BH") %>%
    #     add_significance()
    # }
  }
}

fig3$TRIP13 <- "-"
fig3$PTTG1 <- "-"
for(celltype in unique(fig3$celltype)){
  message(celltype)
  fig3$TRIP13[which(fig3$celltype==celltype & fig3$silencer %in% tf_sil_list[[celltype]]$TRIP13)] <- "+"
  fig3$PTTG1[which(fig3$celltype==celltype & fig3$silencer %in% tf_sil_list[[celltype]]$PTTG1)] <- "+"
  
}

fig3$score <- "NoScore"
for(celltype in unique(fig3$celltype)){
  message(celltype)
  for(enh in enh_all){
    message(enh)
    for(relat in names(oligo_split_list[[celltype]][[enh]])){
      fig3$score[which(fig3$silencer %in% oligo_split_list[[celltype]][[enh]][[relat]] & 
                         fig3$celltype==celltype & fig3$enhancer==enh)] <- relat
    }
  }
}

fig3$score[which(fig3$score=="NoScore" & fig3$project=="motif")] <- "below"



for(celltype in unique(fig3$celltype)){
  message(celltype)
  for(enh in enh_all){
    message(enh)
    for(rel in c("above","below")){
      message(rel)
      plot_temp <- fig3[which(fig3$celltype==celltype & fig3$enhancer==enh & fig3$score %in% c("NoScore",rel)),]
      message(nrow(plot_temp))
      for(inter in c("REST_RCOR_HDAC", "REST_MIER_EHMT2")){
        message(inter)
        skew_levels <- c("-/-/-","+/-/-","+/+/-","+/-/+","+/+/+","-/+/-","-/-/+","-/+/+")
        l2fc_levels <- c("NegCtrl","+/-/-","+/+/-","+/-/+","+/+/+","-/-/-","-/+/-","-/-/+","-/+/+")
        
        if(inter=="REST_RCOR_HDAC"){
          lab_name <- "REST/RCOR/HDAC binding"
        }
        if(inter=="REST_MIER_EHMT2"){
          lab_name <- "REST/MIER1/EHMT2"
        }
        
        red_plot <- plot_temp[,c("ID", "enhancer","silencer", "celltype", "log2FoldChange", "scr_l2fc", "LogSkew", "project",inter)]
        message(nrow(red_plot))
        colnames(red_plot) <- c("ID", "enhancer","silencer", "celltype", "log2FoldChange", "scr_l2fc", "LogSkew", "project","binding")
        skew_temp <- red_plot[complete.cases(red_plot$LogSkew),]
        
        red_plot$binding <- factor(red_plot$binding, levels = l2fc_levels)
        
        skew_temp$binding <- factor(skew_temp$binding, levels = skew_levels)
        
        # l2fc_min <- -2.5
        # l2fc_max <- 10.5
        # skew_min <- -3
        # skew_max <- 5.5
        
        neg_ctrl_med <- median(red_plot$log2FoldChange[which(red_plot$project=="NegCtrl")], na.rm = T)
        
        pdf(paste0("plots/tf_analysis/box_plots/neg_ctrl_motif_comp/fig_S8_above_below/",enh,"_",celltype,"_",inter,"_",rel,".pdf"), width=40/25.4, height=45/25.4)
        print(ggplot(red_plot, aes(x=binding,y=log2FoldChange, fill=binding)) + geom_hline(yintercept = neg_ctrl_med,lwd=0.2) + geom_boxplot(outlier.shape = NA,lwd=0.2) + theme_light(base_size = 6) +
                # coord_cartesian(ylim = c(l2fc_min,l2fc_max)) + 
                scale_fill_manual(values=c("grey",brewer_paired[c(1,3,5,7,2,4,6,8)])) + xlab(lab_name) +
                theme(legend.position = "none",axis.text.x = element_text(angle = 45, hjust = 1)))
        dev.off()
        pdf(paste0("plots/tf_analysis/box_plots/skew/fig_s8_above_below/",enh,"_",celltype,"_",inter,"_",rel,".pdf"), width=40/25.4, height=45/25.4)
        print(ggplot(skew_temp, aes(x=binding,y=LogSkew, fill=binding)) + geom_hline(yintercept = 0,lwd=0.2) + geom_boxplot(outlier.shape = NA, lwd=0.2) + theme_light(base_size = 6) +
                # coord_cartesian(ylim = c(skew_min,skew_max)) + 
                scale_fill_manual(values=c(brewer_paired[c(2,1,3,5,7,4,6,8)])) + xlab(lab_name) +
                theme(legend.position = "none",axis.text.x = element_text(angle = 45, hjust = 1)))
        dev.off()
        # 
      }
    }
  }
}

fig3$ZNF644 <- "-"

for(celltype in unique(fig3$celltype)){
  message(celltype)
  fig3$ZNF644[which(fig3$celltype==celltype & fig3$silencer %in% tf_sil_list[[celltype]]$ZNF644)] <- "+"
  
}
