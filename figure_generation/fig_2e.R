#plot figure 2e see README for explanation of what is happening

oligo_split_list <- list()
for(celltype in names(split_lower)){
  oligo_split_list[[celltype]] <- list()
  for(enh in names(split_lower[[celltype]])){
    oligo_split_list[[celltype]][[enh]] <- list()
    split_val <- 20.86
    if(split_val < 25 & split_val > 15){
      oligo_split_list[[celltype]][[enh]]$below <- emvar_scores_motif_df$silencer[which(emvar_scores_motif_df$celltype==celltype &
                                                                                          emvar_scores_motif_df$enhancer==enh &
                                                                                          emvar_scores_motif_df$score < split_val)]
      oligo_split_list[[celltype]][[enh]]$above <- emvar_scores_motif_df$silencer[which(emvar_scores_motif_df$celltype==celltype &
                                                                                          emvar_scores_motif_df$enhancer==enh &
                                                                                          emvar_scores_motif_df$score >= split_val)]
    }
    if(split_val > 25 | split_val < 15){
      oligo_split_list[[celltype]][[enh]]$below <- emvar_scores_motif_df$silencer[which(emvar_scores_motif_df$celltype==celltype &
                                                                                          emvar_scores_motif_df$enhancer==enh )]
    }
  }
}



scr_skew_stats <- list()
for(enh in enh_all){
  message(enh)
  scale_max <- 9
  if(enh=="En02"){
    scale_min <- -6
  }
  if(enh=="En09"){
    scale_min <- -2
  }
  if(enh=="En11"){
    scale_min <- -5
  }
  if(enh=="En19"){
    scale_min <- -2
    scale_max <- 8
  }
  if(enh=="En21"){
    scale_min <- 0
    scale_max <- 10
  }
  to_plot <- as.data.frame(all_38_enh_res[[enh]])
  to_plot <- to_plot[which(to_plot$proj %in% c("NegCtrl","motif_ChIP_pos","motif_ChIP_neg")),]
  # plot_split <- colsplit(rownames(to_plot),"\\^", c("enhancer","silencer"))
  # sil_split <- colsplit(plot_split$silencer,"\\.", c("silencer","rep"))
  # to_plot$sil <- sil_split$silencer
  to_plot$group <- "NoScore"
  to_plot$proj <- factor(to_plot$proj, levels = c("NegCtrl","motif_ChIP_pos","motif_ChIP_neg"))
  # to_plot$group <- factor(to_plot$group, levels = c("above", "below", "NoScore"))
  for(celltype in unique(to_plot$celltype)){
    message(celltype)
    for(type in names(oligo_split_list[[celltype]][[enh]])){
      message(type)
      to_plot$group[which(to_plot$sil %in% oligo_split_list[[celltype]][[enh]][[type]] & to_plot$celltype==celltype)] <- type
    }
  }
  to_plot$group <- factor(to_plot$group, levels = c("above", "below", "NoScore"))
  to_plot$group[which(to_plot$proj!="NegCtrl" & to_plot$group=="NoScore")] <- "below"
  to_plot$inter <- interaction(to_plot$proj, to_plot$group)
  to_plot$inter <- factor(to_plot$inter, levels = c("NegCtrl.NoScore", 
                                                    "motif_ChIP_pos.below","motif_ChIP_pos.above",
                                                    "motif_ChIP_neg.below","motif_ChIP_neg.above"))
  scr_skew_stats[[enh]] <- list()
  for(celltype in unique(to_plot$celltype)){
    # message(celltype)
    to_stat <- to_plot
    to_stat <- remove.factors(to_stat)
    stat_temp <- wilcox_test(to_stat[which(to_stat$celltype==celltype),], log2FoldChange~inter, ref.group = "NegCtrl.NoScore", p.adjust.method = "BH")
    scr_skew_stats[[enh]][[celltype]] <- as.data.frame(stat_temp)
  }
  
  med_plot <- to_plot[which(to_plot$proj=="NegCtrl"),] %>%
    group_by(celltype) %>%
    summarise(med_line = median(log2FoldChange, na.rm = T))
  
  pdf(paste0("plots/log2foldchange_box_plots/",enh,"_piecewise_box_swarm_alpha.pdf"), width=174/25.4, height = (174/3.25)/25.4)
  print(ggplot(to_plot, aes(x=inter, y=log2FoldChange, fill=inter, group=inter)) + geom_hline(data=med_plot, aes(yintercept=med_line)) + geom_quasirandom(aes(col=inter), outlier.shape=NA, size=0.5, method = "smiley", alpha=0.2) + geom_boxplot(outlier.shape = NA) +
          ylim(scale_min,scale_max) + facet_grid(~celltype) + theme_light(base_size = 6) + xlab("") + ylab("log2FoldChange") +
          theme(axis.text.x = element_blank(), legend.position = "bottom", strip.background = element_blank(), strip.text.x = element_blank(),legend.margin=margin(t = -.5, unit='cm'), panel.grid = element_line(size = 0.24), panel.border = element_rect(size = 0.24), legend.text = element_text(size=4)) +
          guides(fill=guide_legend(nrow = 1), shape=guide_legend(override.aes = list(size = 0.5)), color=guide_legend(override.aes = list(size = 0.5))) +
          scale_fill_manual(values = c("grey","#66C2A5","#FC8D62","#8DA0CB","#E78AC3","#A6D854","#FFD92F")) + scale_color_manual(values = c("grey","#66C2A5","#FC8D62","#8DA0CB","#E78AC3","#A6D854","#FFD92F")))
  dev.off()
}