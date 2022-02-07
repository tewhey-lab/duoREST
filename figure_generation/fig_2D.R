#plot figure 2D see README for explanation of what is happening

model_list <- list()
for(celltype in unique(emvar_scores_motif_df$celltype)){
  message(celltype)
  model_list[[celltype]] <- list()
  for(enh in enh_all){
    message(enh)
    model_list[[celltype]][[enh]]
    cell_enh_score <- emvar_scores_motif_df$score[which(emvar_scores_motif_df$celltype==celltype & emvar_scores_motif_df$enhancer==enh)]
    cell_enh_skew <- emvar_scores_motif_df$LogSkew[which(emvar_scores_motif_df$celltype==celltype & emvar_scores_motif_df$enhancer==enh)]
    lin_mod <- lm(cell_enh_skew ~ cell_enh_score)
    seg_mod <- segmented(lin_mod, seg.Z=~cell_enh_score, psi = c(20))
    model_list[[celltype]][[enh]]$linear <- lin_mod
    model_list[[celltype]][[enh]]$segmented <- seg_mod
  }
}

reg_lines_plot_both <- function(celltype, enh, df){
  group_temp <- df[which(df[,"celltype"]==celltype & df[,"enhancer"]==enh),c("celltype","enhancer","silencer","LogSkew","score","skew_group")]
  
  piecewise_fitted <- fitted(model_list[[celltype]][[enh]]$segmented)
  linear_fitted <- fitted(model_list[[celltype]][[enh]]$linear)
  
  message("Fitted Models Complete")
  
  p.model.temp <- data.frame(score=group_temp$score, LogSkew = piecewise_fitted)
  l.model.temp <- data.frame(score=group_temp$score, LogSkew = linear_fitted)
  
  message(colnames(p.model.temp))
  
  message("DFs set up")
  p <- ggplot(group_temp, aes(x=score, y=LogSkew)) + geom_point(aes(col=skew_group), alpha=0.5, size=0.5) + xlim(0,35) + 
    geom_line(l.model.temp, mapping=aes(x=score, y=LogSkew), col="purple", lwd=0.6) + 
    geom_line(p.model.temp, mapping=aes(x=score, y=LogSkew), col="salmon", lwd=0.6) +
    theme_light(base_size = 6) + theme(legend.position = "none", axis.title.y = element_text(size=4),
                                       axis.title.x = element_text(size = 4), axis.text = element_text(size=3))
  return(p)
}

for(celltype in names(model_list)){
  message(celltype)
  for(enh in names(model_list[[celltype]])){
    message(enh)
    p <- reg_lines_plot_both(celltype, enh, emvar_scores_motif_df)
    
    pdf(paste0("plots/piecewise_v_linear/", celltype, "_", enh, "_segmented.pdf"), width = 45/25.4, height = 45/25.4)
    print(p)
    dev.off()
  }
}
