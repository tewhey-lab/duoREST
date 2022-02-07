### Gene - TF Analysis ###

cell_binding <- fig3[,c(4,2,3,1,8:9,11,13,16:17,19)]

tpm_df <- read.delim("nearest_gene/05282021REST_ref_RNAseq_TPM.allgene.score.txt", stringsAsFactors = F, header = T)
colnames(tpm_df)[5] <- "TPM_SK.N.SH"

mean_breakpoints <- data.frame()
for(celltype in unique(model_comp_sum$celltype)){
  mbp <- mean(model_comp_sum$seg_breakpoints[which(abs(model_comp_sum$seg_breakpoints - breakpoints_mean) < breakpoints_sd*2 & model_comp_sum$celltype==celltype)])
  temp <- data.frame(celltype, mbp)
  mean_breakpoints <- rbind(mean_breakpoints, temp)
}

head(cell_binding)
cell_binding <- merge(cell_binding, emvar_scores_motif_df[,c(19:22)], by=c("celltype","enhancer","silencer"), all=T)

cell_binding <- cell_binding[which(cell_binding$project!="NegCtrl"),]

tpm_df_cp <- tpm_df
tpm_sil_sep <- setDT(tpm_df_cp)[,lapply(.SD, function(x) unlist(tstrsplit(x,",", fixed=T)))]
table(is.na(tpm_sil_sep$gene_name))

tpm_sil_sep$nearest <- "No"
tpm_sil_sep$nearest[which(!is.na(tpm_sil_sep$silencer))] <- "Yes"

tpm_sil_sep <- data.frame(tpm_sil_sep)

tpm_list2 <- list()
for(celltype in c("GM12878","HepG2","K562","SK.N.SH")){
  message(celltype)
  grab_cols <- c("gene_name","silencer","nearest")
  tpm_col <- paste("TPM_",celltype,sep = "")
  grab_cols <- c(grab_cols,tpm_col)
  tpm_list2[[celltype]] <- tpm_sil_sep[,grab_cols]
  tpm_list2[[celltype]]$REST <- NA
  tpm_list2[[celltype]]$REST[which(tpm_list2[[celltype]]$silencer %in% 
                                     cell_binding$silencer[which(cell_binding$REST=="+" & cell_binding$celltype==celltype)])] <- "+"
  tpm_list2[[celltype]]$REST[which(tpm_list2[[celltype]]$silencer %in% 
                                     cell_binding$silencer[which(cell_binding$REST=="-" & cell_binding$celltype==celltype)])] <- "-"
  tpm_list2[[celltype]]$EHMT2 <- NA
  tpm_list2[[celltype]]$EHMT2[which(tpm_list2[[celltype]]$silencer %in% 
                                      cell_binding$silencer[which(cell_binding$EHMT2=="+" & cell_binding$celltype==celltype)])] <- "+"
  tpm_list2[[celltype]]$EHMT2[which(tpm_list2[[celltype]]$silencer %in% 
                                      cell_binding$silencer[which(cell_binding$EHMT2=="-" & cell_binding$celltype==celltype)])] <- "-"
  tpm_list2[[celltype]]$HDAC2 <- NA
  tpm_list2[[celltype]]$HDAC2[which(tpm_list2[[celltype]]$silencer %in% 
                                      cell_binding$silencer[which(cell_binding$HDAC2=="+" & cell_binding$celltype==celltype)])] <- "+"
  tpm_list2[[celltype]]$HDAC2[which(tpm_list2[[celltype]]$silencer %in% 
                                      cell_binding$silencer[which(cell_binding$HDAC2=="-" & cell_binding$celltype==celltype)])] <- "-"
  tpm_list2[[celltype]]$TRIP13 <- NA
  tpm_list2[[celltype]]$TRIP13[which(tpm_list2[[celltype]]$silencer %in% 
                                       cell_binding$silencer[which(cell_binding$TRIP13=="+" & cell_binding$celltype==celltype)])] <- "+"
  tpm_list2[[celltype]]$TRIP13[which(tpm_list2[[celltype]]$silencer %in% 
                                       cell_binding$silencer[which(cell_binding$TRIP13=="-" & cell_binding$celltype==celltype)])] <- "-"
  tpm_list2[[celltype]]$PTTG1 <- NA
  tpm_list2[[celltype]]$PTTG1[which(tpm_list2[[celltype]]$silencer %in% 
                                      cell_binding$silencer[which(cell_binding$PTTG1=="+" & cell_binding$celltype==celltype)])] <- "+"
  tpm_list2[[celltype]]$PTTG1[which(tpm_list2[[celltype]]$silencer %in% 
                                      cell_binding$silencer[which(cell_binding$PTTG1=="-" & cell_binding$celltype==celltype)])] <- "-"
  tpm_list2[[celltype]]$ZNF644 <- NA
  tpm_list2[[celltype]]$ZNF644[which(tpm_list2[[celltype]]$silencer %in% 
                                       cell_binding$silencer[which(cell_binding$ZNF644=="+" & cell_binding$celltype==celltype)])] <- "+"
  tpm_list2[[celltype]]$ZNF644[which(tpm_list2[[celltype]]$silencer %in% 
                                       cell_binding$silencer[which(cell_binding$ZNF644=="-" & cell_binding$celltype==celltype)])] <- "-"
}


for(celltype in names(tpm_list2)){
  message(celltype)
  temp <- tpm_list2[[celltype]]
  temp$relative <- "No Score"
  for(sil in tpm_list2[[celltype]]$silencer){
    temp_bp <-20.86
    sil_score <- cell_binding$score[which(cell_binding$silencer==sil & cell_binding$celltype==celltype)]
    sil_score <- mean(sil_score, na.rm=T)
    if(is.na(sil_score)) next
    if(sil_score < temp_bp){
      temp$relative[which(temp$silencer==sil)] <- "below"
    }
    if(sil_score > temp_bp){
      temp$relative[which(temp$silencer==sil)] <- "above"
    }
  }
  tpm_list2[[celltype]] <- temp
}


tpm_log_df <- data.frame()
for(celltype in names(tpm_list2)){
  temp <- tpm_list2[[celltype]]
  colnames(temp)[4] <- "TPM"
  temp$celltype <- celltype
  tpm_log_df <- rbind(tpm_log_df, temp)
}

tpm_log_df <- tpm_log_df[complete.cases(tpm_log_df$TPM),]
tpm_log_df$TPM <- as.numeric(tpm_log_df$TPM)
tpm_log_df$logTPM <- log(tpm_log_df$TPM)

tpm_log_df_fin <- tpm_log_df[is.finite(tpm_log_df$logTPM),]
summary(tpm_log_df_fin)

tpm_log_df_fin$relative[which(tpm_log_df_fin$relative=="No Score")] <- "below"

tpm_log_all_df <- tpm_log_df_fin
tpm_log_all_df$relative <- "All Genes"

tpm_log_df_fin$relative[which(tpm_log_df_fin$REST=="-")] <- "REST-"

tpm_log_all_df <- rbind(tpm_log_all_df, tpm_log_df_fin)
tpm_log_all_df$EHMT2_group <- interaction(tpm_log_all_df$relative, tpm_log_all_df$EHMT2, sep = "/")
tpm_log_all_df$HDAC2_group <- interaction(tpm_log_all_df$relative, tpm_log_all_df$HDAC2, sep = "/")
tpm_log_all_df$TRIP13_group <- interaction(tpm_log_all_df$relative, tpm_log_all_df$TRIP13, sep = "/")
tpm_log_all_df$PTTG1_group <- interaction(tpm_log_all_df$relative, tpm_log_all_df$PTTG1, sep = "/")
tpm_log_all_df$ZNF644_group <- interaction(tpm_log_all_df$relative, tpm_log_all_df$ZNF644, sep = "/")


tpm_log_df_fin <- remove.factors(tpm_log_df_fin)
tpm_log_df_fin$EHMT2_group <- interaction(tpm_log_df_fin$relative, tpm_log_df_fin$EHMT2, sep = "/")
tpm_log_df_fin$HDAC2_group <- interaction(tpm_log_df_fin$relative, tpm_log_df_fin$HDAC2, sep = "/")
tpm_log_df_fin$TRIP13_group <- interaction(tpm_log_df_fin$relative, tpm_log_df_fin$TRIP13, sep = "/")
tpm_log_df_fin$PTTG1_group <- interaction(tpm_log_df_fin$relative, tpm_log_df_fin$PTTG1, sep = "/")
tpm_log_df_fin$ZNF644_group <- interaction(tpm_log_df_fin$relative, tpm_log_df_fin$ZNF644, sep = "/")

tpm_log_df_fin <- remove.factors(tpm_log_df_fin)

tpm_log_df_fin$EHMT2_group[is.na(tpm_log_df_fin$silencer)] <- "No Nearest"
tpm_log_df_fin$HDAC2_group[is.na(tpm_log_df_fin$silencer)] <- "No Nearest"
tpm_log_df_fin$TRIP13_group[is.na(tpm_log_df_fin$silencer)] <- "No Nearest"
tpm_log_df_fin$PTTG1_group[is.na(tpm_log_df_fin$silencer)] <- "No Nearest"
tpm_log_df_fin$ZNF644_group[is.na(tpm_log_df_fin$silencer)] <- "No Nearest"

tpm_log_df_fin$relative[is.na(tpm_log_df_fin$silencer)] <- "No Nearest"

table(tpm_log_df_fin[,c("EHMT2_group","celltype")])



for(celltype in names(tpm_list2)){
  message(celltype)
  for(tf in c("ZNF644","PTTG1","TRIP13","HDAC2","EHMT2")){
    message(tf)
    plot_temp <- tpm_log_all_df[which(tpm_log_all_df$celltype==celltype),c("REST",tf,"relative","celltype","logTPM",paste0(tf,"_group"))]
    colnames(plot_temp)[2] <- "binding"
    colnames(plot_temp)[6] <- "grouping"
    plot_temp <- remove.factors(plot_temp)
    plot_temp$grouping[which(plot_temp$REST=="-")] <- "REST-"
    plot_temp$grouping[which(plot_temp$relative=="All Genes")] <- "All Genes"
    plot_temp <- plot_temp[complete.cases(plot_temp$grouping),]
    plot_temp$relative <- factor(plot_temp$relative, levels = c("All Genes", "REST-","below","above"))
    plot_temp$grouping <- factor(plot_temp$grouping, levels = c("All Genes", "REST-", "below/-","below/+","above/-","above/+"))
    pdf(paste0("plots/nearest_gene/", celltype,"_",tf,"_violin_mval.pdf"), width=45/25.4, height=45/25.4)
    print(ggplot(plot_temp, aes(x=grouping, y=logTPM)) + geom_hline(yintercept = 0) + geom_violin(aes(fill=grouping)) +
            geom_boxplot(outlier.shape = NA, width=0.25) + theme_light(base_size = 6, base_line_size = .24) + ylab("ln(TPM)") + xlab("") +
            theme(legend.position = "none") + scale_fill_manual(values=tf_grouping_colors) + scale_color_manual(values = tf_grouping_colors))
    dev.off()
  }
}

# geom_boxplot(position=position_dodge2(width = 1, preserve="single", padding=0.3),

tpm_stats <- list()
for(celltype in names(tpm_list)){
  message(celltype)
  tpm_stats[[celltype]] <- list()
  for(tf in c("ZNF644","PTTG1","TRIP13","HDAC2","EHMT2")){
    message(tf)
    plot_temp <- tpm_log_df_fin[which(tpm_log_df_fin$celltype==celltype),c("REST",tf,"relative","celltype","logTPM",paste0(tf,"_group"))]
    colnames(plot_temp)[2] <- "binding"
    colnames(plot_temp)[6] <- "grouping"
    plot_temp <- remove.factors(plot_temp)
    plot_temp$grouping[which(plot_temp$REST=="-")] <- "REST-"
    # plot_temp$grouping[which(plot_temp$relative=="All Genes")] <- "All Genes"
    plot_temp <- plot_temp[complete.cases(plot_temp$grouping),]
    plot_temp$relative <- factor(plot_temp$relative, levels = c("No Nearest", "REST-","No Score","below","above"))
    plot_temp$grouping <- factor(plot_temp$grouping, levels = c("No Nearest", "REST-","below/+","below/-","above/+","above/-"))
    
    to_stat <- plot_temp
    to_stat <- remove.factors(to_stat)
    if(nrow(to_stat[which(to_stat$relative=="REST-"),])>0){
      stat_temp <- wilcox_test(to_stat, logTPM~grouping, ref.group = "REST-", p.adjust.method = "BH")
      tpm_stats[[celltype]][[tf]] <- as.data.frame(stat_temp)
    }
    
    pdf(paste0("plots/nearest_gene/no_nearest_", celltype,"_",tf,"_mval.pdf"), width=174/25.4, height=50/25.4)
    print(ggplot(plot_temp, aes(x=relative, y=logTPM, group=grouping, fill=grouping)) + geom_hline(yintercept = 0) + geom_quasirandom(dodge.width = dodge$width ,aes(col=grouping), outlier.shape=NA, size=0.5, method = "smiley", alpha=0.2) +
            geom_boxplot(position=dodge,outlier.shape = NA) + theme_light(base_size = 6, base_line_size = .24) + ylab("ln(TPM)") + xlab("") +
            theme(legend.position = "bottom", legend.margin=margin(t = -.5, unit='cm'), legend.text = element_text(size=4)) +
            guides(fill=guide_legend(nrow = 1), shape=guide_legend(override.aes = list(size = 0.5)), color=guide_legend(override.aes = list(size = 0.5))))
    dev.off()
  }
}

tpm_stats_df <- data.frame()
for(celltype in names(tpm_stats)){
  message(celltype)
  for(tf in names(tpm_stats[[celltype]])){
    message(tf)
    temp <- as.data.frame(tpm_stats[[celltype]][[tf]])
    temp$celltype <- celltype
    temp$tf <- tf
    tpm_stats_df <- rbind(tpm_stats_df, temp)
  }
}

write.table(tpm_stats_df,"results/TF_bound_unbound_analysis/log_tpm_stats_mval.txt", quote=F, sep = "\t", row.names = F)
