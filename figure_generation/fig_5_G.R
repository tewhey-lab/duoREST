emvar_comp <- data.frame()
for(celltype in names(emVAR_all)){
  keep <- emVAR_all[[celltype]][which(emVAR_all[[celltype]]$SNP %in% SNP_incl & !is.na(emVAR_all[[celltype]]$ref_allele)), c("ID","SNP","LogSkew","Skew_logFDR")]
  keep$celltype <- celltype
  emvar_comp <- rbind(emvar_comp, keep)
}

emvar_comp <- merge(emvar_comp[,1:10], maf_rare_scores, by.x="SNP", by.y="snp", all.x=T)

emvar_comp$DBS <- emvar_comp$score_a - emvar_comp$score_r

emvar_comp$maf <- ">1"
emvar_comp$maf[which(emvar_comp$SNP %in% pop_AF$SNP[which(pop_AF$max_AF < 0.01)])] <- "<1"
emvar_comp$maf[which(emvar_comp$SNP %in% pop_AF$SNP[which(pop_AF$max_AF > 0.05)])] <- ">5"
emvar_comp$maf[which(emvar_comp$SNP %in% pop_AF$SNP[which(pop_AF$max_AF > 0.1)])] <- ">10"

emvar_comp$maf <- factor(emvar_comp$maf, levels = c("<1",">1",">5"))

emvar_comp$window <- "outside"
emvar_comp$window[which(emvar_comp$loc_a=="inside" & emvar_comp$loc_r=="inside")] <- "inside"

comp_sil_split <- colsplit(emvar_comp$silencer,":wP",c("sil_allele","wP"))
# emvar_comp$window[which(emvar_comp$loc_a!=emvar_comp$loc_r)] <- ifelse(comp_sil_split$wP[which(emvar_comp$loc_a!=emvar_comp$loc_r)] < 112 & comp_sil_split$wP[which(emvar_comp$loc_a!=emvar_comp$loc_r)] > 90, "inside","outside")
emvar_comp$window <- ifelse(comp_sil_split$wP < 112 & comp_sil_split$wP > 90, "inside","outside")

emvar_comp$rel <- "below"
emvar_comp$rel[which(emvar_comp$score_a >= 20.86 | emvar_comp$score_r >= 20.86)] <- "above"

#### Figure 5G ####
emvar_comp_stat <- list()
for(celltype in unique(emvar_comp$celltype)){
  emvar_comp_stat[[celltype]] <- list()
  message(celltype)
  for(enh in enh_all){
    emvar_comp_stat[[celltype]][[enh]] <- list()
    message(enh)
    for(rel_pos in c("above","below")){
      emvar_comp_stat[[celltype]][[enh]][[rel_pos]] <- list()
      message(rel_pos)
      for(win in c("inside","outside")){
        message(win)
        to_plot <- emvar_comp[which(emvar_comp$celltype==celltype & emvar_comp$enhancer==enh & emvar_comp$maf!=">1" & emvar_comp$window==win & emvar_comp$rel==rel_pos),]
        
        pdf(paste0("plots/minor_allele_frequency_analysis/density/",celltype,"_",enh,"_maf_density_",rel_pos,"_",win,".pdf"), width=40/25.4, height=40/25.4)
        print(ggplot(to_plot, aes(x=LogSkew, group=maf, fill=maf)) + geom_density(adjust=1.5, alpha=0.2, lwd=0.1) + xlim(-5,5) +
                theme_minimal(base_size = 6) + theme(legend.position = "none"))
        dev.off()
        
        to_stat <- remove.factors(to_plot)
        stat_temp <- wilcox_test(to_stat, LogSkew~maf, ref.group = "<1", p.adjust.method = "BH")
        emvar_comp_stat[[celltype]][[enh]][[rel_pos]][[win]] <- stat_temp
        
        # pdf(paste0("plots/minor_allele_frequency_analysis/density/",celltype,"_",enh,"_maf_density_",rel_pos,"_",win,".pdf"), width=40/25.4, height=40/25.4)
        # print(ggplot(to_plot, aes(x=LogSkew, group=maf3b, fill=maf3b)) + geom_density(adjust=1.5, alpha=0.2, lwd=0.1) +
        #         theme_minimal(base_size = 6) + theme(legend.position = "none"))
        # dev.off()
        
      }
    }
  }
}

emvar_comp_stat_df <- data.frame()
for(celltype in names(emvar_comp_stat)){
  message(celltype)
  for(enhancer in names(emvar_comp_stat[[celltype]])){
    message(enhancer)
    for(relativePosition in names(emvar_comp_stat[[celltype]][[enhancer]])){
      message(relativePosition)
      for(window in names(emvar_comp_stat[[celltype]][[enhancer]][[relativePosition]])){
        message(window)
        stat_temp <- data.frame(emvar_comp_stat[[celltype]][[enhancer]][[relativePosition]][[window]])
        stat_temp$celltype <- celltype
        stat_temp$enhancer <- enhancer
        stat_temp$relativePosition <- relativePosition
        stat_temp$window <- window
        emvar_comp_stat_df <- rbind(emvar_comp_stat_df, stat_temp)
      }
    }
  }
}

emvar_comp_stat_df$ks_p <- NA
for(celltype in unique(emvar_comp_stat_df$celltype)){
  message(celltype)
  for(enh in unique(emvar_comp_stat_df$enhancer)){
    message(enh)
    for(rel in unique(emvar_comp_stat_df$relativePosition)){
      message(rel)
      for(win in unique(emvar_comp_stat_df$window)){
        message(win)
        g1 <- emvar_comp$LogSkew[which(emvar_comp$celltype==celltype & emvar_comp$enhancer==enh & emvar_comp$rel==rel & emvar_comp$window==win & emvar_comp$maf3b=="<1")]
        g2 <- emvar_comp$LogSkew[which(emvar_comp$celltype==celltype & emvar_comp$enhancer==enh & emvar_comp$rel==rel & emvar_comp$window==win & emvar_comp$maf3b==">5")]
        
        emvar_comp_stat_df$ks_p[which(emvar_comp_stat_df$celltype==celltype & emvar_comp_stat_df$enhancer==enh & emvar_comp_stat_df$relativePosition==rel & emvar_comp_stat_df$window==win)] <- ks.test(g1,g2)$p.value
      }
    }
  }
}
