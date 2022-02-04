This is a repository to keep supplementary files and analysis steps in one place for the **Whole genome functional characterization of RE1 silencers using a modified massively parallel reporter assay** paper by Mouri et all.

<details open>
<summary>Initial analysis steps</summary>
<br>
  This step uses the DUOmodel functions in the <a href="https://github.com/tewhey-lab/MPRAduo" target="_top">MPRAduo</a> pipeline
  
    RESTscreen_condition <- as.data.frame(c(rep("DNA",4), rep("GM12878",4), rep("K562",4), rep("HepG2",4), rep("SKNSH", 4)), stringsAsFactors=F)

    colnames(RESTscreen_condition) <- "condition"

    rownames(RESTscreen_condition) <- colnames(RESTscreen_combined_fixed)[4:23]

    RESTscreen_counts <- read.delim("01092020_REST_screen_fixed_oligos.count", stringsAsFactors = F)

    RESTscreen_counts <- duoStats(RESTscreen_counts,RESTscreen_condition)

    RESTscreen_names <- duoNames(RESTscreen_counts,RESTscreen_counts,libExcl1 = c("E","S"), libExcl2 = c("E","S"),duoOnly = T)

    RESTscreen_attributes <- read.delim("attributes/RESTscreen_full.attributes", stringsAsFactors = F)
    
    enh_all <- c("En02", "En09", "En11", "En19", "En21")

    SNC <- RESTscreen_attributes$ID[which(RESTscreen_attributes$project=="NegCtrl")]

    RESTscreen_attr <- duoAttr(RESTscreen_attributes, enhList = enh_all)

    ENC <- "En02"

    RESTscreen_seq <- duoSeq(RESTscreen_attr,RESTscreen_counts, RESTscreen_condition, 1, "20201207_RESTscreen",RESTscreen_names, negListE = ENC, negListS = SNC, libExcl = c("E","S"))
    
    # QC plots
    duoCor(dataCount=as.data.frame(counts(RESTscreen_seq$ES_1)),dataCond=RESTscreen_condition,namesList = RESTscreen_names,filePrefix = "agg_counts_all_reps",run=1,libExcl = c("E","S"), libIncl = "ES_1")

    duoLogCor(RESTscreen_seq, RESTscreen_condition, "20201207_RESTscreen")
</details>

<details open>
<summary>Read results files back into environment</summary>
<br>
  
    emVAR_all <- list()
    emVAR_all$GM12878 <- read.delim("results/20201207_RESTscreen_GM12878_emVAR.out", stringsAsFactors=F)
    emVAR_all$HepG2 <- read.delim("results/20201207_RESTscreen_HepG2_emVAR.out", stringsAsFactors=F)
    emVAR_all$K562 <- read.delim("results/20201207_RESTscreen_K562_emVAR.out", stringsAsFactors=F)
    emVAR_all$SKNSH <- read.delim("results/20201207_RESTscreen_SKNSH_emVAR.out", stringsAsFactors=F)
                       
    standard_res <- list()
    standard_res$GM12878 <- read.delim("results/20201207_RESTscreen_GM12878_results.run1.txt", stringsAsFactors = F)
    standard_res$HepG2 <- read.delim("results/20201207_RESTscreen_HepG2_results.run1.txt", stringsAsFactors = F)
    standard_res$K562 <- read.delim("results/20201207_RESTscreen_K562_results.run1.txt", stringsAsFactors = F)
    standard_res$SKNSH <- read.delim("results/20201207_RESTscreen_SKNSH_results.run1.txt", stringsAsFactors = F)
</details>

<details open>
<summary>Categorize ChIP positivity of oligos</summary>
<br>

    ref_2.0 <- list()
    ref_2.0$GM12878 <- read.delim("cell_specificity/Ref_GM12878.wa.txt", header = F, stringsAsFactors = F)
    ref_2.0$HepG2 <- read.delim("cell_specificity/Ref_HepG2.wa.txt", header = F, stringsAsFactors = F)
    ref_2.0$K562 <- read.delim("cell_specificity/Ref_K562.wa.txt", header = F, stringsAsFactors = F)
    ref_2.0$SK.N.SH <- read.delim("cell_specificity/Ref_SKNSH.wa.txt", header = F, stringsAsFactors = F)
    for(celltype in names(ref_2.0)){
      colnames(ref_2.0[[celltype]]) <- "ID"
    }

    for(celltype in names(ref_2.0)){
      under_split <- colsplit(ref_2.0[[celltype]]$ID, "_", c("REST","info","Ref"))
      snp_split <- colsplit(under_split$info,"r",c("ch","SNP"))
      snp_split$SNP <- paste0(snp_split$SNP,":NA:NA")
      ref_2.0[[celltype]]$SNP <- snp_split$SNP
    }

    all_38_res <- list()
    for(celltype in names(ref_list)){
      message(celltype)
      all_38_res[[celltype]] <- data.frame()
      res_temp <- as.data.frame(results(RESTscreen_seq$ES_1, contrast=c("condition",celltype,"DNA")))
      message("expanding duplicates")
      res_temp2 <- expandDups(res_temp)
      res_temp <- res_temp2[which(rownames(res_temp2) %in% c(negCtrl_ID_all$ID,motif_ID_all$ID,noMotif_ID_all$ID,scr_ID_all$ID)),]
      res_temp <- unique(res_temp)
      oligo_split <- colsplit(rownames(res_temp),"\\^",names = c("Enhancer","Silencer"))
      for(enh in enh_all){
        message(enh)
        enh_res <- res_temp[which(oligo_split$Enhancer==enh),]

        ### subset the negative controls
        message("negative controls")
        neg_ctrl_res <- enh_res[which(rownames(enh_res) %in% neg_ctrls),2,drop=F]
        neg_ctrl_res$proj <- "NegCtrl"
        neg_ctrl_res$enh <- enh

        message(paste0("total negative controls: ", nrow(as.data.frame(neg_ctrl_res))))

        ### subset the noMotif set
        message("no motif ChIP positive")
        nm_chip_p_res <- enh_res[which(rownames(enh_res) %in% nm_chip_p[[celltype]]$V1),2,drop=F]
        if(nrow(nm_chip_p_res)>0){
          message(nrow(nm_chip_p_res))
          nm_chip_p_res$proj <- "noMotif_ChIP_pos"
          nm_chip_p_res$enh <- enh
        }
        message("no motif ChIP negative")
        nm_chip_n_res <- enh_res[which(rownames(enh_res) %in% nm_chip_n[[celltype]]),2,drop=F]
        if(nrow(nm_chip_n_res)>0){
          message(nrow(nm_chip_n_res))
          nm_chip_n_res$proj <- "noMotif_ChIP_neg"
          nm_chip_n_res$enh <- enh
        }

        message(paste0("total noMotif: ", nrow(as.data.frame(nm_chip_p_res))+nrow(as.data.frame(nm_chip_n_res))))

        ### subset the scrambled sets
        message("scrambled ChIP positive")
        scr_chip_p_res <- enh_res[which(rownames(enh_res) %in% scr_chip_p[[celltype]]),2,drop=F]
        if(nrow(scr_chip_p_res)>0){
          message(nrow(scr_chip_p_res))
          scr_chip_p_res$proj <- "scrambled_ChIP_pos"
          scr_chip_p_res$enh <- enh
        }
        message("scrambled ChIP negative")
        scr_chip_n_res <- enh_res[which(rownames(enh_res) %in% scr_chip_n[[celltype]]),2,drop=F]
        if(nrow(scr_chip_n_res)>0){
          message(nrow(scr_chip_n_res))
          scr_chip_n_res$proj <- "scrambled_ChIP_neg"
          scr_chip_n_res$enh <- enh
        }

        message(paste0("total scrambled: ", nrow(as.data.frame(scr_chip_p_res))+nrow(as.data.frame(scr_chip_n_res))))

        ### subset the motif sets
        message("motif ChIP positive")
        motif_chip_p_res <- enh_res[which(rownames(enh_res) %in% motif_chip_p_id[[celltype]]),2,drop=F]
        if(nrow(motif_chip_p_res)>0){
          message(nrow(motif_chip_p_res))
          motif_chip_p_res$proj <- "motif_ChIP_pos"
          motif_chip_p_res$enh <- enh
        }
        message("motif ChIP negative")
        motif_chip_n_res <- enh_res[which(rownames(enh_res) %in% motif_chip_n_id[[celltype]]),2,drop=F]
        if(nrow(motif_chip_n_res)>0){
          message(nrow(motif_chip_n_res))
          motif_chip_n_res$proj <- "motif_ChIP_neg"
          motif_chip_n_res$enh <- enh
        }

        message(paste0("total motif: ", nrow(as.data.frame(motif_chip_p_res))+nrow(as.data.frame(motif_chip_n_res))))

        all_38_res[[celltype]] <- rbind(all_38_res[[celltype]],neg_ctrl_res,nm_chip_p_res,nm_chip_n_res,motif_chip_p_res,scr_chip_p_res,motif_chip_n_res,scr_chip_n_res)
      }
    }
</details>      
  
<details open>
<summary>Set up FIMO scores for ref/scr/noMotif oligos</summary>
<br>
  
    ref_all_scores <- read.delim("meme-suite/Ref_all/fimo.txt", stringsAsFactors = F)
    scr_all_scores <- read.delim("meme-suite/Scr_noSNP/fimo.txt", stringsAsFactors = F)
    nomotif_all_scores <- read.delim("meme-suite/NoMotif/fimo.txt", stringsAsFactors = F)

    ref_all_scores$project <- "motif"
    scr_all_scores$project <- "scrambled"
    
    emvar_all_motif_scores <- data.frame()
    for(celltype in names(emVAR_all)){
      message(celltype)
      emvar_temp <- emVAR_all[[celltype]][is.na(emVAR_all[[celltype]]$ref_allele),]
      emvar_temp$celltype <- celltype
      emvar_all_motif_scores <- rbind(emvar_all_motif_scores, emvar_temp)
    }

    head(emvar_all_motif_scores)
    motif_all_IDs <- colsplit(emvar_all_motif_scores$ID, "\\^", c("enhancer","silencer"))
    emvar_all_motif_scores$enhancer <- motif_all_IDs$enhancer
    emvar_all_motif_scores$silencer <- motif_all_IDs$silencer

    emvar_all_motif_scores <- merge(emvar_all_motif_scores, ref_all_scores[which(ref_all_scores$start > 90 & ref_all_scores$stop < 112 & ref_all_scores$X.pattern.name=="MA0138.2"),], by.x="silencer", by.y="sequence.name", drop=F)
    emvar_all_motif_scores$skew_group <- NA
    emvar_all_motif_scores$skew_group[which(emvar_all_motif_scores$LogSkew < 0)] <- 1
    emvar_all_motif_scores$skew_group[which(emvar_all_motif_scores$LogSkew > 0)] <- 0

    emvar_scores_motif_df <- emvar_all_motif_scores[,c("skew_group","A_Ctrl_Mean","A_Exp_Mean","A_log2FC","A_log2FC_SE","A_logP","A_logPadj_BH","A_logPadj_BF",
                                                   "B_Ctrl_Mean","B_Exp_Mean","B_log2FC","B_log2FC_SE","B_logP","B_logPadj_BH","B_logPadj_BF","LogSkew",
                                                   "Skew_logP","Skew_logFDR","celltype","enhancer","silencer","score")]

    emvar_scores_motif_df <- emvar_scores_motif_df[complete.cases(emvar_scores_motif_df$skew_group),]
    emvar_scores_motif_df$skew_group <- as.numeric(emvar_scores_motif_df$skew_group)
  </details>
