## July 17 2017
## Goal: plot the results from "mv_EWAS_lm_semen_chems.R"


rm(list=ls()) # clear workspace;  # ls() # list objects in the workspace
cat("\014")   # same as ctrl-L
options(max.print=3000) # default 1000, and rstudio set at 5000
options(warn=1) # default = 0. To check loop warnings

# Load data
load("results/mv_EWAS_lm_semen_chems.Rdata")


# ******** -----
# A. Prep plots  -------------------------------------------


# ******** -----
# B. Manahatten plot ------------------------------------------------------
library(ggplot2); library(ggrepel)


# \\ prep ----
title <- sprintf("%s against chemicals", "mv")


p1_e3 <- 1e-3       # 0.001
p5_e3 <- 5e-3       # 0.005
pbonm <- 0.05/128   # 0.000375



pbonmlab <- subset(ret_df_mvlm, p < pbonm) #
p1_e3lab <- subset(ret_df_mvlm, p < 1e-3)
p5_e3lab <- subset(ret_df_mvlm, p < 5e-3)


# FDR

ret_df_mvlm$fdr <- p.adjust(ret_df_mvlm$p, method='fdr')  # min(fdr)
fdrsign <- subset(ret_df_mvlm, fdr < 0.1) 


# create my color scale for grp
    

cl_f = list(  
    #PCBs
    PCBs = c(
        "PCB_28_f", 
        "PCB_44_f", 
        "PCB_49_f", 
        "PCB_52_f", 
        "PCB_66_f", 
        "PCB_74_f", 
        "PCB_87_f", 
        "PCB_99_f", 
        "PCB_101_f", 
        "PCB_105_f", 
        "PCB_110_f", 
        "PCB_118_f", 
        "PCB_114_f", 
        "PCB_128_f", 
        "PCB_138_f", 
        "PCB_146_f", 
        "PCB_149_f", 
        "PCB_151_f", 
        "PCB_153_f", 
        "PCB_156_f", 
        "PCB_157_f", 
        "PCB_167_f", 
        "PCB_170_f", 
        "PCB_172_f", 
        "PCB_177_f", 
        "PCB_178_f", 
        "PCB_180_f", 
        "PCB_183_f", 
        "PCB_187_f", 
        "PCB_189_f", 
        "PCB_194_f", 
        "PCB_195_f", 
        "PCB_196_f", 
        "PCB_201_f", 
        "PCB_206_f", 
        "PCB_209_f" 
    ),
    #OCPs
    OCPs = c(
        "HCB_f",
        "b_HCB_f",
        "g_HCB_f",
        "op_DDT_f",
        "pp_DDE_f",
        "pp_DDT_f",
        "oxychlordane_f",
        "tr_nonachlor_f",
        "mirex_f"),
    #PBC
    PBBs = c(
        "BB_153_f"),
    PBDEs = c(
        "BDE_17_f", 
        "BDE_28_f", 
        "BDE_47_f", 
        "BDE_66_f", 
        "BDE_85_f", 
        "BDE_99_f", 
        "BDE_100_f", 
        "BDE_153_f", 
        "BDE_154_f", 
        "BDE_183_f"),
    #PFASs
    PFASs = c(
        "Et_PFOSA_AcOH_f", 
        "Me_PFOSA_AcOH_f", 
        "PFDeA_f", 
        "PFNA_f", 
        "PFOSA_f", 
        "PFOS_f", 
        "PFOA_f"),
    #blood metals
    Blood_metals = c(
        "blood_Cd_f", 
        "blood_Pb_f", 
        "blood_Hg_f"),
    #cotinine
    Cotinine = c(
        "cotinine_f"),
    #phytoestrogens
    Phytoestrogens = c(
        "genistein_f", 
        "daidzein_f", 
        "O_DMA_f", 
        "equol_f", 
        "enterodiol_f", 
        "enterolactone_f"),
    #phthalates
    Phthalates = c(
        "mMP_f", 
        "mEP_f", 
        "mCPP_f", 
        "mBP_f", 
        "miBP_f", 
        "mECPP_f", 
        "mCMHP_f", 
        "mEHHP_f", 
        "mEOHP_f", 
        "mCHP_f", 
        "mBzP_f", 
        "mEHP_f", 
        "mOP_f", 
        "mNP_f"),
    #phenols
    Bisphenol_A = c(
        "BPA_f"),
    Benzophenones = c(
        "X2_OH_4MeO_BP_m", 
        "X4_OH_BP_f", 
        "X24_OH_BP_f", 
        "X22_OH_4MeO_BP_f", 
        "X2244_OH_BP_f"),
    #anti microbial
    Anti_microbial_cpds = c(
        "MP_f", 
        "EP_f", 
        "PP_f", 
        "BP_f", 
        "BzP_f", 
        "HP_f", 
        "X4_HB_f", 
        "X34_DHB_f", 
        "OH_Me_P_f", 
        "OH_Et_P_f", 
        "TCS_f", 
        "TCC_f"),
    #paracetamol
    Paracetamols = c(
        "paracetamol_f", 
        "X4_aminophenol_f"),
    #urine metals
    Urine_metals = c(
        "manganese_f", 
        "chromium_f", 
        "beryllium_f", 
        "cobalt_f", 
        "molybdenum_f", 
        "cadmium_f", 
        "tin_f", 
        "caesium_f", 
        "barium_f", 
        "nickel_f", 
        "copper_f", 
        "zinc_f", 
        "tungsten_f", 
        "platinum_f", 
        "thallium_f", 
        "lead_f", 
        "uranium_f"),
    #urine metalloids
    Urine_metalloids = c(
        "selenium_f", 
        "arsenic_f", 
        "antimony_f", 
        "tellurium_f")
)

cl_m <- cl_f

# replace the _f in the cl_m to _m
for(i in 1:length(cl_m)){
    cl_m[[i]] <- gsub("_f$", "_m", cl_m[[i]], ignore.case = T)
    # print(cl_m[[i]])
}

# create category
EDC_class <- "n____"

for(i in 1:length(cl_m)) {
    mtchC <- which(ret_df_mvlm$indvar %in% cl_m[[i]])
    EDC_class[mtchC] <- names(cl_m)[i]
}

ret_df_mvlm <- cbind(ret_df_mvlm, EDC_class)

ret_df_mvlm$indvar3 <- factor(ret_df_mvlm$indvar, levels = unlist(cl_m), ordered = TRUE) # manhat X in my order
ret_df_mvlm$EDC_class <- factor(ret_df_mvlm$EDC_class, levels = names(cl_m), ordered = TRUE) # manhat legend in my order

# color code
colorCodes <- c("#ea7b00",
                "#0195fb",
                "#3aae24",
                "#c821a7",
                "#01df98",
                "#da0085",
                "#e5c440",
                "#f18cff",
                "#535800",
                "#972064",
                "#00b2b6",
                "#964400",
                "#5e4882",
                "#ff9288",
                "#b67696") # 15 color

names(colorCodes) <- names(cl_m)

colorCodes_gg <- scale_colour_manual(name = "EDC_class",values = colorCodes)



# \\ plot-w-grp ----

p <- ggplot(ret_df_mvlm, aes(indvar3, -log10(p))) 
p <- p + geom_point(aes(colour = EDC_class))
p <- p + ylab(expression(paste(-log[10], (italic(p))))) #?plotmath
p <- p + colorCodes_gg

# Threshold lines
p <- p + geom_hline(yintercept=-log10(pbonm), col="red") 
p <- p + geom_hline(yintercept=-log10(p1_e3), col="orange")
p <- p + geom_hline(yintercept=-log10(p5_e3), col="green")


# label lines
p <- p + annotate("text", x = 12, y = 3.5, label = "italic(p) == Bonferroni~correction", size = 3, parse = TRUE) 
p <- p + annotate("text", x = 8, y = 3.1, label = "italic(p) == 0.001", size = 3, parse = TRUE)
p <- p + annotate("text", x = 8, y = 2.4, label = "italic(p) == 0.005", size = 3, parse = TRUE)

# label axis
p <- p + xlab("EDCs") + ggtitle(sprintf("%s. Plot -log10(pvalue) vs. chem", title))

p <- p + geom_text(data=pbonmlab, aes(indvar, -log10(p), label=as.character(indvar)), size=3, vjust=0) #pbonmlab has no indvar3

p <- p + geom_text(data=p1_e3lab, aes(indvar, -log10(p), label=as.character(indvar)), size=3, vjust=0)
p <- p + geom_text(data=p5_e3lab, aes(indvar, -log10(p), label=as.character(indvar)), size=3, vjust= -0.8)


# label dot: FDR
p <- p + geom_text(data=fdrsign, aes(indvar, -log10(p), label=as.character(indvar)), hjust=0, vjust = -0.8, color = "red")


# annotate those p < 0.05 with ggrel
tmp1 <- subset(ret_df_mvlm, p <= 0.05)
tmp1$labelfig <- c("PCB 99", "PCB 105", "PCB 114", "PCB 167", "PBDE 17", "Me-PFOSA-AcOH", "PFOSA")

p <- p + geom_text_repel(data = tmp1, aes(label = labelfig), nudge_x = 1)


# adjust y axis padding
p <- p + scale_y_continuous(expand = c(0.02,0)) # adj. padding thickness

# adjust x axis text
p <- p + theme(plot.title = element_text(hjust = 0.5), 
               panel.grid.major.x = element_blank(),
               axis.title.x = element_blank(),
               axis.text.x = element_blank(),
               axis.ticks.x = element_blank())
# p <- p + theme_bw() # use a b/w color theme
p


## for ggplot only: you can specify or it will take default including to save the last ggplot
savingpath <- sprintf("results/man_pval_chem__%s_EWAS.png", "mv")
ggsave(savingpath, scale=1, dpi=400)
# ggsave(savingpath, scale=1, dpi=200)

# ggsave("results/man_pval_chem__mv_EWAS.svg", scale=1, dpi=400)

