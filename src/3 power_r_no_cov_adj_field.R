# Aug 21 2017
## Goal: Power check for r (field-wise). Can't do r^2 comparison as they don't report it with MR model

rm(list=ls()) # clear workspace;  # ls() # list objects in the workspace
cat("\014")   # same as ctrl-L
options(max.print=3000) # default 1000, and rstudio set at 5000
options(warn=1) # default = 0. To check loop warnings



# Load -----------------------------------------------------
load("results/mv_EWAS_lm_semen_chems_in.Rdata") # from "glm_sem_chems__adj_EWAS_non_WHO"


## \\ name-change ----

## These 133 name are not R fault-proof (e.g.4_HB_f). Order of (papername) content is matched with (name(papername))
papername_m <- c("pcb101amt", "pcb105amt", "pcb110amt", "pcb114amt", "pcb118amt", "pcb128amt", "pcb138amt", "pcb146amt", "pcb149amt", "pcb151amt", "pcb153amt", "pcb156amt", 
                 "pcb157amt", "pcb167amt", "pcb170amt", "pcb172amt", "pcb177amt", "pcb178amt", "pcb180amt", "pcb183amt", "pcb187amt", "pcb189amt", "pcb194amt", "pcb195amt", 
                 "pcb196amt", "pcb201amt", "pcb206amt", "pcb209amt", "pcb028amt", "pcb044amt", "pcb049amt", "pcb052amt", "pcb066amt", "pcb074amt", "pcb087amt", "pcb099amt", 
                 "pfcepahamt", "pfcmpahamt", "pfcpfdeamt", "pfcpfnaamt", "pfcpfsaamt", "pfcpfosamt", "pfcpfoaamt", "metbcdamt", "metbpbamt", "metthgamt", "popbb1amt", "popbhcamt", 
                 "popghcamt", "pophcbamt", "popmiramt", "popodtamt", "popoxyamt", "poppdeamt", "poppdtamt", "poptnaamt", "pbdebr1amt", "pbdebr2amt", "pbdebr3amt", "pbdebr4amt", 
                 "pbdebr5amt", "pbdebr6amt", "pbdebr7amt", "pbdebr8amt", "pbdebr9amt", "pbdebr66amt", 
                 "DAZAMOUNT", "DMAAMOUNT", "EQUAMOUNT", "ETDAMOUNT", "ETLAMOUNT", "GNSAMOUNT", "CREAMOUNT", "fcholamt", "cholamt", "trigamt", "phosamt", "cotamt", "Selenium", 
                 "Arsenic", "Manganese", "Chromium", "Beryllium", "Cobalt", "Molybdenum", "Cadmium_Corrected", "Tin", "Antimony", "Tellurium", "Caesium", "Barium", "Nickel", 
                 "Copper", "Zinc", "Tungsten", "Platinum", "Thallium", "Lead", "Uranium", "mMethylPhthalate", "mEthylPhthalate", "mCarboxyPropylPhthalate", "mButylPhthalate", 
                 "mIsobutylPhthalate", "mCarboxyEthylPentylPhthalate", "mCarboxyMethylHexylPhthalate", "mEthylHydroxyHexylPhthalate", "mEthylOxoHexylPhthalate", "mCycloHexylPhthalate", 
                 "mBenzylPhthalate", "mEthylHexylPhthalate", "mOctylPhthalate", "mIsononylPhthalate", "BPA", "HydroxyMethoxyBenzoPhenone", "HydroxyBenzoPhenone", "DiHydroxyBenzoPhenone", 
                 "DiHydroxyMethoxyBenzoPhenone", "TetraHydroxyBenzoPhenone", "MeP", "EtP", "PrP", "BuP", "BzP", "HeP", "X_4_HB", "X_3_4_DHB", "OH_MeP", "OH_EtP", "TCS", "TCC", "PAP", "APAP")

names(papername_m) <- c("PCB_101_m", "PCB_105_m", "PCB_110_m", "PCB_118_m", "PCB_114_m", "PCB_128_m", "PCB_138_m", "PCB_146_m", "PCB_149_m", "PCB_151_m", "PCB_153_m", "PCB_156_m", "PCB_157_m", "PCB_167_m", "PCB_170_m", "PCB_172_m", "PCB_177_m", "PCB_178_m", "PCB_180_m", "PCB_183_m", "PCB_187_m", "PCB_189_m", "PCB_194_m", "PCB_195_m", "PCB_196_m", "PCB_201_m", "PCB_206_m", "PCB_209_m", "PCB_28_m", "PCB_44_m", "PCB_49_m", "PCB_52_m", "PCB_66_m", "PCB_74_m", "PCB_87_m", "PCB_99_m", "Et_PFOSA_AcOH_m", "Me_PFOSA_AcOH_m", "PFDeA_m", "PFNA_m", "PFOSA_m", "PFOS_m", "PFOA_m", "blood_Cd_m", "blood_Pb_m", "blood_Hg_m", "BB_153_m", "b_HCB_m", "g_HCB_m", "HCB_m", "mirex_m", "op_DDT_m", "oxychlordane_m", "pp_DDE_m", "pp_DDT_m", "tr_nonachlor_m", "BDE_17_m", "BDE_28_m", "BDE_47_m", "BDE_85_m", "BDE_99_m", "BDE_100_m", "BDE_153_m", "BDE_154_m", "BDE_183_m", "BDE_66_m", "daidzein_m", "O_DMA_m", "equol_m", "enterodiol_m", "enterolactone_m", "genistein_m", "CREAMOUNT_m", "fcholamt_m", "cholamt_m", "trigamt_m", "phosamt_m", "cotinine_m", "selenium_m", "arsenic_m", "manganese_m", "chromium_m", "beryllium_m", "cobalt_m", "molybdenum_m", "cadmium_m", "tin_m", "antimony_m", "tellurium_m", "caesium_m", "barium_m", "nickel_m", "copper_m", "zinc_m", "tungsten_m", "platinum_m", "thallium_m", "lead_m", "uranium_m", "mMP_m", "mEP_m", "mCPP_m", "mBP_m", "miBP_m", "mECPP_m", "mCMHP_m", "mEHHP_m", "mEOHP_m", "mCHP_m", "mBzP_m", "mEHP_m", "mOP_m", "mNP_m", "BPA_m", "2_OH_4MeO_BP_m", "4_OH_BP_m", "24_OH_BP_m", "22_OH_4MeO_BP_m", "2244_OH_BP_m", "MP_m", "EP_m", "PP_m", "BP_m", "BzP_m", "HP_m", "4_HB_m", "34_DHB_m", "OH_Me_P_m", "OH_Et_P_m", "TCS_m", "TCC_m", "paracetamol_m", "4_aminophenol_m")

# remove _m
no_m <- gsub("_m$", "", colnames(lifemerg), ignore.case = T)

# match and get the index order by female_no_f
index_m <- vector(mode = "numeric")
for(i in 1:length(no_m)) {
    mtchM <- which(papername_m %in% no_m[i])
    print(mtchM) # error check
    index_m <- append(index_m, mtchM)
}

#
colnames(lifemerg)[which(no_m %in% papername_m)] <- names(papername_m[index_m])

# prevent loop error because of the name
colnames(lifemerg) <- make.names(colnames(lifemerg))

## \\ mice ----
library(mice)
# read in master df for MICE
lifemice <- as.mids(data = lifemerg, .imp=1, .id = 2)


# ******** -----
# A. Y-full -----------------------------------------------------

# \\ models ----

## Male
mvlm_ml <- function(indvar, dat) {
    setform <- sprintf("log(PERMOT + 1e-10) ~ I(scale(log10(%s+1))) + lipids_m", indvar)
    mod <- with(data = dat, lm(as.formula(setform)))
}


mvlm_mc <- function(indvar, dat) {
    setform <- sprintf("log(PERMOT + 1e-10) ~ I(scale(log10(%s+1))) + CREAMOUNT_m", indvar)
    mod <- with(data = dat, lm(as.formula(setform)))
}

mvlm_mn <- function(indvar, dat) {
    setform <- sprintf("log(PERMOT + 1e-10) ~ I(scale(log10(%s+1)))", indvar)
    mod <- with(data = dat, lm(as.formula(setform)))
}



# \\ chem-list ----

# 3.1. Create list for matching in the loop to choose model


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
    Polybrominated_cpds = c(
        "BB_153_f", 
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
    Phenols = c(
        "BPA_f", 
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

# Operators 
# not 15 groups in germaine paper OK
list_m_oper <- list(
    M_creat = c("Phytoestrogens", "Phthalates", "Phenols", "Urine_metalloids", "Urine_metals", "Anti_microbial_cpds", "Paracetamols"),
    M_lipid = c("PCBs", "OCPs", "Polybrominated_cpds"),
    M_null = c("PFASs", "Cotinine", "Blood_metals"))


# \\ loop ----
ret_list <- list()

for(i in 1:length(unlist(cl_m))) {
    tmpm_var <- unlist(cl_m)[i]
    tmpm_oper <- names(cl_m[sapply(cl_m, "%in%", x = tmpm_var)])
    final_oper <- names(list_m_oper[sapply(list_m_oper, "%in%", x = tmpm_oper)])
    
    if (final_oper == "M_null"){
        frm <- mvlm_mn(tmpm_var, lifemice)
    } else if (final_oper == "M_creat"){
        frm <- mvlm_mc(tmpm_var, lifemice)
    } else {
        frm <- mvlm_ml(tmpm_var, lifemice)
    }
    print(tmpm_var) # after the warnings
    ret_list[[tmpm_var]] <- frm 
}

## \\ pool ----

# run
library(miceadds)

ret_df_full_r2 <- data.frame()

for(i in 1:length(unlist(cl_m))) {
    tmpm_var <- unlist(cl_m)[i]
    input_m <- ret_list[[tmpm_var]]
    frm <- as.data.frame(pool.r.squared(input_m)[1])
    colnames(frm) <- "r2est"
    frm$indvar <- tmpm_var
    print(tmpm_var)
    ret_df_full_r2 <- rbind(ret_df_full_r2, frm) 
}

# r^2 to pearson r
ret_df_full_r2$r <- sqrt(ret_df_full_r2$r2est)

## summary
hist(ret_df_full_r2$r2est)
hist(ret_df_full_r2$r)

hist(ret_df_full_r2$r, prob = TRUE)
hist(ret_df_full_r2$r, prob = TRUE, xlim = range(0,1))
lines(density(ret_df_full_r2$r, na.rm = TRUE))

quantile(ret_df_full_r2$r)


## \\ literature r ----

## all r, without spearman
lit_r <- c(0.14, 0.06, 0.05, -0.01, 0.04, -0.05, 0.05, 0.05, -0.0248, -0.095, -0.102, -0.128, -0.094, -0.201, -0.229, 0.165, 0.152, 0.043, 0.18, -0.041, -0.1, 0.057, 0.229, -0.092, 0.023, -0.071, 0.123, -0.045, 0.223, 0.057, -0.16, -0.43, -0.29, -0.39, -0.38, -0.25, 0.1, -0.28, -0.06, -0.13, -0.19, -0.15, 0.02, -0.01, 0.01, 0.01, -0.09, 0.11, -0.13, 0.11, 0.02, 0.09, -0.04, 0.05, 0.25, -0.25, -0.14, -0.682, -0.022, -0.403, -0.477, 0.124, -0.111, -0.791, -0.754, -0.076, 0.564, -0.198, -0.221, -0.046, -0.142, -0.769, -0.436, -0.125, -0.165, 0.04, 0.855)

lit_r <- abs(lit_r)


## summary
# sample size
mean(c(197, 116, 50, 191, 300, 21))
# 145.83

hist(lit_r, prob = TRUE)

hist(lit_r, prob = TRUE, xlim = range(0,1))
lines(density(lit_r, na.rm = TRUE))

quantile(lit_r)


lit_r_2 <- c(-0.19, -0.15, 0.02, -0.01, 0.01, 0.01, -0.09, 0.11, -0.13, 0.11, 0.02, 0.09, -0.04, 0.05, 0.25, -0.25, -0.14, 0.14, 0.06, 0.05, -0.01, 0.04, -0.05, 0.05, 0.05, -0.0248, -0.095, -0.102, -0.128, -0.094, -0.201, -0.229, 0.165, 0.152, 0.043, 0.18, -0.041, -0.1, 0.057, 0.229, -0.092, 0.023, -0.071, 0.123, -0.045, 0.223, 0.057)


lit_r_2 <- abs(lit_r_2)


## summary
# average sample size
mean(c(197, 116, 191, 300))
# 201

# 145.83
hist(lit_r_2, prob = TRUE)

hist(lit_r_2, prob = TRUE, xlim = range(0,1))
lines(density(lit_r_2, na.rm = TRUE))


quantile(lit_r_2)



# ******** -----
## C. Saving Rdata ----

outdirectory <- "results"
outfilename <- "mv_EWAS_power_r_no_cov_semen.Rdata"
save(file=file.path(outdirectory,outfilename), 
     ret_df_full_r2, lit_r, lit_r_2)

# Load data
# load("results/mv_EWAS_power_r_no_cov_semen.Rdata") 



# ******** -----
# C. Effect-size -----------------------------------------------------

# create LIFE_r df in quantile
esr <- ret_df_full_r2[order(ret_df_full_r2[,3]),] 
esr$r <- abs(esr$r) 
tmp1 <- quantile(esr$r)

quantile(esr$r, probs = 0.95)


lifeq1 <- findInterval(tmp1[2], esr$r) 
lifeq2 <- findInterval(tmp1[3], esr$r) 
lifeq3 <- findInterval(tmp1[4], esr$r)
life95 <- findInterval(quantile(esr$r, probs = 0.95), esr$r)

esr2 <- esr[c(lifeq1, lifeq2, lifeq3, life95), ]
esr2$quan <- c("LIFE_Q1", "LIFE_Q2", "LIFE_Q3", "LIFE_95")
esr2 <- esr2[, c(3,4)]

# create lit_r list in quantile
litr <- sort(lit_r_2) 

tmp2 <- quantile(litr)

quantile(litr, probs = 0.95)


litq1 <- findInterval(tmp2[2], litr) 
litq2 <- findInterval(tmp2[3], litr) 
litq3 <- findInterval(tmp2[4], litr) 
lit95 <- findInterval(quantile(litr, probs = 0.95), litr)

litr2 <- as.data.frame(litr[c(litq1, litq2, litq3, lit95)])
litr2$quan <- c("Literature_Q1", "Literature_Q2", "Literature_Q3", "Literature_95")
colnames(litr2) <- c("r", "quan")


# \\ cp-r-power-1 ----
# power vs sample size, comparing LIFE and literature, all 0.05 p value

library(pwr)

# LIFE_r_power
cp_p_life <- list()
cpseq <- seq(10, 560, by = 50) 

for(i in 1:nrow(esr2)) {
    input_1 <- esr2$r[i]
    power_life <- sapply(cpseq, function(x){
        pwr.r.test(r = input_1, n = x, sig.level = 0.05, alternative = "two.sided")$power # fixed mulitplicity = 100
    })
    cp_p_life[[i]] <- data.frame(rep(esr2$quan[i], length(cpseq)), rep(input_1, length(cpseq)), cpseq, power_life)
}

# LIterature_r_power
cp_p_lit <- list()
cpseq <- seq(10, 560, by = 50) 

for(i in 1:nrow(litr2)) {
    input_1 <- litr2$r[i]
    power_lit <- sapply(cpseq, function(x){
        pwr.r.test(r = input_1, n = x, sig.level = 0.05, alternative = "two.sided")$power # fixed mulitplicity = 100
    })
    cp_p_lit[[i]] <- data.frame(rep(litr2$quan[i], length(cpseq)), rep(input_1, length(cpseq)), cpseq, power_lit)
}


# create df
powerplot_life <- do.call("rbind", cp_p_life)
colnames(powerplot_life) <- c("quan", "effect_size_r", "sample_size", "power")
powerplot_life$grp <- "LIFE_study"

powerplot_lit <- do.call("rbind", cp_p_lit)
colnames(powerplot_lit) <- c("quan", "effect_size_r", "sample_size", "power")
powerplot_lit$grp <- "Literature"

powerplot <- rbind(powerplot_life, powerplot_lit)


# plot
library(ggplot2); library(ggrepel)

p <- ggplot(powerplot, aes(sample_size, power, group = quan))
p <- p + geom_line(aes(colour = grp)) + xlim(0, 620)
p <- p + annotate("text", x = 565, y = 0.43, label = "Literature_Q3", size = 3, hjust = 0)
p <- p + annotate("text", x = 565, y = 0.095, label = "Literature_Q2", size = 3, hjust = 0)
p <- p + annotate("text", x = 565, y = 0.075, label = "LIFE_Q3", size = 3, hjust = 0)
p <- p + annotate("text", x = 565, y = 0.04, label = "LIFE_Q2", size = 3, hjust = 0)
p <- p + annotate("text", x = 565, y = 0.015, label = "LIFE_Q1", size = 3, hjust = 0)
p <- p + annotate("text", x = 565, y = 0, label = "Literature_Q1", size = 3, hjust = 0)
p <- p + ggtitle("Power vs sample size for the effects of EDCs on semen quality") +
    xlab("Sample size") + ylab("Power")
p

# ggsave("results/power_n_lit_LIFE_005.png", scale=1, dpi=400)



# \\ cp-r-power-2 ----
# power vs sample size, literature comparison, multiplicity 12 vs 100

library(pwr)

litr2$quan <- c("25th P", "50th P", "75th P", "95th P")

# Literature_r_power_12
cp_p_lit_12 <- list()
cpseq <- seq(10, 560, by = 50) 

for(i in 1:nrow(litr2)) { 
    input_1 <- litr2$r[i]
    power_lit <- sapply(cpseq, function(x){
        pwr.r.test(r = input_1, n = x, sig.level = 0.05/12, alternative = "two.sided")$power # fixed mulitplicity 
    })
    cp_p_lit_12[[i]] <- data.frame(rep(litr2$quan[i], length(cpseq)), rep(input_1, length(cpseq)), cpseq, power_lit)
}


# LIterature_r_power_100

litr2$quan <- c("25th P ", "50th P ", "75th P ", "95th P ")

cp_p_lit_100 <- list()
cpseq <- seq(10, 560, by = 50) 

for(i in 1:nrow(litr2)) { 
    input_1 <- litr2$r[i]
    power_lit <- sapply(cpseq, function(x){
        pwr.r.test(r = input_1, n = x, sig.level = 0.05/100, alternative = "two.sided")$power # fixed mulitplicity
    })
    cp_p_lit_100[[i]] <- data.frame(rep(litr2$quan[i], length(cpseq)), rep(input_1, length(cpseq)), cpseq, power_lit)
}


# create df
powerplot_lit_12 <- do.call("rbind", cp_p_lit_12)
colnames(powerplot_lit_12) <- c("quan", "effect_size_r", "sample_size", "power")
powerplot_lit_12$grp <- "12 comparisons"

powerplot_lit_100 <- do.call("rbind", cp_p_lit_100)
colnames(powerplot_lit_100) <- c("quan", "effect_size_r", "sample_size", "power")
powerplot_lit_100$grp <- "100 comparisons"

powerplot <- rbind(powerplot_lit_12, powerplot_lit_100)


# plot
library(ggplot2); library(ggrepel)

# vline data
vldata <- data.frame(n = c(197, 116, 191, 300))

p <- ggplot(powerplot, aes(sample_size, power, group = quan))
p <- p + geom_line(aes(colour = grp)) + xlim(0, 620)
p <- p + geom_text_repel(data = subset(powerplot, sample_size == 560), aes(label = quan), nudge_x = 20)
p <- p + geom_vline(data = vldata, aes(xintercept = n), linetype = 4, colour = "blue", alpha = 0.3)
p <- p + ggtitle("Power vs sample size for the effects of EDCs on semen quality") +
    xlab("Sample size") + ylab("Power")
p

# ggsave("results/power_n_lit_12_100.png", scale=1, dpi=400)




subset(powerplot, sample_size == 210)



# \\ cp-r-power-3 ----
# power vs sample size, literature comparison, multiplicity 12 vs 100 vs no adjust

library(pwr)


# Literature_r_power_0_adj

litr2$quan <- c("25th P", "50th P", "75th P", "95th P")

cp_p_lit_0 <- list()
cpseq <- seq(10, 560, by = 50) 

for(i in 1:nrow(litr2)) { 
    input_1 <- litr2$r[i]
    power_lit <- sapply(cpseq, function(x){
        pwr.r.test(r = input_1, n = x, sig.level = 0.05, alternative = "two.sided")$power # fixed mulitplicity 
    })
    cp_p_lit_0[[i]] <- data.frame(rep(litr2$quan[i], length(cpseq)), rep(input_1, length(cpseq)), cpseq, power_lit)
}



# Literature_r_power_12

litr2$quan <- c("25th P ", "50th P ", "75th P ", "95th P ") 

cp_p_lit_12 <- list()
cpseq <- seq(10, 560, by = 50) 

for(i in 1:nrow(litr2)) { 
    input_1 <- litr2$r[i]
    power_lit <- sapply(cpseq, function(x){
        pwr.r.test(r = input_1, n = x, sig.level = 0.05/12, alternative = "two.sided")$power # fixed mulitplicity 
    })
    cp_p_lit_12[[i]] <- data.frame(rep(litr2$quan[i], length(cpseq)), rep(input_1, length(cpseq)), cpseq, power_lit)
}


# LIterature_r_power_100

litr2$quan <- c("25th P  ", "50th P  ", "75th P  ", "95th P  ") 

cp_p_lit_100 <- list()
cpseq <- seq(10, 560, by = 50) 

for(i in 1:nrow(litr2)) { 
    input_1 <- litr2$r[i]
    power_lit <- sapply(cpseq, function(x){
        pwr.r.test(r = input_1, n = x, sig.level = 0.05/100, alternative = "two.sided")$power # fixed mulitplicity
    })
    cp_p_lit_100[[i]] <- data.frame(rep(litr2$quan[i], length(cpseq)), rep(input_1, length(cpseq)), cpseq, power_lit)
}


# create df
powerplot_lit_0 <- do.call("rbind", cp_p_lit_0)
colnames(powerplot_lit_0) <- c("quan", "effect_size_r", "sample_size", "power")
powerplot_lit_0$group <- "No comparison"

powerplot_lit_0 <- subset(powerplot_lit_0, quan == "95th P") 

powerplot_lit_12 <- do.call("rbind", cp_p_lit_12)
colnames(powerplot_lit_12) <- c("quan", "effect_size_r", "sample_size", "power")
powerplot_lit_12$group <- "12 comparisons"

powerplot_lit_100 <- do.call("rbind", cp_p_lit_100)
colnames(powerplot_lit_100) <- c("quan", "effect_size_r", "sample_size", "power")
powerplot_lit_100$group <- "100 comparisons"

powerplot <- rbind(powerplot_lit_0, powerplot_lit_12, powerplot_lit_100)




# plot
library(ggplot2); library(ggrepel)

# vline data
vldata <- data.frame(n = c(197, 116, 191, 300))

p <- ggplot(powerplot, aes(sample_size, power, group = quan))
p <- p + geom_line(aes(colour = group)) + xlim(0, 620)
p <- p + geom_text_repel(data = subset(powerplot, sample_size == 560), aes(label = quan), nudge_x = 20)
p <- p + geom_vline(data = vldata, aes(xintercept = n), linetype = 4, colour = "blue", alpha = 0.3)
p <- p + ggtitle("Power vs sample size for the effects of EDCs on semen quality") +
    xlab("Sample size") + ylab("Statistical power")
p

# ggsave("results/power_n_lit_0_12_100_0_95_only.png", scale=1, dpi=400)
# ggsave("results/power_n_lit_0_12_100_0_95_only.svg", scale=1, dpi=400)



subset(powerplot, sample_size == 210)

# \\ n-for-95th-r ----


pwr.r.test(r=0.23,sig.level=0.05/100, power=0.8, alternative="greater")





