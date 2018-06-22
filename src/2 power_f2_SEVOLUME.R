# May 2 2017
## Goal: Power check for multiple regression (NOT MV)
# selected x1 morphology outcome, as germaine demand in the revision

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

lifemice <- as.mids(lifemerg, .imp=1, .id = 2)


# ******** -----
# A. Y-full -----------------------------------------------------



# \\ models ----

## Male
mvlm_ml <- function(indvar, dat) {
    setform <- sprintf("log(SEVOLUME + 1e-10) ~ I(scale(log10(%s+1))) + lipids_m + Age_m + catbmi_m + LMSMKNOW + LMEXERCS + parity_m", indvar)
    mod <- with(data = dat, lm(as.formula(setform)))
}



mvlm_mc <- function(indvar, dat) {
    setform <- sprintf("log(SEVOLUME + 1e-10) ~ I(scale(log10(%s+1))) + CREAMOUNT_m + Age_m + catbmi_m + LMSMKNOW + LMEXERCS + parity_m", indvar)
    mod <- with(data = dat, lm(as.formula(setform)))
}

mvlm_mn <- function(indvar, dat) {
    setform <- sprintf("log(SEVOLUME + 1e-10) ~ I(scale(log10(%s+1))) + Age_m + catbmi_m + LMSMKNOW + LMEXERCS + parity_m", indvar)
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
    frm <- data.frame(r2est = pool.r.squared(input_m)[1], pvaledc = summary(pool(input_m))[2,5])
    frm$indvar <- tmpm_var
    print(tmpm_var)
    ret_df_full_r2 <- rbind(ret_df_full_r2, frm) 
}


# ******** -----
# A.1. For-simulate -----------------------------------------------------
# 060118 FDR simulation

ret_list_full <- list()

for(i in 1:length(unlist(cl_m))) {
    tmpm_var <- unlist(cl_m)[i]
    input_m <- ret_list[[tmpm_var]]
    out1 <-  summary(pool(input_m))[, 1]
    ret_list_full[[tmpm_var]] <- out1
}



ret_df_full_r2_adj <- data.frame()

for(i in 1:length(unlist(cl_m))) {
    tmpm_var <- unlist(cl_m)[i]
    input_m <- ret_list[[tmpm_var]]
    frm <- data.frame(r2est = pool.r.squared(input_m)[1], ad_r2est = pool.r.squared(input_m, adjusted = TRUE)[1], pvaledc = summary(pool(input_m))[2,5])
    
    # get y SD (sigma)
    ysigma <- sapply(input_m$analyses, function(x){
        tmp1 <- summary(x)
        tmp1$sigma
    })
    frm$ysigma_sd <- mean(ysigma)
    frm$indvar <- tmpm_var
    print(tmpm_var)
    ret_df_full_r2_adj <- rbind(ret_df_full_r2_adj, frm) 
}

# save(ret_list_full, file = "results/sim_volume_SEVOLUME.RData", compress = FALSE)
# save(ret_df_full_r2_adj, file = "results/sim_r2_volume_SEVOLUME.RData", compress = FALSE)


# ******** -----
# B. Y-redu -----------------------------------------------------



# \\ models ----

## Male
mvlm_ml <- function(indvar, dat) {
    setform <- sprintf("log(SEVOLUME + 1e-10) ~ lipids_m + Age_m + catbmi_m + LMSMKNOW + LMEXERCS + parity_m", indvar)
    mod <- with(data = dat, lm(as.formula(setform)))
}


mvlm_mc <- function(indvar, dat) {
    setform <- sprintf("log(SEVOLUME + 1e-10) ~ CREAMOUNT_m + Age_m + catbmi_m + LMSMKNOW + LMEXERCS + parity_m", indvar)
    mod <- with(data = dat, lm(as.formula(setform)))
}

mvlm_mn <- function(indvar, dat) {
    setform <- sprintf("log(SEVOLUME + 1e-10) ~ Age_m + catbmi_m + LMSMKNOW + LMEXERCS + parity_m", indvar)
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

ret_df_redu_r2 <- data.frame()

for(i in 1:length(unlist(cl_m))) {
    tmpm_var <- unlist(cl_m)[i]
    input_m <- ret_list[[tmpm_var]]
    frm <- as.data.frame(pool.r.squared(input_m)[1])
    # frm <- as.data.frame(pool.r.squared(input_m, adjusted = TRUE)[1])
    colnames(frm) <- "r2est"
    frm$indvar <- tmpm_var
    print(tmpm_var)
    ret_df_redu_r2 <- rbind(ret_df_redu_r2, frm) 
}



# ******** -----
## X. Saving Rdata ----

## merging

es <- cbind(ret_df_full_r2, ret_df_redu_r2)[, c(1, 4, 3, 2)]
colnames(es) <- c("r2full", "r2redu", "indvar", "pvaledc")
es$f2 <- (es$r2full-es$r2redu)/(1-es$r2full)

# remove -ve f2
es <- es[(es$f2 > 0),] 

# dir
outdirectory <- "results"
outfilename <- "mv_EWAS_power_semen_chems_SEVOLUME.Rdata"
# outfilename <- sprintf("%s_reg_7.Rdata", depVariable)

save(file=file.path(outdirectory,outfilename),
     lifemerg, impute_f, impute_m, life, es)

# load("results/mv_EWAS_power_semen_chems_SEVOLUME.Rdata") 

# ******** -----
# C. Effect-size -----------------------------------------------------

library(pwr)

# \\ n_80%

v_list <- list()
for(i in 1:nrow(es)) {
    input_1 <- es$f2[i]
    if (input_1 > 0){
        v <- pwr.f2.test(u = 1, f2 = input_1, sig.level = 0.05/128, power = 0.8)$v  # *Bonf p* value
        v_list[[i]] <- v 
    } else {
        v_list[[i]] <- 0 # optional, removed f2 <0
    }
}

es$n_80 <- ceiling(unlist(v_list)) + 7 + 1 


# \\ power 


p_list <- list()
for(i in 1:nrow(es)) { # 126 rows
    input_1 <- es$f2[i]
    if (input_1 > 0){
        p <- pwr.f2.test(u = 1, v = 473 - 7 - 1, f2 = input_1, sig.level = 0.05/128)$p # *Bonf p* value; # full model residual/error df
        p_list[[i]] <- p
    } else {
        p_list[[i]] <- 0 # optional, removed f2 <0
    }
}

es$power <- round(unlist(p_list),2)
es$power2 <- round(unlist(p_list),4)


# sample size for 95% f2
quantile(es$f2, 0.95)
pwr.f2.test(u = 1, f2 = quantile(es$f2, 0.95), sig.level = 0.05/128, power = 0.8) 



# write
# write.csv(es, file = "results/power_bon_adj_p_morphology.csv") 
# write.csv(es, file = "results/power_nil_adj_p_morphology.csv") 





# \\ cp-power-plot----
# figure in the paper

# power vs sample size
cp_p_list <- list()
cpseq <- seq(300, 700, by = 40) # step in sample size

for(i in 1:nrow(es)) { 
    input_1 <- es$f2[i]
    power_life <- sapply(cpseq, function(x){
        pwr.f2.test(u = 1, v = x - 7 - 1, f2 = input_1, sig.level = 0.05/128)$p
    })
    cp_p_list[[i]] <- data.frame(rep(es$indvar[i], length(cpseq)), rep(input_1, length(cpseq)), cpseq, power_life)
    
}


# create df
powerplot_df <- do.call("rbind", cp_p_list)
colnames(powerplot_df) <- c("EDC", "effect_size", "sample_size", "power")



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
    mtchC <- which(powerplot_df$EDC %in% cl_m[[i]])
    EDC_class[mtchC] <- names(cl_m)[i]
}

powerplot_df <- cbind(powerplot_df, EDC_class)

# create var for ordering x axis (thru factor lv)
powerplot_df$EDC <- factor(powerplot_df$EDC, levels = unlist(cl_m), ordered = TRUE) # manhat X in my order
powerplot_df$EDC_class <- factor(powerplot_df$EDC_class, levels = names(cl_m), ordered = TRUE) # manhat legend in my order

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
# names(colorCodes) <- levels(factor(get(tempname)$categ)) # by factor name to match the plot


# ggplot
library(ggplot2); library(ggrepel)

colorCodes_gg <- scale_colour_manual(name = "ECD_class",values = colorCodes)

# txt df
tmp1 <- subset(powerplot_df, sample_size == 700)
tmp1 <- tmp1[order(tmp1[,4]),]
tmp1 <- tail(tmp1, 5)
# dput(as.character(tmp1$EDC))
tmp1$EDC <- c("nickel", "2-OH-4-MeO-BP", "2,4-OH-BP", "Me-PFOSA-AcOH", "manganese") # manual

# plot
p <- ggplot(powerplot_df, aes(sample_size, power, group = EDC))
p <- p + geom_line(aes(colour = factor(EDC_class))) + xlim(300, 750)
p <- p + colorCodes_gg
p <- p + ggtitle("Power vs sample size for the endpoint seminal volume") +
    xlab("Sample size") + 
    ylab("Statistical power")
p <- p + geom_text_repel(data = tmp1, aes(label = EDC), 
                         size = 3,

                         box.padding = 0.2, # Add extra padding around each text label.
                         point.padding = 0.3, # Add extra padding around each data point.
                         segment.color = '#000000', # Color of the line segments.
                         segment.size = 0.5, # Width of the line segments.
                         arrow = arrow(length = unit(0.01, 'npc')), # Draw an arrow from the label to the data point.
                         force = 10, # Strength of the repulsion force.
                         nudge_x = 30, 
                         nudge_y = 0.025)
p

ggsave("results/power_n_SELVOLUME.png", scale=1, dpi=400)
# ggsave("results/power_n_SELVOLUME.svg", scale=1, dpi=400)


# library(plotly)
# ggplotly(p)

# quantile
library(dplyr); library(magrittr)
powerquan <- powerplot_df %>% filter(sample_size == 500)

quantile(powerquan$power)
# 0%          25%          50%          75%         100% 
# 0.0003948346 0.0008584900 0.0021847938 0.0116988981 0.4028237221 



