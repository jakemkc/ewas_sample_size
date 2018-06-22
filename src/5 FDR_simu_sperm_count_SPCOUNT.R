###' Jun 1 2018
###' Goal: FDR simulation

rm(list=ls()) # clear workspace;  # ls() # list objects in the workspace
# cmd-shit-f10: restart R. reset loaded settings & library etc
cat("\014")   # same as ctrl-L
options(max.print = 3000) # default 1000, and rstudio set at 5000
options(warn=1) # default = 0. To check loop warnings


# Load
load("results/mv_EWAS_lm_semen_chems_in.Rdata") # from "glm_sem_chems__adj_EWAS_non_WHO"

## \\ name-change 

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
no_m <- gsub("_m$", "", colnames(life), ignore.case = T)

# match and get the index order by female_no_f
index_m <- vector(mode = "numeric")
for(i in 1:length(no_m)) {
    mtchM <- which(papername_m %in% no_m[i])
    print(mtchM) # error check
    index_m <- append(index_m, mtchM)
}

#
colnames(life)[which(no_m %in% papername_m)] <- names(papername_m[index_m])

# prevent loop error because of the name
colnames(life) <- make.names(colnames(life))


## crete lipid

life$lipids_m <- (1.494*life$cholamt_m) + life$trigamt_m + life$phosamt_m


## clear object env
rm(index_m, mtchM, no_m, papername_m)




### .\ 1 var-vectors -----------------------------------------------------------
library(tidyverse)

yvar <- c("WHONORM", "SEVOLUME", "SPCOUNT", "TOTCNT", "PERMOT", "SCSADFI", "SCSAHDS")

xedc <- c("PCB_101_m", "PCB_105_m", "PCB_110_m", "PCB_118_m", "PCB_114_m", "PCB_128_m", "PCB_138_m", "PCB_146_m", "PCB_149_m", "PCB_151_m", "PCB_153_m", "PCB_156_m", "PCB_157_m", "PCB_167_m", "PCB_170_m", "PCB_172_m", "PCB_177_m", "PCB_178_m", "PCB_180_m", "PCB_183_m", "PCB_187_m", "PCB_189_m", "PCB_194_m", "PCB_195_m", "PCB_196_m", "PCB_201_m", "PCB_206_m", "PCB_209_m", "PCB_28_m", "PCB_44_m", "PCB_49_m", "PCB_52_m", "PCB_66_m", "PCB_74_m", "PCB_87_m", "PCB_99_m", "Et_PFOSA_AcOH_m", "Me_PFOSA_AcOH_m", "PFDeA_m", "PFNA_m", "PFOSA_m", "PFOS_m", "PFOA_m", "blood_Cd_m", "blood_Pb_m", "blood_Hg_m", "BB_153_m", "b_HCB_m", "g_HCB_m", "HCB_m", "mirex_m", "op_DDT_m", "oxychlordane_m", "pp_DDE_m", "pp_DDT_m", "tr_nonachlor_m", "BDE_17_m", "BDE_28_m", "BDE_47_m", "BDE_85_m", "BDE_99_m", "BDE_100_m", "BDE_153_m", "BDE_154_m", "BDE_183_m", "BDE_66_m", "daidzein_m", "O_DMA_m", "equol_m", "enterodiol_m", "enterolactone_m", "genistein_m", "CREAMOUNT_m", "fcholamt_m", "cholamt_m", "trigamt_m", "phosamt_m", "cotinine_m", "selenium_m", "arsenic_m", "manganese_m", "chromium_m", "beryllium_m", "cobalt_m", "molybdenum_m", "cadmium_m", "tin_m", "antimony_m", "tellurium_m", "caesium_m", "barium_m", "nickel_m", "copper_m", "zinc_m", "tungsten_m", "platinum_m", "thallium_m", "lead_m", "uranium_m", "mMP_m", "mEP_m", "mCPP_m", "mBP_m", "miBP_m", "mECPP_m", "mCMHP_m", "mEHHP_m", "mEOHP_m", "mCHP_m", "mBzP_m", "mEHP_m", "mOP_m", "mNP_m", "BPA_m", "2_OH_4MeO_BP_m", "4_OH_BP_m", "24_OH_BP_m", "22_OH_4MeO_BP_m", "2244_OH_BP_m", "MP_m", "EP_m", "PP_m", "BP_m", "BzP_m", "HP_m", "4_HB_m", "34_DHB_m", "OH_Me_P_m", "OH_Et_P_m", "TCS_m", "TCC_m", "paracetamol_m", "4_aminophenol_m")

xedc <- make.names(xedc)

xcov <- c("Age_m", "catbmi_m", "LMSMKNOW", "LMEXERCS", "parity_m", "CREAMOUNT_m", "lipids_m")

xcov <- make.names(xcov)





### ........ ---------------------------------------------------------------
### B.  model-para ------------------------------------------------------------

### .\ 1 extract-X-para -----------------------------------------------------------



## xedc parameters

tmp1_xedc_pm <- skimr::skim(life[, xedc]) 

xedc_pm <- as.data.frame(xedc)
xedc_pm$mean <- tmp1_xedc_pm %>% filter(stat == "mean") %>% pull(value) 
xedc_pm$sd <- tmp1_xedc_pm %>% filter(stat == "sd") %>% pull(value) 


all(xedc_pm$mean >= 0) # are all vales True? 

xedc_pm$mean[xedc_pm$mean <= 0] <- xedc_pm$sd[xedc_pm$mean <= 0] 






### .\ 1 extract-Y-para -----------------------------------------------------------

#' Y = motbility
#' from power_f2_motility.R

load("results/sim_sperm_count_SPCOUNT.RData")  # ret_list_full
load("results/sim_r2_sperm_count_SPCOUNT.RData")  # ret_df_full_r2_adj (with adjusted r, with y sigma (sd)

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




# ## Find the 95P ES (beta)
emp_p <- data.frame()[1:128, ] # empry df

for(i in 1:length(unlist(cl_m))) {
    tmpm_var <- unlist(cl_m)[i]
    ext1 <- ret_list_full[[tmpm_var]] %>% .[2] %>% as.data.frame()
    ext1$edc <- tmpm_var
    emp_p <- rbind(emp_p, ext1)
}

colnames(emp_p) <- c("est", "edc")
emp_p$abs <- abs(emp_p$est)
# emp_p$edc <- rownames(emp_p)
emp_p <- emp_p %>% arrange(abs)
quantile(emp_p$abs, probs = 0.95)
# 95%
# 0.7991797


findInterval(quantile(emp_p$abs, probs = 0.95), emp_p$abs) # (closest without going over)
# 121
es_beta <- emp_p[121,2] # mCMHP

## compare r2, adj r2 and sigma sd of Y
y_sd <- ret_df_full_r2_adj$ysigma_sd %>% mean #8.192147



### .\ 1 gen-XY-data -----------------------------------------------------------



## function to gen XY data
xydata <- function (NSUB = 100, x = "PCB_101_m") {
    
    ## X data
    # 1) xedc, e.g. pcb101 data, rlnorm
    m <- xedc_pm %>% filter(xedc == x) %>% .$mean 
    s <- xedc_pm %>% filter(xedc == x) %>% .$sd 
    
    location <- log(m^2 / sqrt(s^2 + m^2))
    shape <- sqrt(log(1 + (s^2 / m^2)))
    edc_data <- rlnorm(n = NSUB, location, shape)
    # mean(edc_data)
    # sd(edc_data)
    
    
    ## Y data

    
    # Y betas
    # pcb lipid: betas
    int_beta <- ret_list_full[[x]][1]
    
    model_edc_beta <- ret_list_full[[x]][2]
    
    

    y_mean = int_beta + model_edc_beta*scale(log10(edc_data+1))
    

    ydata <- rnorm(n = NSUB, mean = y_mean, sd = y_sd) 
    
    xydf <- cbind(ydata, edc_data) %>% as.data.frame()
    colnames(xydf) <- c("ydata", x)
    return(xydf)
}    



### ........ ---------------------------------------------------------------
### C.  code-sim ------------------------------------------------------------

# \\ models ----

## Male
mvlm_ml <- function(indvar, dat) {
    setform <- sprintf("ydata ~ I(scale(log10(%s+1)))", indvar)
    mod <- with(data = dat, lm(as.formula(setform)))
}


# \\ loop ----

# ret_list <- list()

NSUB = 2000
n.sims = 5
sign <- rep(NA, n.sims)

for(j in 1:n.sims) {
    
    if(j %% 10==0) {
        # Print every 10 j loop passed
        cat(paste0("iteration: ", j, "\n"))
    }
    
    ret_df <- data.frame()[1:length(unlist(cl_m)), ] #' empy df with 0 columns and 128 rows
    for(i in 1:length(unlist(cl_m))) {
        tmpm_var <- unlist(cl_m)[i]
        
        # print(tmpm_var) # before the warnings
        frm <- mvlm_ml(indvar = tmpm_var, dat = xydata(NSUB, x = tmpm_var))
        
        # ret_list[[tmpm_var]] <- frm 
        ext1 <- summary(frm) %>% coef %>% as.data.frame() %>% .[2, ]
        colnames(ext1) <- c("est", "se", "t", "p")
        rownames(ext1) <- tmpm_var
        ext1$rsq <- summary(frm)$r.squared
        ext1$ad_rsq <- summary(frm)$adj.r.squared
        ret_df <- rbind(ret_df, ext1) 
    }
    ret_df$var <- rownames(ret_df)
    ret_df$fdr <- p.adjust(ret_df$p, method='fdr')  # min(fdr)
    targetedc <- ret_df %>% filter(var == es_beta)
    sign[j] <- targetedc$fdr <= 0.05
    
}

power_ret <- mean(sign) #power, max = 1





## function

powerfun <- function(NSUB, n.sims) {
    sign <- rep(NA, n.sims)
    
    for(j in 1:n.sims) {
        
        if(j %% 10==0) {
            # Print every 10 j loop passed
            cat(paste0("iteration: ", j, "\n"))
        }
        
        ret_df <- data.frame()[1:length(unlist(cl_m)), ] #' empy df with 0 columns and 128 rows
        for(i in 1:length(unlist(cl_m))) {
            tmpm_var <- unlist(cl_m)[i]
            
            # print(tmpm_var) # before the warnings
            frm <- mvlm_ml(indvar = tmpm_var, dat = xydata(NSUB, x = tmpm_var))
            
            # ret_list[[tmpm_var]] <- frm 
            ext1 <- summary(frm) %>% coef %>% as.data.frame() %>% .[2, ]
            colnames(ext1) <- c("est", "se", "t", "p")
            rownames(ext1) <- tmpm_var
            ext1$rsq <- summary(frm)$r.squared
            ext1$ad_rsq <- summary(frm)$adj.r.squared
            ret_df <- rbind(ret_df, ext1) 
        }
        ret_df$var <- rownames(ret_df)
        ret_df$fdr <- p.adjust(ret_df$p, method='fdr')  # min(fdr)
        targetedc <- ret_df %>% filter(var == es_beta)
        sign[j] <- targetedc$fdr <= 0.05
        
    }
    
    return(mean(sign))
}


## linear reg to get nsize at power = 0.8

## parameters
NSUB <- c(1800, 2000, 2200)
n.sims = 100


ret_power <- data.frame(NSUB, rep(NA,3))
colnames(ret_power) <- c("nsize", "power")

for(i in 1:length(NSUB)) {
    print(i)
    tmpsize <- NSUB[i]
    tmp1 <- powerfun(tmpsize, n.sims) 
    ret_power[i, 2] <- tmp1
}

coeflm <- lm(data = ret_power, nsize ~ power) %>% coef
power = 0.8
coeflm[1] + coeflm[2]*power 








