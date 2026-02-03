##############################################################################
##
## Load packages
##
##############################################################################
library(doParallel)
library(BakerMetabolomics)
library(dplyr)

registerDoParallel(cores = 10)

##############################################################################
##
## Setup working environment
##
##############################################################################
setwd("Q:/Metabolomics/Habtamu B/PhD work/The Ausdiab Baseline lipidomics/AusDiab ratios/")

##############################################################################
##
## Load AusDiab data
##
##############################################################################
loadAusdiab()

## Derive HOMA-IR
Ausdiab_covar$HOMA_IR <- (Ausdiab_covar$insulinu_00 * Ausdiab_covar$fbg_00) / 22.5

## Select covariates
Ausdiab_covar <- Ausdiab_covar[, c(
  "id","bmi_00","drage_00","drsex_00_n",
  "chol_00","hdl_00","trig_00","plg_00",
  "insulinu_00","fbg_00","hba1c_00",
  "Homa_b_00","HOMA_IR"
)]

# write.csv(Ausdiab_covar, "Ausdiab_covar2.csv", row.names = FALSE)
# Ausdiab_covar <- read.csv("Ausdiab_covar2.csv", row.names = "ID")

##############################################################################
##
## Load lipid matrix (Cer subset)
##
##############################################################################
Ausdiab <- as.matrix(read.csv("Ausdiab2.csv",
                              check.names = FALSE,
                              row.names = "ID"))

Ausdiab_lipids_names <- colnames(Ausdiab)
colnames(Ausdiab) <- make.names(Ausdiab_lipids_names)

## Match covariates order
Ausdiab_covar <- Ausdiab_covar[row.names(Ausdiab), ]

##############################################################################
##
## Analysis setup
##
##############################################################################
covariates <- c("bmi_00","drage_00","drsex_00_n",
                "chol_00","hdl_00","trig_00")

covar <- scale(as.matrix(Ausdiab_covar[, covariates]))
predictor <- "plg_00"
pred <- scale(Ausdiab_covar[, predictor])

dir.create(predictor, showWarnings = FALSE)

##############################################################################
##
## Function: log10 p-values
##
##############################################################################
get_log_p <- function(t_val, df){
  log.p <- pt(abs(t_val), df, lower.tail = FALSE, log.p = TRUE) + log(2)
  log.p / log(10)
}

##############################################################################
##
## Univariate lipid associations
##
##############################################################################
df.used <- lm(scale(log(Ausdiab[,1, drop = FALSE])) ~ pred + covar)$df.residual

scaled_lipids <- scale(log(Ausdiab))

uniVar <- apply(scaled_lipids, 2, function(x){
  coef(summary(lm(x ~ pred + covar)))["pred", ]
})

if(is.null(dim(uniVar))){
  uniVar <- matrix(uniVar, nrow = 1)
  rownames(uniVar) <- colnames(scaled_lipids)
}

uniVar <- t(uniVar)

uniVar <- data.frame(
  lipid   = rownames(uniVar),
  Beta    = uniVar[,1],
  SE      = uniVar[,2],
  t       = uniVar[,3],
  p       = uniVar[,4],
  p.log10 = sapply(uniVar[,3], get_log_p, df = df.used),
  BH      = p.adjust(uniVar[,4], method = "BH"),
  stringsAsFactors = FALSE
)

write.csv(uniVar, paste0(predictor,"/univariate.csv"), row.names = FALSE)

##############################################################################
##
## Ratio analysis (FULL directional ratios – FIXED)
##
##############################################################################
lipid_names <- colnames(Ausdiab)

results <- foreach(i = seq_along(lipid_names),
                   .combine = rbind,
                   .errorhandling = "pass") %dopar% {
                     
                     lipid2 <- lipid_names[i]   # denominator
                     
                     ratios <- apply(Ausdiab[, -i, drop = FALSE], 2, function(x)
                       scale(log(x) - log(Ausdiab[, lipid2]))
                     )
                     
                     ratios <- scale(ratios)
                     
                     res <- t(apply(ratios, 2, function(x){
                       coef(summary(lm(x ~ pred + covar)))["pred", ]
                     }))
                     
                     res.p.log10 <- sapply(res[,3], get_log_p, df = df.used)
                     
                     df <- data.frame(
                       lipid1 = rownames(res),   # numerator
                       lipid2 = lipid2,
                       Beta   = res[,1],
                       SE     = res[,2],
                       t      = res[,3],
                       p      = res[,4],
                       p.log10 = res.p.log10,
                       stringsAsFactors = FALSE
                     )
                     
                     ## Add univariate stats
                     idx1 <- match(df$lipid1, uniVar$lipid)
                     idx2 <- match(df$lipid2, uniVar$lipid)
                     
                     df$lipid1_Beta <- uniVar$Beta[idx1]
                     df$lipid1_SE   <- uniVar$SE[idx1]
                     df$lipid1_t    <- uniVar$t[idx1]
                     df$lipid1_p.log10 <- uniVar$p.log10[idx1]
                     
                     df$lipid2_Beta <- uniVar$Beta[idx2]
                     df$lipid2_SE   <- uniVar$SE[idx2]
                     df$lipid2_t    <- uniVar$t[idx2]
                     df$lipid2_p.log10 <- uniVar$p.log10[idx2]
                     
                     df$univar_minP.log10 <- pmin(df$lipid1_p.log10, df$lipid2_p.log10)
                     df$pGain.log10 <- df$univar_minP.log10 - df$p.log10
                     
                     df
                   }

##############################################################################
##
## Global BH and significance
##
##############################################################################
p.adjust.BH.log <- function(log.p){
  n <- length(log.p)
  i <- n:1
  o <- order(log.p, decreasing = TRUE)
  ro <- order(o)
  pmin(0, cummin(log10(n/i) + log.p[o]))[ro]
}

results$Global_BH <- p.adjust(results$p, method = "BH")
results$Global_BH_log10 <- p.adjust.BH.log(results$p.log10)

crit <- nrow(results) / (2 * 0.05)
results$significant <- results$pGain.log10 > log10(crit)

##############################################################################
##
## Map lipid names back
##
##############################################################################
results$lipid1 <- Ausdiab_lipids_names[
  match(results$lipid1, make.names(Ausdiab_lipids_names))
]
results$lipid2 <- Ausdiab_lipids_names[
  match(results$lipid2, make.names(Ausdiab_lipids_names))
]

##############################################################################
##
## Write outputs
##
##############################################################################
write.csv(results, paste0(predictor,"/Complete.results.csv"), row.names = FALSE)
write.csv(results[results$significant == 1,],
          paste0(predictor,"/Significant.results.csv"), row.names = FALSE)







#####SOME PLOTS
# Install if necessary
# install.packages("readxl")

library(readxl)

# File path
file_path <- "\\\\fs1\\labs\\Metabolomics\\Habtamu B\\Postdoc\\Collaborations\\Barbara Khan\\Complete.results_2hPLG_CerRatios.xlsx"

# Read the sheet named "d182_d181"
df <- read_excel(path = file_path, sheet = "d182_d181")

# Inspect first rows
head(df)


library(ggplot2)
library(dplyr)
library(tidyr)
library(ggrepel)

# Prepare ratio name
df <- df %>% 
  mutate(ratio_name = paste(lipid1, "/", lipid2))

# Reshape for plotting
plot_df <- df %>%
  select(ratio_name, 
         Beta, p.log10, 
         lipid1, lipid1_Beta, lipid1_p.log10,
         lipid2, lipid2_Beta, lipid2_p.log10) %>%
  pivot_longer(
    cols = c("Beta", "lipid1_Beta", "lipid2_Beta"),
    names_to = "Type",
    values_to = "Beta_value"
  ) %>%
  mutate(
    log10p_value = case_when(
      Type == "Beta" ~ p.log10,
      Type == "lipid1_Beta" ~ lipid1_p.log10,
      Type == "lipid2_Beta" ~ lipid2_p.log10
    ),
    label = case_when(
      Type == "Beta"       ~ ratio_name,
      Type == "lipid1_Beta" ~ lipid1,
      Type == "lipid2_Beta" ~ lipid2
    )
  )

# Plot
PLOT<-ggplot(plot_df, aes(x = Beta_value, y = -log10p_value, color = Type)) +
  geom_point(size = 4) +
  geom_text_repel(aes(label = label), size = 3, nudge_y = 0.5) +
  geom_text(data = annot_df, aes(x = x, y = y, label = label),
            inherit.aes = FALSE, hjust = -0.1, vjust = 1.2, size = 4, color = "black") +
  facet_wrap(~ ratio_name, scales = "free") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  labs(
    x = "Coefficient (Beta)",
    y = "-log10(p-value)",
    color = "Type"
  ) +
  scale_color_manual(values = c("Beta" = "red", "lipid1_Beta" = "blue", "lipid2_Beta" = "green"),
                     labels = c("Ratio", "Lipid 1", "Lipid 2")) +
  theme_bw(base_size = 14) +
  theme(
    strip.text = element_blank(),     # <-- hides facet titles
    panel.grid.minor = element_blank()
  )



PLOT

####SAVE PLOT
library(ggplot2)
library(eoffice)
topptx(PLOT, file = "d182_181_ratio_COMPLETE.pptx", height = 7, width = 10)
