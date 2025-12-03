library(radEmu)
library(TreeSummarizedExperiment)
library(ANCOMBC)
library(ALDEx2)
library(DESeq2)
library(fastEmu)
library(tidyverse)
library(MicrobiomeStat)
library(maaslin3)

# get seed
args <- commandArgs(trailingOnly = FALSE)
if (length(args) == 0) {
  seed <- 1
} else {
  arg <- args[length(args)]
  seed <- abs(readr::parse_number(arg))
}

sim_dir_name <- "stamps/new_t1e_sims/"
source("radEmu/scripts/generate_test_data.R")

# generate b's (constant across simulation, use radEmu setup, make mean and median equal to 0)
n <- 50
J <- 250
dist <- "ZINB"
B_list <- radEmu:::get_sim_bs(250)
b0 <- B_list$b0
b1 <- B_list$b1
X <- cbind(intercept = 1, 
           cov1 = rep(0:1, each = n / 2))

# generate ZINB data with ~60% 0's 
Y <- generate_test_data_include_delta(n = n,
                                      J = J,
                                      B = rbind(b0, b1), 
                                      X = X,
                                      distn = "ZINB",
                                      zinb_size = 5,
                                      zinb_zero_prop = 0.6,
                                      mean_count_before_ZI = 50)$Y
colnames(Y) <- paste0("category", 1:J)
Y <- Y[, colSums(Y) > 0]
new_J <- ncol(Y)

# fit methods
res <- data.frame(
  seed = seed, 
  est = NA,
  p = NA,
  mse = NA,
  na_prop = NA, 
  method = c("ALDEx2", "ANCOM-BC2", "ANCOM-BC2 with SA", "CLR t-test", "DESeq2", "linda", "maaslin", "radEmu")
)

# ALDEx2 
clr_data <- aldex.clr(t(Y), X) 
aldex_model <- try({
  suppressMessages(
    aldex.glm(clr = clr_data, X)
  )
})
if (!inherits(aldex_model, "try-error")) {
  res$mse[1] <- mean((B_list$b1[parse_number(rownames(aldex_model))] - aldex_model$`cov1:Est` * log(2))^2, na.rm = TRUE)
  if ("category125" %in% rownames(aldex_model)) {
    ind <- which(rownames(aldex_model) == "category125")
    res$est[1] <- aldex_model$`cov1:Est`[ind] * log(2)
    res$p[1] <- aldex_model$`cov1:pval`[ind]
    res$na_prop[1] <- (new_J - nrow(aldex_model)) / new_J
  }
}

# ANCOM-BC2
rownames(Y) <- paste0("sample", 1:nrow(Y))
tse_data <- TreeSummarizedExperiment(assays = list("counts" = t(Y)),
                                     colData = data.frame(cov1 = X[, 2]))
ancom_model <- try({
  suppressMessages(
    ancombc2(data = tse_data,
             assay_name = "counts", 
             p_adj_method = "none",
             fix_formula = "cov1",
             prv_cut = 0,
             alpha = 0.05,
             pseudo_sens = TRUE)
  )
})
if (!inherits(ancom_model, "try-error")) {
  ind <- which(ancom_model$res$taxon == "category125") 
  res$mse[2:3] <- mean((B_list$b1[parse_number(ancom_model$res$taxon)] - 
                          ancom_model$res$lfc_cov1)^2, na.rm = TRUE)
  if (length(ind) > 0) {
    res$est[2:3] <- ancom_model$res$lfc_cov1[ind]
    res$p[2] <- ancom_model$res$p_cov1[ind]
    res$p[3] <- ifelse(ancom_model$res$passed_ss_cov1[ind], res$p[2], NA)
  }
  res$na_prop[2:3] <- mean(is.na(ancom_model$res$lfc_cov1))
}

# CLR t-test
Y_pseudo <- Y + 1
log_Y <- log(Y_pseudo)
mean_log_Y <- rowMeans(log_Y)
clr_data <- log_Y - matrix(mean_log_Y, nrow = n, ncol = new_J, byrow = FALSE)
coef_vec <- rep(NA, new_J)
names(coef_vec) <- colnames(Y)
ind <- which(names(coef_vec) == "category125")
for (j in 1:new_J) {
  mod <- lm(clr_Y ~ cov1, 
             data = data.frame(clr_Y = clr_data[, j], 
                               cov1 = X[, 2]))
  coef_vec[j] <- mod$coefficients[2]
  if (j == ind) {
    res$est[4] <- mod$coefficients[2]
    res$p[4] <- summary(mod)$coefficients[2, 4]
  }
}
res$mse[4] <- mean((B_list$b1[parse_number(names(coef_vec))] - coef_vec)^2, na.rm = T)
res$na_prop[4] <- 0

# DESeq2
dds <- DESeqDataSetFromMatrix(countData = t(Y),
                              colData = data.frame(X),
                              design = ~ cov1)
dds <- try({DESeq(dds, sfType = "poscounts")})
if (!inherits(dds, "try-error")) {
  deseq_res <- results(dds, name="cov1")
  ind <- which(rownames(deseq_res) == "category125")
  res$mse[5] <- mean((B_list$b1[parse_number(rownames(deseq_res))] - 
                        deseq_res$log2FoldChange * log(2))^2, na.rm = TRUE)
  res$est[5] <- deseq_res$log2FoldChange[ind] * log(2)
  res$p[5] <- deseq_res$pvalue[ind]
  res$na_prop[5] <- mean(is.na(deseq_res$log2FoldChange))
}

# linda
linda_res <- try({linda(feature.dat = t(Y), meta.dat = data.frame(X),
                   formula = "~ cov1", feature.dat.type = "count", 
                   is.winsor = FALSE)})
if (!inherits(linda_res, "try-error")) {
  ind <- which(rownames(linda_res$output$cov1) == "category125")
  res$mse[6] <- mean((B_list$b1[parse_number(rownames(linda_res$output$cov1))] - 
                        linda_res$output$cov1$log2FoldChange * log(2))^2, na.rm = T)
  res$est[6] <- linda_res$output$cov1$log2FoldChange[ind] * log(2)
  res$p[6] <- linda_res$output$cov1$pvalue[ind]
  res$na_prop[6] <- mean(is.na(linda_res$output$cov1$log2FoldChange))
}

# maaslin
file <- paste0("maaslin_temp", seed)
df <- data.frame(X)
rownames(df) <- paste0("sample", 1:nrow(df))
rownames(Y) <- rownames(df)
dir.create(file, recursive = TRUE, showWarnings = FALSE)
maaslin_run <- maaslin3(input_data = Y,
                        input_metadata = df,
                        formula = ~ cov1,
                        output = file, 
                        plot_summary_plot = FALSE,
                        plot_associations = FALSE,
                        save_models = FALSE)
maaslin_model <- read.delim(paste0(file, "/all_results.tsv")) %>%
  filter(model == "abundance", metadata == "cov1")
feat <- maaslin_model %>% filter(feature == "category125")
if (nrow(feat) == 1) {
  res$est[7] <- feat$coef * log(2)
  res$p[7] <- feat$pval_individual
  res$mse[7] <- mean((B_list$b1[parse_number(maaslin_model$feature)] - 
                        maaslin_model$coef * log(2))^2, na.rm = TRUE)
  res$na_prop[7] <- (new_J - nrow(maaslin_model) + sum(is.na(maaslin_model$coef))) / new_J
}

# radEmu
ind <- which(colnames(Y) == "category125")
emuRes <- try({
  emuFit(Y = Y, X = X, test_kj = data.frame(k = 2, j = ind), tolerance = 1e-3)
})
if (!inherits(emuRes, "try-error")) {
  ind <- which(colnames(emuRes$B) == "category125")
  res$mse[8] <- mean((B_list$b1[parse_number(colnames(emuRes$B))] - emuRes$B[2, ])^2, na.rm = T)
  res$est[8] <- emuRes$coef$estimate[125]
  res$p[8] <- emuRes$coef$pval[125]
  res$na_prop[8] <- mean(is.na(emuRes$coef$estimate))
}

saveRDS(res, paste0(sim_dir_name, "new_t1e_res", seed, ".rds"))
