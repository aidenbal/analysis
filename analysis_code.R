# =====================================================================
# analysis_code.R
# White Matter Tractometry, Stress, and Cardiovascular Health (LEMON)
#
# HOW TO RUN
#   Open analysis.Rproj in RStudio, then run:  source("analysis_code.R")
#   Takes a few minutes. Progress is printed as it goes.
#
# SECTIONS
#   0. Settings and packages
#   1. Behavioural and physiological data (IDs, age, sex, BP -> MAP, PSQ, TICS)
#   2. PyAFQ tract profiles
#   3. Node-wise regressions: FA ~ outcome + Age + Sex, FDR-corrected
#   4. Whole-brain screening summary and ROI node clusters
#   5. Elastic Net prediction with NESTED 10-fold cross-validation
#  5b. Baseline comparison (age + sex only), standardized betas, OLS
#   6. Figures
#   7. Checks against previously saved results (only if .RData exists)
#
# OUTPUTS
#   results/  CSV tables      figures/  PNG figures
# =====================================================================


# ---- 0. Settings and packages ----------------------------------------
pkgs <- c("dplyr", "tidyr", "ggplot2", "glmnet")
to_install <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(to_install)) install.packages(to_install)
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(glmnet)
})

SEED      <- 2026   # makes cross-validation folds reproducible
ALPHA     <- 0.5    # Elastic Net mixing: 0 = ridge, 1 = lasso
N_FOLDS   <- 10     # outer and inner cross-validation folds
N_REPEATS <- 20     # outer CV is repeated with different splits for stability
METRIC    <- "dti_fa"
P_THRESH  <- 0.05   # node-wise threshold used to define ROI clusters
MIN_RUN   <- 3      # minimum consecutive nodes to count as a cluster

# Tracts used in the ROI and machine-learning analyses.
# Three systems implicated in stress regulation in prior literature:
# corpus callosum (CCMid), bilateral uncinate fasciculus, bilateral SLF I-III.
ROI_TRACTS <- c(CCMid    = "CCMid",
                UF_L     = "UF_L",
                UF_R     = "UF_R",
                LSLF_I   = "Left Superior Longitudinal I",
                LSLF_II  = "Left Superior Longitudinal II",
                LSLF_III = "Left Superior Longitudinal III",
                RSLF_I   = "Right Superior Longitudinal I",
                RSLF_II  = "Right Superior Longitudinal II",
                RSLF_III = "Right Superior Longitudinal III")

OUTCOMES <- c(MAP = "MAP", PSQ = "PSQ_OverallScore", TICS = "TICS_ScreeningScale")

FILES <- list(
  profiles_main = "pyafq_bundle_profiles.csv",           # PyAFQ bundles (AF, CCMid, UF, CST, ...)
  profiles_slf  = "pyafq_slf_callosum_profiles.csv",     # PyAFQ SLF I-III + callosum segments
  name_match    = "name_match.csv",
  participants  = "Participants_LEMON.csv",
  bp            = "Blood_Pressure_LEMON.csv",
  psq           = "PSQ.csv",
  tics          = "TICS.csv"
)

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

# One-time step: the two PyAFQ profile sets currently exist only inside .RData
# (node_data_clean = main bundles, afq_all = SLF/callosum segments).
# Export them to CSV so the pipeline runs from files.
exports <- list(c(FILES$profiles_main, "node_data_clean"),
                c(FILES$profiles_slf,  "afq_all"))
for (ex in exports) {
  if (file.exists(ex[1])) next
  if (!file.exists(".RData")) stop("Cannot find ", ex[1], " or .RData")
  tmp <- new.env(); load(".RData", envir = tmp)
  if (!exists(ex[2], envir = tmp)) stop(ex[2], " not found in .RData")
  write.csv(as.data.frame(get(ex[2], envir = tmp))[, c("tractID", "nodeID", "dti_fa", "dti_md", "subjectID")],
            ex[1], row.names = FALSE)
  rm(tmp)
  message("Exported ", ex[2], " to ", ex[1])
}

num <- function(x) suppressWarnings(as.numeric(x))
step <- function(txt) cat("\n==========", txt, "==========\n")


# ---- 1. Behavioural and physiological data ---------------------------
step("1. Loading behavioural and BP data")

read_id_csv <- function(path) {
  d <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  names(d)[1] <- "INDI_ID"   # first (unnamed) column holds the LEMON INDI ID
  d
}

ids <- read.csv(FILES$name_match, stringsAsFactors = FALSE)  # Initial_ID <-> INDI_ID

part <- read_id_csv(FILES$participants)
names(part)[grep("Gender", names(part))] <- "Sex"            # 1 = female, 2 = male
part$Sex     <- num(part$Sex)
part$Age_mid <- vapply(strsplit(as.character(part$Age), "-"),
                       function(z) mean(num(z)), numeric(1))  # "20-25" -> 22.5

bp_raw   <- read_id_csv(FILES$bp)
sys_cols <- grep("systole",  names(bp_raw), value = TRUE)
dia_cols <- grep("diastole", names(bp_raw), value = TRUE)
bp <- data.frame(
  INDI_ID        = bp_raw$INDI_ID,
  mean_systolic  = rowMeans(sapply(bp_raw[sys_cols], num), na.rm = TRUE),
  mean_diastolic = rowMeans(sapply(bp_raw[dia_cols], num), na.rm = TRUE)
)
bp[is.na(bp)] <- NA                                          # NaN -> NA
bp$MAP <- bp$mean_diastolic + (bp$mean_systolic - bp$mean_diastolic) / 3

psq  <- read_id_csv(FILES$psq)
tics <- read_id_csv(FILES$tics)
psq  <- data.frame(INDI_ID = psq$INDI_ID,  PSQ_OverallScore    = num(psq$PSQ_OverallScore))
tics <- data.frame(INDI_ID = tics$INDI_ID, TICS_ScreeningScale = num(tics$TICS_ScreeningScale))

subj <- ids %>%
  left_join(part[, c("INDI_ID", "Sex", "Age", "Age_mid")], by = "INDI_ID") %>%
  left_join(bp,   by = "INDI_ID") %>%
  left_join(psq,  by = "INDI_ID") %>%
  left_join(tics, by = "INDI_ID") %>%
  rename(subjectID = Initial_ID) %>%
  distinct(subjectID, .keep_all = TRUE)

cat("Participants with IDs:", nrow(subj),
    "| with MAP:", sum(!is.na(subj$MAP)),
    "| with PSQ:", sum(!is.na(subj$PSQ_OverallScore)),
    "| with TICS:", sum(!is.na(subj$TICS_ScreeningScale)), "\n")


# ---- 2. PyAFQ tract profiles ------------------------------------------
step("2. Loading PyAFQ tract profiles")

read_profiles <- function(path) {
  d <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  d <- d[, c("tractID", "nodeID", METRIC, "subjectID")]
  names(d)[3] <- "value"
  d
}

profiles <- bind_rows(read_profiles(FILES$profiles_main),
                      read_profiles(FILES$profiles_slf)) %>%
  group_by(subjectID, tractID, nodeID) %>%
  summarise(value = mean(value, na.rm = TRUE), .groups = "drop")

missing_roi <- setdiff(ROI_TRACTS, unique(profiles$tractID))
if (length(missing_roi)) {
  cat("Tracts found:\n"); print(sort(unique(profiles$tractID)))
  stop("ROI tract(s) not found in the profile files: ", paste(missing_roi, collapse = ", "))
}

node_df <- inner_join(profiles, subj, by = "subjectID")
cat("Tracts:", length(unique(profiles$tractID)),
    "| Subjects with profiles:", length(unique(profiles$subjectID)),
    "| Matched to behavioural data:", length(unique(node_df$subjectID)), "\n")
if (length(unique(node_df$subjectID)) == 0) stop("No subject IDs matched. Check name_match.csv.")


# ---- 3. Node-wise regressions ----------------------------------------
step("3. Node-wise regressions (this step takes a minute or two)")

run_nodewise <- function(dat, outcome) {
  dat %>%
    group_by(tractID, nodeID) %>%
    group_modify(function(d, key) {
      d <- d[complete.cases(d[, c("value", outcome, "Age_mid", "Sex")]), ]
      empty <- tibble(beta = NA_real_, p = NA_real_, n = nrow(d))
      if (nrow(d) < 10) return(empty)
      fit <- lm(reformulate(c(outcome, "Age_mid", "Sex"), response = "value"), data = d)
      ct  <- coef(summary(fit))
      if (!outcome %in% rownames(ct)) return(empty)
      tibble(beta = ct[outcome, 1], p = ct[outcome, 4], n = nrow(d))
    }) %>%
    ungroup() %>%
    mutate(fdr = p.adjust(p, method = "BH"))   # FDR across all nodes, per outcome
}

node_res <- list()
for (o in names(OUTCOMES)) {
  node_res[[o]] <- run_nodewise(node_df, OUTCOMES[[o]])
  write.csv(node_res[[o]], file.path("results", paste0("nodewise_", o, ".csv")), row.names = FALSE)
  cat(o, "done:", sum(node_res[[o]]$p < P_THRESH, na.rm = TRUE), "nodes with p <", P_THRESH,
      "|", sum(node_res[[o]]$fdr < 0.05, na.rm = TRUE), "nodes with FDR < .05\n")
}


# ---- 4. Screening summary and ROI clusters ---------------------------
step("4. Whole-brain screening and ROI clusters")

screening <- bind_rows(lapply(names(node_res), function(o) {
  node_res[[o]] %>%
    filter(!is.na(p)) %>%
    group_by(tractID) %>%
    summarise(min_p = min(p), node_at_min = nodeID[which.min(p)],
              min_fdr = min(fdr), .groups = "drop") %>%
    mutate(outcome = o)
}))
write.csv(screening, "results/wholebrain_screening.csv", row.names = FALSE)

screening_wide <- screening %>%
  select(tractID, outcome, min_p) %>%
  pivot_wider(names_from = outcome, values_from = min_p, names_prefix = "min_p_") %>%
  arrange(pmin(min_p_MAP, min_p_PSQ, min_p_TICS, na.rm = TRUE))
print(as.data.frame(screening_wide), digits = 3)

find_clusters <- function(nodes, p) {
  o <- order(nodes); nodes <- nodes[o]
  sig <- !is.na(p[o]) & p[o] < P_THRESH
  r <- rle(sig); ends <- cumsum(r$lengths); starts <- ends - r$lengths + 1
  keep <- r$values & r$lengths >= MIN_RUN
  if (!any(keep)) return(NULL)
  data.frame(start_node = nodes[starts[keep]], end_node = nodes[ends[keep]],
             n_nodes = r$lengths[keep])
}

roi_clusters <- bind_rows(lapply(names(node_res), function(o) {
  bind_rows(lapply(ROI_TRACTS, function(t) {
    d <- filter(node_res[[o]], tractID == t)
    cl <- find_clusters(d$nodeID, d$p)
    if (is.null(cl)) return(NULL)
    cbind(outcome = o, tractID = t, cl)
  }))
}))
write.csv(roi_clusters, "results/roi_clusters.csv", row.names = FALSE)
cat("\nROI clusters (runs of >=", MIN_RUN, "nodes with p <", P_THRESH, "):\n")
print(roi_clusters)


# ---- 5. Elastic Net with nested cross-validation ---------------------
step("5. Elastic Net prediction (nested CV, takes a few minutes)")
# Outer 10-fold CV estimates performance. Inside each training split,
# cv.glmnet picks lambda with its own 10-fold CV. Every prediction is made
# for a person the model never saw, so R-squared is truly out-of-sample.

roi_means <- profiles %>%
  filter(tractID %in% ROI_TRACTS) %>%
  group_by(subjectID, tractID) %>%
  summarise(mean_value = mean(value, na.rm = TRUE), .groups = "drop") %>%
  mutate(tract = names(ROI_TRACTS)[match(tractID, ROI_TRACTS)]) %>%
  select(subjectID, tract, mean_value) %>%
  pivot_wider(names_from = tract, values_from = mean_value)

ml_data    <- inner_join(roi_means, subj, by = "subjectID")
PREDICTORS <- c(names(ROI_TRACTS), "Age_mid", "Sex")

r2_oos <- function(y, p) 1 - sum((y - p)^2) / sum((y - mean(y))^2)
rmse   <- function(y, p) sqrt(mean((y - p)^2))

nested_cv <- function(X, y, seed) {
  set.seed(seed)
  folds <- sample(rep(seq_len(N_FOLDS), length.out = length(y)))
  pred  <- numeric(length(y))
  for (k in seq_len(N_FOLDS)) {
    tr    <- folds != k
    inner <- cv.glmnet(X[tr, , drop = FALSE], y[tr], alpha = ALPHA, nfolds = N_FOLDS)
    pred[!tr] <- as.numeric(predict(inner, newx = X[!tr, , drop = FALSE], s = "lambda.min"))
  }
  pred
}

ml <- list()
for (o in names(OUTCOMES)) {
  d <- ml_data[complete.cases(ml_data[, c(PREDICTORS, OUTCOMES[[o]])]), ]
  X <- as.matrix(d[, PREDICTORS]); y <- d[[OUTCOMES[[o]]]]

  preds <- sapply(seq_len(N_REPEATS), function(r) nested_cv(X, y, SEED + r))
  r2s   <- apply(preds, 2, function(p) r2_oos(y, p))
  rmses <- apply(preds, 2, function(p) rmse(y, p))

  set.seed(SEED)
  full    <- cv.glmnet(X, y, alpha = ALPHA, nfolds = N_FOLDS)   # final model, all data
  pred_in <- as.numeric(predict(full, newx = X, s = "lambda.min"))

  ml[[o]] <- list(
    n = length(y), y = y, pred_oof = rowMeans(preds), model = full,
    summary = data.frame(
      outcome = o, N = length(y),
      cv_R2_mean = mean(r2s), cv_R2_min = min(r2s), cv_R2_max = max(r2s),
      cv_RMSE_mean = mean(rmses),
      insample_R2 = r2_oos(y, pred_in),
      insample_cor2 = cor(y, pred_in)^2
    )
  )
  cat(o, "done (N =", length(y), ")\n")
}

ml_summary <- bind_rows(lapply(ml, `[[`, "summary"))
write.csv(ml_summary, "results/ml_performance.csv", row.names = FALSE)
cat("\nPrediction performance (cv_ = nested cross-validated; insample_ = NOT cross-validated):\n")
print(ml_summary, digits = 3)

coef_map <- as.matrix(coef(ml$MAP$model, s = "lambda.min"))
coef_df  <- data.frame(Predictor = rownames(coef_map), Beta = coef_map[, 1], row.names = NULL) %>%
  filter(Predictor != "(Intercept)")

# ---- 5b. Does white matter add anything beyond age and sex? ----------
step("5b. Baseline comparison (age + sex only vs. tracts + age + sex)")
cv_summary <- function(X, y) {
  preds <- sapply(seq_len(N_REPEATS), function(r) nested_cv(X, y, SEED + r))
  r2s   <- apply(preds, 2, function(p) r2_oos(y, p))
  data.frame(cv_R2_mean = mean(r2s), cv_R2_min = min(r2s), cv_R2_max = max(r2s),
             cv_RMSE_mean = mean(apply(preds, 2, function(p) rmse(y, p))))
}

# Both models are evaluated on the same people, with the same CV folds
baseline <- do.call(rbind, lapply(names(OUTCOMES), function(o) {
  dd <- ml_data[complete.cases(ml_data[, c(PREDICTORS, OUTCOMES[[o]])]), ]
  y  <- dd[[OUTCOMES[[o]]]]
  rbind(cbind(outcome = o, model = "Age + Sex",          N = nrow(dd),
              cv_summary(as.matrix(dd[, c("Age_mid", "Sex")]), y)),
        cbind(outcome = o, model = "Tracts + Age + Sex", N = nrow(dd),
              cv_summary(as.matrix(dd[, PREDICTORS]), y)))
}))
write.csv(baseline, "results/ml_baseline_comparison.csv", row.names = FALSE)
print(baseline, digits = 3, row.names = FALSE)

# Standardized betas for MAP: mmHg per 1 SD of each predictor
d_map  <- ml_data[complete.cases(ml_data[, c(PREDICTORS, "MAP")]), ]
X_full <- as.matrix(d_map[, PREDICTORS])
coef_df$Beta_std <- coef_df$Beta * apply(X_full, 2, sd)[coef_df$Predictor]
write.csv(coef_df, "results/ml_MAP_coefficients.csv", row.names = FALSE)
cat("\nMAP Elastic Net coefficients (lambda.min, full sample):\n"); print(coef_df, digits = 3)

# Ordinary regression for p-values, and a joint test of all tracts
ols  <- lm(reformulate(PREDICTORS, response = "MAP"), data = d_map)
base <- lm(MAP ~ Age_mid + Sex, data = d_map)
write.csv(coef(summary(ols)), "results/ols_MAP_coefficients.csv")
cat("\nOrdinary regression (MAP):\n"); print(coef(summary(ols)), digits = 3)
cat("\nJoint test: do the tracts add to age + sex?\n"); print(anova(base, ols))


# ---- 6. Figures --------------------------------------------------------
step("6. Saving figures")

s <- ml$MAP$summary
p1 <- ggplot(data.frame(Actual = ml$MAP$y, Predicted = ml$MAP$pred_oof),
             aes(Actual, Predicted)) +
  geom_point(color = "#0072B2", alpha = 0.7, size = 3) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "#D55E00", linewidth = 1.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
  labs(title = "Observed vs Predicted Mean Arterial Pressure",
       subtitle = paste0("Nested cross-validated R² = ", round(s$cv_R2_mean, 3),
                         " | RMSE = ", round(s$cv_RMSE_mean, 2), " mmHg | N = ", s$N),
       x = "Observed MAP (mmHg)", y = "Predicted MAP (mmHg, out-of-fold)") +
  theme_classic(base_size = 16) +
  theme(plot.title = element_text(face = "bold", size = 18),
        plot.subtitle = element_text(size = 13))
ggsave("figures/map_observed_vs_predicted.png", p1, width = 8, height = 7, dpi = 300)

p2 <- ggplot(coef_df, aes(reorder(Predictor, abs(Beta_std)), Beta_std, fill = Beta_std > 0)) +
  geom_col() + coord_flip() +
  scale_fill_manual(values = c("#D55E00", "#0072B2"), guide = "none") +
  labs(title = "Elastic Net Coefficients for MAP Prediction",
       subtitle = "Standardized: mmHg per 1 SD of each predictor", x = NULL,
       y = "Standardized coefficient (mmHg per SD)") +
  theme_classic(base_size = 16) + theme(plot.title = element_text(face = "bold"))
ggsave("figures/map_coefficients.png", p2, width = 8, height = 6, dpi = 300)

cols <- c(MAP = "#D55E00", PSQ = "#0072B2", TICS = "#009E73")
for (o in names(OUTCOMES)) {
  d <- filter(screening, outcome == o)
  p <- ggplot(d, aes(reorder(tractID, -min_p), -log10(min_p))) +
    geom_col(fill = cols[[o]]) + coord_flip() +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    labs(title = paste0("Whole-Brain Screening (", o, ")"),
         subtitle = "Lowest node-wise p-value within each tract",
         x = NULL, y = expression(-log[10](p))) +
    theme_classic(base_size = 16)
  ggsave(file.path("figures", paste0("screening_", o, ".png")), p, width = 8, height = 7, dpi = 300)
}
cat("Figures saved to figures/\n")


# ---- 7. Checks against previously saved results -----------------------
if (file.exists(".RData")) {
  step("7. CHECKS against your previously saved results")
  old <- new.env(); load(".RData", envir = old)
  verdict <- function(x) if (is.finite(x) && x < 1e-6) "MATCH" else "DIFFERENT"

  if (exists("bp_mapped", envir = old)) {
    chk <- inner_join(
      as.data.frame(old$bp_mapped) %>% select(subjectID = Initial_ID, old_sys = mean_systolic),
      subj %>% select(subjectID, mean_systolic), by = "subjectID")
    dmax <- max(abs(chk$old_sys - chk$mean_systolic), na.rm = TRUE)
    cat("Blood pressure (", nrow(chk), "people): max difference =", signif(dmax, 3), "->", verdict(dmax), "\n")
  }

  if (exists("ccmid", envir = old)) {
    chk <- inner_join(
      as.data.frame(old$ccmid) %>% select(subjectID, old_fa = mean_FA),
      roi_means %>% select(subjectID, CCMid), by = "subjectID")
    dmax <- max(abs(chk$old_fa - chk$CCMid), na.rm = TRUE)
    cat("CCMid tract means (", nrow(chk), "people): max difference =", signif(dmax, 3), "->", verdict(dmax), "\n")
  }

  for (pair in list(c("PSQ", "psq_node_results"), c("TICS", "tics_node_results"))) {
    if (!exists(pair[2], envir = old)) next
    saved <- as.data.frame(ungroup(get(pair[2], envir = old))) %>%
      select(tractID, nodeID, p_old = p, beta_old = beta)
    cmp <- inner_join(saved, node_res[[pair[1]]], by = c("tractID", "nodeID"))
    dmax <- max(abs(cmp$p - cmp$p_old), na.rm = TRUE)
    cat(pair[1], "node-wise p-values (", nrow(cmp), "nodes): max difference =", signif(dmax, 3),
        "| correlation =", round(cor(cmp$p, cmp$p_old, use = "complete.obs"), 4),
        "->", verdict(dmax), "\n")
    cat("   beta ratio new/old (median):", signif(median(cmp$beta / cmp$beta_old, na.rm = TRUE), 4),
        "(1 = same scale)\n")
  }

  # Diagnostic: which age coding reproduces the saved top PSQ node (AF_L node 18)?
  d <- node_df %>% filter(tractID == "AF_L", nodeID == 18) %>%
    filter(complete.cases(value, PSQ_OverallScore, Age_mid, Sex))
  if (nrow(d) > 10) {
    p_mid <- coef(summary(lm(value ~ PSQ_OverallScore + Age_mid + Sex, d)))[2, 4]
    p_fac <- coef(summary(lm(value ~ PSQ_OverallScore + factor(Age) + Sex, d)))[2, 4]
    cat("AF_L node 18 PSQ p-value: saved = 0.001702 | Age as number =", signif(p_mid, 4),
        "| Age as categories =", signif(p_fac, 4), "\n")
  }
  rm(old)
}

writeLines(capture.output(sessionInfo()), "results/session_info.txt")
step("DONE")
cat("Tables in results/, figures in figures/.\n")
