# ===============================================================
#  Tethers ~ Contacts  (DTW-only)  v3
#  - Goal: do contacts (count/volume) co-vary over time with tether expression
#  - Uses DTW to compare shapes across time (z-scored within pair)
#  - Reads exact columns from contact CSVs: dataset, image_name, object, count_volume, sum_volume, ...
#  - Replaces NA -> 0 in count_volume and sum_volume (pre-filter)
#  - Filters pairs by number of 'X' in contact name (default 1 => pairs)
#  - Optional pairing file (object,gene) to restrict comparisons
#  - Exports: summaries, merged, complete pairs, z-scored series, DTW with ranks (+ optional Spearman)
#  - Saves DTW ribbon plots for each (object,gene) for count & volume
#  - Puts key data.frames into Environment; writes DEBUG parsing issues per timepoint
# ===============================================================

suppressPackageStartupMessages({
  library(tidyverse); library(readr); library(readxl); library(openxlsx)
  library(janitor);  library(dtw);   library(ggrepel)
})

# ----------------- USER PARAMS -----------------
tether_xlsx   <- "Tethers expression_05192025.xlsx"
contact_csvs  <- c(
  iPSCs  = "iPSCs_08102025_SumStatper_contact_summarystats_m.csv",
  iN_D7  = "iNday7_soma_08102025_SumStatper_contact_summarystats_m.csv",
  iN_D14 = "12302025_iNday14SegmmMIX_soma_SumStatper_contact_summarystats_m.csv",
  iN_D21 = "12282025_iNday21D14_soma_soma_SumStatper_contact_summarystats_m.csv"
)
pairing_file   <- "pairing_tether.xlsx"   # set NULL to compare all-against-all
out_dir        <- "pairing_analysis_outputs_NOTweighted"

num_X_required <- 1    # 1=pairs; 2=3-way; 0=no filter
do_spearman    <- FALSE  # set TRUE to also compute Spearman rho/p/FDR (optional)
do_perm_pval   <- TRUE   # DTW empirical p-value via permutation of expression series
perm_B         <- 120    # number of permutations (max 120 for 5 timepoints)

tp_levels <- c("iPSCs","iN_D7","iN_D14","iN_D21")

# ----------------- HELPERS -----------------
sem <- function(x) sd(x, na.rm=TRUE) / sqrt(sum(!is.na(x)))

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
debug_dir <- file.path(out_dir, "DEBUG")
dir.create(debug_dir, showWarnings = FALSE, recursive = TRUE)

# ----------------- 1) LOAD CONTACTS (strict columns, NA->0) -----------------
load_contacts <- function(path, tp) {
  raw <- readr::read_csv(
    file = path,
    col_types = cols(
      dataset      = col_character(),
      image_name   = col_character(),
      object       = col_character(),
      count_volume = col_double(),
      sum_volume   = col_double(),
      .default     = col_guess()
    ),
    guess_max = 100000,
    show_col_types = FALSE
  )
  names(raw) <- make.unique(names(raw), sep = "__dup")

  required <- c("dataset","image_name","object","count_volume","sum_volume")
  missing  <- setdiff(required, names(raw))
  if (length(missing) > 0) {
    stop(sprintf("[load_contacts] In %s missing columns: %s\nAvailable: %s",
                 basename(path), paste(missing, collapse=", "),
                 paste(names(raw), collapse=", ")))
  }

  assign(paste0("contacts_raw_", tp), raw, envir = .GlobalEnv)
  prb <- readr::problems(raw)
  assign(paste0("problems_", tp), prb, envir = .GlobalEnv)
  readr::write_csv(prb, file.path(debug_dir, paste0("problems_", tp, ".csv")))

  df <- raw %>%
    transmute(
      image_name,
      object,
      count_volume = tidyr::replace_na(count_volume, 0),
      sum_volume   = tidyr::replace_na(sum_volume, 0),
      TimePoint    = tp
    )
  return(df)
}

contacts_all <- purrr::imap_dfr(contact_csvs, ~load_contacts(.x, .y))

contacts_clean <- contacts_all %>%
  mutate(n_x = stringr::str_count(object, "X")) %>%
  filter(if (num_X_required>0) n_x==num_X_required else TRUE) %>%
  mutate(TimePoint = factor(TimePoint, levels = tp_levels))

contacts_summary <- contacts_clean %>%
  group_by(object, TimePoint) %>%
  summarise(
    count_mean  = mean(count_volume, na.rm=TRUE),
    count_sd    = sd(count_volume,   na.rm=TRUE),
    count_sem   = sem(count_volume),
    volume_mean = mean(sum_volume,   na.rm=TRUE),
    volume_sd   = sd(sum_volume,     na.rm=TRUE),
    volume_sem  = sem(sum_volume),
    n_cells     = sum(!is.na(count_volume)),
    .groups="drop"
  ) %>%
  mutate(
    count_sd   = tidyr::replace_na(count_sd, 0),
    count_sem  = tidyr::replace_na(count_sem, 0),
    volume_sd  = tidyr::replace_na(volume_sd, 0),
    volume_sem = tidyr::replace_na(volume_sem, 0)
  )

openxlsx::write.xlsx(contacts_summary, file.path(out_dir, "01_contacts_summary.xlsx"))

assign("contacts_all", contacts_all, envir = .GlobalEnv)
assign("contacts_clean", contacts_clean, envir = .GlobalEnv)
assign("contacts_summary", contacts_summary, envir = .GlobalEnv)

# ----------------- 2) LOAD TETHERS -----------------
suppressPackageStartupMessages({ library(readxl); library(dplyr); library(tidyr); library(stringr); library(openxlsx) })

# upload Excel
tether_raw <- readxl::read_xlsx(tether_xlsx)

# choose the colum of the gene of interest
gene_col <- dplyr::case_when(
  "PG.Genes" %in% names(tether_raw) ~ "PG.Genes",
  "Gene"     %in% names(tether_raw) ~ "Gene",
  TRUE ~ names(tether_raw)[1]
)

# Alias on the column names (es. iPSC_1 -> iPSCs_1)
nm <- names(tether_raw)

# iPSCs variant -> iPSCs_
nm <- sub("^iPSC[s]?[_\\s\\.-]?", "iPSCs_", nm, ignore.case = TRUE, perl = TRUE)

# iNday7/14/21/28  -> iN_Dx_
nm <- sub("^iN[_\\s\\.-]*day[_\\s\\.-]*7[_\\s\\.-]?",  "iN_D7_",  nm, ignore.case = TRUE, perl = TRUE)
nm <- sub("^iN[_\\s\\.-]*day[_\\s\\.-]*14[_\\s\\.-]?", "iN_D14_", nm, ignore.case = TRUE, perl = TRUE)
nm <- sub("^iN[_\\s\\.-]*day[_\\s\\.-]*21[_\\s\\.-]?", "iN_D21_", nm, ignore.case = TRUE, perl = TRUE)
nm <- sub("^iN[_\\s\\.-]*day[_\\s\\.-]*28[_\\s\\.-]?", "iN_D28_", nm, ignore.case = TRUE, perl = TRUE)

# (optional: if you have iNeurons* normalise it in here)
# nm <- sub("^iNeurons?[_\\s\\.-]*D[_\\s\\.-]*7[_\\s\\.-]?",  "iN_D7_",  nm, ignore.case = TRUE, perl = TRUE)
# nm <- sub("^iNeurons?[_\\s\\.-]*D[_\\s\\.-]*14[_\\s\\.-]?", "iN_D14_", nm, ignore.case = TRUE, perl = TRUE)
# nm <- sub("^iNeurons?[_\\s\\.-]*D[_\\s\\.-]*21[_\\s\\.-]?", "iN_D21_", nm, ignore.case = TRUE, perl = TRUE)
# nm <- sub("^iNeurons?[_\\s\\.-]*D[_\\s\\.-]*28[_\\s\\.-]?", "iN_D28_", nm, ignore.case = TRUE, perl = TRUE)

names(tether_raw) <- nm

# final pattern 
pat <- "^(iPSCs|iN_D7|iN_D14|iN_D21|iN_D28)_(\\d+)$"

expr_cols <- names(tether_raw)[ stringr::str_detect(names(tether_raw), pat) ]

if (length(expr_cols) == 0) {
  message("No column match the following format iPSCs_# / iN_D7_# / ... after alias.")
  print(head(names(tether_raw), 40))
  stop("Update alias here above based on raw names.")
}

# Pivot in long e find timepoint/rep
tether_long <- tether_raw %>%
  dplyr::select(gene = dplyr::all_of(gene_col), dplyr::all_of(expr_cols)) %>%
  tidyr::pivot_longer(-gene, names_to = "cond", values_to = "expression") %>%
  dplyr::mutate(
    timepoint = stringr::str_replace(cond, pat, "\\1"),
    rep       = as.integer(stringr::str_replace(cond, pat, "\\2")),
    timepoint = factor(timepoint, levels = c("iPSCs","iN_D7","iN_D14","iN_D21","iN_D28"))
  ) %>%
  tidyr::drop_na(timepoint)

# Summary per gene × timepoint (mean, median, SD, SEM, n)
tether_summary <- tether_long %>%
  dplyr::group_by(gene, timepoint) %>%
  dplyr::summarise(
    mean_expr   = mean(expression, na.rm = TRUE),
    median_expr = median(expression, na.rm = TRUE),
    sd_expr     = sd(expression,   na.rm = TRUE),
    sem_expr    = sd(expression,   na.rm = TRUE) / sqrt(sum(!is.na(expression))),
    n_reps      = sum(!is.na(expression)),
    .groups     = "drop"
  )

# Sanity check: all TP, including iPSCs
print(tether_long %>% dplyr::count(timepoint))

# Save
openxlsx::write.xlsx(tether_summary, file.path(out_dir, "02_tether_summary.xlsx"))

# Show in Environment 
assign("tether_long",    tether_long,    envir = .GlobalEnv)
assign("tether_summary", tether_summary, envir = .GlobalEnv)


# ----------------- 3) PAIRING: object,gene-----------------
suppressPackageStartupMessages({ library(readxl); library(janitor) })

pairing_path <- "pairing_tether.xlsx"     # <— Your table file
use_pairing  <- file.exists(pairing_path) # automatic toggle 

if (use_pairing) {
  pairing <- readxl::read_xlsx(pairing_path) %>%
    janitor::clean_names() %>%
    dplyr::rename(object = dplyr::any_of(c("object","contact","pair","contact_object")),
                  gene   = dplyr::any_of(c("gene","tether","tether_gene"))) %>%
    dplyr::select(object, gene) %>%
    dplyr::mutate(
      object = stringr::str_squish(as.character(object)),
      gene   = stringr::str_squish(as.character(gene))
    ) %>%
    dplyr::distinct()
  
  # Diagnostic of match
  objects_contacts <- contacts_summary %>% dplyr::distinct(object) %>% dplyr::mutate(in_contacts = TRUE)
  genes_tethers    <- tether_summary   %>% dplyr::distinct(gene)   %>% dplyr::mutate(in_tethers  = TRUE)
  
  pairing_diag <- pairing %>%
    dplyr::left_join(objects_contacts, by = "object") %>%
    dplyr::left_join(genes_tethers,    by = "gene") %>%
    dplyr::mutate(
      status = dplyr::case_when(
        is.na(in_contacts) & is.na(in_tethers) ~ "MISSING object & gene",
        is.na(in_contacts)                     ~ "MISSING object",
        is.na(in_tethers)                      ~ "MISSING gene",
        TRUE                                   ~ "OK"
      )
    )
  
  # salva report pairing
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  openxlsx::write.xlsx(pairing_diag, file.path(out_dir, "03_pairing_diagnostics.xlsx"))
  
  # tieni solo le coppie OK
  pairing_ok <- pairing_diag %>% dplyr::filter(status == "OK") %>% dplyr::select(object, gene) %>% dplyr::distinct()
  
  if (nrow(pairing_ok) == 0) {
    stop("No pair (object,gene) match with data. Check '03_pairing_diagnostics.xlsx'.")
  }
} else {
  pairing <- NULL
}

# ----------------- 4) MERGE (many-to-many per timepoint, then filter with pairing if present) -----------------
merged <- contacts_summary %>%
  dplyr::rename(timepoint = TimePoint) %>%
  dplyr::inner_join(tether_summary, by = "timepoint", relationship = "many-to-many")

# sif you have pairing, analysis will be limited to only the desired pairing
if (use_pairing) {
  merged <- merged %>% dplyr::semi_join(pairing_ok, by = c("object","gene"))
}

openxlsx::write.xlsx(merged, file.path(out_dir, "03_merged_long.xlsx"))
assign("merged", merged, envir = .GlobalEnv)

# --------
complete_df <- merged %>%
  dplyr::group_by(object, gene) %>%
  dplyr::filter(dplyr::n_distinct(timepoint) == length(tp_levels)) %>%
  dplyr::ungroup()

# save
openxlsx::write.xlsx(complete_df, file.path(out_dir, "04_complete_pairs.xlsx"))

# ---- Z-SCORED SERIES ----
z_series <- complete_df %>%
  dplyr::arrange(object, gene, timepoint) %>%
  dplyr::group_by(object, gene) %>%
  dplyr::mutate(
    z_expr   = as.numeric(scale(mean_expr)),
    z_count  = as.numeric(scale(count_mean)),
    z_volume = as.numeric(scale(volume_mean))
  ) %>%
  dplyr::ungroup() %>%
  dplyr::select(object, gene, timepoint,
                mean_expr, count_mean, volume_mean,
                z_expr, z_count, z_volume,
                sd_expr, count_sd, volume_sd,
                n_cells, n_reps)

openxlsx::write.xlsx(z_series, file.path(out_dir, "05_series_zscored.xlsx"))

# =========================
# complete DTW: distance, p-value (null permutation), uncertainty (bootstrap SEM)
# =========================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble); library(openxlsx)
})

# 0) generate complete_df
stopifnot(exists("complete_df"))
if ("TimePoint" %in% names(complete_df) && !"timepoint" %in% names(complete_df)) {
  complete_df <- dplyr::rename(complete_df, timepoint = TimePoint)
}
complete_df <- complete_df %>%
  mutate(
    object = as.character(object),
    gene   = as.character(gene),
    timepoint = as.character(timepoint)
  )

# 1) z-score  (if sd=0 -> all 0)
zsafe <- function(v) {
  v <- as.numeric(v); s <- stats::sd(v, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(0, length(v)))
  as.numeric(scale(v))
}

# 2) Distance DTW helper
dtw_dist <- function(x, y) dtw::dtw(x, y, step.pattern = dtw::symmetric2)$distance

# 3) General Runner :
#    - ycol: "count_mean" or "volume_mean"
#    - ysem: "count_sem"  or "volume_sem"
#    - p-value: null = Temporal permutation of X; otherwise Monte Carlo
#    - uncertainty: bootstrap SEM (for each TP N(mean, SEM))
run_dtw_full <- function(df, ycol, ysem, metric_label,
                         use_weights = FALSE,            # TRUE: down-weight with SEM (more SEM -> less weight)
                         do_perm = TRUE, B_null = 5000, # perm/null
                         do_boot = TRUE, B_boot = 2000, # uncertainty bootstrap 
                         eps = 1e-9) {
  
  need <- c("object","gene","timepoint","mean_expr","sem_expr", ycol, ysem)
  miss <- setdiff(need, names(df))
  if (length(miss)) stop("Missing columns in df: ", paste(miss, collapse=", "))
  
  # Measuring distance (weighted optional) given a frame (pair)
  dist_one <- function(x_mean, x_sem, y_mean, y_sem) 
    # z-score 
    xe <- zsafe(x_mean)
    ye <- zsafe(y_mean)
    
    if (use_weights) {
      # weight = 1/SEM for down-weight of uncertain points; normalised on stability
      wx <- 1 / (x_sem + eps)
      wy <- 1 / (y_sem + eps)
      wx <- wx / median(wx[is.finite(wx) & wx > 0], na.rm=TRUE)
      wy <- wy / median(wy[is.finite(wy) & wy > 0], na.rm=TRUE)
      # weighted series
      xw <- xe * wx
      yw <- ye * wy
    } else {
      xw <- xe; yw <- ye
    }
    
    keep <- is.finite(xw) & is.finite(yw)
    xw <- xw[keep]; yw <- yw[keep]
    if (length(xw) < 3L) return(NA_real_)
    dtw_dist(xw, yw)
  }
  
  out <- df %>%
    group_by(object, gene) %>%
    arrange(timepoint, .by_group = TRUE) %>%
    reframe({
      x_mean <- mean_expr; x_sem <- sem_expr
      y_mean <- .data[[ycol]]; y_sem <- .data[[ysem]]
      
      # Observed distance
      d_obs <- dist_one(x_mean, x_sem, y_mean, y_sem)
      
      if (!isTRUE(do_perm) || !is.finite(d_obs)) {
        p_emp <- NA_real_; null_mean <- NA_real_; null_sd <- NA_real_; dist_z <- NA_real_
      } else {
        n <- length(x_mean)
        # Builds weighted series after null
        xe <- zsafe(x_mean); ye <- zsafe(y_mean)
        if (use_weights) {
          wx <- 1 / (x_sem + eps); wy <- 1 / (y_sem + eps)
          wx <- wx / median(wx[is.finite(wx) & wx > 0], na.rm=TRUE)
          wy <- wy / median(wy[is.finite(wy) & wy > 0], na.rm=TRUE)
          xw <- xe * wx; yw <- ye * wy
        } else {
          xw <- xe; yw <- ye
        }
        keep <- is.finite(xw) & is.finite(yw); xw <- xw[keep]; yw <- yw[keep]; n <- length(xw)
        
        if (n < 3L) {
          p_emp <- NA_real_; null_mean <- NA_real_; null_sd <- NA_real_; dist_z <- NA_real_
        } else {
          # null: permutation order xw
          if (n == 5 && requireNamespace("gtools", quietly = TRUE)) {
            ords <- gtools::permutations(n, n)
            dnull <- apply(ords, 1, function(o) dtw_dist(xw[o], yw))
          } else {
            dnull <- replicate(B_null, { o <- sample.int(n, n); dtw_dist(xw[o], yw) })
          }
          null_mean <- mean(dnull)
          null_sd   <- stats::sd(dnull)
          p_emp     <- (1 + sum(dnull <= d_obs)) / (length(dnull) + 1)  # add-one smoothing
          dist_z    <- if (is.finite(null_sd) && null_sd > 0) (d_obs - null_mean)/null_sd else NA_real_
        }
      }
      
      # bootstrap (uncertainty based on SEM/time point)
      if (isTRUE(do_boot) && is.finite(d_obs)) {
        # generate means from N(mean, SEM) for each TP and calculate distance
        db <- replicate(B_boot, {
          xb <- rnorm(length(x_mean), mean = x_mean, sd = pmax(x_sem, eps))
          yb <- rnorm(length(y_mean), mean = y_mean, sd = pmax(y_sem, eps))
          dist_one(xb, x_sem, yb, y_sem)
        })
        sd_boot   <- stats::sd(db, na.rm=TRUE)
        ci_low    <- stats::quantile(db, 0.025, na.rm=TRUE)
        ci_high   <- stats::quantile(db, 0.975, na.rm=TRUE)
      } else {
        sd_boot <- NA_real_; ci_low <- NA_real_; ci_high <- NA_real_
      }
      
      tibble(
        dist = as.numeric(d_obs),
        pval = as.numeric(p_emp),
        null_mean = as.numeric(null_mean),
        null_sd   = as.numeric(null_sd),
        dist_z    = as.numeric(dist_z),
        sd_boot   = as.numeric(sd_boot),
        ci95_low  = as.numeric(ci_low),
        ci95_high = as.numeric(ci_high)
      )
    }) %>%
    ungroup() %>%
    mutate(metric = metric_label)
  
  out

# 4) Analysis per COUNT and VOLUME
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

dtw_count_full  <- run_dtw_full(complete_df, "count_mean",  "count_sem",  "count",
                                use_weights = TRUE, do_perm = TRUE, B_null = 5000, do_boot = TRUE, B_boot = 2000)
dtw_volume_full <- run_dtw_full(complete_df, "volume_mean", "volume_sem", "volume",
                                use_weights = TRUE, do_perm = TRUE, B_null = 5000, do_boot = TRUE, B_boot = 2000)

# 5) Ranking + FDR e export
rank_and_fdr <- function(tab){
  tab %>%
    arrange(metric, dist) %>%
    group_by(metric) %>%
    mutate(
      rank  = row_number(),
      p_fdr = if_else(is.na(pval), NA_real_, p.adjust(pval, method = "BH"))
    ) %>% ungroup()
}

dtw_full_ranked <- rank_and_fdr(bind_rows(dtw_count_full, dtw_volume_full))

openxlsx::write.xlsx(
  list(
    DTW_full_raw     = bind_rows(dtw_count_full, dtw_volume_full),
    DTW_full_ranked  = dtw_full_ranked
  ),
  file.path(out_dir, "06_DTW_full_distance_p_uncertainty.xlsx")
)

# check 
dtw_full_ranked %>% select(object, gene, metric, dist, pval, p_fdr, sd_boot, ci95_low, ci95_high, dist_z) %>% head()


# ----------------- 7) Optional: Spearman -----------------
if (isTRUE(do_spearman)) {
  run_spearman <- function(df, ycol, metric_label){
    df %>%
      dplyr::group_by(object, gene) %>%
      dplyr::arrange(TimePoint, .by_group = TRUE) %>%
      dplyr::summarise(
        rho = suppressWarnings(cor(mean_expr, .data[[ycol]], method = "spearman", use = "pairwise.complete.obs")),
        p   = suppressWarnings(cor.test(mean_expr, .data[[ycol]], method = "spearman")$p.value),
        .groups = "drop"
      ) %>%
      dplyr::mutate(metric = metric_label)
  }
  sp_count  <- run_spearman(complete, "count_mean",  "count")  %>% dplyr::mutate(FDR = p.adjust(p, method = "BH"))
  sp_volume <- run_spearman(complete, "volume_mean", "volume") %>% dplyr::mutate(FDR = p.adjust(p, method = "BH"))
  openxlsx::write.xlsx(list(
    Spearman_count  = sp_count,
    Spearman_volume = sp_volume
  ), file.path(out_dir, "07_Correlations_Spearman.xlsx"))
  assign("spearman_count", sp_count, envir = .GlobalEnv)
  assign("spearman_volume", sp_volume, envir = .GlobalEnv)
}

# =========================
# PREP DATI PER PLOT 
# =========================
suppressPackageStartupMessages({library(dplyr); library(tidyr); library(ggplot2)})

stopifnot(exists("complete_df"))

# name from TP
dfp <- complete_df
if ("TimePoint" %in% names(dfp) && !"timepoint" %in% names(dfp)) {
  dfp <- dplyr::rename(dfp, timepoint = TimePoint)
}

# SEM (if missing, derived; or put 0)
if (!"sem_expr" %in% names(dfp)) {
  if (all(c("sd_expr","n_reps") %in% names(dfp))) {
    dfp <- dfp %>% mutate(sem_expr = sd_expr / sqrt(pmax(n_reps, 1)))
  } else dfp <- dfp %>% mutate(sem_expr = 0)
}
if (!"count_sem" %in% names(dfp)) {
  if (all(c("count_sd","n_cells") %in% names(dfp))) {
    dfp <- dfp %>% mutate(count_sem = count_sd / sqrt(pmax(n_cells, 1)))
  } else if (all(c("count_sd","n_reps") %in% names(dfp))) {
    dfp <- dfp %>% mutate(count_sem = count_sd / sqrt(pmax(n_reps, 1)))
  } else dfp <- dfp %>% mutate(count_sem = 0)
}
if (!"volume_sem" %in% names(dfp)) {
  if (all(c("volume_sd","n_cells") %in% names(dfp))) {
    dfp <- dfp %>% mutate(volume_sem = volume_sd / sqrt(pmax(n_cells, 1)))
  } else if (all(c("volume_sd","n_reps") %in% names(dfp))) {
    dfp <- dfp %>% mutate(volume_sem = volume_sd / sqrt(pmax(n_reps, 1)))
  } else dfp <- dfp %>% mutate(volume_sem = 0)
}

# keep columns
base_df <- dfp %>%
  select(object, gene, timepoint,
         mean_expr, sd_expr, sem_expr,
         count_mean, count_sd, count_sem,
         volume_mean, volume_sd, volume_sem)

# optional: fixed time point onX axis
tp_levels <- c("iPSCs","iN_D7","iN_D14","iN_D21","iN_D28")
base_df <- base_df %>% mutate(timepoint = factor(timepoint, levels = tp_levels))

# =========================
# RIBBON PLOTS: build base_df, define plot fn, save ALL plots
# Requires: complete_df with columns:
#  object, gene, timepoint (or TimePoint), mean_expr, sd_expr (or sem_expr),
#  count_mean, count_sd (or count_sem), volume_mean, volume_sd (or volume_sem)
# =========================
suppressPackageStartupMessages({library(dplyr); library(tidyr); library(ggplot2)})

# 0) output dir
if (!exists("out_dir")) out_dir <- "pairing_analysis_outputs"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
plots_dir <- file.path(out_dir, "plots_ribbon")
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

# 1) base_df from complete_df (handles TimePoint -> timepoint, and SEM if missing)
stopifnot(exists("complete_df"))
dfp <- complete_df
if ("TimePoint" %in% names(dfp) && !"timepoint" %in% names(dfp)) {
  dfp <- dplyr::rename(dfp, timepoint = TimePoint)
}

# --- EXPORT BASE PER I PLOT ---
if (!exists("out_dir")) out_dir <- "pairing_analysis_outputs"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

openxlsx::write.xlsx(base_df, file.path(out_dir, "07_base_df_used_for_plots.xlsx"))
readr::write_csv(base_df,       file.path(out_dir, "07_base_df_used_for_plots.csv"))

if (exists("z_series")) {
  openxlsx::write.xlsx(z_series, file.path(out_dir, "07_z_series_used_for_plots.xlsx"))
  readr::write_csv(z_series,     file.path(out_dir, "07_z_series_used_for_plots.csv"))
}

# derive SEM where missing (fallback 0)
if (!"sem_expr" %in% names(dfp)) {
  if (all(c("sd_expr","n_reps") %in% names(dfp))) {
    dfp <- dfp %>% mutate(sem_expr = sd_expr / sqrt(pmax(n_reps, 1)))
  } else dfp <- dfp %>% mutate(sem_expr = 0)
}
if (!"count_sem" %in% names(dfp)) {
  if (all(c("count_sd","n_cells") %in% names(dfp))) {
    dfp <- dfp %>% mutate(count_sem = count_sd / sqrt(pmax(n_cells, 1)))
  } else if (all(c("count_sd","n_reps") %in% names(dfp))) {
    dfp <- dfp %>% mutate(count_sem = count_sd / sqrt(pmax(n_reps, 1)))
  } else dfp <- dfp %>% mutate(count_sem = 0)
}
if (!"volume_sem" %in% names(dfp)) {
  if (all(c("volume_sd","n_cells") %in% names(dfp))) {
    dfp <- dfp %>% mutate(volume_sem = volume_sd / sqrt(pmax(n_cells, 1)))
  } else if (all(c("volume_sd","n_reps") %in% names(dfp))) {
    dfp <- dfp %>% mutate(volume_sem = volume_sd / sqrt(pmax(n_reps, 1)))
  } else dfp <- dfp %>% mutate(volume_sem = 0)
}

tp_levels <- c("iPSCs","iN_D7","iN_D14","iN_D21","iN_D28")
base_df <- dfp %>%
  select(object, gene, timepoint,
         mean_expr, sd_expr, sem_expr,
         count_mean, count_sd, count_sem,
         volume_mean, volume_sd, volume_sem) %>%
  mutate(timepoint = factor(as.character(timepoint), levels = tp_levels))

# 2) plotting function (prints to Viewer + optionally saves PNG)
plot_ribbon_pair <- function(object_id, gene_id,
                             metric = c("count","volume"),
                             spread = c("sem","sd"),
                             out_file = NULL,
                             col_expr = "#1f77b4", fill_expr = "#1f77b4",
                             col_y    = "#ff7f0e", fill_y    = "#ff7f0e") {
  metric <- match.arg(metric); spread <- match.arg(spread)
  dat <- base_df %>% dplyr::filter(object == object_id, gene == gene_id) %>% dplyr::arrange(timepoint)
  if (nrow(dat) == 0) stop("No rows for requested object/gene.")
  
  x_mean <- dat$mean_expr
  y_mean <- if (metric == "count") dat$count_mean else dat$volume_mean
  x_sp <- if (spread == "sem") dat$sem_expr else dat$sd_expr
  y_sp <- if (metric == "count") { if (spread == "sem") dat$count_sem else dat$count_sd
  } else                     { if (spread == "sem") dat$volume_sem else dat$volume_sd }
  
  sx <- stats::sd(x_mean, na.rm=TRUE); sy <- stats::sd(y_mean, na.rm=TRUE)
  zx <- if (is.finite(sx) && sx>0) as.numeric(scale(x_mean)) else rep(0, length(x_mean))
  zy <- if (is.finite(sy) && sy>0) as.numeric(scale(y_mean)) else rep(0, length(y_mean))
  zx_sp <- if (is.finite(sx) && sx>0) x_sp / sx else rep(0, length(x_sp))
  zy_sp <- if (is.finite(sy) && sy>0) y_sp / sy else rep(0, length(y_sp))
  
  plot_df <- tibble::tibble(
    timepoint = dat$timepoint,
    z_expr = zx,    z_expr_lo = zx - zx_sp, z_expr_hi = zx + zx_sp,
    z_y    = zy,    z_y_lo    = zy - zy_sp, z_y_hi    = zy + zy_sp
  )
  
  p <- ggplot(plot_df, aes(x = timepoint)) +
    geom_ribbon(aes(ymin = z_expr_lo, ymax = z_expr_hi, group = 1),
                fill = scales::alpha(fill_expr, 0.20)) +
    geom_line(aes(y = z_expr, group = 1), linewidth = 1.2, color = col_expr) +
    geom_point(aes(y = z_expr), size = 2.4, color = col_expr) +
    geom_ribbon(aes(ymin = z_y_lo, ymax = z_y_hi, group = 2),
                fill = scales::alpha(fill_y, 0.18)) +
    geom_line(aes(y = z_y, group = 2), linewidth = 1.2, color = col_y, linetype = "dashed") +
    geom_point(aes(y = z_y), size = 2.4, color = col_y, shape = 17) +
    labs(
      title = sprintf("%s × %s  (%s, ribbon=%s)", object_id, gene_id, metric, toupper(spread)),
      y = "z-score", x = "Time point",
      caption = "Solid: Tether expression (z)\nDashed: Contacts (z)"
    ) +
    theme_minimal(base_size = 13) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold"),
          axis.title.x = element_text(margin = margin(t = 6)))
  
  if (!is.null(out_file)) ggplot2::ggsave(out_file, p, width = 8.5, height = 5, dpi = 300)
  print(p) # ensure it shows in RStudio Viewer
  invisible(p)
}

# 3) choose which pairs to plot:
#    Option A: use ranked DTW (preferred: top pairs first)
rank_tab <- get0("dtw_full_ranked", ifnotfound = get0("dtw_ranked", ifnotfound = NULL))
if (!is.null(rank_tab)) {
  pairs_tbl <- rank_tab %>%
    arrange(metric, dist) %>%
    select(object, gene, metric) %>%
    distinct()
} else {
  # Option B: plot every object×gene present in base_df (will be many!)
  pairs_tbl <- base_df %>% distinct(object, gene) %>%
    tidyr::crossing(metric = c("count","volume"))
}

# --- HELPER: exact data for plotting (linee + ribbon) ---
build_plotdata <- function(object_id, gene_id,
                           metric = c("count","volume"),
                           spread = c("sem","sd")) {
  metric <- match.arg(metric); spread <- match.arg(spread)
  
  dat <- base_df %>%
    dplyr::filter(object == object_id, gene == gene_id) %>%
    dplyr::arrange(timepoint)
  if (nrow(dat) == 0) stop("no row per object/gene.")
  
  x_mean <- dat$mean_expr
  y_mean <- if (metric == "count") dat$count_mean else dat$volume_mean
  
  x_sp <- if (spread == "sem") dat$sem_expr else dat$sd_expr
  y_sp <- if (metric == "count") {
    if (spread == "sem") dat$count_sem else dat$count_sd
  } else {
    if (spread == "sem") dat$volume_sem else dat$volume_sd
  }
  
  sx <- stats::sd(x_mean, na.rm=TRUE); sy <- stats::sd(y_mean, na.rm=TRUE)
  zx <- if (is.finite(sx) && sx>0) as.numeric(scale(x_mean)) else rep(0, length(x_mean))
  zy <- if (is.finite(sy) && sy>0) as.numeric(scale(y_mean)) else rep(0, length(y_mean))
  zx_sp <- if (is.finite(sx) && sx>0) x_sp / sx else rep(0, length(x_sp))
  zy_sp <- if (is.finite(sy) && sy>0) y_sp / sy else rep(0, length(y_sp))
  
  tibble::tibble(
    object    = object_id,
    gene      = gene_id,
    metric    = metric,
    spread    = spread,
    timepoint = dat$timepoint,
    # raw
    expr_mean = x_mean, expr_spread = x_sp,
    y_mean    = y_mean, y_spread    = y_sp,
    # z-scale (used for plots)
    z_expr    = zx,  z_expr_lo = zx - zx_sp, z_expr_hi = zx + zx_sp,
    z_y       = zy,  z_y_lo    = zy - zy_sp, z_y_hi    = zy + zy_sp
  )
}

# 4) SAVE ALL plots (SEM ribbons by default). Also preview the first one in Viewer.
options(device = "RStudioGD")
n <- nrow(pairs_tbl)
message("Generating ribbon plots for ", n, " pairs...")
for (i in seq_len(n)) {
  rr <- pairs_tbl[i,]
  f <- file.path(plots_dir, sprintf("ribbon_%s_%s_%s_SEM.png", rr$object, rr$gene, rr$metric))
  try(plot_ribbon_pair(rr$object, rr$gene, metric = rr$metric, spread = "sem", out_file = f), silent = TRUE)
  if (i == 1L) { # preview the very first
    plot_ribbon_pair(rr$object, rr$gene, metric = rr$metric, spread = "sem")
  }
}
message("Done. PNGs saved in: ", normalizePath(plots_dir, winslash = "/"))

#---------------------------------------------------------------------------------
# --- Export all csv per each plot ---
#plots_data_dir <- file.path(out_dir, "plot_data")
#dir.create(plots_data_dir, showWarnings = FALSE, recursive = TRUE)

#rank_tab <- get0("dtw_full_ranked", ifnotfound = get0("dtw_ranked", ifnotfound = NULL))
#if (!is.null(rank_tab)) {
#  pairs_tbl <- rank_tab %>%
#    dplyr::arrange(metric, dist) %>%
#    dplyr::select(object, gene, metric) %>%
#    dplyr::distinct()
#} else {
#  pairs_tbl <- base_df %>% dplyr::distinct(object, gene) %>%
#    tidyr::crossing(metric = c("count","volume"))
#}
#
#for (i in seq_len(nrow(pairs_tbl))) {
#  rr <- pairs_tbl[i,]
#  pd <- build_plotdata(rr$object, rr$gene, metric = rr$metric, spread = "sem")
#  fname <- sprintf("plotdata_%s_%s_%s_SEM.csv",
#                   gsub("[^A-Za-z0-9._-]","_", rr$object),
#                   gsub("[^A-Za-z0-9._-]","_", rr$gene),
#                   rr$metric)
#  readr::write_csv(pd, file.path(plots_data_dir, fname))
#}
#
#cat("Export completato.\n",
#    "Base DF: ", normalizePath(file.path(out_dir, "07_base_df_used_for_plots.xlsx")), "\n",
#    "Per-plot CSV in: ", normalizePath(plots_data_dir), "\n")
