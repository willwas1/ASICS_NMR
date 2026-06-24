# =============================================================================
#  NMR METABOLOMICS PIPELINE (ASICS, built-in library)  -- v25
#  800 MHz Bruker spectra  ->  absolute metabolite concentrations (mM)
#  Developed : Will Smith - william.smith-10@manchester.ac.uk
# =============================================================================
#  An R script to automate metabolite detection using ASICS R package:
#  Tardivel, P.J.C., Canlet, C., Lefort, G. et al. ASICS: an automatic method 
#  for identification and quantification of metabolites in complex 1D 1H NMR spectra. 
#  Metabolomics 13, 109 (2017). https://doi.org/10.1007/s11306-017-1244-5
#
#  https://www.bioconductor.org/packages/release/bioc/html/ASICS.html
#
#  Developed on :
#  R version 4.6.0
#  ASICS v2.28.0
#
#  WHAT THIS SCRIPT DOES
#  ---------------------
#  Takes a folder of TopSpin-processed Bruker spectra and produces, per sample:
#    * A panel of 30 metabolites in mM  (manual qHNMR,
#      TSP-anchored, dilution-corrected)  <-- this is the number you report.
#    * A broad ASICS deconvolution table (193-metabolite library) for screening
#      and pathway-level trends only.
#    * QC plots, a run manifest, and a focused CHO-target table.
#
#  #  ---------------------------------------------------------------------------
#  BEFORE YOU RUN (one-time setup)
#  ---------------------------------------------------------------------------
#   1. Install R (>= 4.2) and RStudio.
#   2. it auto-installs the packages Listed in Section 0 
#      (ASICS comes from Bioconductor; can take 10-20 min).
#   3. Your spectra must be processed (phase etc.) in TopSpin (each sample folder needs
#      a "1r" file). 
#   4. Have a dilution CSV ready (columns: sample, dilution_factor). See below.
#
#  ---------------------------------------------------------------------------
#  HOW TO RUN
#  ---------------------------------------------------------------------------
#   1. Edit ONLY the "SETTINGS YOU MUST CHECK" block in Section 1.
#   2. In RStudio: Code > Source  (or click "Source", top-right of the editor).
#   3. Watch the Console. A quick path check runs first; if a path is wrong it
#      stops immediately with a plain-English message telling you what to fix.
#   4. Results land in OUT_DIR. Nothing in your data folder is modified.
#
#  ---------------------------------------------------------------------------
#  Output
#  ---------------------------------------------------------------------------
#   manual_panel_real_mM.csv .... mM concentrations, 30 metabolites.
#   manual_panel_tube_mM.csv .... Same, before dilution correction (in-tube mM).
#   qc_panel_windows.pdf ........ Shows every integration window on
#                                 the median spectrum (black=high confidence,
#                                 orange=moderate, red=upper-bound). If a window
#                                 sits off-peak, nudge it in Section 7.
#   absolute_ASICS_crosscalibrated_real_mM.csv .. Screening table (use for trends).
#   CHO_target_metabolites_*.csv  Focused CHO panel pulled from the ASICS table.
#   run_manifest.txt / sessionInfo.txt .. Exactly how this run was configured.
#
#  ---------------------------------------------------------------------------
#  IF SOMETHING GOES WRONG
#  ---------------------------------------------------------------------------
#   * "SPECTRA_DIR not found"  -> fix the path in Section 1. Use FORWARD slashes
#     "/" even on Windows (copy the Explorer path and flip the backslashes).
#   * Run seems frozen at quantification -> almost always an unprocessed sample
#     (missing "1r"). Reprocess that sample in TopSpin; it is NOT a core/RAM issue.
#   * "DILUTION_FILE matched 0 samples" -> the 'sample' column must match the
#     Bruker FOLDER NUMBERS (e.g. 10, 11, 12 ...), not sample names.
#   * Slow file reads -> if data sits in OneDrive, copy it to a local disk first.
#
# # =============================================================================


# -- 0. Packages (auto-installed on first run) --------------------------------
required <- c("ASICS", "ggplot2", "dplyr", "tidyr", "pheatmap")
for (p in required) if (!requireNamespace(p, quietly = TRUE)) {
  message("First run: installing '", p, "' (this can take a while) ...")
  if (p == "ASICS") { if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager"); BiocManager::install("ASICS", ask = FALSE) }
  else install.packages(p, repos = "https://cloud.r-project.org")
}
library(ASICS); library(ggplot2); library(dplyr); library(tidyr); library(pheatmap)


# =============================================================================
# -- 1. SETTINGS --------------------------------------------------------------
# =============================================================================

# #############################################################################
# ##   >>>>>>>>>>   SETTINGS YOU MUST CHECK FOR EACH DATASET   <<<<<<<<<<     ##
# ##   Use FORWARD slashes "/" in every path, even on Windows.               ##
# #############################################################################

# 1) Folder of TopSpin-processed Bruker spectra (the parent of the numbered
#    sample folders). Each sample folder must contain a processed "1r" file.
SPECTRA_DIR <- "C:/data/Topspin_processed"

# 2) Where to write results. Created automatically if it does not exist.
#    TIP: give each run its own folder so outputs never overwrite each other.
OUT_DIR     <- "C:/data/results_name"

# 3) Dilution factors, as a CSV with EXACTLY two columns: sample, dilution_factor
#       sample           = Bruker folder number (10, 11, 12, ...)
#       dilution_factor  = total dilution for that tube
#                          (this study: 5.5 for 1:4.5 sample:PBS/D2O/TSP,
#                                       1.1 for 500:50 preparations)
#    Set DILUTION_FILE <- NULL to apply DEFAULT_DF to every sample.
DILUTION_FILE <- "C:/data/asics_dilution.csv"
DEFAULT_DF    <- 1.1     # used for any sample not listed in the dilution file

# 4) Optional sample metadata CSV (sample, condition, time_h, biorep).
#    Leave as NULL if you do not have one.
META_FILE     <- NULL

# ##   end of the block most users need to edit.                            ##



# ---- EXPERIMENT SETTINGS (check once per study, then leave alone) -----------
FIELD_MHZ      <- 800           # spectrometer field
ALPHA_FRACTION <- 0.36          # alpha-D-glucose anomeric fraction (~0.36 @ 25 C)

# TSP internal standard (final in-tube concentration is the absolute anchor)
TSP_PCT_WV     <- 0.01          # % w/v TSP actually in the NMR tube
TSP_IS_D4      <- TRUE          # TRUE = TSP-d4 (MW 172.27); FALSE = TSP (162.20)
TSP_PROTONS    <- 9             # equivalent protons giving the 0 ppm singlet
TSP_WINDOW     <- c(-0.10, 0.10)# ppm integration window for the TSP peak
TSP_IN_TUBE_MM <- (TSP_PCT_WV * 10) / (if (TSP_IS_D4) 172.27 else 162.20) * 1000


# ---- ADVANCED (sensible defaults; change only if you know why) --------------
STOP_IF_DILUTION_BAD <- TRUE    # halt if the dilution file matches no samples
N_CORES        <- min(11L, max(1L, parallel::detectCores() - 1L))  # parallelism
EXCLUSION_AREAS <- matrix(c(4.7, 5.0), ncol = 2, byrow = TRUE)      # water region
NOISE_THRES <- 0.02; ADD_NOISE <- 0.05; MULT_NOISE <- 0.172
MAX_SHIFT   <- 0.05; CLEAN_THRES <- 10; DO_ALIGN <- TRUE; NORM_METHOD <- "pqn"
DESIGN_FILE  <- NULL            # reserved (not required by this pipeline)
GROUP_COLUMN <- "condition"     # reserved (used only if META_FILE has groups)

# Focused CHO target panel (post-hoc extraction from the full ASICS table)
TARGET_PATTERNS <- c(
  glucose="glucose", lactate="lactate|lactic", alanine="alanine",
  glutamine="glutamine", glutamate="glutamate|glutamic",
  aspartate="aspartate|aspartic", asparagine="asparagine", glycine="glycine",
  serine="serine", threonine="threonine", valine="valine", leucine="leucine",
  isoleucine="isoleucine", lysine="lysine", histidine="histidine",
  arginine="arginine", methionine="methionine", phenylalanine="phenylalanine",
  tyrosine="tyrosine", tryptophan="tryptophan", proline="proline",
  ornithine="ornithine", citrulline="citrulline", acetate="acetate|acetic",
  formate="formate|formic", pyruvate="pyruvate|pyruvic", citrate="citrate|citric",
  succinate="succinate|succinic", fumarate="fumarate|fumaric", malate="malate|malic",
  choline="choline", phosphocholine="phosphocholine",
  glycerophosphocholine="glycerophosphocholine",
  myo_inositol="myo.?inositol|inositol", taurine="taurine", betaine="betaine",
  creatine="creatine", creatinine="creatinine", uridine="uridine",
  hypoxanthine="hypoxanthine", adenosine="adenosine")


# ---- Quick path check: fail fast & clearly BEFORE the long run --------------
.check_inputs <- function() {
  ok <- TRUE
  if (!dir.exists(SPECTRA_DIR)) {
    message("X  SPECTRA_DIR not found:\n     ", SPECTRA_DIR,
            "\n   Fix the path in Section 1 (use forward slashes '/').")
    ok <- FALSE
  } else message("OK  SPECTRA_DIR found.")
  if (!is.null(DILUTION_FILE)) {
    if (!file.exists(DILUTION_FILE)) {
      message("X  DILUTION_FILE not found:\n     ", DILUTION_FILE,
              "\n   Fix the path, or set DILUTION_FILE <- NULL to use DEFAULT_DF.")
      ok <- FALSE
    } else message("OK  DILUTION_FILE found.")
  } else message("--  DILUTION_FILE is NULL: DEFAULT_DF (", DEFAULT_DF,
                 ") applied to all samples.")
  if (!is.null(META_FILE) && !file.exists(META_FILE)) {
    message("X  META_FILE set but not found:\n     ", META_FILE,
            "\n   Fix the path, or set META_FILE <- NULL."); ok <- FALSE
  }
  dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
  if (!dir.exists(OUT_DIR)) {
    message("X  Could not create OUT_DIR:\n     ", OUT_DIR); ok <- FALSE
  } else message("OK  OUT_DIR ready:  ", OUT_DIR)
  if (!ok) stop("Please fix the path(s) above in Section 1, then Source again.",
                call. = FALSE)
  message("All input paths OK -- starting the run.\n")
}
.check_inputs()
message("Using ", N_CORES, " of ", parallel::detectCores(), " cores.")
message(sprintf("TSP in-tube anchor: %.3f mM", TSP_IN_TUBE_MM))


# -- 2. Helper functions -------------------------------------------------------
integrate_region <- function(ppm, y, from, to) {
  lo <- min(from, to); hi <- max(from, to)
  idx <- which(ppm >= lo & ppm <= hi)
  if (length(idx) < 2) return(NA_real_)
  x <- ppm[idx]; yy <- y[idx]; o <- order(x); x <- x[o]; yy <- yy[o]
  sum(diff(x) * (head(yy, -1) + tail(yy, -1)) / 2)
}
# anchored grep: pattern not preceded by a letter (alanine!=phenylalanine etc.)
grep_anchored <- function(pat, names)
  grep(paste0("(?<![[:alpha:]])(?:", pat, ")"), names, ignore.case = TRUE, perl = TRUE)

read_dilution <- function(path, samples, default_df) {
  df <- setNames(rep(default_df, length(samples)), samples)
  if (!is.null(path) && file.exists(path)) {
    d <- read.csv(path, stringsAsFactors = FALSE, colClasses = "character")
    if (!all(c("sample","dilution_factor") %in% names(d)))
      stop("DILUTION_FILE must have columns: sample, dilution_factor")
    d$sample <- trimws(d$sample); m <- match(d$sample, samples); ok <- !is.na(m)
    df[m[ok]] <- as.numeric(d$dilution_factor[ok])
    attr(df, "n_matched") <- sum(ok); attr(df, "n_in_file") <- nrow(d)
  } else { attr(df, "n_matched") <- 0L; attr(df, "n_in_file") <- NA_integer_ }
  df
}
load_metadata_flexible <- function(meta_file) {
  if (is.null(meta_file) || !file.exists(meta_file)) return(NULL)
  meta <- read.csv(meta_file, stringsAsFactors = FALSE)
  nm <- gsub("\\.", "_", names(meta)); nm <- gsub("time_H","time_h",nm,fixed=TRUE)
  nm <- gsub("diltution_factor","dilution_factor",nm,fixed=TRUE); names(meta) <- nm
  if (!("sample" %in% names(meta))) return(NULL)
  if (!("condition" %in% names(meta))) meta$condition <- "Unknown"
  if (!("time_h" %in% names(meta)))    meta$time_h <- NA_real_
  if (!("biorep" %in% names(meta)))    meta$biorep <- NA
  meta$sample <- trimws(as.character(meta$sample)); meta
}


# -- 3. Built-in library -------------------------------------------------------
message("Loading ASICS built-in pure library ...")
builtin_library <- ASICS::pure_library
message("  ", length(getSampleName(builtin_library)), " metabolites (curated).")


# -- 4. Import + pre-flight validation + preprocess ----------------------------
message("Importing all spectra ...")
spectra_data <- importSpectraBruker(SPECTRA_DIR, which.spectra = NULL)
spectra_ppm  <- as.numeric(rownames(spectra_data))
spectra_raw  <- spectra_data
sample_ids   <- colnames(spectra_data)
message("Imported: ", ncol(spectra_data), " spectra | ", nrow(spectra_data), " ppm points")

message("\n=== Pre-flight file validation (before the long run) ===")
message("  Sample IDs: ", paste(utils::head(sample_ids, 8), collapse = ", "),
        if (length(sample_ids) > 8) ", ..." else "")
dilution_df <- read_dilution(DILUTION_FILE, sample_ids, DEFAULT_DF)
if (!is.null(DILUTION_FILE)) {
  if (!file.exists(DILUTION_FILE)) {
    m <- paste0("DILUTION_FILE not found: ", DILUTION_FILE)
    if (STOP_IF_DILUTION_BAD) stop(m) else warning(m)
  } else {
    nmatch <- attr(dilution_df, "n_matched")
    message(sprintf("  DILUTION_FILE: %d rows, %d matched sample IDs.",
                    attr(dilution_df, "n_in_file"), nmatch))
    if (nmatch == 0) {
      m <- paste0("DILUTION_FILE matched 0 samples -- 'sample' must be folder ",
                  "numbers like ", sample_ids[1], ".")
      if (STOP_IF_DILUTION_BAD) stop(m) else warning(m)
    }
  }
}
meta <- load_metadata_flexible(META_FILE)
if (!is.null(meta)) message(sprintf("  METADATA: %d/%d sample IDs matched.",
                                    sum(sample_ids %in% meta$sample), length(sample_ids)))
write.csv(data.frame(sample = names(dilution_df), dilution_factor = as.numeric(dilution_df)),
          file.path(OUT_DIR, "dilution_factors_used.csv"), row.names = FALSE)
message("=== Pre-flight OK ===\n")

message("Normalising (", NORM_METHOD, ") ...")
spectra_data <- normaliseSpectra(spectra_data, type.norm = NORM_METHOD)
if (DO_ALIGN) { message("Aligning spectra ..."); spectra_data <- alignSpectra(spectra_data) }
spectra_obj <- createSpectra(spectra_data)


# -- 5. ASICS quantification ---------------------------------------------------
n_spectra <- ncol(spectra_data)
log_file  <- file.path(OUT_DIR, "progress_log.txt")
sentinel  <- file.path(OUT_DIR, ".timer_stop")
if (file.exists(sentinel)) file.remove(sentinel)
message("Running ASICS quantification ...  Spectra: ", n_spectra,
        " | Library: ", length(getSampleName(builtin_library)), " | Cores: ", N_CORES)
t_start <- Sys.time()
writeLines(c("ASICS v25 (built-in) -- RUNNING",
             paste(" Started:", format(t_start, "%Y-%m-%d %H:%M:%S")),
             paste(" Spectra:", n_spectra, "| Cores:", N_CORES)), log_file)

timer_script <- file.path(OUT_DIR, "timer.R")
writeLines(c(
  sprintf('t_start  <- as.POSIXct("%s")', format(t_start, "%Y-%m-%d %H:%M:%S")),
  sprintf('log_file <- "%s"', gsub("\\\\", "/", log_file)),
  sprintf('sentinel <- "%s"', gsub("\\\\", "/", sentinel)),
  sprintf('est_total_m <- (%d * 3) / %d', n_spectra, N_CORES),
  'repeat { if (file.exists(sentinel)) break',
  '  em <- as.numeric(Sys.time() - t_start, units = "mins")',
  '  writeLines(c("ASICS v25 (built-in) -- RUNNING",',
  '    paste(" Elapsed:", round(em,1), "min"),',
  '    paste(" Progress: ~", round(min(em/est_total_m*100,99)), "%"),',
  '    paste(" Updated:", format(Sys.time(), "%H:%M:%S"))), log_file); Sys.sleep(60) }'
), timer_script)
tryCatch(shell(paste("start /B Rscript", shQuote(gsub("/", "\\\\", timer_script))),
               wait = FALSE, intern = FALSE),
         error = function(e) message("  Background timer unavailable: ", e$message))

ASICS_results <- ASICS(
  spectra_obj, exclusion.areas = EXCLUSION_AREAS, pure.library = builtin_library,
  joint.align = TRUE, quantif.method = "both", max.shift = MAX_SHIFT,
  noise.thres = NOISE_THRES, add.noise = ADD_NOISE, mult.noise = MULT_NOISE,
  clean.thres = CLEAN_THRES, ncores = N_CORES, seed = 1234, verbose = TRUE)

elapsed <- round(as.numeric(Sys.time() - t_start, units = "mins"), 1)
writeLines("stop", sentinel); Sys.sleep(2)
writeLines(c("ASICS v25 (built-in) -- COMPLETE", paste(" Total:", elapsed, "min")), log_file)
tryCatch(unlink(timer_script), error = function(e) NULL)
message("Quantification complete in ", elapsed, " minutes."); print(ASICS_results)


# -- 6. Save results + reconstruction ------------------------------------------
saveRDS(ASICS_results, file.path(OUT_DIR, "ASICS_results.rds"))
quant <- getQuantification(ASICS_results)
write.csv(quant, file.path(OUT_DIR, "quantifications_all_samples.csv"))
write.csv(data.frame(Sample = colnames(quant), Metabolites_detected = colSums(quant > 0)),
          file.path(OUT_DIR, "detections_per_sample.csv"), row.names = FALSE)
write.csv(data.frame(Metabolite = rownames(quant), Samples_detected_in = rowSums(quant > 0)),
          file.path(OUT_DIR, "metabolite_prevalence.csv"), row.names = FALSE)
for (i in seq_len(min(3L, ncol(spectra_data)))) {
  pdf(file.path(OUT_DIR, paste0("reconstruction_sample", i, ".pdf")), width = 16, height = 5)
  plot(ASICS_results, idx = i); dev.off()
}
message("Results saved to: ", OUT_DIR)


# -- 7. Manual qHNMR panel (tiered, library-independent absolutes) ------------
message("\n=== Section 7: manual qHNMR panel ===")
ppm_axis <- as.numeric(rownames(spectra_raw))
tsp_integral <- vapply(colnames(spectra_raw), function(s)
  integrate_region(ppm_axis, spectra_raw[[s]], TSP_WINDOW[1], TSP_WINDOW[2]), numeric(1))
tsp_cv <- sd(tsp_integral, na.rm = TRUE) / mean(tsp_integral, na.rm = TRUE) * 100
message(sprintf("  TSP integral CV: %.1f%%", tsp_cv))
write.csv(data.frame(Sample = names(tsp_integral), TSP_integral = tsp_integral),
          file.path(OUT_DIR, "tsp_qc.csv"), row.names = FALSE)

# correction = 1 except glucose (anomeric fraction). ppm windows are STARTING
# points at ~pH 7 -- verify on qc_panel_windows.pdf and adjust before trusting.
reporters <- data.frame(
  metabolite = c(
    "Glucose","Lactate","Alanine","Acetate","Pyruvate",
    "Succinate","Formate","Phenylalanine","Tyrosine","Choline",
    "Betaine","Glycine","Methionine","Fumarate","Tryptophan",
    "Histidine","Valine","Leucine","Isoleucine","Creatine",
    "Citrate","Taurine","Asparagine","Aspartate","Arginine",
    "Lysine","Pyroglutamate","Ethanol","Glutamate","Glutamine"),
  ppm_from = c(
    5.21,1.31,1.46,1.9,2.35,
    2.39,8.44,7.31,6.87,3.19,
    3.25,3.54,2.12,6.5,7.71,
    7.86,0.97,0.94,0.92,3.02,
    2.5,3.4,2.94,2.78,3.22,
    3,4.16,1.16,2.32,2.43),
  ppm_to = c(
    5.25,1.35,1.5,1.94,2.39,
    2.42,8.48,7.44,6.92,3.21,
    3.27,3.56,2.15,6.54,7.75,
    7.92,1.01,0.97,0.95,3.04,
    2.7,3.44,2.99,2.84,3.26,
    3.04,4.19,1.19,2.37,2.47),
  n_protons = c(
    1,3,3,3,3,
    4,1,5,2,9,
    9,2,3,2,1,
    1,3,6,3,3,
    4,2,1,1,2,
    2,1,3,2,2),
  correction = c(
    ALPHA_FRACTION,1,1,1,1,
    1,1,1,1,1,
    1,1,1,1,1,
    1,1,1,1,1,
    1,1,1,1,1,
    1,1,1,1,1),
  confidence = c(
    "high","high","high","high","high",
    "high","high","high","high","high",
    "high","high","high","high","high",
    "moderate","moderate","moderate","moderate","moderate",
    "moderate","moderate","moderate","moderate","moderate",
    "moderate","moderate","moderate","upper_bound","upper_bound"),
  note = c(
    "alpha-anomeric H1 / alpha-fraction","CH3 doublet","CH3 doublet","CH3 singlet","CH3 singlet",
    "(CH2)2 singlet","CH singlet","aromatic ring (5H)","ring 3,5-H doublet","N(CH3)3 singlet",
    "N(CH3)3 singlet","CH2 singlet","S-CH3 singlet","vinyl singlet (if present)","indole H (clean aromatic)",
    "imidazole H2 (pH-mobile)","CH3 doublet (BCAA overlap)","2x CH3 (BCAA overlap)","delta-CH3 triplet (BCAA overlap)","N-CH3 singlet (Cr/creatinine/Lys overlap)",
    "AB quartet (pH-mobile)","N-CH2 triplet","beta-CH dd (Asp overlap)","beta-CH dd (Asn overlap)","delta-CH2 triplet (choline tail)",
    "epsilon-CH2 triplet (Cr overlap)","alpha-CH dd (Gln degradation marker)","CH3 triplet (carryover)","C4 CH2 (overlaps Gln)","C4 CH2 (overlaps Glu)"),
  stringsAsFactors = FALSE)

manual_tube_mM <- matrix(NA_real_, nrow(reporters), ncol(spectra_raw),
                         dimnames = list(reporters$metabolite, colnames(spectra_raw)))
for (s in colnames(spectra_raw)) {
  I_tsp <- tsp_integral[[s]]
  for (k in seq_len(nrow(reporters))) {
    I_reg <- integrate_region(ppm_axis, spectra_raw[[s]],
                              reporters$ppm_from[k], reporters$ppm_to[k])
    manual_tube_mM[k, s] <- (I_reg / I_tsp) * (TSP_PROTONS / reporters$n_protons[k]) *
      TSP_IN_TUBE_MM / reporters$correction[k]
  }
}
manual_real_mM <- sweep(manual_tube_mM, 2, dilution_df[colnames(manual_tube_mM)], "*")
write.csv(data.frame(metabolite = rownames(manual_tube_mM), confidence = reporters$confidence,
                     manual_tube_mM, check.names = FALSE),
          file.path(OUT_DIR, "manual_panel_tube_mM.csv"), row.names = FALSE)
write.csv(data.frame(metabolite = rownames(manual_real_mM), confidence = reporters$confidence,
                     manual_real_mM, check.names = FALSE),
          file.path(OUT_DIR, "manual_panel_real_mM.csv"), row.names = FALSE)
write.csv(reporters, file.path(OUT_DIR, "manual_panel_windows.csv"), row.names = FALSE)
message(sprintf("  Manual panel: %d metabolites (%d high / %d moderate / %d upper-bound).",
                nrow(reporters), sum(reporters$confidence=="high"),
                sum(reporters$confidence=="moderate"), sum(reporters$confidence=="upper_bound")))


# -- 8. Cross-calibrated full ASICS table to mM (Route B) ----------------------
anchor_patterns <- c(Lactate="lactate|lactic", Alanine="alanine", Acetate="acetate|acetic",
                     Pyruvate="pyruvate|pyruvic", Formate="formate|formic")
anchor_rows <- sapply(anchor_patterns, function(p) {
  h <- grep_anchored(p, rownames(quant)); if (length(h)) h[1] else NA_integer_ })
matched <- names(anchor_rows)[!is.na(anchor_rows)]
message("  Anchors matched: ", if (length(matched)) paste(matched, collapse=", ") else "NONE")
if (length(matched) >= 2) {
  cs <- intersect(colnames(quant), colnames(manual_tube_mM))
  sf <- vapply(cs, function(s) {
    num <- manual_tube_mM[matched, s]; den <- quant[anchor_rows[matched], s]
    ok <- is.finite(num) & is.finite(den) & den > 0
    if (sum(ok) >= 2) median(num[ok]/den[ok]) else NA_real_ }, numeric(1))
  abs_asics_real_mM <- sweep(sweep(quant[, cs, drop=FALSE], 2, sf, "*"), 2, dilution_df[cs], "*")
  write.csv(abs_asics_real_mM, file.path(OUT_DIR, "absolute_ASICS_crosscalibrated_real_mM.csv"))
  message("  Route B written (screening-grade; trust the manual panel where they overlap).")
}


# -- 9. QC / reproducibility ---------------------------------------------------
message("\n=== Section 9: QC / reproducibility ===")
param_names <- c("SPECTRA_DIR","OUT_DIR","FIELD_MHZ","NORM_METHOD","DO_ALIGN","NOISE_THRES",
                 "ADD_NOISE","MULT_NOISE","MAX_SHIFT","CLEAN_THRES","N_CORES","TSP_IS_D4",
                 "TSP_PCT_WV","TSP_IN_TUBE_MM","DEFAULT_DF","ALPHA_FRACTION")
params <- mget(param_names, ifnotfound = list(NA), inherits = TRUE)
writeLines(c("=== v25 BUILT-IN RUN MANIFEST ===",
             paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
             paste("R:", R.version.string), paste("ASICS:", as.character(utils::packageVersion("ASICS"))),
             paste("Spectra:", ncol(spectra_raw)),
             paste("Library:", length(getSampleName(ASICS_results)), "metabolites (built-in)"),
             "", "--- parameters ---",
             paste(names(params), vapply(params, function(x) paste(x, collapse=","), character(1)), sep=" = ")),
           file.path(OUT_DIR, "run_manifest.txt"))
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))

# tiered window QC plot (all panel windows + TSP on the median spectrum)
med_spec <- apply(as.matrix(spectra_raw), 1, median)
qc <- rbind(data.frame(metabolite="TSP", ppm_from=TSP_WINDOW[1], ppm_to=TSP_WINDOW[2],
                       confidence="ref"),
            reporters[, c("metabolite","ppm_from","ppm_to","confidence")])
n <- nrow(qc); ncg <- 4; nrg <- ceiling(n/ncg)
pdf(file.path(OUT_DIR, "qc_panel_windows.pdf"), width = 3*ncg, height = 2.4*nrg)
op <- par(mfrow = c(nrg, ncg), mar = c(3.5,3.5,2,1), mgp = c(2,0.6,0))
for (i in seq_len(n)) {
  w <- qc[i,]; pad <- (w$ppm_to - w$ppm_from)*3 + 0.04
  sel <- which(ppm_axis >= (w$ppm_from-pad) & ppm_axis <= (w$ppm_to+pad))
  if (length(sel) < 2) { plot.new(); title(main=paste(w$metabolite,"(no data)")); next }
  cm <- switch(w$confidence, high="black", moderate="darkorange3",
               upper_bound="red3", "blue3")
  plot(ppm_axis[sel], med_spec[sel], type="l", xlim=rev(range(ppm_axis[sel])),
       xlab="ppm", ylab="median int.", main=w$metabolite, col.main=cm)
  abline(v=c(w$ppm_from, w$ppm_to), col="red", lty=2)
}
par(op); dev.off()
message("  qc_panel_windows.pdf written (black=high, orange=moderate, red=upper-bound).")


# -- 10. Focused CHO metabolite extraction -------------------------------------
message("\n=== Section 10: Focused CHO metabolite extraction ===")
target_hits <- lapply(TARGET_PATTERNS, grep_anchored, names = rownames(quant))
thits <- do.call(rbind, lapply(names(target_hits), function(tg) {
  h <- target_hits[[tg]]; if (!length(h)) return(NULL)
  data.frame(target = tg, asics_row = rownames(quant)[h], stringsAsFactors = FALSE) }))
if (!is.null(thits) && nrow(thits)) {
  stats <- data.frame(asics_row = rownames(quant), prevalence = rowSums(quant > 0),
                      detection_pct = round(rowMeans(quant > 0)*100, 1),
                      mean_in_detected = apply(quant, 1, function(r){r<-r[r>0]; if(length(r)) mean(r) else 0}),
                      stringsAsFactors = FALSE)
  thits <- merge(thits, stats, by = "asics_row", all.x = TRUE)
  write.csv(thits, file.path(OUT_DIR, "CHO_target_matches_all.csv"), row.names = FALSE)
  best <- thits %>% arrange(target, desc(prevalence), desc(mean_in_detected)) %>%
    group_by(target) %>% slice(1) %>% ungroup()
  tq <- quant[best$asics_row, , drop = FALSE]; rownames(tq) <- best$target
  write.csv(tq, file.path(OUT_DIR, "CHO_target_metabolites_relative_intensity.csv"))
  write.csv(best, file.path(OUT_DIR, "CHO_target_metabolite_summary.csv"), row.names = FALSE)
  message("  Matched ", nrow(best), " of ", length(TARGET_PATTERNS), " targets.")
  if (nrow(tq) >= 2 && ncol(tq) >= 2) {
    pdf(file.path(OUT_DIR, "CHO_target_heatmap.pdf"), width = 10, height = 10)
    pheatmap(log10(tq + 1e-6), scale = "row",
             main = "CHO target metabolites (relative ASICS intensity)"); dev.off()
  }
}

message("\n=== Pipeline v25 (built-in, consolidated) complete. Results: ", OUT_DIR, " ===")
message("Report from: manual_panel_real_mM.csv (defensible). Screen from: ",
        "absolute_ASICS_crosscalibrated_real_mM.csv. Verify windows: qc_panel_windows.pdf.")