check_testing_main_db_top1 <- function(
    testing_main_db,
    organism = "hsa",
    kegg_lib = "current",
    method = "hyperg",
    background = "rbc",
    id_col = "kegg_id",
    run_on = c("sample", "tolerances"),
    out_csv = "metaboanalyst_top1_checks.csv",
    out_rds = "metaboanalyst_top1_checks.rds",
    progress_every = 25,
    verbose = TRUE
) {
  run_on <- match.arg(run_on, several.ok = TRUE)
  stopifnot(requireNamespace("MetaboAnalystR", quietly = TRUE))
  
  # ---- count total runs for progress ----
  total_runs <- 0L
  for (pathway_id in names(testing_main_db)) {
    frac_level <- testing_main_db[[pathway_id]]
    if (!is.list(frac_level)) next
    
    for (frac_key in names(frac_level)) {
      tol_level <- frac_level[[frac_key]]
      if (!is.list(tol_level)) next
      
      if ("sample" %in% run_on && is.data.frame(tol_level[["sample"]]) &&
          (id_col %in% colnames(tol_level[["sample"]]))) {
        total_runs <- total_runs + 1L
      }
      if ("tolerances" %in% run_on) {
        tol_keys <- setdiff(names(tol_level), "sample")
        for (tol_key in tol_keys) {
          df <- tol_level[[tol_key]]
          if (is.data.frame(df) && (id_col %in% colnames(df))) {
            total_runs <- total_runs + 1L
          }
        }
      }
    }
  }
  
  if (verbose) cat("Total MetaboAnalystR runs to perform:", total_runs, "\n")
  
  # ---- runner ----
  rows <- list()
  idx <- 0L
  done <- 0L
  t0 <- Sys.time()
  
  fmt_eta <- function(done, total, t0) {
    if (done == 0L) return(NA_character_)
    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    rate <- done / elapsed
    remaining <- (total - done) / rate
    sprintf("%0.1fs elapsed, ETA %0.1fs", elapsed, remaining)
  }
  
  run_one <- function(pathway_id, frac_key, compound_ids, list_type, tol_key = NA_character_) {
    compound_ids <- unique(compound_ids)
    compound_ids <- compound_ids[!is.na(compound_ids) & compound_ids != ""]
    n_input <- length(compound_ids)
    
    top_pathway <- NA_character_
    top_raw_p <- NA_real_
    top_fdr <- NA_real_
    ok_run <- FALSE
    err <- NA_character_
    
    if (n_input == 0) {
      return(data.frame(
        pathway = pathway_id,
        frac = frac_key,
        list_type = list_type,
        tolerance_ppm = tol_key,
        n_input = 0L,
        top_pathway = NA_character_,
        correct_top1 = FALSE,
        top_raw_p = NA_real_,
        top_fdr = NA_real_,
        ok_run = FALSE,
        error = "Empty compound list",
        stringsAsFactors = FALSE
      ))
    }
    
    tryCatch({
      nSet <- MetaboAnalystR::InitDataObjects("conc", "pathora", FALSE, default.dpi = 300)
      nSet <- MetaboAnalystR::Setup.MapData(nSet, compound_ids)
      nSet <- MetaboAnalystR::CrossReferencing(nSet, "kegg")
      nSet <- MetaboAnalystR::CreateMappingResultTable(nSet)
      
      nSet <- MetaboAnalystR::SetKEGG.PathLib(nSet, organism, kegg_lib)
      nSet <- MetaboAnalystR::SetMetabolomeFilter(nSet, FALSE)
      nSet <- MetaboAnalystR::CalculateOraScore(nSet, background, method)
      
      ora <- nSet$analSet$ora.mat
      
      if (!is.null(ora) && nrow(ora) > 0) {
        top_pathway <- rownames(ora)[1]
        if ("Raw p" %in% colnames(ora)) top_raw_p <- as.numeric(ora[1, "Raw p"])
        if ("FDR"   %in% colnames(ora)) top_fdr   <- as.numeric(ora[1, "FDR"])
      }
      
      ok_run <- TRUE
    }, error = function(e) {
      err <<- conditionMessage(e)
    })
    
    data.frame(
      pathway = pathway_id,
      frac = frac_key,
      list_type = list_type,
      tolerance_ppm = tol_key,
      n_input = n_input,
      top_pathway = top_pathway,
      correct_top1 = isTRUE(ok_run) && !is.na(top_pathway) && identical(top_pathway, pathway_id),
      top_raw_p = top_raw_p,
      top_fdr = top_fdr,
      ok_run = ok_run,
      error = err,
      stringsAsFactors = FALSE
    )
  }
  
  # ---- main loop w/ progress ----
  for (pathway_id in names(testing_main_db)) {
    frac_level <- testing_main_db[[pathway_id]]
    if (!is.list(frac_level)) next
    
    for (frac_key in names(frac_level)) {
      tol_level <- frac_level[[frac_key]]
      if (!is.list(tol_level)) next
      
      if ("sample" %in% run_on && !is.null(tol_level[["sample"]])) {
        df <- tol_level[["sample"]]
        if (is.data.frame(df) && (id_col %in% colnames(df))) {
          done <- done + 1L
          idx <- idx + 1L
          rows[[idx]] <- run_one(pathway_id, frac_key, df[[id_col]], "sample", NA_character_)
          
          if (verbose && (done %% progress_every == 0L || done == total_runs)) {
            cat(sprintf("[%d/%d] %s\n", done, total_runs, fmt_eta(done, total_runs, t0)))
          }
        }
      }
      
      if ("tolerances" %in% run_on) {
        tol_keys <- setdiff(names(tol_level), "sample")
        for (tol_key in tol_keys) {
          df <- tol_level[[tol_key]]
          if (!is.data.frame(df) || !(id_col %in% colnames(df))) next
          
          done <- done + 1L
          idx <- idx + 1L
          rows[[idx]] <- run_one(pathway_id, frac_key, df[[id_col]], "tol", tol_key)
          
          if (verbose && (done %% progress_every == 0L || done == total_runs)) {
            cat(sprintf("[%d/%d] %s\n", done, total_runs, fmt_eta(done, total_runs, t0)))
          }
        }
      }
    }
  }
  
  res <- if (length(rows)) do.call(rbind, rows) else data.frame()
  
  utils::write.csv(res, out_csv, row.names = FALSE)
  saveRDS(res, out_rds)
  
  if (verbose) cat("Done. Wrote:", out_csv, "and", out_rds, "\n")
  res
}

res <- check_testing_main_db_top1(
  testing_main_db,
  progress_every = 20,   # print every 20 runs
  verbose = TRUE
)

