run_pathway_min_compounds <- function(pathway_id,
                                      all_pathway_lists,
                                      fracs = seq(1, 0.1, by = -0.1),
                                      n_reps = 1,
                                      seed = 1,
                                      organism = "hsa",
                                      kegg_lib = "current",
                                      method = "hyperg",
                                      background = "rbc") {
  set.seed(seed)
  
  # Compounds available for this pathway
  comp_ids <- all_pathway_lists[[pathway_id]]$kegg_id
  comp_ids <- unique(comp_ids[!is.na(comp_ids) & comp_ids != ""])
  n_total <- length(comp_ids)
  
  if (n_total == 0) return(data.frame())
  
  out <- vector("list", length(fracs) * n_reps)
  idx <- 1
  
  for (frac in fracs) {
    size <- max(1L, floor(n_total * frac))
    
    for (rep in seq_len(n_reps)) {
      sampled <- sample(comp_ids, size = size, replace = FALSE)
      
      # Run MetaboAnalystR ORA using KEGG IDs as input
      nSet <- InitDataObjects("conc", "pathora", FALSE, default.dpi = 300)
      nSet <- Setup.MapData(nSet, sampled)
      nSet <- CrossReferencing(nSet, "kegg")
      nSet <- CreateMappingResultTable(nSet)
      
      nSet <- SetKEGG.PathLib(nSet, organism, kegg_lib)
      nSet <- SetMetabolomeFilter(nSet, FALSE)
      nSet <- CalculateOraScore(nSet, background, method)
      
      ora <- nSet$analSet$ora.mat
      
      top_pathway <- if (!is.null(ora) && nrow(ora) > 0) rownames(ora)[1] else NA_character_
      hit <- identical(top_pathway, pathway_id)
      
      # record some metrics if present
      top_raw_p <- if (!is.null(ora) && nrow(ora) > 0) ora[1, "Raw p"] else NA_real_
      top_fdr   <- if (!is.null(ora) && nrow(ora) > 0) ora[1, "FDR"] else NA_real_
      
      out[[idx]] <- data.frame(
        pathway = pathway_id,
        n_total = n_total,
        frac = frac,
        n_sample = size,
        rep = rep,
        top_pathway = top_pathway,
        correct_top1 = hit,
        top_raw_p = top_raw_p,
        top_fdr = top_fdr,
        stringsAsFactors = FALSE
      )
      idx <- idx + 1
    }
  }
  
  do.call(rbind, out)
}

# Example settings:
fracs <- seq(1, 0.1, by = -0.1)
n_reps <- 1

all_pids <- names(all_pathway_lists)
ran_pids <- sample(all_pids, size = 5, replace = FALSE)

all_results <- do.call(rbind, lapply(ran_pids, function(pid) {
  cat("Running:", pid, "\n")
  run_pathway_min_compounds(pid, all_pathway_lists, fracs = fracs, n_reps = n_reps, seed = 123)
}))

# Save for later plotting
saveRDS(all_results, file = "pathway_subsampling_results.rds")


# 
# for (pathway in names(all_pathway_lists)){
#   print(pathway)
#   i = 1
#   
#   while (i > 0) {
#     new_compounds <- all_pathway_lists[[pathway]]$kegg_id
#     new_compounds <- sample(new_compounds, size = NROW(new_compounds)*i, replace = FALSE)
#     nSet<-InitDataObjects("conc", "pathora", FALSE, default.dpi = 300)
#     nSet<-Setup.MapData(nSet, new_compounds);
#     nSet<-CrossReferencing(nSet, "kegg");
#     nSet<-CreateMappingResultTable(nSet);
#     
#     if (rownames(nSet$analSet$ora.mat)[1] == pathway){
#       # TRUE
#     }
#     
#     
#     i=i-0.1
#   }
# }