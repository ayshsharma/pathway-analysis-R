.kegg_cache <- new.env(parent = emptyenv())
cache_key <- function(mass, tol) paste0(format(mass, digits = 12), "_", tol)



FindMolEqualWeight <- function(query_weight, tolerance = 1000){
  key <- cache_key(query_weight, tolerance)
  if (exists(key, envir = .kegg_cache, inherits = FALSE)) return(get(key, envir = .kegg_cache))
  
  tolerance_da = (query_weight/1e6)*tolerance # Calculates tolerance in daltons
  print(tolerance_da)
  lower_limit = query_weight - tolerance_da
  upper_limit = query_weight + tolerance_da
  similar <- keggFind("compound", floor(lower_limit):ceiling(upper_limit), "exact_mass")
  
  # Defines dataframe for result
  df <- data.frame( 
    kegg_id = names(similar),
    mol_mass = as.numeric(as.vector(similar)),
    stringsAsFactors = FALSE
  )
  
  df$delta_da <- abs(df$mol_mass - query_weight)
  
  df <- df[df$delta_da <= tolerance_da, , drop = FALSE]

  df <- df[order(df$delta_da), , drop = FALSE]
  # if (nrow(df) > max_return) df <- df[seq_len(max_return), , drop = FALSE]

  df$query_weight <- rep(query_weight, nrow(df))
  df$tolerance_ppm <- rep(tolerance, nrow(df))
  df$tolerance_da <- rep(tolerance_da, nrow(df))
  rownames(df) <- NULL
  
  assign(key, df, envir = .kegg_cache)
  return(df)
}





# Removes the compound that is already in the pathway from the list
getAdditionalMol <- function(query_kegg, query_weight, tolerance = 1000, max_return = 5){
  df <- FindMolEqualWeight(query_weight, tolerance = tolerance)
  df <- df[df$kegg_id != query_kegg, , drop = FALSE]
  
  return(df)
}






# Creating a dataframe with all combinations of tolerances and percentage of pathway compounds
fracs <- list(0.01, 0.05, 0.1, 0.2, 0.5, 0.75, 1)
tolerances <- list(1, 3, 5, 10, 30, 60, 100)









sampleCompounds <- function(pathway_id, frac, pathway_db = all_pathway_lists) {
  df <- pathway_db[[pathway_id]][, c("kegg_id", "mono_mass")]
  
  # remove missing/blank ids
  df <- df[!is.na(df$kegg_id) & df$kegg_id != "", , drop = FALSE]
  
  # ensure unique compounds (keep first row per kegg_id)
  df <- df[!duplicated(df$kegg_id), , drop = FALSE]
  
  n_total <- nrow(df)
  if (n_total == 0) return(df)
  
  size <- max(1L, floor(n_total * frac))
  df <- df[sample.int(n_total, size), , drop = FALSE]
  
  data.frame(
    kegg_id = df$kegg_id,
    mol_mass = as.numeric(df$mono_mass),
    stringsAsFactors = FALSE
  )
}

sampleCompounds("hsa00260", 0.5)




testing_main_db <- list()

fracs <- c(1, 0.75, 0.5, 0.2, 0.1, 0.05, 0.01)
tolerances <- c(1, 3, 5, 10, 30, 60, 100)
for (pathway in names(all_pathway_lists)) {
  testing_main_db[[pathway]] <- list()
  
  for (frac in fracs) {
    frac_key <- as.character(frac)
    testing_main_db[[pathway]][[frac_key]] <- list()
    
    tempSampleCompounds <- sampleCompounds(pathway, frac = frac)
    testing_main_db[[pathway]][[frac_key]][["sample"]] <- tempSampleCompounds
    
    for (tol in tolerances) {
      tol_key <- as.character(tol)
      
      # start with the sampled compounds (or use [0,] if you want only additions)
      out <- tempSampleCompounds
      
      for (i in seq_len(nrow(tempSampleCompounds))) {
        pathway_compound <- tempSampleCompounds$kegg_id[i]
        comp_mass <- tempSampleCompounds$mol_mass[i]
        if (!is.finite(comp_mass)) next
        
        add_mol_max <- getAdditionalMol(pathway_compound, comp_mass, tolerance = 100)
        tol_da <- (comp_mass / 1e6) * tol
        
        add_mol <- add_mol_max[add_mol_max$delta_da <= tol_da, , drop = FALSE]
        add_mol$delta_da <- NULL
        add_mol$query_weight <- NULL
        add_mol$tolerance_ppm <- NULL
        add_mol$tolerance_da <- NULL
        
        out <- rbind(out, add_mol)
      }
      
      out <- out[!duplicated(out$kegg_id), , drop = FALSE]
      testing_main_db[[pathway]][[frac_key]][[tol_key]] <- out
    }
  }
}

saveRDS(testing_main_db, file="Main_Testing.Rda")


# Get a graph of number of molecules against T for a set of 50 random query molecules


# Adding an in_pathway column
for (pathway in names(testing_main_db)) {
  for (frac in names(testing_main_db[[pathway]])){
    for (tol in names(testing_main_db[[pathway]][[as.character(frac)]])){
      
      df <- testing_main_db[[pathway]][[frac]][[tol]]
      
      df$in_pathway <- df$kegg_id %in% testing_main_db[[pathway]][[frac]]$sample$kegg_id
      
      testing_main_db[[pathway]][[frac]][[tol]] <- df
    }
  }
}


sampled_pathways = sample(names(testing_main_db), 5, replace = FALSE)
sampled_tol = as.character(c(1, 5))
sampled_fracs = as.character(c(1, 0.5, 0.01))

all_pathway_lists <- vector("list", length(sampled_pathways))
names(all_pathway_lists) <- sampled_pathways

for (pathway in sampled_pathways) {
  for (frac in sampled_fracs){
    for (tol in sampled_tol){
      pathway_compound_temp <- as.character(testing_main_db[[pathway]][[frac]][[tol]]$kegg_id)
      print(pathway_compound_temp)
      oSet<-InitDataObjects("conc", "pathora", FALSE, default.dpi = 300)
      oSet<-Setup.MapData(oSet, pathway_compound_temp);
      oSet<-CrossReferencing(oSet, "kegg");
      oSet<-CreateMappingResultTable(oSet);

      oSet<-SetKEGG.PathLib(oSet, "hsa", "current")
      oSet<-SetMetabolomeFilter(oSet, FALSE);
      oSet<-CalculateOraScore(oSet, "rbc", "hyperg")
      
      # res <- tryCatch({
      #   oSet <- InitDataObjects("conc", "pathora", FALSE, default.dpi = 300)
      #   oSet <- Setup.MapData(oSet, pathway_compound_temp)
      #   oSet <- CrossReferencing(oSet, "kegg")
      #   oSet <- CreateMappingResultTable(oSet)
      #   oSet <- SetKEGG.PathLib(oSet, "hsa", "current")
      #   oSet <- SetMetabolomeFilter(oSet, FALSE)
      #   oSet <- CalculateOraScore(oSet, "rbc", "hyperg")
      #   oSet$analSet$ora.mat
      # }, error = function(e) e)
      # 
      # if (inherits(res, "error")) {
      #   all_pathway_lists[[pathway]][[frac]][[tol]]$error <- conditionMessage(res)
      # } else {
      #   all_pathway_lists[[pathway]][[frac]][[tol]]$pathway_anal <- res
      # }
      
      if (is.null(all_pathway_lists[[pathway]]) || !is.list(all_pathway_lists[[pathway]]))
        all_pathway_lists[[pathway]] <- list()
      
      if (is.null(all_pathway_lists[[pathway]][[frac]]) || !is.list(all_pathway_lists[[pathway]][[frac]]))
        all_pathway_lists[[pathway]][[frac]] <- list()
      
      if (is.null(all_pathway_lists[[pathway]][[frac]][[tol]]) || !is.list(all_pathway_lists[[pathway]][[frac]][[tol]]))
        all_pathway_lists[[pathway]][[frac]][[tol]] <- list()
      
      all_pathway_lists[[pathway]][[frac]][[tol]]$pathway_anal <- oSet$analSet$ora.mat
    }
  }
} # The issue that I think that I am running into is that once there are less than 3 compounds left, metaboanalyst throws an error
