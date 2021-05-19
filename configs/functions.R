plot_below_notch_per_prot <-  function (notch_per_protein) {p_notch_per_protein <- notch_per_protein %>% group_by(.data$n_below, 
                                                                                                                  sample = remove_x(.data$sample)) %>% tally() %>% ggplot(aes(sample, n, fill =cutr::cutf2(.data$n_below, cuts = c(0, 1, 5,  10, 15, 20, 30, 50)))) + geom_bar(stat = "identity",  position = "fill") + scale_fill_manual(name = "# PSMs below notch",  values = c("grey", get_cat_palette(7))) + xlab("Tag") + scale_y_continuous(name = "Fraction of proteins", expand = c(0,0)) + theme_camprot(base_size = 12, border = FALSE)
return(p_notch_per_protein)
}

makeMSNSet <- function(obj, samples_inf, ab_col_ix=3, level="peptide", quant_name="Abundance"){
  # make dataframes for MSnset object
  rownames(obj) <- seq(1, length(obj[,1]))
  
  meta_columns <- colnames(obj)
  meta_columns <- meta_columns[grep("Found.*", meta_columns, invert=TRUE)]
  meta_columns <- meta_columns[grep(sprintf('%s.*', quant_name), meta_columns, invert=TRUE)]
  
  if(level=="PSM"){
    abundance_columns <- colnames(obj)[grep(sprintf('%s.*', quant_name), colnames(obj))]
    renamed_abundance_columns <- sapply(strsplit(abundance_columns, "\\."), "[[", ab_col_ix)
  }
  else if(level=="peptide"){
    abundance_columns <- colnames(obj)[grep(sprintf('%s.*', quant_name), colnames(obj))]
    abundance_columns <- abundance_columns[grep(sprintf('%s.Count', quant_name), abundance_columns, invert=TRUE)]
    renamed_abundance_columns <- sapply(strsplit(abundance_columns, "\\."), "[[", ab_col_ix)
  }
  else{
    stop("level must be PSM or peptide")
  }
  
  exprsCsv <- obj[,abundance_columns]
  colnames(exprsCsv) <- renamed_abundance_columns
  
  exprsCsv[exprsCsv==""] <- NA
  
  fdataCsv <- obj[,meta_columns]
  
  pdataCsv <- read.table(samples_inf, sep="\t", header=T, row.names = 1, colClasses="character")
  
  exprsCsv <- exprsCsv[,rownames(pdataCsv)]
  
  res <- MSnSet(as.matrix(sapply(exprsCsv, as.numeric)), fdataCsv, pdataCsv)
  
  summariseMissing(res)
  
  cat(sprintf("\n%s peptides do not have any quantification values\n\n", sum(rowSums(is.na(exprs(res)))==ncol(exprs(res)))))
  res <- res[rowSums(is.na(exprs(res)))!=ncol(exprs(res)),] # exclude peptides without any quantification
  return(res)
}

plot_quant_density_tg <- function(obj){
  
  e_data <- obj %>% exprs() %>% data.frame()
  e_data[e_data == ""] <- NA
  e_data <- e_data %>%
    gather(key = "sample", value = "intensity") %>% 
    mutate(sample = factor(remove_x(sample), levels = remove_x(colnames(e_data)))) %>%
    merge(pData(obj), by.x='sample', by.y='row.names') %>%
    mutate(fraction=factor(fraction, levels=unique(pData(obj)$fraction)))
  
  p <- ggplot(e_data) +
    theme_camprot() +
    geom_density(aes(.data$intensity, col = fraction, linetype=condition)) + 
    xlab("Feature intensity") +
    ylab("Density")
  
  return(p)
}



# document and move to camprotR?
detect_tmt_psm_outliers <- function(obj, master_prot_col='Master.Protein.Accessions'){
  message('Identifying outlier PSMs')
  
  camprotR:::message_parse(fData(obj), master_prot_col, 'Input')
  obj <- obj %>% filterNA(pNA=0) %>% normalise('sum')
  camprotR:::message_parse(fData(obj), master_prot_col, 'Filtering PSMs with NAs')
  
  retain_proteins <- obj %>% fData() %>%
    group_by(Master.Protein.Accessions) %>%
    tally() %>%
    filter(n>=2) %>%
    pull(Master.Protein.Accessions)
  
  obj <- obj[fData(obj)[[master_prot_col]] %in% retain_proteins]
  camprotR:::message_parse(fData(obj), master_prot_col, 'Filtering Proteins with <2 PSMs')
  
  prot_ids <- unique(fData(obj)[[master_prot_col]])
  distances <- vector('list', length(prot_ids))
  names(distances) <- prot_ids
  
  pb <-txtProgressBar(min = 0, max = length(prot_ids), style = 1)
  
  for(prot_ix in seq_along(prot_ids)){
    prot <- prot_ids[prot_ix]
    prot_psms <- obj[fData(obj)[[master_prot_col]]==prot,]
    
    distances[[prot]] <- prot_psms %>%
      exprs() %>%
      dist() %>%
      as.matrix() %>% data.frame() %>%
      tibble::rownames_to_column('PSM2') %>%
      pivot_longer(-PSM2, names_to='PSM1', values_to='distance') %>%
      mutate(PSM1=remove_x(PSM1)) %>%
      filter(PSM1!=PSM2) %>%
      group_by(PSM1) %>%
      summarise(median_distance=median(distance), .groups='drop') %>%
      mutate(protein=prot)
    
    setTxtProgressBar(pb, prot_ix)
  }
  
  distances <- do.call('rbind', distances)
  
  distances_annotated <- distances %>% merge(fData(obj), by.x='PSM1', by.y='row.names')
  
  return(distances_annotated)
}



detect_tmt_psm_outliers <- function(obj, master_prot_col='Master.Protein.Accessions'){
  message('Identifying outlier PSMs')
  
  camprotR:::message_parse(fData(obj), master_prot_col, 'Input')
  obj <- obj %>% filterNA(pNA=0) %>% normalise('sum')
  camprotR:::message_parse(fData(obj), master_prot_col, 'Filtering PSMs with NAs')
  
  retain_proteins <- obj %>% fData() %>%
    group_by(Master.Protein.Accessions) %>%
    tally() %>%
    filter(n>=2) %>%
    pull(Master.Protein.Accessions)
  
  obj <- obj[fData(obj)[[master_prot_col]] %in% retain_proteins]
  camprotR:::message_parse(fData(obj), master_prot_col, 'Filtering Proteins with <2 PSMs')
  
  prot_ids <- unique(fData(obj)[[master_prot_col]])
  distances <- vector('list', length(prot_ids))
  names(distances) <- prot_ids
  
  pb <-txtProgressBar(min = 0, max = length(prot_ids), style = 1)
  
  for(prot_ix in seq_along(prot_ids)){
    prot <- prot_ids[prot_ix]
    prot_psms <- obj[fData(obj)[[master_prot_col]]==prot,]
    
    distances[[prot]] <- prot_psms %>%
      exprs() %>%
      dist() %>%
      as.matrix() %>% data.frame() %>%
      tibble::rownames_to_column('PSM2') %>%
      pivot_longer(-PSM2, names_to='PSM1', values_to='distance') %>%
      mutate(PSM1=remove_x(PSM1)) %>%
      filter(PSM1!=PSM2) %>%
      group_by(PSM1) %>%
      summarise(median_distance=median(distance), .groups='drop') %>%
      mutate(protein=prot)
    
    setTxtProgressBar(pb, prot_ix)
  }
  
  distances <- do.call('rbind', distances)
  
  distances_annotated <- distances %>% merge(fData(obj), by.x='PSM1', by.y='row.names')
  
  return(distances_annotated)
}
