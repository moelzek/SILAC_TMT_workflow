#library(tcltk)

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
