# A project-specific update for camprotR:plot_quant(method='density')
# with appropriate aesthetics for project
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
