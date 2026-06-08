####################
# Custom Functions #
####################


#################
# Plotting funs #
#################

## Make a pdf -------------------------
makepdf <- function(x, name, is_pheatmap=F, ncol=3, nrow=2, figNum=NULL, width=20, height=12) {
  #figs per page
  if (is.null(figNum)) figNum = ncol*nrow
  if (length(x) == 0) return(NULL)
  # draw plots
  if(is_pheatmap==T){
    ## pheatmap
    require(grid)
    pdf(name, width=width, height=height)
    for(i in seq(1,length(x))){
      grid::grid.newpage()
      grid::grid.draw(x[[i]]$gtable)
    }
    dev.off() %>% invisible #close pdf
  }else{
    ## normal ggplots
    require(gridExtra)
    pdf(name, width=width, height=height)
    for (i in seq(1, length(x), by = figNum)) {
      j  <- min(i + figNum - 1, length(x))
      do.call(gridExtra::grid.arrange, c(x[i:j], ncol = ncol, nrow = nrow)) 
    }
    dev.off() %>% invisible #close pdf
  }
}

## get ready-plot PCA data  -------------------------
get_PCAdata = function(mat, transformation=c('DESeq2::vst', 'base::scale'), 
                       center=T, scale=T, ## only relevant when transformation base::scale is used
                       ntop=500, PCs=1:10, obsName='patient_ID'){
  if (transformation=='DESeq2::vst'){
    require(DESeq2)
    mat=DESeq2::vst(mat, blind=T, nsub=nrow(mat)) 
  } else{
    mat=base::scale(mat, center=center, scale=scale)
  }
  ### select ntop features by var
  rv=rowVars(mat)
  ordered_features=order(rv, decreasing=T)
  mat=head(mat[ordered_features,], ntop)
  ### PCA
  require(tidyverse)
  #### compute PCA
  res_pca=prcomp(t(mat),center=T,scale.=F)
  #### get % var per PC
  percentVar = head(round((res_pca$sdev^2 / sum(res_pca$sdev^2))*100,2),length(PCs))
  percentVar = data.frame(PC_name=paste0('PC',PCs), percentVar=percentVar)
  df=data.frame(res_pca$x[,PCs]) %>% tibble::rownames_to_column(var=obsName) %>%
    tidyr::pivot_longer(cols=-any_of(obsName), names_to = 'PC_name', values_to='PC_coord') %>%
    mutate(percentVar=percentVar[match(PC_name,percentVar$PC_name),]$percentVar)
  
  return(df)
}


###################
# Plot aesthetics #
###################

## plot aesthetics -------------------------
p_theme = theme_bw() + 
  theme(plot.title = element_text(hjust=0.5, size=11, face='bold'), 
        plot.subtitle = element_text(hjust=0,size=9, face='plain'),
        axis.title = element_text(size=9), axis.text = element_text(size=8), 
        legend.title = element_text(size=9), legend.text = element_text(size=8),
        panel.grid = element_blank(), panel.border = element_rect(linewidth=1))







