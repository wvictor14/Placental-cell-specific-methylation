# Taken from https://github.com/hhhh5/ewastools/blob/master/R/check_snp_agreement.R on May 21, 2019
# Edited to return all agreement scores

check_snp_agreement2 = function(genotypes,donor_ids,sample_ids, threshold = F, 
                                return = 'within-donor'){
  
  J = ncol(genotypes$snps)
  weights = 1-genotypes$outliers
  gamma = genotypes$gamma
  
  stopifnot(anyDuplicated(sample_ids)==0)
  stopifnot(length(sample_ids) == J)
  
  sample_ids = as.character(sample_ids)
  
  conflicts = cbind(data.table::CJ(donor_ids,donor_ids,sorted=FALSE),
                    data.table::CJ(sample_ids,sample_ids,sorted=FALSE))
  data.table::setcolorder(conflicts,c(1,3,2,4))
  data.table::setnames(conflicts,1:4,c('donor1','sample1','donor2','sample2'))
  
  d = matrix(NA_real_,nrow=J,ncol=J)
  for(j in 1:J){
    tmp = 
      (weights * gamma[[1]]) * (weights[,j] * gamma[[1]][,j]) +
      (weights * gamma[[2]]) * (weights[,j] * gamma[[2]][,j]) +
      (weights * gamma[[3]]) * (weights[,j] * gamma[[3]][,j])
    
    d[,j] = colSums(tmp,na.rm=TRUE) / colSums(weights*weights[,j],na.rm=TRUE)
    
  }
  
  # first drop the duplicate entries and the diagonal
  ## EDIT: do not drop the diagonal
  d[upper.tri(d,diag=F)] = NA_real_
  conflicts$agreement = as.numeric(d)
  rm(tmp,d)
  conflicts = conflicts[!is.na(agreement)]
  
  # drop cases of high agreement
  if (threshold) {
    conflicts = conflicts[!(donor1!=donor2 & agreement<0.90)]
    conflicts = conflicts[!(donor1==donor2 & agreement>0.90)]
    
    if(nrow(conflicts)==0) return(NULL)
    
    ### find the weakly connected components of the graph (consider conflicts as edges in a graph)
    e = rep(NA,times=2*nrow(conflicts))
    e[c(TRUE,FALSE)] = conflicts$sample1
    e[c(FALSE,TRUE)] = conflicts$sample2
    
    g = igraph::make_graph(edges=e,directed=FALSE)
    g = igraph::components(g,mode="weak")$membership
    
    conflicts$group = g[conflicts$sample1]
    conflicts = split(conflicts,by='group',keep.by=FALSE)
  }
  
  if (return == 'within-donor') {
    conflicts <- conflicts %>% dplyr::filter(donor1 == donor2)
  }
  
  return(conflicts)
    
}
