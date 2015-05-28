

#' get_dist_between_calls
#'
#' @param ranges GRanges object 
#' @param max_dist to consider
#' @return returns tabulated distance between calls
#' @seealso \code{\link{distRanges}} 
#' @export
get_dist_between_calls<- function(ranges, max_dist = 1000){
  sep = start(ranges)[2:(length(ranges))] - end(ranges)[1:(length(ranges)-1)] 
  chrom_same = which(seqnames(ranges)[1:(length(ranges)-1)] == seqnames(ranges)[2:(length(ranges))])
  close = intersect(which(sep < max_dist),chrom_same)
  return(tabulate(sep[close]))
}


#' distRanges
#'
#' @param ranges1 GRanges object 
#' @param ranges2 GRanges object
#' @param strand either 1,2, or 0 if dist is relative to strand of ranges1, ranges2, or neither
#' @return returns either subset of nucs that are +1 or annotated version of sites
#' @seealso \code{\link{get_m1_nucs}}  \code{\link{get_p1_nucs}}  \code{\link{get_dist_between_calls}} 
#' @export
distRanges <- function(ranges1,ranges2, strand = 0){
  ranges1_mod = ranges1
  strand(ranges1_mod)="*"
  neighbor = nearest(ranges1_mod,ranges2)  
  if (strand == 2){
    dists = rep(NA,length(ranges1))
    plus = which(strand(ranges2)[neighbor]!="-")
    minus = which(strand(ranges2)[neighbor]=="-")
    dists[minus] = ifelse(overlapsAny(ranges1[minus],ranges2[neighbor[minus]]),
                          0,
                          ifelse(start(ranges2[neighbor[minus]])>start(ranges1[minus]),
                                 start(ranges2)[neighbor[minus]] - end(ranges1)[minus],
                                 end(ranges2)[neighbor[minus]] - start(ranges1)[minus]))
    dists[plus] = ifelse(overlapsAny(ranges1[plus],ranges2[neighbor[plus]]),
                         0,
                         ifelse(start(ranges2[neighbor[plus]])>start(ranges1[plus]),
                                end(ranges1)[plus] - start(ranges2)[neighbor[plus]],
                                start(ranges1)[plus] - end(ranges2)[neighbor[plus]]))                    
  } 
  else if (strand == 0){
    dists = ifelse(overlapsAny(ranges1,ranges2[neighbor]),
                         0,
                         ifelse(start(ranges2[neighbor])>start(ranges1),
                                end(ranges1) - start(ranges2)[neighbor],
                                start(ranges1) - end(ranges2)[neighbor]))                    
  }
  else if (strand ==1){
    dists = rep(NA,length(ranges1))
    plus = which(strand(ranges1)!="-")
    minus = which(strand(ranges1)=="-")
    dists[plus] = ifelse(overlapsAny(ranges1[plus],ranges2[neighbor[plus]]),
                          0,
                          ifelse(start(ranges2[neighbor[plus]])>start(ranges1[plus]),
                                 start(ranges2)[neighbor[plus]] - end(ranges1)[neighbor[plus]],
                                 end(ranges2)[neighbor[plus]] - start(ranges1)[neighbor[plus]]))
    dists[minus] = ifelse(overlapsAny(ranges1[minus],ranges2[neighbor[minus]]),
                         0,
                         ifelse(start(ranges2[neighbor[minus]])>start(ranges1[minus]),
                                end(ranges1)[minus] - start(ranges2)[neighbor[minus]],
                                start(ranges1)[minus] - end(ranges2)[neighbor[minus]]))                    
    return(dists)
  }
  else{
    stop("strand must equal 0, 1, or 2.")
  }
  return(dists)
}





#' get_p1_nucs
#'
#' @param nuc.ranges GRanges object with nucleosome positions
#' @param sites GRanges object with sites to get +1 nucs from.
#' @param annotate_sites what to return-- sites with +1 nuc annotation or +1 nucs?
#' @param max_dist what is maximum distance to consider
#' @return returns either subset of nucs that are +1 or annotated version of sites
#' @seealso \code{\link{get_m1_nucs}} 
#' @export
get_p1_nuc<-function(nucs.ranges,sites, annotate_sites = F, max_dist = 250){
  plus=which(strand(sites)!="-")
  minus=which(strand(sites)=="-")
  
  p1 = precede(sites,nucs.ranges)
  p1_nna = which(!is.na(p1))
  p1_nna = p1_nna[which(seqnames(nucs.ranges[p1[p1_nna]])==seqnames(sites[p1_nna]))]
  if (annotate_sites){
    sites$p1 = NA
    sites$p1[p1_nna] = start(nucs.ranges[p1[p1_nna]])
    sites$p1_ix = NA
    sites$p1_ix[p1_nna] = p1[p1_nna]  
    sites$p1_dist = NA
    p1_dist = rep(NA,length(sites))
    p1_dist[intersect(p1_nna,plus)] = start(nucs.ranges[p1[intersect(p1_nna,plus)]])-start(sites[intersect(p1_nna,plus)])
    p1_dist[intersect(p1_nna,minus)] = start(sites[intersect(p1_nna,minus)])-start(nucs.ranges[p1[intersect(p1_nna,minus)]])
    close = which(p1_dist[p1_nna] <= max_dist)
    sites$p1_dist[p1_nna[close]]=p1_dist[p1_nna[close]]
    return(sites)
  } else{
    p1_dist = rep(NA,length(sites))
    p1_dist[intersect(p1_nna,plus)] = start(nucs.ranges[p1[intersect(p1_nna,plus)]])-start(sites[intersect(p1_nna,plus)])
    p1_dist[intersect(p1_nna,minus)] = start(sites[intersect(p1_nna,minus)])-start(nucs.ranges[p1[intersect(p1_nna,minus)]])
    close = which(p1_dist[p1_nna] <= max_dist)
    p1_nucs = nucs.ranges[p1[p1_nna[close]]]
    strand(p1_nucs)=strand(sites)[p1_nna[close]]
    p1_nucs$dist = p1_dist[p1_nna[close]]
    return(p1_nucs)
  }
}

#' get_m1_nucs
#'
#' @param nuc.ranges GRanges object with nucleosome positions
#' @param sites GRanges object with sites to get -1 nucs from.
#' @param annotate_sites what to return-- sites with -1 nuc annotation or +1 nucs?
#' @param max_dist what is maximum distance to consider
#' @return returns either subset of nucs that are -1 or annotated version of sites
#' @seealso \code{\link{get_p1_nucs}} 
#' @export
get_m1_nuc<-function(nucs.ranges,sites, annotate_sites = F, max_dist = 350){
  plus=which(strand(sites)!="-")
  minus=which(strand(sites)=="-")
  
  m1 = follow(sites,nucs.ranges)
  m1_nna = which(!is.na(m1))
  m1_nna = m1_nna[which(seqnames(nucs.ranges[m1[m1_nna]])==seqnames(sites[m1_nna]))]
  if (annotate_sites){
    sites$m1 = NA
    sites$m1[m1_nna] = start(nucs.ranges[m1[m1_nna]])
    sites$m1_ix = NA
    sites$m1_ix[m1_nna] = m1[m1_nna]
    sites$m1_dist = NA
    m1_dist= rep(NA,length(sites))
    m1_dist[intersect(m1_nna,plus)] = start(nucs.ranges[m1[intersect(m1_nna,plus)]])-start(sites[intersect(m1_nna,plus)])
    m1_dist[intersect(m1_nna,minus)] = start(sites[intersect(m1_nna,minus)])-start(nucs.ranges[m1[intersect(m1_nna,minus)]])    
    close = which(m1_dist[m1_nna] >= -1*max_dist)
    sites$m1_dist[m1_nna[close]] = m1_dist[m1_nna[close]]
    return(sites)
  } else{
    m1_dist= rep(NA,length(sites))
    m1_dist[intersect(m1_nna,plus)] = start(nucs.ranges[m1[intersect(m1_nna,plus)]])-start(sites[intersect(m1_nna,plus)])
    m1_dist[intersect(m1_nna,minus)] = start(sites[intersect(m1_nna,minus)])-start(nucs.ranges[m1[intersect(m1_nna,minus)]])
    close = which(m1_dist[m1_nna] >= -1*max_dist)
    m1_nucs = nucs.ranges[m1[m1_nna[close]]]
    strand(m1_nucs)=strand(sites)[m1_nna[close]]
    m1_nucs$dist = m1_dist[m1_nna[close]]
    return(m1_nucs)
  }
}