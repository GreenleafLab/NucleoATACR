

#' get_dist_between_calls
#'
#' @param ranges GRanges object 
#' @param max_dist to consider
#' @return returns tabulated distance between calls
#' @seealso \code{\link{distRanges}} 
#' @export
get_dist_between_calls<- function(ranges, max_dist = 1000){
  sep = BiocGenerics::start(ranges)[2:(length(ranges))] - BiocGenerics::end(ranges)[1:(length(ranges)-1)] 
  chrom_same = which(GenomicRanges::seqnames(ranges)[1:(length(ranges)-1)] == GenomicRanges::seqnames(ranges)[2:(length(ranges))])
  close = intersect(which(sep < max_dist),chrom_same)
  return(tabulate(sep[close]))
}

#' getstrand
#'
#' @param ranges GRanges object 
#' @return returns strand as character (rather than factor Rle)
#' @seealso \code{\link{shiftRanges}} 
#' @export
getstrand<- function(ranges){
  return(as.character(BiocGenerics::strand(ranges)))
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
  BiocGenerics::strand(ranges1_mod)="*"
  neighbor = nearest(ranges1_mod,ranges2)  
  if (strand == 2){
    dists = rep(NA,length(ranges1))
    plus = which(getstrand(ranges2)[neighbor]!="-")
    minus = which(getstrand(ranges2)[neighbor]=="-")
    dists[minus] = ifelse(overlapsAny(ranges1[minus],ranges2[neighbor[minus]]),
                          0,
                          ifelse(BiocGenerics::start(ranges2[neighbor[minus]])>BiocGenerics::start(ranges1[minus]),
                                 BiocGenerics::start(ranges2)[neighbor[minus]] - BiocGenerics::end(ranges1)[minus],
                                 BiocGenerics::end(ranges2)[neighbor[minus]] - BiocGenerics::start(ranges1)[minus]))
    dists[plus] = ifelse(overlapsAny(ranges1[plus],ranges2[neighbor[plus]]),
                         0,
                         ifelse(BiocGenerics::start(ranges2[neighbor[plus]])>BiocGenerics::start(ranges1[plus]),
                                BiocGenerics::end(ranges1)[plus] - BiocGenerics::start(ranges2)[neighbor[plus]],
                                BiocGenerics::start(ranges1)[plus] - BiocGenerics::end(ranges2)[neighbor[plus]]))                    
  } 
  else if (strand == 0){
    dists = ifelse(overlapsAny(ranges1,ranges2[neighbor]),
                         0,
                         ifelse(BiocGenerics::start(ranges2[neighbor])>BiocGenerics::start(ranges1),
                                BiocGenerics::end(ranges1) - BiocGenerics::start(ranges2)[neighbor],
                                BiocGenerics::start(ranges1) - BiocGenerics::end(ranges2)[neighbor]))                    
  }
  else if (strand ==1){
    dists = rep(NA,length(ranges1))
    plus = which(getstrand(ranges1)!="-")
    minus = which(getstrand(ranges1)=="-")
    dists[plus] = ifelse(overlapsAny(ranges1[plus],ranges2[neighbor[plus]]),
                          0,
                          ifelse(BiocGenerics::start(ranges2[neighbor[plus]])>BiocGenerics::start(ranges1[plus]),
                                 BiocGenerics::start(ranges2)[neighbor[plus]] - BiocGenerics::end(ranges1)[plus],
                                 BiocGenerics::end(ranges2)[neighbor[plus]] - BiocGenerics::start(ranges1)[plus]))
    dists[minus] = ifelse(overlapsAny(ranges1[minus],ranges2[neighbor[minus]]),
                         0,
                         ifelse(BiocGenerics::start(ranges2[neighbor[minus]])>BiocGenerics::start(ranges1[minus]),
                                BiocGenerics::end(ranges1)[minus] - BiocGenerics::start(ranges2)[neighbor[minus]],
                                BiocGenerics::start(ranges1)[minus] - BiocGenerics::end(ranges2)[neighbor[minus]]))                    
    return(dists)
  }
  else{
    stop("strand must equal 0, 1, or 2.")
  }
  return(dists)
}



#' shiftRanges
#'
#' @param ranges GRanges object 
#' @param shift amount by which to shift ranges (negative means shift upstream, positive downstream)
#' @return returns GRanges object for which each range has been shifted appropriately
#' @seealso \code{\link{get_m1_nucs}} \code{\link{get_p1_nucs}} \code{\link{distRanges}} 
#' @export
shiftRanges <- function(ranges, shift = 0){
  minus = which(getstrand(ranges)=="-")
  plus = which(getstrand(ranges)!="-")
  if (shift < 0){
    BiocGenerics::end(ranges)[minus] = BiocGenerics::end(ranges)[minus] - shift
    BiocGenerics::start(ranges)[minus] = BiocGenerics::start(ranges)[minus] - shift
    BiocGenerics::start(ranges)[plus] = BiocGenerics::start(ranges)[plus] + shift
    BiocGenerics::end(ranges)[plus] = BiocGenerics::end(ranges)[plus] + shift
  }
  else if (shift > 0){
    BiocGenerics::start(ranges)[minus] = BiocGenerics::start(ranges)[minus] - shift
    BiocGenerics::end(ranges)[minus] = BiocGenerics::end(ranges)[minus] - shift
    BiocGenerics::end(ranges)[plus] = BiocGenerics::end(ranges)[plus] + shift
    BiocGenerics::start(ranges)[plus] = BiocGenerics::start(ranges)[plus] + shift
  }
  return(ranges)
}


#' get_p1_nuc
#'
#' @param nuc.ranges GRanges object with nucleosome positions
#' @param sites GRanges object with sites to get +1 nucs from.
#' @param annotate_sites what to return-- sites with +1 nuc annotation or +1 nucs?
#' @param max_dist what is maximum distance to consider
#' @param shift integer shift to apply to sites before finding -1
#' @return returns either subset of nucs that are +1 or annotated version of sites
#' @seealso \code{\link{get_m1_nucs}} 
#' @export
get_p1_nuc<-function(nucs.ranges, sites, annotate_sites = F, max_dist = 250, shift = 0){
  plus=which(getstrand(sites) != "-")
  minus=which(getstrand(sites) == "-")
  if (shift!=0){
     p1 = precede(shiftRanges(sites,shift), nucs.ranges)
  }
  else{
      p1 = precede(sites,nucs.ranges)  
  }  
  p1_nna = which(!is.na(p1))
  if (length(p1_nna)>0){
      p1_nna = p1_nna[which(as.character(GenomicRanges::seqnames(nucs.ranges[p1[p1_nna]]))==as.character(GenomicRanges::seqnames(sites[p1_nna])))]
  }
  if (annotate_sites){
    sites$p1 = NA
    sites$p1[p1_nna] = BiocGenerics::start(nucs.ranges[p1[p1_nna]])
    sites$p1_ix = NA
    sites$p1_ix[p1_nna] = p1[p1_nna]  
    sites$p1_dist = NA
    p1_dist = rep(NA,length(sites))
    p1_dist[intersect(p1_nna,plus)] = BiocGenerics::start(nucs.ranges[p1[intersect(p1_nna,plus)]])-BiocGenerics::end(sites[intersect(p1_nna,plus)])
    p1_dist[intersect(p1_nna,minus)] = BiocGenerics::start(sites[intersect(p1_nna,minus)])-BiocGenerics::end(nucs.ranges[p1[intersect(p1_nna,minus)]])
    close = which(p1_dist[p1_nna] <= max_dist)
    sites$p1_dist[p1_nna[close]]=p1_dist[p1_nna[close]]
    return(sites)
  } else{
    p1_dist = rep(NA,length(sites))
    p1_dist[intersect(p1_nna,plus)] = BiocGenerics::start(nucs.ranges[p1[intersect(p1_nna,plus)]])-BiocGenerics::end(sites[intersect(p1_nna,plus)])
    p1_dist[intersect(p1_nna,minus)] = BiocGenerics::start(sites[intersect(p1_nna,minus)])-BiocGenerics::end(nucs.ranges[p1[intersect(p1_nna,minus)]])
    close = which(p1_dist[p1_nna] <= max_dist)
    p1_nucs = nucs.ranges[p1[p1_nna[close]]]
    BiocGenerics::strand(p1_nucs)=BiocGenerics::strand(sites)[p1_nna[close]]
    p1_nucs$dist = p1_dist[p1_nna[close]]
    return(p1_nucs)
  }
}

#' get_m1_nuc
#'
#' @param nuc.ranges GRanges object with nucleosome positions
#' @param sites GRanges object with sites to get -1 nucs from.
#' @param annotate_sites what to return-- sites with -1 nuc annotation or +1 nucs?
#' @param max_dist what is maximum distance to consider
#' @param shift integer shift to apply to sites before finding -1
#' @return returns either subset of nucs that are -1 or annotated version of sites
#' @seealso \code{\link{get_p1_nucs}} 
#' @export
get_m1_nuc<-function(nucs.ranges,sites, annotate_sites = F, max_dist = 350, shift = 0){
  plus=which(getstrand(sites)!="-")
  minus=which(getstrand(sites)=="-")
  if (shift!=0){
    m1 = follow(shiftRanges(sites,shift), nucs.ranges)
  }
  else{
    m1 = follow(sites,nucs.ranges)  
  }    
  m1 = follow(sites,nucs.ranges)
  m1_nna = which(!is.na(m1))
  m1_nna = m1_nna[which(as.character(GenomicRanges::seqnames(nucs.ranges[m1[m1_nna]]))==as.character(GenomicRanges::seqnames(sites[m1_nna])))]
  if (annotate_sites){
    sites$m1 = NA
    sites$m1[m1_nna] = BiocGenerics::start(nucs.ranges[m1[m1_nna]])
    sites$m1_ix = NA
    sites$m1_ix[m1_nna] = m1[m1_nna]
    sites$m1_dist = NA
    m1_dist= rep(NA,length(sites))
    m1_dist[intersect(m1_nna,plus)] = BiocGenerics::end(nucs.ranges[m1[intersect(m1_nna,plus)]])-BiocGenerics::start(sites[intersect(m1_nna,plus)])
    m1_dist[intersect(m1_nna,minus)] = BiocGenerics::end(sites[intersect(m1_nna,minus)])-BiocGenerics::start(nucs.ranges[m1[intersect(m1_nna,minus)]])
    close = which(m1_dist[m1_nna] >= -1*max_dist)
    sites$m1_dist[m1_nna[close]] = m1_dist[m1_nna[close]]
    return(sites)
  } else{
    m1_dist= rep(NA,length(sites))
    m1_dist[intersect(m1_nna,plus)] = BiocGenerics::end(nucs.ranges[m1[intersect(m1_nna,plus)]])-BiocGenerics::start(sites[intersect(m1_nna,plus)])
    m1_dist[intersect(m1_nna,minus)] = BiocGenerics::end(sites[intersect(m1_nna,minus)])-BiocGenerics::start(nucs.ranges[m1[intersect(m1_nna,minus)]])
    close = which(m1_dist[m1_nna] >= -1*max_dist)
    m1_nucs = nucs.ranges[m1[m1_nna[close]]]
    BiocGenerics::strand(m1_nucs)=BiocGenerics::strand(sites)[m1_nna[close]]
    m1_nucs$dist = m1_dist[m1_nna[close]]
    return(m1_nucs)
  }
}