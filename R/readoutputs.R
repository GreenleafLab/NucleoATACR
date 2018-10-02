
#' readNucs
#'
#' @param nucfile file of nucleosome positions from NucleoATAC
#' @param out desired object type, either data.frame or GRanges
#' @return returns object of type specified by format, either data.frame or GRanges
#' @seealso \code{\link{readNFRs}}
#' @importFrom readr cols col_character
#' @export
readNucs<-function(nucfile, out = "GRanges"){
  out = as.character(out)
  if (grepl("nucpos.bed",nucfile)){
    nucDF=as.data.frame(readr::read_tsv(nucfile,col_names = c('chr','start','end',"z","occ","occ_lower","occ_upper","lr","nuc_signal","raw_signal","reads","nfr","fuzz"), col_types=cols(chr = col_character())))
  }
  else if (grepl("nucmap_combined.bed",nucfile)){
    nucDF=as.data.frame(readr::read_tsv(nucfile,col_names = c('chr','start','end',"occ","occ_lower","occ_upper","reads","type"), col_types=cols(chr = col_character())))
  }
  else if (grepl("occpeaks.bed",nucfile)){
    nucDF=as.data.frame(readr::read_tsv(nucfile,col_names = c('chr','start','end',"occ","occ_lower","occ_upper","reads"), col_types=cols(chr = col_character())))       
  }
  else{
    stop("File name doesn't seem to include stantard format of NucleoATAC output")
  }
  if (out=="data.frame"){
    return(nucDF)
  }
  else if (out=="GRanges"){
    if (grepl("nucpos.bed",nucfile)){
      nucGR = with(nucDF,GenomicRanges::GRanges(chr,IRanges::IRanges(start,start),
                        z=z,
                        occ = occ,
                        occ_lower = occ_lower,
                        occ_upper = occ_upper,
                        lr = lr,
                        nuc_signal=nuc_signal,
                        raw_signal = raw_signal,
                        reads=reads,
                        nfr=nfr,
                        fuzz=fuzz))
    }
    else if (grepl("nucmap_combined.bed",nucfile)){
      nucGR = with(nucDF,GenomicRanges::GRanges(chr,IRanges::IRanges(start,start),
                             occ = occ,
                             occ_lower = occ_lower,
                             occ_upper = occ_upper,
                             type = type))
    }
    else if (grepl("occpeaks.bed",nucfile)){
      nucGR = with(nucDF,GenomicRanges::GRanges(chr,IRanges::IRanges(start,start),
                                 occ = occ,
                                 occ_lower = occ_lower,
                                 occ_upper = occ_upper))
    }
  }
  else{
    stop("out must equal either GRanges or data.frame")
  }
}

#' readNFRs
#'
#' @param nfrfile file of nucleosome positions from NucleoATAC
#' @param format desired object type, either data.frame or GRanges
#' @return returns object of type specified by format, either data.frame or GRanges
#' @seealso \code{\link{readNucs}}
#' @export
readNFRs<-function(nfrfile, format = "GRanges"){
  format = as.character(format)
  nfrDF=as.data.frame(readr::read_tsv(nfrfile, col_names = FALSE))
  if (ncol(nfrDF) == 7){
     colnames(nfrDF) = c('chr','start','end',"occ","min_occ_upper","ins_density","bias_density")
  } else if (ncol(nfrDF) == 8){
    colnames(nfrDF) = c('chr','start','end',"occ","min_occ_upper","ins_density","bias_density","nfr_signal")
  }
  if (format=="data.frame"){
    return(nfrDF)
  }
  else if (format=="GRanges"){
    if (ncol(nfrDF) == 7){
      return(with(nfrDF,GenomicRanges::GRanges(chr,IRanges::IRanges(start,end-1),
                          occ = occ,
                          min_occ_upper = min_occ_upper,
                          ins = ins_density,
                          bias = bias_density)))
    }else if (ncol(nfrDF) == 8){
      return(with(nfrDF,GenomicRanges::GRanges(chr,IRanges::IRanges(start,end-1),
                                occ = occ,
                                min_occ_upper = min_occ_upper,
                                ins = ins_density,
                                bias = bias_density,
                                nfr_signal = nfr_signal)))
    }
  }
  else{
    stop("Format must equal either GRanges or data.frame")
  }
}

#' readBedgraph
#'
#' @param bgfile tabix-indexed bedgraph file
#' @param chrom chromosome
#' @param start start (0-based)
#' @param end end (1-based)
#' @param empty value for positions not in bedgraph to take (default NA)
#' @return returns vector of bedgraph values for interval
#' @seealso \code{\link{readNucs}}  \code{\link{readNFRs}}
#' @export
readBedgraph<-function(bgfile, chrom, start, end, empty = NA){
  tmp = Rsamtools::scanTabix(bgfile, param = GenomicRanges::GRanges(chrom, IRanges::IRanges(start,end-1)))
  out = rep(empty, end - start)
  for (rec in tmp[[1]]){
    vals = as.numeric(strsplit(rec,"\t")[[1]][2:4])
    out[(vals[1]-start+1):(vals[2]-start+1)] = vals[3]
  }
  return(out)
}
