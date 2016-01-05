# NucleoATACR
R package for reading in and working with [NucleoATAC](https://github.com/GreenleafLab/NucleoATAC) outputs.

Also includes functions for finding +1 and -1 nucleosomes, finding distances between nucleosomes, and plotting V-plots.

Install using devtools:
```R
devtools::install_github("GreenleafLab/NucleoATACR")
```

Load:
```R
library(NucleoATACR)
```

Read in nucleosome and nfr positions:
```R
nucs <- readNucs("test.nucmap_combined.bed.gz")
nfrs <- readNFRs("test.nfrpos.bed.gz")
```

Read in nucleoatac signal track for particular locus:
```R
signal <- readBedgraph("test.nucleoatac_signal.bedgraph.gz", "chrII", 706551, 707705)
```

Read in vplot and plot:
```R
v <- read_vplot("test.VMat")
plotV(v)
```

Get +1 and -1 nucleosomes:
```R
fake_tss <- GenomicRanges::GRanges("chrII", IRanges::IRanges(start = c(707200,707500), width = 1), 
              strand = c("+","-")) 
p1 <- get_p1_nuc(nucs.ranges = nucs, sites = fake_tss)
m1 <- get_m1_nuc(nucs.ranges = nucs, sites = fake_tss)
```







