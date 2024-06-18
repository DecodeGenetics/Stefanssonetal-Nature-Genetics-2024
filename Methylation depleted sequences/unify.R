# load necessary packages
suppressPackageStartupMessages(library(GenomicRanges)); suppressPackageStartupMessages(library(data.table))

# set path to output files from "collect_sequential_github.R"
dpath <- "/.../representative-sequential-unmethylHap-per-cluster_20x/"; d <- dir(dpath)

# set "dt0" as the file index
dt0 <- data.table(path2file=paste0(dpath,d), fileName=d)

# set class
dt0[,Class:=unlist(tstrsplit(fileName,"_",keep=1))]

# set chromosome
dt0[,Chrom:=unlist(tstrsplit(fileName,"_",keep=4))]

# collect files into object "f"
  f <- NULL; for(i in 1:nrow(dt0)) { 
    
    print(paste0(i, " out of ", nrow(dt0)))
    
      if(file.size(as.character(dt0$path2file[i]))>1){
      
        f0 <- fread(as.character(dt0$path2file[i])); f0$Class <- dt0[i,Class]; f0$fileName <- dt0[i,fileName]; f <- rbindlist(list(f, f0)) 
      
      }
    
  }

# set pos
f$pos <- unlist(tstrsplit(f$methylHap.id,":",keep=2))

# set chromosome
f$chromosome <- unlist(tstrsplit(f$methylHap.id,":",keep=1))

# set sex
f$sex <- unlist(tstrsplit(f$fileName,"_",keep=5)); f$sex[(f$sex%in%c("males","females"))==F] <- "auto"; f$sex[f$chromosome%in%"chrY"] <- "males"; f$sex[f$chromosome%in%"chrX" & f$sex%in%"auto"] <- "females"

# show sex by chromosome
f[,.N,c("sex","chromosome")]

# ensure Start is numeric
f$Start <- as.numeric(unlist(tstrsplit(f$pos,"-",keep=1)))

# ensure End is numeric
f$End <- as.numeric(unlist(tstrsplit(f$pos,"-",keep=2)))

# set width
f$width <- f$End+1-f$Start

# show summary of width
summary(f$width)

# total width
sum(f$width)

# retain relevant information
fOut <- f[,.(Chrom=chromosome, 
             Start=as.numeric(unlist(tstrsplit(pos,"-",keep=1))), 
             End=as.numeric(unlist(tstrsplit(pos,"-",keep=2))), 
             methylHap.id, 
             Class, 
             sex,
             cluster, 
             N_withinCluster,
             SN_depleted,
             N_measured_in_PNSN_median,
             N_depleted_in_PNSN_median,
             N_total_hq,
             Dist_between_CpGs_mean,
             flanking_measured_pre,flanking_measured_pre_methylRate,
             flanking_measured_next,flanking_measured_next_methylRate,
             most_upstream_position, most_downstream_position)]

# set width
fOut$width <- fOut$End+1-fOut$Start

options(scipen=999)

# write to disk

tabFile <- "/.../SequentialUnmethylatedSequences.tab"

write.table(fOut, file=tabFile, quote=FALSE, sep="\t", row.names=FALSE)

# End!
