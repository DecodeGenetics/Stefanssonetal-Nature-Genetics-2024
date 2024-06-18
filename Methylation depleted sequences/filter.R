args <- commandArgs(TRUE)

options(scipen=999)

i <- as.numeric(as.character(args[1]))

suppressPackageStartupMessages(library(data.table)); suppressPackageStartupMessages(library(GenomicRanges)); suppressPackageStartupMessages(library(dplyr)); suppressPackageStartupMessages(library(tictoc))

dt0 <- fread("/.../SequentialUnmethylHaps_PNSN.idx")

# Show data format of "dt0": 

# head(dt0,3)
# 
# path2file
# <char>
# 1: /.../SequentialUnmethylHaps_by_PNSN-20x/cleaned_mat_chr1_0.hq
# 2: /.../SequentialUnmethylHaps_by_PNSN-20x/cleaned_mat_chr1_1.hq
# 3: /.../SequentialUnmethylHaps_by_PNSN-20x/cleaned_mat_chr1_10.hq
# fileName        Class
# <char>       <char>
# 1: cleaned_mat_chr1_0.hq unmethylated
# 2: cleaned_mat_chr1_1.hq unmethylated
# 3: cleaned_mat_chr1_10.hq unmethylated

f <- NULL

  print(paste0("i=", i))
  
  # open
  f0 <- fread(as.character(dt0$path2file[i]))

  # set as logical
  f0$nextPos.ok <- as.logical(f0$nextPos.ok); f0$prePos.ok <- as.logical(f0$prePos.ok)
  
  # only use haplotypes that are "unbroken", i.e. where the CpGs preceding and following the start and end positions, respectively, of a methylation depleted haplotype are known!
  f0 <- f0[prePos.ok==TRUE & nextPos.ok==TRUE, ]
  
  # determine whether any haplotypes were found that satisfied the criteria from above:   
  if(nrow(f0)>0){

    # order 
    f0 <- f0[order(Start)]

    # make GRanges object    
    gr0 <- makeGRangesFromDataFrame(f0, keep.extra.columns = T)
    
    gr <- reduce(gr0, with.revmap=T)
    
    n0 <- NULL
    
    for(j in 1:length(gr)){
      
      # get overlapping methylHaplotypes
      grHaps <- gr0[gr$revmap[[j]]]
      
      tmp <- as.data.table(grHaps)
      
      # keep info on the "overarching" cluster containing methylHaplotypes that overlap in one way or another
      tmp0 <- gr[j]; clusterOverlap <- paste0(seqnames(tmp0),":", start(tmp0), "-", end(tmp0))
      
      # find most common methylHaplotype
      n <- tmp[,.N,c("methylHap.id","width")][order(-N, -width)]
      
      # re-define "tmp" and "n" to include only those methylHapltypes that overlap with the most frequent methylHaplotype
      tmp <- as.data.table( grHaps[ grHaps%over%  grHaps[grHaps$methylHap.id%in%n[1]][1] ] )
      
      n <- tmp[,.N,c("methylHap.id","width")][order(-N, -width)]

      # keep results on the most frequent methylHaplotype
      n0 <- rbindlist(list(n0, 
                           data.table(cluster=clusterOverlap, 
                                      methylHap.id=n$methylHap.id[1], 
                                      N_withinCluster=paste0(n$N, collapse=","), 
                                      SN_depleted=sum(n$N))))
      
      # see if there are haplotypes that do not overlap with the most common haplotype selected as a representative for this cluster
      while(length(grHaps)>0){
        
        grHaps <- grHaps[(grHaps %over% grHaps[grHaps$methylHap.id %in% n$methylHap.id[1]])==FALSE]
        
        if(length(grHaps)>0){
          
          # redefine "tmp" 
          tmp <- as.data.table(grHaps)
          
          n <- tmp[,.N,c("methylHap.id","width")][order(-N, -width)]

          n0 <- rbindlist(list( n0, 
                                data.table(cluster=clusterOverlap, 
                                          methylHap.id=n$methylHap.id[1], 
                                          N_withinCluster=paste0(n$N, collapse=","), 
                                          SN_depleted=sum(n$N)) ))
        }
        
      }
      
    }
    
    f <- rbindlist(list(f,n0))
    
  }

# see if there is any output to report
if(is.null(f)==F){

if(nrow(f)>0){

# ensure no duplicated "methylHap.id" (just in case!)
f <- f[duplicated(methylHap.id)==F, ]
    
# add info on total number of hq CpG units relevant to each methylHap.id
f <- left_join(f, f0[,.(N_total_hq=median(N_total_hq)),by="methylHap.id"], by="methylHap.id")

# add info on number of measured CpG units in each methylHap.id
f <- left_join(f, f0[,.(N_measured_in_PNSN_median=median(N_measured_in_PNSN)),by="methylHap.id"], by="methylHap.id")

# add info on the number of CpG units depleted of 5-mCpGs in each methylHap.id
f <- left_join(f, f0[,.(N_depleted_in_PNSN_median=median(N_depleted)),by="methylHap.id"], by="methylHap.id")

# add info on the mean distance between CpG units in each methylHap.id
f <- left_join(f, f0[,.(Dist_between_CpGs_mean =signif(median(mean_dist_between_CpGs ), digits=3)),by="methylHap.id"], by="methylHap.id")

# add info on measured CpG units flanking up- and downstream of unmethylated sequences
f0[, flanking_measured_pre:=as.numeric(unlist(tstrsplit(flanking_measured_pos,",",keep=1)))]
f0[, flanking_measured_next:=as.numeric(unlist(tstrsplit(flanking_measured_pos,",",keep=2)))]

# add info on the 5-mCpG rate of CpG units flanking up and downstream of unmethylated sequences
f0[,flanking_measured_pre_methylRate:=as.numeric(unlist(tstrsplit(flanking_measured_rate,",",keep=1)))]
f0[,flanking_measured_next_methylRate:=as.numeric(unlist(tstrsplit(flanking_measured_rate,",",keep=2)))]

# set the grouping column
setkey(f0,"methylHap.id")

# find the methylRate of the CpG unit that is closest to the start position of each methylHap.id
f <- left_join(f, f0[,.SD[which.max(flanking_measured_pre)], by=key(f0)][,.(methylHap.id, flanking_measured_pre,flanking_measured_pre_methylRate)], by="methylHap.id")

# find the methylRate of the CpG unit that is closest to the end position of each methylHap.id
f <- left_join(f, f0[,.SD[which.min(flanking_measured_next)], by=key(f0)][,.(methylHap.id, flanking_measured_next,flanking_measured_next_methylRate)], by="methylHap.id")

# store information on the most upstream position observed for each methylHap.id
f <- left_join(f, f0[,.SD[which.min(flanking_measured_pre)], by=key(f0)][,.(methylHap.id, most_upstream_position=flanking_measured_pre)], by="methylHap.id")

# store information on the most downstream position observed for each methylHap.id
f <- left_join(f, f0[,.SD[which.max(flanking_measured_next)], by=key(f0)][,.(methylHap.id, most_downstream_position=flanking_measured_next)], by="methylHap.id")

# set path and filename of output
fOut <- paste0("/.../representative-sequential-unmethylHap-per-cluster_20x/",as.character(dt0$Class[i]),"_", as.character(dt0$fileName[i]))

# write out the output
write.table(f, fOut, quote=FALSE, sep="\t", row.names=FALSE)

}

}

# End!