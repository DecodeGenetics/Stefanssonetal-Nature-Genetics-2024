# specify argument "l", as any single value in {1,2,3,... nrow(idx.m)}
args <- commandArgs(TRUE); l <- as.numeric(as.character(args[1]))

# load necessary libraries
suppressPackageStartupMessages(library(data.table)); suppressPackageStartupMessages(library(dplyr)); suppressPackageStartupMessages(library(GenomicRanges))

# load methylation-index file (idx.m): #
idx.m <- fread("/.../list-of-phasedMeth-mats-230817.gor", header=T); print(idx.m[l])

# Show data format of "idx.m": 

# head(idx.m,3)
# Chrom   Start     End     path2file
# <char>   <int>   <num>    <char>
# 1:   chr1  891453  999977 /.../cleaned_mat_chr1_0.hq
# 2:   chr1 1000094 1999991 /.../cleaned_mat_chr1_1.hq
# 3:   chr1 2000117 2999969 /.../cleaned_mat_chr1_2.hq
# fileName    sex     i
# <char> <char> <int>
# 1: cleaned_mat_chr1_0.hq   auto  1998
# 2: cleaned_mat_chr1_1.hq   auto 16034
# 3: cleaned_mat_chr1_2.hq   auto 14221

# set the .hq filename
filename.x <-idx.m$fileName[l]

# set the minium number of lines (CpG units) to analyse in each .hq data file
min.CpG <- 100

# open the relevant .hq file specified as the row number "l" from data.table "idx.m"
f <- fread(as.character(idx.m$path2file[l]), header=T)

# Show data format of .hq files referred to in path2file column of "idx.m" from above: 

# fread("/.../cleaned_mat_chr1_2.hq")[1:3,1:10]
# chr     pos PN1_SN1_P_0 PN1_SN1_P_1 PN1_SN1_M_0 PN1_SN1_M_1 PN2_SN1_P_0
# <char>   <int>              <int>              <int>              <int>              <int>              <int>
# 1:   chr1 2000117                  1                  5                  2                  9                  0
# 2:   chr1 2000141                  0                 12                  0                 15                  0
# 3:   chr1 2000299                  0                  7                  2                  9                  1
# PN2_SN1_P_1 PN2_SN1_M_0 PN2_SN1_M_1
# <int>              <int>              <int>
# 1:                  1                  0                  2
# 2:                  4                  0                  7
# 3:                  3                  0                  4

# Notes on data format of .hq files: 
#
# Column names are composed as follows: PN_SN_Haplotype_Status, where: 
#
# PN is the individual (censored) and SN is the sample identifier (censored), e.g. PN1_SN1 = whole blood sample nr. 1 obtained from individual nr.1.  
#
# The ending of the column names (except the first two columns) represent the parental haplotype (P=paternal or M=maternal) and the methylation status (0=unmethylated, 1=methylated)

# stop if "f" is empty
if((nrow(f)>2)==F) stop("no data to process!")

# ensure column nr. 2 is named "pos"
colnames(f)[2] <- "pos"

# set the maximum position
pos.max <- f[,max(pos)]

# set "sex" and "chromo"
sex <- as.character(idx.m$sex[l]); chromo <- as.character(idx.m$Chrom[l])

# add some CpGs from the next file
if(l<nrow(idx.m) & chromo !="chrX"){
  
  # open file to add
  f.add <- fread(as.character(idx.m$path2file[l+1]), nrow=min.CpG)
  
  # ensure column nr. 2 is called "pos"
  colnames(f.add)[2] <- "pos"
  
  if(length(colnames(f))==length(colnames(f.add))){
  
    # add data 
    if(all(colnames(f)==colnames(f.add)) & idx.m$Chrom[l]==idx.m$Chrom[l+1]) f <- rbindlist(list(f,f.add))
    
  }
  
}

# address the special case when chromo="chrX"
if(l<nrow(idx.m) & chromo=="chrX"){

fileName.hq <- as.character(idx.m$path2file[l])

idx.m0 <- as.data.frame(idx.m); idx.m0 <- idx.m0[idx.m0$sex %in% sex & idx.m0$Chrom %in% chromo, ]

w0 <- which(idx.m0$path2file %in% fileName.hq)
  
  if(w0+1<nrow(idx.m0)){

    # open file to add
    f.add <- fread(as.character(idx.m0$path2file[w0+1]), nrow=min.CpG)
    
    # ensure column nr. 2 is called "pos"
    colnames(f.add)[2] <- "pos"
    
      if(length(colnames(f))==length(colnames(f.add))){
        
        # add data
        if(all(colnames(f)==colnames(f.add))) f <- rbindlist(list(f,f.add))
        
      }
    
  }
  
}

# keep information on measured CpGs in "f.gr"
f.gr <- GRanges(seqnames=chromo, IRanges(start=f$pos, end=f$pos))

# load covariates
coVar <- fread("/.../covariates.tab")

# Show data format of "coVar":

# head(coVar[,.(pn,dna_sn,coverage)],3)
# pn dna_sn coverage
# <char>  <int>    <num>
# 1: PN1 SN1      6.6
# 2: PN2 SN1     29.7
# 3: PN3 SN1      2.2

# set SN ids
coVar[,PN_SN:=paste0(pn,"_",dna_sn)]

# find PN_SNs where coverage is higher than, or equal to, 20x
pn_sn_20x <- coVar[coverage>=20, unique(PN_SN)]

# initiate "m", to be used for collecting relevant data
m <- NULL

# transpose matrix, filter and compute ratios #
print("transpose to long format")

for(j in 1:nrow(f)){
  
  # select CpG in line "j"
  dat <- data.table(transpose(f[j,-(1:2)]), id=colnames(f)[-(1:2)])
  
  # column name
  colnames(dat)[1] <- "methyl_counts"
  
  # set IDs
  dat[,PN:=unlist(tstrsplit(id,"_",keep=1))]
  dat[,PN_SN:=paste0(unlist(tstrsplit(id,"_",keep=1)), "_", unlist(tstrsplit(id,"_",keep=2)))]
  
  # restrict to "PN_SN" where average coverage is >20x
  dat <- dat[PN_SN %in% pn_sn_20x]
  
  # set state
  if(sex!="males"){
    dat[,hap:=unlist(tstrsplit(id,"_",keep=3))]
    dat[,state:=unlist(tstrsplit(id,"_",keep=4))]
  } else {
    if(chromo=="chrY") dat[,hap:="P"]
    if(chromo=="chrX") dat[,hap:="M"]
    dat[,state:=unlist(tstrsplit(id,"_",keep=3))]
  }
  
  # set haplotype IDs as: PN_SN_phase
  dat[,PN_SN_phase:=paste0(PN_SN,"_",hap)]
  
  # remove invalid lines (NA and -1)
  dat <- dat[is.na(methyl_counts)==F, ]; dat <- dat[!(methyl_counts==(-1)), ]
  
  # compute sum of methylation counts (total number of measured CpGs per parental haplotype)
  datSum <- dat[,lapply(.SD,sum), .SDcol="methyl_counts", by="PN_SN_phase"]; colnames(datSum)[2] <- c("methyl_sum")
  
  # join sum to "dat"
  dat <- left_join(dat,datSum,by="PN_SN_phase")
  
  # restrict to state==1
  dat <- dat[state==1, ]
  
  # compute ratio
  dat$methylRatio <- dat$methyl_counts/dat$methyl_sum
  
  # set minimum coverage per each CpG unit
  dat <- dat[methyl_sum<10, methylRatio:=NA]
  
  # extract relevant column into temporary data.table called "mTmp"
  mTmp <- na.omit(data.table(dat[,.(Pos=as.numeric(f$pos[j]), PN_SN_phase, methylRatio)]))
  
  # collect relevant data into data.table called "m"
  m <- rbindlist(list(m, mTmp))
  
}

##########################################################
# main routine goes through each PN_SN_phase, one-by-one #
##########################################################

# see whether "m" is still null
if(is.null(m)==F){
  
# run if "m" contains haplotype specific measurements of CpG units
if(nrow(m)>0){

# find unique PN_SN_phase
u <- m[,.N,c("PN_SN_phase")][N>1, unique(PN_SN_phase)]

# intitalize "r" which will be used to collect information on umethylated intervals/haplotypes
r <- data.table(Chrom=NA, Start=NA, End=NA, methylHap.id=NA, 
                N_depleted=NA, N_measured_in_PNSN=NA, N_total_hq=NA, 
                nextPos.ok=NA, prePos.ok=NA,
                mean_dist_between_CpGs=NA, 
                meanRate=NA, 
                PN_SN_phase=NA,
                flanking_hq_pos=NA,
                flanking_hq_rate=NA,
                flanking_measured_pos=NA,
                flanking_measured_rate=NA) 

# go through each PN_SN_phase in "u"
for(i in 1:length(u)){
  
  # show progress
  print(paste0("PN_SN_phase: i=",i," | ", u[i]))
  
  # select a unique PN_SN to analyse
  m.sub0 <- m[PN_SN_phase%in%u[i] & is.na(methylRatio)==F & !(methylRatio==(-1)), ][order(Pos)]
  
  # assign states (0=unmethylated units; 1=methylated units)
  m.sub0$State <- 1*( (m.sub0$methylRatio>0.15) )
  
  # compute distance to next CpG
  m.sub0[, PosDiff := c(0, diff(Pos))]  # Initialize with 0 for the first difference
  
  # quickly find out the line numbers in m.sub0 where there is a CpG unit with State==0
  w <- which(m.sub0$State==0)
  
  # see if there are any such CpG units, i.e. where State==0
  if(length(w)>0){

  # if so, then go to that line number
  for(j in 1:length(w)){
    
    # see data from that CpG unit where State==0
    m.sub00 <- m.sub0[w[j]:nrow(m.sub0), ]
    
    # initialize "pos"
    pos <- NULL
    
    # see if there is any data to work with
    if(nrow(m.sub00)>1){
    
    # go to the next CpGs
    for(z in 2:nrow(m.sub00)){
      
      # see if the next CpG unit is also State==0, and also determine whether it is located less than 500bp from the CpG unit in "z"
      l <- (m.sub00[z, State]==0 & m.sub00[z, PosDiff]<500)
      
      # if l==TRUE, then keep information on CpG unit position
      if(l){ pos <- c(pos,z)}
    
      # if "l" is no longer TRUE, then stop the for-loop
      if(l==FALSE) break()
        
    }
    
    # see if there is data to keep
    if(is.null(pos)==F){
    
      # set the start position of the unmethylated interval/haplotype
      startPosition <- m.sub00$Pos[1]
      
      # set the end position of the unmethylated interval/haplotype
      endPosition <- m.sub00$Pos[pos[length(pos)]]
      
      # set "methylHap.id"
      methylHap.id <- paste0(chromo,":",startPosition,"-",endPosition)
      
      # compute mean of methylRatio for the unmethylated interval/haplotype
      meanRate <- signif(mean(m.sub00$methylRatio[1:pos[length(pos)]], na.rm=T), digits=4)
      
      # determine the number of CpG units in the unmethylated interval/haplotype
      N_depleted <- length(pos) + 1
      
      # by definition, "N_measured_in_PNSN" is equal to "N_depleted"
      N_measured_in_PNSN <- N_depleted
      
      # determine the mean distance between CpG units for the unmetylated interval/haplotype
      mean_dist_between_CpGs <- signif(mean(m.sub00$PosDiff[2:pos[length(pos)]], na.rm=T), digits=4)
      
      # determine whether the flanking CpG units were actually measured in this PN_SN_phase
      pre.hq.cpg <- f$pos[which(f$pos %in% startPosition)-1]
      next.hq.cpg <- f$pos[which(f$pos %in% endPosition)+1]

      prePos.ok <- pre.hq.cpg %in% m.sub0$Pos
      nextPos.ok <- next.hq.cpg %in% m.sub0$Pos

      # determine the rate of the previous and next hq CpG units
      preState <- m.sub0[Pos %in% c(pre.hq.cpg), methylRatio]; if(length(preState)==0) preState <- NA
      nextState <- m.sub0[Pos %in% c(next.hq.cpg), methylRatio]; if(length(nextState)==0) nextState <- NA

      # determine the rate of the previous and next measured CpG units
      pre.measured.cpg <- m.sub0$Pos[which(m.sub0$Pos %in% startPosition)-1]
      next.measured.cpg <- m.sub0$Pos[which(m.sub0$Pos %in% endPosition)+1]
      preMeasuredState <- m.sub0[Pos %in% c(pre.measured.cpg), methylRatio]; if(length(preMeasuredState)==0) preMeasuredState <- NA
      nextMeasuredState <- m.sub0[Pos %in% c(next.measured.cpg), methylRatio]; if(length(nextMeasuredState)==0) nextMeasuredState <- NA
      
      # set the total number of hq CpG units within the unmethylated interval/haplotype
      N_total_hq <- length(f$pos[f$pos>=startPosition & f$pos<=endPosition])
      
        # determine whether interval has already been identified in the PN_SN_phase being analysed (i.e. u[i])
        if( (endPosition %in% r$End[r$PN_SN_phase %in% u[i]])==FALSE ){
        
        # if not already identified then keep relevant information about the interval in "r": 
        r <- rbindlist(list( r, data.table(Chrom=chromo, Start=startPosition, End=endPosition, methylHap.id, 
                                           N_depleted, N_measured_in_PNSN, N_total_hq, 
                                           nextPos.ok, prePos.ok,
                                           mean_dist_between_CpGs, 
                                           meanRate, 
                                           PN_SN_phase=u[i],
                                           flanking_hq_pos=paste0(pre.hq.cpg,",",next.hq.cpg),
                                           flanking_hq_rate=paste0(signif(preState,digits=3),",",signif(nextState,digits=3)),
                                           flanking_measured_pos=paste0(pre.measured.cpg,",",next.measured.cpg),
                                           flanking_measured_rate=paste0( signif(preMeasuredState,digits=3),",",signif(nextMeasuredState,digits=3)) )))
        
        }
      
    }
    
  }
    
  }
  
}
   
}

}

}

if(exists("r")==F) stop("no unmethylated haplotypes were found in this data chunk!")

options(scipen=999)

outpath <- ".../SequentialUnmethylHaps_by_PNSN-20x/"

if(file.exists(outpath)==F) system(paste0("mkdir ", outpath))

if(nrow(r)>1){

# write unmethylated haplotypes to disk
fileOut <- paste0(outpath, filename.x)

r$width <- (r$End + 1) - r$Start

write.table(r[-1,], file=fileOut, quote=FALSE, sep="\t", row.names=FALSE)

}

# End!
