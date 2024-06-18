arg <- commandArgs(TRUE)

# arg1: path to file containing markers; tab-separated file format: Chrom,Pos,Name ... these markers should be well-separated from each other e.g. spaced by 1Mb. 
# arg2: path to file containing the genomic attributes tested for association e.g. MDS or CpG coordinates for ASM-QTLs; tab-separated file format: Chrom,Start,End
# arg3: length (bp) from each moltrait element (e.g. unmethylated region, or gene) to consider for selecting random markers, e.g. 10000.   
# arg4: set the number of sequence variants in final output e.g. 200000
# arg5: path to output file.
# arg6: marker file from which to draw (Chrom,Pos,Name,ImpMAF)

# Example command: $ Rscript /.../make-random-marker-lists.R /.../markers-lead_spaced-by-1Mb.tab /.../MDSs.gor 10000 200000 /.../bag-of-random-markers.tab /.../markers_v3-varjoin.gor

markerlist <- as.character(arg[1])
moltrait <- as.character(arg[2])
e <- as.numeric(as.character(arg[3]))
n <- as.numeric(as.character(arg[4]))
outputfile <- as.character(arg[5])
markerfile <- as.character(arg[6])

print(markerlist); print(moltrait); print(e); print(n); print(outputfile); print(markerfile)

library(data.table); library(GenomicRanges);library(dplyr)

print("loading marker set, as set in argument nr. 6")

hq <- fread(markerfile)
hq.gr <- GRanges(seqnames=hq$Chrom, IRanges(start=hq$Pos, end=hq$Pos))
hq.maf <- hq[,.(Name,ImpMAF)]

mcols(hq.gr) <- hq.maf

print("loading markerlist, as set in argument nr. 1")

m <- fread(markerlist); if(colnames(m)[2]=="Start") colnames(m)[2] <- "Pos"

m <-left_join(m,hq.maf,by="Name"); m <- m[is.na(ImpMAF)==FALSE,]

amplifyBy <- round((n/nrow(m))+.5)

print("loading genomic attributes file, as set in argument nr. 2")
mol <- fread(moltrait)

newstart <- mol$Start - e; newstart[newstart<0] <- 0; newend <- mol$End + e
gr_extended <- GRanges(seqnames = mol$Chrom, IRanges(start=newstart, end=newend))
gr_extended$id <- paste0(seqnames(gr_extended), ":", start(gr_extended), "-", end(gr_extended))

f <- mergeByOverlaps(hq.gr,gr_extended)

ff <- data.table(Name=f$Name,ImpMAF=f$ImpMAF)
ff <- ff[duplicated(Name)==FALSE, ]

mafbreaks <- c(0,0.005,0.01,0.05,0.10,0.20,0.30,0.4,0.5)

m[,MAFbin:=cut(ImpMAF, breaks=mafbreaks)]
ff[,MAFbin:=cut(ImpMAF, breaks=mafbreaks)]

mafBins <- m[,.N,MAFbin]; mafBins[,N:=N*amplifyBy]

ff.c <- NULL

for(i in 1:nrow(mafBins)){
  
  howmany <- mafBins[i, N]; binid <- mafBins[i, MAFbin]
  
  fff <- ff[MAFbin%in%binid, ]
  s <- sample(1:nrow(fff), howmany, replace=FALSE)
  
  ff.c <- rbindlist(list(ff.c, fff[s, ]))
  
}

ff.c <- ff.c[1:n, ]

o <- data.table(Chrom=unlist(tstrsplit(ff.c$Name,":",keep=1)), Pos=unlist(tstrsplit(ff.c$Name,":",keep=2)), Name=ff.c$Name)

options(scipen=999)
write.table(o, file=outputfile, quote=FALSE, sep="\t", row.names=FALSE)
