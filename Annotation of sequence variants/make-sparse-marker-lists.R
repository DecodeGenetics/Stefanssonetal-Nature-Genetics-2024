args <- commandArgs(TRUE)

# args1: marker file (Chrom, Pos, Name)
# args2: minimum distance between markers (e.g. 1e+6)
# args3: logical; TRUE/FALSE: is column #4 = -log10(P-value) ? ordering according to significance of markers.  
# args4: path to output file

# example run:
# $ Rscript /.../make-sparse-marker-lists.R /.../markerfile.gor 1e+6 FALSE /.../markerfile-per1Mb.gor

library(data.table); options(scipen=999)

markerlist <- as.character(args[1])

minDist <- as.numeric(as.character(args[2]))

pval <- as.logical(as.character(args[3]))

outputname <- as.character(args[4])

p <- fread(markerlist)

if(pval) colnames(p)[4] <- "P"

u <- unique(p$Chrom)

pin.c <- NULL

for(i in 1:length(u)){
  
  pp <- p[Chrom%in%u[i], ]
  
  if(pval) pp <- pp[order(-P),]  

  pin <- pp[1, ]
  
  if(nrow(pp)>1){  
    
  for(w in 2:nrow(pp)){
    
    pc <- pp[w, ]
    
    len <- abs(pc$Pos-pin[,Pos])
    
    if(all(len>minDist)) pin <- rbindlist(list(pin,pc))
    
  }
    
  }
  
  pin.c <- rbindlist(list(pin.c, pin))
  
}

write.table(pin.c, file=outputname, quote=FALSE, sep="\t", row.names=FALSE)
