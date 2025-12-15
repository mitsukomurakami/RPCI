#!/usr/bin/sh
# bvlist.sh

PREFIX=/projects/rpci/songliu/mkorobkin
PID=RQ023532-Paragh

JOBDIR=$PREFIX/projects/$PID
SRCDIR=$JOBDIR/src

# Write an R script
cat > $SRCDIR/bvlist.R << EOF
pids <- dir("$JOBDIR", "^A", recursive =FALSE)
pids

for(j in 1:length(pids)){

## bam files
bdir <- paste0("$JOBDIR/", pids[j], "/output/BAM")
bams <- list.files(bdir, ".bam$", recursive = TRUE, full.names = TRUE)
sids <- basename(dirname(bams))

## sj files
vdir <- paste0("$JOBDIR/", pids[j], "/output/variants")
sjs <- list.files(vdir, "*strelka2_mutect2_muse_vardict.sj", recursive = TRUE, full.names = TRUE)  
vdir2 <- paste0("$JOBDIR/", pids[j], "/output/variants2")
sjs2 <- list.files(vdir2, "*bambino.sj", recursive = TRUE, full.names = TRUE)  
sjs <- append(sjs,sjs2)
sjs

sjs.ids <- do.call(rbind, lapply(strsplit(basename(sjs), split = "[_.]"), function(x)x[1:2]))
tids <- sjs.ids[,1] 
nids <- sjs.ids[,2]

for(t in c("strelka2", "mutect2", "muse", "vardict", "bambino")){
    if(t=="bambino") { 
      prefix <- vdir2 
    } else {
      prefix <- vdir
    }
    bv <- cbind(Sample=as.character(tids),
                Reference="hg19",
                dbam=normalizePath(bams[match(tids, sids)]),
                gbam=normalizePath(bams[match(nids, sids)]),
                GID=as.character(nids),
                PID=as.character(pids[j]),
                Variants=paste0(prefix, "/",tids,"/", tids,"_",nids, "_", t, ".sj"))
                )

    ex <- rbind(apply(rbind(bv), 2, file.exists))
    ex[,c(1,2,5,6)] <- TRUE
    bv[!ex] <- NA

  ofile <- paste0("$JOBDIR/${PID}_bam_trimmed_variant_list_", t)  ## trimmed
  if(j==1) {
    write.table(bv, file=ofile, append = FALSE, quote=FALSE, row.names=FALSE, sep = "\t")
  } else {
    write.table(bv, file=ofile, append = TRUE, quote=FALSE, col.names = FALSE, row.names=FALSE, sep = "\t")
  }

}

}
# Check results
list.files("$JOBDIR", "*.bam_trimmed_variant_list_*", recursive = FALSE, full.names = TRUE)
EOF

#Rscript2 $SRCDIR/bvlist.R &
/user/mkorobkin/rpci/R/bin/Rscript $SRCDIR/bvlist.R &
