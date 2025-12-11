#!/bin/bash
# summary2.sh

#module use /projects/rpci/shared/modulefiles
#module load gcc mosdepth
#module load openmpi
#module load samtools        ## for bgzip

PREFIX=/projects/rpci/songliu/mkorobkin
PID=RQ026467-Huss

JOBDIR=$PREFIX/projects/$PID
SRCDIR=$JOBDIR/src

xxx=CMm 
#ct=10000,20000,50000,100000  ## coverage thresholds
ct=3000,5000,8000,10000 
#ct=10,100,200,300

PLIST=`ls -d *${xxx}*`
for pid in $PLIST; do
  OUTDIR=$JOBDIR/$pid/output
  SLIST=`ls $OUTDIR/BAM`
  for sid in $SLIST; do
    echo $sid:
    ofile=$OUTDIR/coverage2/$sid/$sid.sample_summary
    echo sampleID total mean $ct | sed 's/,/X /g' | sed 's/$/X/' > $ofile
    msfile=$OUTDIR/coverage2/$sid/$sid.mosdepth.summary.txt
    total=`tail -1 $msfile | awk '{print $3}'`
    mean=`tail -1 $msfile | awk '{print $4}'`
    bfile=$OUTDIR/coverage2/$sid/$sid.thresholds.bed
    if [ ! -f $bfile ]; then
      if [ -f $bfile.gz ]; then
        bgzip -d $bfile.gz  ## de-compress bed.gz file
      else
        echo Error: $bfile.gz does not exist 
        exit 1
      fi
    fi
    sed 1d $bfile | awk -v sid=$sid -v total=$total -v mean=$mean '{sum+=$3-$2;sum1+=$5;sum2+=$6;sum3+=$7;sum4+=$8;}END{print sid,total,mean,sum1/sum*100,sum2/sum*100,sum3/sum*100,sum4/sum*100;}' >> $ofile
  done
done

ofile=$JOBDIR/${PID}_reads_summary2.csv
if [ -f $ofile ]; then
  rm $ofile
fi

# Write an R script
cat > $SRCDIR/summary2.R << EOF
pids <- dir("$JOBDIR", "^.$xxx", recursive =FALSE)
pids

for(j in 1:length(pids)){

prefix <- paste0("$JOBDIR/", pids[j], "/output/BAM")
bams <- list.files(prefix, ".bam$", recursive = TRUE, full.names = TRUE) 
sids <- basename(dirname(bams))
sids

## Write summary
reads <- c()
for(i in 1:length(sids)){
    r1 <- readLines(list.files(paste0("$JOBDIR/",pids[j],"/output/BAM/", sids[i]), "*.flagstat.txt", full.names = T))
    r1 <- as.numeric(sub("\\\s.*", "", r1))
    tt <- r1[10]*2
    r1 <- c(tt, r1[7]/r1[1], r1[5]/r1[1])
    c1 <- read.csv(paste0("$JOBDIR/",pids[j],"/output/coverage2/",sids[i],"/", sids[i], ".sample_summary"), check.names = F, header = TRUE, sep = "")[1,3:7]
    r1 <- cbind(rbind(r1), c1)
    reads <- rbind(reads, r1)
}
rownames(reads) <- sids
if(j==1) {
  colnames(reads)[1:3] <- c("total_reads", "map_rate", "dup_rate")
  write.table(reads, file="$ofile", col.names = NA, row.names = TRUE, sep = ",")
} else {
  write.table(reads, file="$ofile", append = TRUE, col.names = FALSE, row.names = TRUE, sep = ",")
}

}
message(paste0("Results written in: ", "$ofile"))
EOF

Rscript2 $SRCDIR/summary2.R &
