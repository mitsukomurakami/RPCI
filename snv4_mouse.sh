#!/usr/bin/sh
# snv4_mouse.sh

if [[ $# -eq 0 ]] ; then
    echo "Usage: ./snv4_mouse.sh <mouseID>"
    exit 1
fi
pid=$1

PREFIX=/projects/rpci/songliu/mkorobkin
PID=RQ024165-Paragh_1

JOBDIR=$PREFIX/projects/$PID/$pid
OUTDIR=$JOBDIR/output
SRCDIR=$JOBDIR/src
if [ ! -d "$SRCDIR" ]; then
  mkdir -p $SRCDIR
fi

rfile=/projects/rpci/shared/references/GRCm38/Mus_musculus.GRCm38.dna.toplevel.fa ## mouse
rgfile=$PREFIX/projects/$PID/PGD843.coveredRegion_Mouse_clean.bed

## output file for tumor-normal pairs
pfile=$OUTDIR/TNpair.txt
if [ -f "$pfile" ]; then
    rm $pfile
fi

## normal sample
Npath=`ls -d $OUTDIR/BAM/${pid}-M001-*`
Npub=`basename $Npath`
nfile=`ls $Npath/*bam`

## tumor samples
TLIST=`ls -d $OUTDIR/BAM/${pid}-Sp*`
for Tpath in $TLIST; do
  Tpub=`basename $Tpath`
  tfile=`ls $Tpath/*bam`
  echo $Tpub $tfile $Npub $nfile >> $pfile
done

# Write an R script
cat > $SRCDIR/snv4_mouse.R << EOF
library(RcwlPipelines)
library(BiocParallel)
library(basilisk)
Sys.setenv(SINGULARITY_TMPDIR="$PREFIX/tmp")
Sys.setenv(SINGULARITY_CACHEDIR="$PREFIX/singularity_cache")
Sys.setenv(SINGULARITY_LOCALCACHEDIR="$PREFIX/tmp/runtime")
Sys.setenv(DEBUGME = "BiocParallel")

# for cwltool
env_Rcwl = Rcwl:::env_Rcwl
if(!file.exists(Sys.which("cwltool"))){
        cl <- basiliskStart(env_Rcwl)
        on.exit(basiliskStop(cl))
}
binPath <- Sys.which("cwltool")
binPath

# Load pipeline
cwlInstall("pl_SomaticCaller_mouse")
inputs(SomaticCaller_mouse)

# Prepare inputs
df <- data.frame(read.table("$pfile"))
tids <- as.list(df\$V1)
tbam <- as.list(df\$V2)  
nids <- as.list(df\$V3)
nbam <- as.list(df\$V4)  
ref <- "$rfile"
names(tbam) <- as.list(df\$V1)

# Read .bed file
regions <- "$rgfile"

idx=TRUE
inputList <- list(tbam = tbam[idx],
                  nbam = nbam[idx],
                  tumor = tids[idx],
                  normal = nids[idx])
paramList <- list(Ref = ref,
                  interval = regions,
                  dbsnp = "/projects/rpci/songliu/qhu/annotation/Mouse/00-All.vcf.gz",
                  threads = 16)

# Submit a job using runCWLBatch() and BatchtoolsParam() functions
runCWLBatch(SomaticCaller_mouse, outdir = "$OUTDIR/variants", inputList, paramList,
            BPPARAM = BatchtoolsParam(workers = lengths(inputList)[1],
                                      cluster = "slurm",
                                      #template = "$PREFIX/rcwl/slurm_rpci.tmpl",
                                      template = "$PREFIX/rcwl/slurm_nih.tmpl",
                                      resources = list(ncpus = 16,
                                                       jobname = "snv4_mouse",
                                                       walltime = 60*60*36,
                                                       memory = 32000),
                                      log = TRUE, logdir = "$JOBDIR",
                                      progressbar = FALSE,
                                      saveregistry = TRUE),
            stderr = "", docker="singularity")
EOF

Rscript2 $SRCDIR/snv4_mouse.R &
