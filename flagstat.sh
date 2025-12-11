#!/bin/bash
# flagstat.sh
# - produces .flagstat.txt from .bam files

if [[ $# -eq 0 ]] ; then
    echo "Usage: ./flagstat.sh <patientID>"
    exit 1
fi
pid=$1

PREFIX=/projects/rpci/songliu/mkorobkin
PID=RQ026467-Huss

JOBDIR=$PREFIX/projects/$PID/$pid
OUTDIR=$JOBDIR/output

SRCDIR=$JOBDIR/src
if [ ! -d "$SRCDIR" ]; then
  mkdir -p $SRCDIR
fi

FLIST=`ls $OUTDIR/BAM/*/*flagstat.txt`
for file in $FLIST; do
  if [ ! -s $file ]; then
    echo $file: empty
    rm $file
  fi
done

# Write an R script
cat > $SRCDIR/flagstat.R << EOF
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
cwlInstall("tl_samtools_flagstat")
inputs(samtools_flagstat)

# Prepare inputs
bams <- list.files("$OUTDIR/BAM", ".bam$", recursive = TRUE, full.names = TRUE) 
bn <- basename(dirname(bams))
names(bams) <- as.list(bn)
done <- list.files("$OUTDIR/BAM", "*.flagstat.txt", recursive = TRUE)
idx <- !bn %in% basename(dirname(done))
#idx <- TRUE
inputList <- list(bam = bams[idx])
paramList <- list()

inputList
paramList

# Submit a job using runCWLBatch() and BatchtoolsParam() functions
runCWLBatch(samtools_flagstat, outdir = "$OUTDIR/BAM", 
            inputList, paramList,
            BPPARAM = BatchtoolsParam(workers = lengths(inputList)[1],
                                      cluster = "slurm",
                                      #template = "/projects/rpci/songliu/mkorobkin/rcwl/slurm_rpci.tmpl",
                                      #template = "/projects/rpci/songliu/mkorobkin/rcwl/slurm_nih.tmpl",
                                      template = "/projects/rpci/songliu/mkorobkin/rcwl/slurm_nih_2.tmpl",
                                      resources = list(ncpus = 1,
                                                       jobname = "flagstat",
                                                       walltime = 60*60*12,
                                                       memory = 32000),
                                      log = TRUE, logdir = "$JOBDIR",
                                      progressbar = FALSE,
                                      saveregistry = FALSE),
            stderr = "", docker="singularity")
EOF

#/user/mkorobkin/rpci/R/bin/Rscript $SRCDIR/flagstat.R & 
Rscript2 $SRCDIR/flagstat.R &
