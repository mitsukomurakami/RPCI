#!/bin/bash
# coverage2.sh

if [[ $# -eq 0 ]] ; then
    echo "Usage: ./coverage2.sh <patientID>"
    exit 1
fi
pid=$1

module load gcc mosdepth

PREFIX=/projects/rpci/songliu/mkorobkin
PID=RQ026467-Huss

JOBDIR=$PREFIX/projects/$PID/$pid
DATADIR=$JOBDIR/data
OUTDIR=$JOBDIR/output

rgfile=$PREFIX/projects/$PID/PGD943.coveredRegion_clean.bed

SRCDIR=$JOBDIR/src
if [ ! -d "$SRCDIR" ]; then
  mkdir -p $SRCDIR
fi

# Write an R script
cat > $SRCDIR/coverage2.R << EOF
library(RcwlPipelines)
library(BiocParallel)
Sys.setenv(SINGULARITY_TMPDIR="$PREFIX/tmp")
Sys.setenv(SINGULARITY_CACHEDIR="$PREFIX/singularity_cache")
Sys.setenv(SINGULARITY_LOCALCACHEDIR="$PREFIX/tmp/runtime")
Sys.setenv(DEBUGME = "BiocParallel")

## Define "mosdepth" command
i1<-InputParam(id="bedfile", type="File", prefix="--by")
i2<-InputParam(id="ct", type="string", prefix="--thresholds")
i3<-InputParam(id="fileID", type="string", position=1)
i4<-InputParam(id="bamfile", type="File", position=2, secondaryFile=".bai")
o1 <- OutputParam(id = "out", type="File[]", glob="\$(inputs.fileID)*")  
#req1 <- requireDocker("quay.io/biocontainers/mosdepth:0.3.3--h37c5b7d_2")
mosdepth <- cwlProcess(baseCommand="mosdepth",  
                       #requirements = list(req1),
                       inputs=InputParamList(i1,i2,i3,i4),
                       outputs = OutputParamList(o1))
mosdepth

arguments(mosdepth, step = NULL) <- list("--fast-mode", "--no-per-base")  ## neglects mate-pair overlap and cigar operations

## Prepare inputs
bams <- list.files("$OUTDIR/BAM", ".bam$", recursive = TRUE, full.names = TRUE) 
sids <- basename(dirname(bams))
bams <- as.list(bams)
names(bams) <- sids
#coverageThreshold <- "10000,20000,50000,100000"  ## no spaces in between
coverageThreshold <- "3000,5000,8000,10000"  ## no spaces in between
#coverageThreshold <- "1,10,20,30"  ## 

cov_done <- list.files("$OUTDIR/coverage2", "*.thresholds.bed.gz", recursive = TRUE)
idx <- !sids %in% basename(dirname(cov_done))
#idx <- TRUE

inputList <- list(bamfile = bams[idx],
                  fileID = as.list(names(bams))[idx])
paramList <- list(bedfile = "$rgfile",
                  ct = coverageThreshold)
inputList
paramList

## Submit a job using runCWLBatch() and BatchtoolsParam() functions
runCWLBatch(mosdepth, outdir = "$OUTDIR/coverage2",
            inputList, paramList,
            BPPARAM = BatchtoolsParam(workers = lengths(inputList)[1],
                                      cluster = "slurm",
                                      #template = "$PREFIX/rcwl/slurm_rpci.tmpl",
                                      #template = "$PREFIX/rcwl/slurm_nih.tmpl",
                                      template = "$PREFIX/rcwl/slurm_nih_2.tmpl",  ## scavenger
                                      resources = list(ncpus = 1,
                                                       jobname = "cov2",
                                                       walltime = 60*60*3,
                                                       memory = 16000),
                                      log = TRUE, logdir = "$JOBDIR",
                                      progressbar = FALSE,
                                      saveregistry = TRUE),
                                      stderr = "", docker="singularity")
EOF

#/user/mkorobkin/rpci/R/bin/Rscript $SRCDIR/coverage2.R & 
Rscript2 $SRCDIR/coverage2.R & 
