#!/bin/bash
# align.sh
# - for target sequencing 

if [[ $# -eq 0 ]] ; then
    echo "Usage: ./align.sh <patientID>"
    exit 1
fi
pid=$1

PREFIX=/projects/rpci/songliu/mkorobkin
PID=RQ026467-Huss

JOBDIR=$PREFIX/projects/$PID/$pid
DATADIR=$JOBDIR/data
OUTDIR=$JOBDIR/output

## Reference genome
#rfile=/projects/rpci/shared/reference/human_g1k_v37/human_g1k_v37.fa  ## human
rfile=/projects/rpci/shared/references/GRCm38/Mus_musculus.GRCm38.dna.toplevel.fa   ## mouse

SRCDIR=$JOBDIR/src
mkdir -p $SRCDIR

# Write an R script
cat > $SRCDIR/align.R << EOF
library(RcwlPipelines)
library(BiocParallel)
Sys.setenv(SINGULARITY_TMPDIR="$PREFIX/tmp")
Sys.setenv(SINGULARITY_CACHEDIR="$PREFIX/singularity_cache")
Sys.setenv(SINGULARITY_LOCALCACHEDIR="$PREFIX/tmp/runtime")
Sys.setenv(DEBUGME = "BiocParallel")

# Load pipeline
cwlInstall("pl_bwa_align")
inputs(bwa_align)

# Prepare inputs
fqs <- list.files("$DATADIR", ".fastq.gz", full.names = TRUE, recursive = TRUE)
fq1 <- fqs[grep("R1", fqs)]
fq2 <- fqs[grep("R2", fqs)]
fqs.ids <- do.call(rbind, lapply(strsplit(basename(fq1), split = "[_.]"), function(x)x[1:6]))
fqs.ids[,2] <- gsub("-[A-Z][a-z]*$","",fqs.ids[,2])
sids <- basename(dirname(fq1))

RGs <- paste("@RG",
             paste0("ID:", fqs.ids[,1], ".", fqs.ids[,3], ".", fqs.ids[,2]),
             paste0("LB:", fqs.ids[,1], ".", fqs.ids[,3], ".", fqs.ids[,2]),
             paste0("DT:", Sys.Date()),
             paste0("PL:", "Illumina"),
             "CN:RPCCC",
             paste0("SM:", sids), sep = "\\\t")
RGs

FID <- paste0(sids, "_", fqs.ids[,1], "_", fqs.ids[,3], "_", fqs.ids[,2], ".bam")  ## [PubID]_[RSIDs].bam
#FID

fq1L <- tapply(fq1, sids, as.list)
fq2L <- tapply(fq2, sids, as.list)
RGL <- tapply(RGs, sids, as.list)
idBam <- tapply(FID, sids, function(x)as.list(unique(x)))

done <- list.files("$OUTDIR/BAM", ".bam$", recursive = TRUE)
idx <- !sids %in% basename(dirname(done))
#idx <- TRUE
inputList <- list(RG = RGL[idx],
                  outBam = idBam[idx],
                  FQ1 = fq1L[idx],
                  FQ2 = fq2L[idx])

paramList <- list(threads = 16,
                  Ref ="$rfile")

inputList

# Submit a job using runCWLBatch() and BatchtoolsParam() functions
runCWLBatch(bwa_align, outdir = "$OUTDIR/BAM", inputList, paramList,
            BPPARAM = BatchtoolsParam(workers = lengths(inputList)[1],
                                      cluster = "slurm",
                                      #template = "$PREFIX/rcwl/slurm_rpci.tmpl", ## For RPCI cluster
                                      template = "$PREFIX/rcwl/slurm_nih.tmpl", ## For up-hpc:general-compute
                                      resources = list(ncpus = 16,
                                                       jobname = "align",
                                                       walltime = 60*60*6,
                                                       memory = 128000),
                                      log = TRUE, logdir = "$JOBDIR",
                                      progressbar = FALSE,
                                      saveregistry = FALSE),
                                      stderr = "", docker="singularity")

EOF

#/user/mkorobkin/rpci/R/bin/Rscript $SRCDIR/align.R &
Rscript2 $SRCDIR/align.R &
