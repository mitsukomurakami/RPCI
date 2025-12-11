#!/usr/bin/sh
# alignDup.sh
# - For WES (whole exome sequencing) data

if [[ $# -eq 0 ]] ; then
    echo "Usage: ./alignDup.sh <patientID>"
    exit 1
fi
pid=$1

#pid=`basename $PWD`

PREFIX=/projects/rpci/songliu/mkorobkin
PID=RQ023608-Wei

JOBDIR=$PREFIX/projects/$PID/$pid
DATADIR=$JOBDIR/data2  ## trimmed
OUTDIR=$JOBDIR/output2

SRCDIR=$JOBDIR/src
mkdir -p $SRCDIR

## Reference genome
rfile=/projects/rpci/shared/reference/human_g1k_v37/human_g1k_v37.fa  ## human
#rfile=/projects/rpci/shared/references/GRCm38/Mus_musculus.GRCm38.dna.toplevel.fa   ## mouse

# Write an R script
cat > $SRCDIR/alignDup.R << EOF
library(RcwlPipelines)
library(BiocParallel)
Sys.setenv(SINGULARITY_TMPDIR="$PREFIX/tmp")
Sys.setenv(SINGULARITY_CACHEDIR="$PREFIX/singularity_cache")
Sys.setenv(SINGULARITY_LOCALCACHEDIR="$PREFIX/tmp/runtime")
Sys.setenv(DEBUGME = "BiocParallel")

# Load pipeline
cwlInstall("pl_bwaDup")
inputs(bwaDup)
arguments(bwaDup, "markdup") <- list("-Xmx32g")  ## for java heapspace with Picard MarkDuplicates

# Prepare inputs
fqs <- list.files("$DATADIR", ".fastq.gz", full.names = TRUE, recursive = TRUE)
fq1 <- fqs[grep("R1", fqs)]
fq2 <- fqs[grep("R2", fqs)]
fqs.ids <- do.call(rbind, lapply(strsplit(basename(fq1), split = "[_.]"), function(x)x[1:6]))
rs <- fqs.ids[grep("RS",fqs.ids)]
dim(rs) <- c(nrow(fqs.ids),2)
sids <- basename(dirname(fq1))

fqs.ids
rs

RGs <- paste("@RG",
             paste0("ID:", rs[,1], ".", rs[,2]),
             paste0("LB:", rs[,1], ".", rs[,2]),
             paste0("DT:", Sys.Date()),
             paste0("PL:", "Illumina"),
             "CN:RPCCC",
             paste0("SM:", sids), sep = "\\\t")
RGs

FID <- paste0(sids, "_", rs[,1], "_", rs[,2], ".bam")
FID

fq1L <- tapply(fq1, sids, as.list, simplify = FALSE)
fq2L <- tapply(fq2, sids, as.list, simplify = FALSE)
RGL <- tapply(RGs, sids, as.list, simplify = FALSE)
idBam <- tapply(FID, sids, function(x)as.list(unique(x)))

done <- list.files("$OUTDIR/BAM", ".bam$", recursive = TRUE)
idx <- !sids %in% dirname(done)
#idx <- TRUE
inputList <- list(RG = RGL[idx],
                  outBam = idBam[idx],
                  FQ1s = fq1L[idx],
                  FQ2s = fq2L[idx])

paramList <- list(threads = 16,
                  Ref ="$rfile")

# Submit a job using runCWLBatch() and BatchtoolsParam() functions
runCWLBatch(bwaDup, outdir = "$OUTDIR/BAM", inputList, paramList,
            BPPARAM = BatchtoolsParam(workers = lengths(inputList)[1],
                                      cluster = "slurm",
                                      template = "$PREFIX/rcwl/slurm_rpci.tmpl", ## for RPCI clusters
                                      #template = "$PREFIX/rcwl/slurm_nih.tmpl", ## for ub-hpc:general-compute
                                      resources = list(ncpus = 16,
                                                       jobname = "alignDup",
                                                       walltime = 60*60*24,
                                                       memory = 128000),
                                      log = TRUE, logdir = "$JOBDIR",
                                      progressbar = FALSE,
                                      saveregistry = TRUE),
                                      stderr = "", docker="singularity")

EOF

#/user/mkorobkin/rpci/R/bin/Rscript $SRCDIR/align.R &
Rscript2 $SRCDIR/alignDup.R &
