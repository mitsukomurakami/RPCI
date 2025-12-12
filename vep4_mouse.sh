#!/bin/bash
# vep4_mouse.sh
# - produces .vep.vcf for strelka2_mutect2_muse_vardict output vcf file

#pid=`basename $PWD`

if [[ $# -eq 0 ]] ; then
    echo "Usage: ./vep4_mouse.sh <mouseID>"
    exit 1
fi
pid=$1

PREFIX=/projects/rpci/songliu/mkorobkin
PID=RQ026467-Huss

JOBDIR=$PREFIX/projects/$PID/$pid
OUTDIR=$JOBDIR/output/variants
SRCDIR=$JOBDIR/src
if [ ! -d "$OUTDIR" ]; then
  mkdir -p $OUTDIR
fi
if [ ! -d "$SRCDIR" ]; then
  mkdir -p $SRCDIR
fi

rfile=/projects/rpci/shared/references/GRCm38/Mus_musculus.GRCm38.dna.toplevel.fa  ## mouse

# Write an R script
cat > $SRCDIR/vep4_mouse.R << EOF
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
cwlInstall("tl_vep")
inputs(vep)
devtools::load_all("/projects/rpci/songliu/qhu/workspace/projects/VariantCombiner/")  

# Prepare inputs
vcfs <- list.files("$OUTDIR", "*strelka2_mutect2_muse_vardict.vcf", recursive = TRUE, full.names = TRUE)
ids <- basename(dirname(vcfs))
ovcf <- sub(".vcf", ".vep.vcf", basename(vcfs))
vcfs <- as.list(vcfs)
ovcf <- as.list(ovcf)
names(vcfs) <- names(ovcf) <- ids
vcfs_done <- list.files("$OUTDIR", "*strelka2_mutect2_muse_vardict.vep.vcf", recursive = TRUE, full.names = TRUE)
idx <- !ids %in% basename(dirname(vcfs_done))

## more inputs
ref <- "$rfile"

inputList <- list(ivcf = vcfs[idx],
                  ovcf = ovcf[idx])
paramList <- list(ref = ref,
                  cacheDir = "/projects/rpci/songliu/qhu/annotation/vep",
                  species = "mus_musculus",
                  version = "102")

inputList
paramList

# Submit a job using runCWLBatch() and BatchtoolsParam() functions
runCWLBatch(vep, outdir = "$OUTDIR", inputList, paramList,
            BPPARAM = BatchtoolsParam(workers = lengths(inputList)[1],
                                      cluster = "slurm",
                                      #template = "$PREFIX/rcwl/slurm_rpci.tmpl",  ## RPCI cluster
                                      #template = "$PREFIX/rcwl/slurm_nih.tmpl", ## ub-hpc:general-compute
                                      template = "$PREFIX/rcwl/slurm_nih_2.tmpl", ## ub-hpc:scavenger
                                      resources = list(ncpus = 1,
                                                       jobname = "vep4_mouse",
                                                       walltime = 60*60*3,
                                                       memory = 16000),
                                      log = TRUE, logdir = "$JOBDIR",
                                      progressbar = FALSE,
                                      saveregistry = FALSE),
            stderr = "", docker="singularity")
EOF

#/user/mkorobkin/rpci/R/bin/Rscript $SRCDIR/vep4_mouse.R & 
Rscript2 $SRCDIR/vep4_mouse.R &

