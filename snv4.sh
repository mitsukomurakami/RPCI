#!/usr/bin/sh
# snv4.sh

if [[ $# -eq 0 ]] ; then
    echo "Usage: ./snv4.sh <patientID>"
    exit 1
fi
pid=$1

PREFIX=/projects/rpci/songliu/mkorobkin
PID=RQ023532-Paragh

JOBDIR=$PREFIX/projects/$PID/$pid
OUTDIR=$JOBDIR/output
#OUTDIR=$JOBDIR/output2      ## if trimmed
SRCDIR=$JOBDIR/src
if [ ! -d "$SRCDIR" ]; then
  mkdir -p $SRCDIR
fi

rfile=/projects/rpci/shared/reference/human_g1k_v37/human_g1k_v37.fa  ## human
rgfile=$PREFIX/projects/$PID/PGD770.coveredRegion_HumanWide_hg19_clean_2.bed

## output file for tumor-normal pairs
pfile=$OUTDIR/TNpair.txt
if [ -f "$pfile" ]; then
    rm $pfile
fi

## normal sample
Npath=`ls -d $OUTDIR/BAM/${pid}-Sd*`
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
cat > $SRCDIR/snv4.R << EOF
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
cwlInstall("pl_SomaticCaller4")
inputs(SomaticCaller4)
devtools::load_all("/projects/rpci/songliu/qhu/workspace/projects/VariantCombiner/")

# Prepare inputs
ref <- "$rfile"
df <- data.frame(read.table("$pfile"))
tids <- as.list(df\$V1)
tbam <- as.list(df\$V2)  
nids <- as.list(df\$V3)
nbam <- as.list(df\$V4)  
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
                  dbsnp = "/projects/rpci/songliu/qhu/annotation/GATK_bundle/b37/dbsnp_138.b37.vcf.gz",
                  gresource = "/projects/rpci/songliu/qhu/annotation/Mutect2/af-only-gnomad.raw.sites.b37.vcf",
                  comvcf = "/projects/rpci/songliu/qhu/annotation/Mutect2/GetPileupSummaries/small_exac_common_3_b37.vcf",
                  pon = "/projects/rpci/songliu/qhu/annotation/Mutect2/Mutect2-exome-panel.vcf",
                  threads = 16)

# Submit a job using runCWLBatch() and BatchtoolsParam() functions
runCWLBatch(SomaticCaller4, outdir = "$OUTDIR/variants", inputList, paramList,
            BPPARAM = BatchtoolsParam(workers = lengths(inputList)[1],
                                      cluster = "slurm",
                                      #template = "$PREFIX/rcwl/slurm_rpci.tmpl", ## For RPCI cluster
                                      template = "$PREFIX/rcwl/slurm_nih.tmpl", ## For up-hpc:general-compute
                                      resources = list(ncpus = 16,
                                                       jobname = "snv4",
                                                       walltime = 60*60*36,
                                                       memory = 32000),
                                      log = TRUE, logdir = "$JOBDIR",
                                      progressbar = FALSE,
                                      saveregistry = TRUE),
            stderr = "", docker="singularity")
EOF

#/user/mkorobkin/rpci/R/bin/Rscript $SRCDIR/snv4.R &
Rscript2 $SRCDIR/snv4.R &
