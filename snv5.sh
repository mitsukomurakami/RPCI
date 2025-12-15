#!/usr/bin/sh
# snv5.sh
# - Bambino: https://github.com/NCIP/cgr-bambino

if [[ $# -eq 0 ]] ; then
    echo "Usage: ./snv5.sh <patientID>"
    exit 1
fi
pid=$1

PREFIX=/projects/rpci/songliu/mkorobkin
PID=RQ023532-Paragh

JOBDIR=$PREFIX/projects/$PID/$pid
OUTDIR=$JOBDIR/output
SRCDIR=$JOBDIR/src
if [ ! -d "$SRCDIR" ]; then
  mkdir -p $SRCDIR
fi

rfile=/projects/rpci/shared/reference/human_g1k_v37/human_g1k_v37.fa  ## human
#rfile=/projects/rpci/shared/references/GRCm38/Mus_musculus.GRCm38.dna.toplevel.fa  ## mouse

## output file for tumor-normal pairs
pfile=$OUTDIR/TNpairBambino.txt
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
  mkdir -p $OUTDIR/variants2/$Tpub
  tfile=`ls $Tpath/*bam`
  ofile=${Tpub}_${Npub}_bambino
  echo $Tpub -- $Npub
  echo $Tpub $tfile $nfile $ofile >> $pfile
done

# Write an R script
cat > $SRCDIR/snv5.R << EOF
library(RcwlPipelines)
library(BiocParallel)
Sys.setenv(SINGULARITY_TMPDIR="$PREFIX/tmp")
Sys.setenv(SINGULARITY_CACHEDIR="$PREFIX/singularity_cache")
Sys.setenv(SINGULARITY_LOCALCACHEDIR="$PREFIX/tmp/runtime")
Sys.setenv(DEBUGME = "BiocParallel")

#######################
# Define a cwlProcess for bambino
#####
req1 <- requireDocker("hubentu/bambino")
p1 <- InputParam(id = "dbam", type = "File", position=1, secondaryFile=".bai")
p2 <- InputParam(id = "gbam", type = "File", position=2, secondaryFile=".bai")
p3 <- InputParam(id = "out", type = "string", position=3)
p4 <- InputParam(id = "ref", type = "File", position=4, secondaryFile=".fai")
o1 <- OutputParam(id = "vout", type = "File", glob = "\$(inputs.out)")
bambino <- cwlProcess(baseCommand = "/opt/run.sh",
                      requirements = list(req1),
                      inputs = InputParamList(p1, p2, p3, p4),
                      outputs = OutputParamList(o1))
#######################
inputs(bambino)

# Prepare inputs
df <- data.frame(read.table("$pfile"))
dbam <- as.list(df\$V2)  
gbam <- as.list(df\$V3)  
names(dbam) <- as.list(df\$V1)
out <- as.list(df\$V4)

idx=TRUE
inputList <- list(dbam = dbam[idx],
                  gbam = gbam[idx],
                  out = out[idx])
paramList <- list(ref = "$rfile")

# Submit a job using runCWLBatch() and BatchtoolsParam() functions
runCWLBatch(bambino, outdir = "$OUTDIR/variants2", inputList, paramList,
            BPPARAM = BatchtoolsParam(workers = lengths(inputList)[1],
                                      cluster = "slurm",
                                      #template = "$PREFIX/rcwl/slurm_rpci.tmpl",  ## RPCI cluster
                                      #template = "$PREFIX/rcwl/slurm_nih.tmpl",  ## ub-hpc:general-compute
                                      template = "$PREFIX/rcwl/slurm_nih_2.tmpl", ## up-hpc:scavenger
                                      resources = list(ncpus = 1,
                                                       jobname = "snv5",
                                                       walltime = 60*60*48,
                                                       memory = 64000),
                                      log = TRUE, logdir = "$JOBDIR",
                                      progressbar = FALSE,
                                      saveregistry = FALSE),
            stderr = "", docker="singularity")

# Check results
list.files("$OUTDIR/variants2/", recursive = TRUE, pattern = "bambino$")
EOF

#/user/mkorobkin/rpci/R/bin/Rscript $SRCDIR/snv5.R &
Rscript2 $SRCDIR/snv5.R &
