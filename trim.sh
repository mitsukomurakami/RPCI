#!/usr/bin/sh
# trim.sh

if [[ $# -eq 0 ]] ; then
    echo "Usage: ./trim.sh <patientID>"
    exit 1
fi
pid=$1

PREFIX=/projects/rpci/songliu/mkorobkin  ## Your project directory
PID=RQ026467-Huss  ## Project ID

JOBDIR=$PREFIX/projects/$PID/$pid
DATADIR=$JOBDIR/data
OUTDIR=$JOBDIR/data2 

SRCDIR=$JOBDIR/src
if [ ! -d "$SRCDIR" ]; then
  mkdir -p $SRCDIR
fi

# Write an R script
cat > $SRCDIR/trim.R << EOF
library(RcwlPipelines)
library(BiocParallel)
Sys.setenv(SINGULARITY_TMPDIR="$PREFIX/tmp")
Sys.setenv(SINGULARITY_CACHEDIR="$PREFIX/singularity_cache")
Sys.setenv(SINGULARITY_LOCALCACHEDIR="$PREFIX/tmp/runtime")
Sys.setenv(DEBUGME = "BiocParallel")

####################
## Cutadapt
## Work for paired-end reads
## Input files need to be in one of these formats: 
##   - FASTA with extensions: .fasta, .fa or .fna
##   - FASTQ with extensions: .fastq or .fq
##   - Any of the above, but compressed as .gz, .bz2 or .xz

## req1 <- list(class = "DockerRequirement", 
##              dockerPull = "kfdrc/cutadapt")
req1 <- requireDocker("quay.io/biocontainers/cutadapt:4.2--py310h1425a21_0")

p1 <- InputParam(id = "threadN", type = "int?", prefix = "-j", position = 1, default = 1L)
p2a <- InputParam(id = "adapter1a", type = "string?", prefix = "-a", position = 2)
p2b <- InputParam(id = "adapter2a", type = "string?", prefix = "-A", position = 3)
p2c <- InputParam(id = "adapter1g", type = "string?", prefix = "-g", position = 4)
p2d <- InputParam(id = "adapter2g", type = "string?", prefix = "-G", position = 5)
p2e <- InputParam(id = "adapter1b", type = "string?", prefix = "-b", position = 6)
p2f <- InputParam(id = "adapter2b", type = "string?", prefix = "-B", position = 7)
p3 <- InputParam(id = "out1prefix", type = "string", prefix = "-o", position = 8)
p4 <- InputParam(id = "out2prefix", type = "string?", prefix = "-p", position = 9)
p5 <- InputParam(id = "in1", type = "File", position = 99)
p6 <- InputParam(id = "in2", type = "File?", position = 100)
o1 <- OutputParam(id = "out1", type = "File", glob = "\$(inputs.out1prefix)")
o2 <- OutputParam(id = "out2", type = "File?", glob = "\$(inputs.out2prefix)")

cutadapt <- cwlProcess(baseCommand = "cutadapt", 
                         requirements = list(req1),
                         inputs = InputParamList(p1, p2a, p2b, p2c, p2d, p2e, p2f, p3, p4, p5, p6), 
                         outputs = OutputParamList(o1, o2))
#######################
#cutadapt

## Prepare inputs
fqs <- list.files("$DATADIR", ".fastq.gz", full.names = TRUE, recursive = TRUE)
fq1 <- fqs[grep("R1", fqs)]
fq2 <- fqs[grep("R2", fqs)]
sids <- basename(dirname(fq1))
fq1.trim <- sub('.fastq.gz','_trimmed.fastq.gz', basename(fq1))
fq2.trim <- sub('.fastq.gz','_trimmed.fastq.gz', basename(fq2))
ADAPTER_FWD <- "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  ## from UG1001-08 user guide
ADAPTER_REV <- "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"   ##

#cutadapt\$in1 <- fq1
#cutadapt\$in2 <- fq2
#cutadapt\$out1prefix <- fq1.trim
#cutadapt\$out2prefix <- fq2.trim
#cutadapt\$threadN <- 1 
#cutadapt\$adapter1a <- ADAPTER_FWD
#cutadapt\$adapter2a <- ADAPTER_REV
#cutadapt

#runCWL(cutadapt, outdir = "$OUTDIR",
#       docker="singularity", showLog=TRUE)

fq1L <- tapply(fq1, sids, as.list)
fq2L <- tapply(fq2, sids, as.list)
fq1L.trim <- tapply(fq1.trim, sids, as.list)
fq2L.trim <- tapply(fq2.trim, sids, as.list)

done <- list.files("$OUTDIR", "*_trimmed.fastq.gz", recursive = TRUE)
idx <- !sids %in% basename(dirname(done))
#idx=TRUE
inputList <- list(in1 = fq1L[idx], 
                  in2 = fq2L[idx],
                  out1prefix = fq1L.trim[idx],
                  out2prefix = fq2L.trim[idx])
paramList <- list(threadN = 16,
                  adapter1a = ADAPTER_FWD,
                  adapter2a = ADAPTER_REV)
inputList
paramList

# Submit a job using runCWLBatch() and BatchtoolsParam() functions
runCWLBatch(cutadapt, outdir = "$OUTDIR",
            inputList, paramList,
            BPPARAM = BatchtoolsParam(workers = lengths(inputList)[1],
                                      cluster = "slurm",
                                      #template = "$PREFIX/rcwl/slurm_rpci.tmpl", ## For RPCI cluster
                                      template = "$PREFIX/rcwl/slurm_nih.tmpl",   ## For ub-hpc:general-compute
                                      #template = "$PREFIX/rcwl/slurm_nih_2.tmpl", ## For ub-hpc:scavenger
                                      resources = list(ncpus = 16,
                                                       jobname = "trim",
                                                       walltime = 60*60*3,
                                                       memory = 128000),
                                      log = TRUE, logdir = "$JOBDIR",
                                      progressbar = FALSE,
                                      saveregistry = FALSE),
                                      stderr = "", docker="singularity")
EOF

#/user/mkorobkin/rpci/R/bin/Rscript $SRCDIR/align.R &
Rscript2 $SRCDIR/trim.R & 
