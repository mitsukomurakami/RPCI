#!/usr/bin/sh
# vcf2sj4.sh

if [[ $# -eq 0 ]] ; then
    echo "Usage: ./vcf2sj4.sh <patientID>"
    exit 1
fi
pid=$1

PREFIX=/projects/rpci/songliu/mkorobkin
PID=RQ023532-Paragh

JOBDIR=$PREFIX/projects/$PID/$pid
OUTDIR=$JOBDIR/output
SNVDIR=$OUTDIR/variants

SRCDIR=$JOBDIR/src
if [ ! -d "$SRCDIR" ]; then
  mkdir $SRCDIR
fi

# Write an R script
cat > $SRCDIR/vcf2sj4.R << EOF
devtools::load_all("/projects/rpci/songliu/qhu/workspace/projects/VariantCombiner/")  ## for vcf2sj() function

# Prepare inputs
vcfs_done <- list.files("$SNVDIR", "*strelka2_mutect2_muse_vardict.vep.vcf", recursive = TRUE, full.names = TRUE)

vcfs_done

# Write sj files
for(i in 1:length(vcfs_done)){
    sj1 <- vcf2sj(vcfs_done[i])

    vcfs.ids <- do.call(rbind, lapply(strsplit(basename(vcfs_done[i]), split = "[_.]"), function(x)x[1:2]))
    id_t <- vcfs.ids[1] 
    id_n <- vcfs.ids[2]
    vd <- sj1[grep("vardict", sj1\$SJQuality),]
    ms <- sj1[grep("muse", sj1\$SJQuality),]
    mt <- sj1[grep("mutect2", sj1\$SJQuality),]
    st <- sj1[grep("strelka2", sj1\$SJQuality),]

    write.table(sj1, paste0("${SNVDIR}/", id_t, "/", id_t, "_", id_n, "_strelka2_mutect2_muse_vardict.sj"), row.names=FALSE, quote = FALSE, sep = "\t")
    write.table(vd, paste0("${SNVDIR}/", id_t, "/", id_t, "_", id_n, "_vardict.sj"), row.names=FALSE, quote = FALSE, sep = "\t")
    write.table(ms, paste0("${SNVDIR}/", id_t, "/", id_t, "_", id_n, "_muse.sj"), row.names=FALSE, quote = FALSE, sep = "\t")
    write.table(mt, paste0("${SNVDIR}/", id_t, "/", id_t, "_", id_n, "_mutect2.sj"), row.names=FALSE, quote = FALSE, sep = "\t")
    write.table(st, paste0("${SNVDIR}/", id_t, "/", id_t, "_", id_n, "_strelka2.sj"), row.names=FALSE, quote = FALSE, sep = "\t")
}

# Check results
list.files("$SNVDIR", "*.sj", recursive = TRUE, full.names = TRUE)
EOF

#/user/mkorobkin/rpci/R/bin/Rscript $SRCDIR/align.R &
Rscript2 $SRCDIR/vcf2sj4.R &

