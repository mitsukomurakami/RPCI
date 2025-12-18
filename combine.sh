#!/usr/bin/sh
# combine.sh

if [[ $# -eq 0 ]] ; then
    echo "Usage: ./combine.sh <patientID>"
    exit 1
fi
pid=$1

PREFIX=/projects/rpci/songliu/mkorobkin
PID=RQ023532-Paragh

JOBDIR=$PREFIX/projects/$PID/$pid
OUTDIR=$JOBDIR/output
SNVDIR=$OUTDIR/variants

## normal sample
Npath=`ls -d $OUTDIR/BAM/${pid}-Sd*` 
nid=`basename $Npath`

## tumor samples
#TLIST=`ls -d $OUTDIR/BAM/${pid}-Sp*`
#for Tpath in $TLIST; do
#  tid=`basename $Tpath`
#done
TLIST="SCMh0036-Spn007"  ## <-------------- Modify here!!!
for tid in $TLIST; do

SRCDIR=$JOBDIR/src/$tid
mkdir -p $SRCDIR

echo $tid -- $nid

# Write an R script
cat > $SRCDIR/combine.R << EOF
#devtools::load_all("/projects/rpci/songliu/qhu/workspace/projects/VariantCombiner/")
devtools::load_all("/projects/rpci/songliu/mkorobkin/rcwl/VariantCombiner/")  ## for debugging

varcombiner <- function(ss, si, m2, mu, vd, id_t, id_n){
    
    library(VariantCombiner)

    ## strelka2
    message("strelka2")
    v1a <- readVcf(ss)
    v1a <- v1a[fixed(v1a)\$FILTER == "PASS"]
    v1b <- readVcf(si)
    v1b <- v1b[fixed(v1b)\$FILTER == "PASS"]
    s1a <- strelka_snv(v1a)
    s1b <- strelka_indel(v1b)
    v_s <- SomaticCombiner(s1a, s1b, sources = c("strelka2", "strelka2"),
                        GENO = c(GT = 1, DP = 1, AD = 1),
                        id_t = id_t, id_n = id_n)

    ## mutect2
    message("mutect2")
    m2v <- readVcf(m2)
    m2v <- m2v[fixed(m2v)\$FILTER %in% c("PASS", "multiallelic")]
    v_m <- SomaticCombiner(m2v, v_s, source = c("mutect2", "strelka2"),
                        GENO = c(GT = 1, DP = 1, AD = 1),
                        id_t = id_t, id_n = id_n)
    
    ## muse
    message("muse")
    mu1 <- readVcf(mu)
    mu1 <- mu1[fixed(mu1)\$FILTER == "PASS"]
    v_m <- SomaticCombiner(v_m, mu1, source = c("", "muse"),
                        GENO = c(GT = 1, DP = 1, AD = 1),
                        id_t = id_t, id_n = id_n)

    ## vardict
    message("verdict")
    vd1 <- readVcf(vd)
    vd1 <- vd1[info(vd1)\$STATUS == "StrongSomatic" & fixed(vd1)\$FILTER == "PASS"]
    vd1 <- vd1[!info(vd1)\$TYPE %in% c("DEL", "DUP", "INV")]
    v_m <- SomaticCombiner(v_m, vd1, source = c("", "vardict"),
                        GENO = c(GT = 1, DP = 1, AD = 1),
                        id_t = id_t, id_n = id_n)

    writeVcf(v_m, paste0("$SNVDIR/$tid/", id_t, "_", id_n, "_strelka2_mutect2_muse_vardict.vcf"))
}

varcombiner( "$SNVDIR/$tid/somatic.snvs.vcf.gz", 
             "$SNVDIR/$tid/somatic.indels.vcf.gz",
             "$SNVDIR/$tid/${nid}.${tid}.filtered.PASS.vcf",
             "$SNVDIR/$tid/${tid}_${nid}_MuSE.vcf",
             "$SNVDIR/$tid/${tid}_${nid}_VarDict.vcf",
             "$tid",
             "$nid"
            )

EOF

Rscript2 $SRCDIR/combine.R

done
