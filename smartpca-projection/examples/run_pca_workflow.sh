#!/usr/bin/env bash
#
# Updated version of the KKA PCA workflow, preserving the original skeleton.
#
# The only structural change vs. the original loop is the split into two
# stages and the addition of a few .par parameters + model/projection paths:
#
#   Stage A : calculate the reference PCA ONCE and save a projection model.
#   Stage B : load that model and project NEW samples only (reference PCA is
#             not recomputed).
#
# You still build a .par per population set and call:
#     smartpca -p ${of2}.par
#
# Set SMARTPCA to the feature build (e.g. .../smartpca-projection/build/smartpca).
set -euo pipefail

SMARTPCA="${SMARTPCA:-smartpca}"
date="260731"
pt1="$(pwd)/"

# reference genotype prefix (the panel the PCA is built on)
fn_ref="/home/projects1/AADRv54_EastAsian/PCA_250805/reference.${date}"
# new-sample genotype prefix (the individuals to be projected later)
fn_new="/home/projects1/mongol_period/analysis/PCA/KARCEM3/${date}/KKA.PCA.1240KHO.${date}"

of1="PCA_${date}"

# population-set lists, exactly as before
fn2="/home/projects1/AADRv54_EastAsian/PCA_250805"
nset=$(ls ${fn2}_*.pops | wc -l)
for K in $(seq 1 "$nset"); do cp "${fn2}_${K}.pops" "${of1}_${K}.pops"; done

nset=$(ls ${of1}_*.pops | wc -l)

########################################################################
## Stage A : reference PCA + save model  (run once per population set) ##
########################################################################
for K in $(seq 1 "$nset"); do
    of2="${of1}_${K}"
    parf="${of2}.refA.par"

    {
      echo "genotypename: ${fn_ref}.geno"
      echo "snpname: ${fn_ref}.snp"
      echo "indivname: ${fn_ref}.ind"

      echo "evecoutname: ${pt1}${of2}.evec"
      echo "evaloutname: ${pt1}${of2}.eval"
      echo "poplistname: ${pt1}${of2}.pops"

      echo "altnormstype: NO"
      echo "numoutevec: 20"
      echo "numoutlieriter: 0"
      echo "numoutlierevec: 0"
      echo "outliersigmathresh: 6.0"
      echo "numthreads: 8"
      echo "qtmode: 0"
      echo "lsqproject: YES"

      # NEW: persist the projection model for this population set
      echo "modeloutname: ${pt1}${of2}.smartpca.model"
    } > "$parf"

    # same call shape as before
    sbatch -c 8 --mem 15000 --wrap="${SMARTPCA} -p ${parf} > ${of2}.refA.log"
done

########################################################################
## Stage B : project NEW samples only, using the saved model          ##
##           (no reference PCA recomputation)                          ##
########################################################################
for K in $(seq 1 "$nset"); do
    of2="${of1}_${K}"
    parf="${of2}.projB.par"

    {
      echo "genotypename: ${fn_new}.geno"
      echo "snpname: ${fn_new}.snp"
      echo "indivname: ${fn_new}.ind"

      # NEW: load the model saved in Stage A and project only
      echo "modelname: ${pt1}${of2}.smartpca.model"
      echo "projectionoutname: ${pt1}${of2}.projected.evec"
      echo "projectiononly: YES"
      # echo "modelstrict: NO"   # uncomment to tolerate model SNPs absent from new data
    } > "$parf"

    # same call shape as before
    sbatch -c 8 --mem 15000 --wrap="${SMARTPCA} -p ${parf} > ${of2}.projB.log"
done
