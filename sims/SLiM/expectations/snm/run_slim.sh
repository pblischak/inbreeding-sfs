#!/bin/bash

printf "****  Running simulations with SLiM...  ****\n\n"
for i in {1..50}
do
	slim -d "Fis=$1" -d "Rep=${i}" snm.slim
done

printf "****  Converting from VCF to frequency spectra...  ****\n\n"
for f in SLiM_F${1}_snm*.vcf
do
	./vcf2sfs.py $f 25
done

printf "****  Compressing VCF files...  ****\n\n"

tar czf SLiM_F${1}_snm_VCFs.tar.gz SLiM_F${1}_snm_*.vcf
rm -f SLiM_F${1}_snm_*.vcf
