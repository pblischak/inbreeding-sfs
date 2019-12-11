#!/bin/bash

for i in {1..50}
do
	slim -d "Fis=${1}" -d "Rep=${i}" domestication.slim
done

for r in {1..50}
do
  for i in {1..100}
  do
    printf "$r\t$i\n"
    ./vcf2jsfs.py SLiM_F${1}_domestication_pop1_T0.2_M0.1_rep${r}_${i}.vcf SLiM_F${1}_domestication_pop2_T0.2_M0.1_rep${r}_${i}.vcf SLiM_F${1}_domestication_rep${r}_${i}.fs
  done
done
