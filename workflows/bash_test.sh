#!/bin/bash

while getopts r:f:o:v: option
do
case "${option}"
in
r) REF=${OPTARG};;
f) FILE1=${OPTARG};;
o) OUTDIR=${OPTARG};;
v) VAR=${OPTARG};;
esac
done


graphtyper genotype $REF --sam=$FILE1 --region=22:42522300-42528400 --output=$OUTDIR --vcf=$VAR
