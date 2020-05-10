#!/bin/bash

while getopts f:o: option
do
case "${option}"
in
f) FILE1=${OPTARG};;
o) FILE2=${OPTARG};;
esac
done


python ../bin/get_snv_def_alleles.py $FILE1 > $FILE2
