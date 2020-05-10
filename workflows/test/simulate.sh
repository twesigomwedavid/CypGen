#!/bin/bash -l
#SBATCH -J test-run
#SBATCH -e errors_%j.txt
#SBATCH -c 4
#SBATCH --mem=10096  


cd /path/to/work_dir

for i in $(cat diplotypes.txt); do
     grep ${i} genotypes3_mod.dbs | awk '{print $2}' > ${i}.dip
     python simulate_vcfs.py ${i}.dip ${i}.vcf ${i}
     bgzip ${i}.vcf
     tabix ${i}.vcf.gz
     bcftools isec ${i}.vcf.gz /home/david/asterigraph/resources/allele_def_var.vcf.gz -p ${i}
     bcftools query -f'[%POS~%REF>%ALT~%GT\n]' ${i}/0002.vcf > ${i}_core_snvs.dip
     sed 's/~0\/1//g' ${i}.dip | sed 's/~1\/1//g' > ${i}_full.dip
     grep ${i} genotypes3_mod.dbs | awk '{print $2}' | sed 's/\;/\n/g' > ${i}_gt.dip
    python /path/to/workflows/bin/test2.py /path/to/database/diplo_db_debugged2.dbs ${i}_core_snvs.dip ${i}_full.dip ${i}_gt.dip /path/to/database/genotypes4.dbs > results/${i}.alleles

    echo ${i} >> results/sim_samples.txt

    if [[ $(grep "*" results/${i}.alleles) ]]; then
    grep "*" results/${i}.alleles >> results/allele_calls.list
    else
    echo check >> results/allele_calls.list
    fi
done

paste results/sim_samples.txt results/allele_calls.list > results/testrun.alleles
