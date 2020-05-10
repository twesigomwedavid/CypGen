#!/usr/bin/env nextflow

in_ch = Channel.fromFilePairs(params.in_bam, type: 'file') { file -> file.name.replaceAll(/.bam|.bai$/,'') }

in_vars = file(params.in_variants, type: 'file')
in_core_var = file(params.core_variants, type: 'file')
hg19 = file(params.hg19_ref, type: 'file')
db_core_snvs = file(params.db_core_file, type: 'file')
db_tiebreak = file(params.db_tb_file, type: 'file')

//in_ch.subscribe {println "$it"}


process call_variants {
//   maxForks 5

   input:
      set val(name), file(bam) from aln_ch

   output:	        
      set val(name), path("${name}_var") into var_ch
      
   script:
    """ 
   	graphtyper genotype $hg19 --sam=${bam} --region=22:42522300-42528400 --output=${name}_var --prior_vcf=$in_vars       
    """

}

process analyse_vcf {
//   maxForks 5
   
   publishDir params.out_dir, mode: 'copy', overwrite: 'true'

   input:
      set val(name), path("${name}_vars") from var_ch

   output:
      set val(name), file("${name}_snv_def.alleles") into dip_candidates_ch

   script:
 
    """
      bcftools norm -m - ${name}_vars/22/042522300-042528400.vcf.gz | bcftools view -e 'GT="1/0"' | bcftools view -e 'GT="0/0"' | bgzip -c > ${name}_vars/22/${name}_all_norm.vcf.gz
      tabix ${name}_vars/22/${name}_all_norm.vcf.gz
      bcftools isec ${name}_vars/22/${name}_all_norm.vcf.gz ${in_core_var} -p ${name}_int
      bgzip ${name}_int/0002.vcf
      tabix ${name}_int/0002.vcf.gz
      bcftools norm -m - ${name}_int/0002.vcf.gz | bcftools view -e 'GT="1/0"' | bcftools view -e 'GT="0/0"' | bgzip -c > ${name}_vars/22/${name}_core.vcf.gz
      tabix ${name}_vars/22/${name}_core.vcf.gz
      bcftools query -f'[%POS~%REF>%ALT~%GT\n]' ${name}_vars/22/${name}_core.vcf.gz  > ${name}_vars/${name}_core_snvs.dip
      bcftools query -f '%POS~%REF>%ALT\n' ${name}_vars/22/${name}_all_norm.vcf.gz > ${name}_vars/${name}_full.dip
      bcftools query -f'[%POS~%REF>%ALT~%GT\n]' ${name}_vars/22/${name}_all_norm.vcf.gz > ${name}_vars/${name}_gt.dip
      python /home/david/asterigraph/workflows/bin/test2.py ${db_core_snvs} ${name}_vars/${name}_core_snvs.dip ${name}_vars/${name}_full.dip ${name}_vars/${name}_gt.dip ${db_tiebreak} > ${name}_snv_def.alleles
  
    """

}

