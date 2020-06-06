#!/usr/bin/env nextflow

in_ch = Channel.fromFilePairs(params.in_bam, type: 'file') { file -> file.name.replaceAll(/.bam|.bai$/,'') }

in_vars = file(params.in_variants, type: 'file')

in_core_var = file(params.core_variants, type: 'file')

hg19 = file(params.hg19_ref, type: 'file')

db_core_snvs = file(params.db_core_file, type: 'file')

db_tiebreak = file(params.db_tb_file, type: 'file')

coded_star_5 = file(params.sv_del_file, type: 'file')

coded_full_gendup = file(params.sv_dup_file, type: 'file')

bed_file = file(params.bed, type: 'file')

db_haps = file(params.db_haplo, type: 'file')

cypgen_mod = file(params.call_script, type: 'file')

//in_ch.subscribe {println "$it"}


process call_variants {
//   maxForks 5

   input:
      set val(name), file(bam) from in_ch

   output:	        
      set val(name), path("${name}_var"), path("${name}_sv_del"), path("${name}_sv_dup"), file("${name}_2d6_ctrl.depth") into var_ch
      
   script:
    """
   	graphtyper genotype $hg19 --sam=${bam.get(0)} --region=22:42522000-42542000 --output=${name}_var_1 --prior_vcf=${in_vars}
	graphtyper genotype $hg19 --sam=${bam.get(0)} --region=22:42522000-42542000 --output=${name}_var_2
	bcftools isec -p ${name}_var -Oz ${name}_var_1/22/042522000-042542000.vcf.gz ${name}_var_2/22/042522000-042542000.vcf.gz 	
	bcftools concat -a -D -r 22:42522300-42528400 ${name}_var/0000.vcf.gz ${name}_var/0001.vcf.gz ${name}_var/0002.vcf.gz -Oz -o ${name}_var/${name}_042522300-042528400.vcf.gz
	tabix ${name}_var/${name}_042522300-042528400.vcf.gz  
	graphtyper genotype_sv $hg19 --sam=${bam.get(0)} --region=22:42522000-42540000 --output=${name}_sv_del ${coded_star_5}
	graphtyper genotype_sv $hg19 --sam=${bam.get(0)} --region=22:42522000-42540000 --output=${name}_sv_dup ${coded_full_gendup}
	samtools bedcov ${bed_file} ${bam.get(0)} > ${name}_2d6_ctrl.depth      

    """

}

process analyse {
//   maxForks 5
   
   publishDir params.out_dir, mode: 'copy', overwrite: 'true'

   errorStrategy 'ignore'
   tag "${name}"   

   input:
      set val(name), path("${name}_vars"), path("${name}_gene_del"), path("${name}_gene_dup"), file("${name}_2d6_dp") from var_ch

   output:
      set val(name), file("${name}_snv_def.alleles") into dip_candidates_ch

   script:
 
    """
      bcftools norm -m - ${name}_vars/${name}_042522300-042528400.vcf.gz | bcftools view -e 'GT="1/0"' | bcftools view -e 'GT="0/0"' | bgzip -c > ${name}_vars/${name}_all_norm.vcf.gz
      tabix ${name}_vars/${name}_all_norm.vcf.gz
      bcftools isec ${name}_vars/${name}_all_norm.vcf.gz ${in_core_var} -p ${name}_int -Oz
 
      bcftools norm -m - ${name}_int/0002.vcf.gz | bcftools view -e 'GT="1/0"' | bcftools view -e 'GT="0/0"' | bgzip -c > ${name}_vars/${name}_core.vcf.gz
      tabix ${name}_vars/${name}_core.vcf.gz
      
      bcftools query -f'%ID\t%ALT\t[\t%GT\t%DP]\t%INFO/ABHet\t%INFO/ABHom\n' ${name}_gene_del/22/042522000-042540000.vcf.gz > ${name}_gene_del/${name}_gene_del_summary.txt
      bcftools norm -m - ${name}_gene_dup/22/042522000-042540000.vcf.gz | bcftools view -e 'GT="1/0"' | bcftools view -e 'GT="0/0"' | bgzip -c > ${name}_gene_dup/22/${name}_dup_norm.vcf.gz
      tabix ${name}_gene_dup/22/${name}_dup_norm.vcf.gz
      bcftools query -f'%POS~%REF>%ALT\t[\t%GT\t%DP]\t%INFO/ABHet\t%INFO/ABHom\n' -i'GT="alt"' ${name}_gene_dup/22/${name}_dup_norm.vcf.gz > ${name}_gene_dup/${name}_gene_dup_summary.txt
      bcftools query -f'%POS~%REF>%ALT\t[\t%GT\t%DP]\t%INFO/ABHet\t%INFO/ABHom\n' -i'GT="alt"' ${name}_vars/${name}_core.vcf.gz >> ${name}_gene_dup/${name}_gene_dup_summary.txt
      bcftools query -f'[%POS~%REF>%ALT~%GT\n]' ${name}_vars/${name}_core.vcf.gz  > ${name}_vars/${name}_core_snvs.dip
      bcftools query -f '%POS~%REF>%ALT\n' ${name}_vars/${name}_all_norm.vcf.gz > ${name}_vars/${name}_full.dip
      bcftools query -f'[%POS~%REF>%ALT~%GT\n]' ${name}_vars/${name}_all_norm.vcf.gz > ${name}_vars/${name}_gt.dip

      python ${cypgen_mod} ${db_core_snvs} ${name}_vars/${name}_core_snvs.dip ${name}_vars/${name}_full.dip ${name}_vars/${name}_gt.dip ${db_tiebreak} ${name}_gene_del/${name}_gene_del_summary.txt ${name}_gene_dup/${name}_gene_dup_summary.txt ${name}_2d6_dp ${db_haps} > ${name}_snv_def.alleles
  
    """

}

