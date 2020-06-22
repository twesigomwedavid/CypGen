#!/usr/bin/env nextflow


params.build='hg38'

if (params.build=='hg19') {
    db = params.build19db
    res_dir = params.build19res
    caller_dir = params.build19call
    region_a1 = "22:42522000-42542000"
    region_a2 = "042522000-042542000"
    region_b1 = "22:42522300-42528400"
    region_b2 = "042522300-042528400"

} else {
    db = params.build38db
    res_dir = params.build38res
    caller_dir = params.build38call
    region_a1 = "22:42126000-42146000"
    region_a2 = "042126000-042146000"
    region_b1 = "22:42126300-42132400"
    region_b2 = "042126300-042132400" 

}


align_file = Channel.fromFilePairs(params.in_bam, type: 'file') { file -> file.simpleName }

ref_ch = Channel.fromFilePairs(params.ref_genome, type: 'file') { file -> file.simpleName }

//ref_ch.subscribe {println "$it"}

// switch (params.ref_genome) {
//     case [null]:
//         ref_genome = "Unspecified!"
//         break
//     default:
//         ref_genome = Channel.fromFilePairs(params.ref_genome, type: 'file', checkIfExists: true)
// }


// switch (params.in_bam) {
//     case [null]:
//         in_bam = "Unspecified!"
//         break
//     default:
//         in_bam = Channel.fromFilePairs(params.in_bam, type: 'file', checkIfExists: true)
// }



align_file.into { data1; data2; data3; data4; data5 }

ref_ch.into {ref_ch1; ref_ch2; ref_ch3; ref_ch4 }


process call_snvs1 {
//   maxForks 5

   input:
      set val(name), file(bam) from data1
      set val(ref_name), file(ref_genome) from ref_ch1
      path res_dir

   output:	        
      set val(name), path("${name}_var_1") into var_ch1
      
   script:
   
    """
   	graphtyper genotype ${ref_genome.get(0)} --sam=${bam.get(0)} --region=${region_a1} --output=${name}_var_1 --prior_vcf=${res_dir}/common_plus_core_var.vcf.gz
    """

}


process call_snvs2 {
//   maxForks 5

   input:
      set val(name), file(bam) from data2
      set val(ref_name), file(ref_genome) from ref_ch2

   output: 
      set val(name), path("${name}_var_2") into var_ch2

   script:
    """  
	graphtyper genotype ${ref_genome.get(0)} --sam=${bam.get(0)} --region=${region_a1} --output=${name}_var_2
    """

}


process call_sv_del {
//   maxForks 5

   input:
      set val(name), file(bam) from data3
      set val(ref_name), file(ref_genome) from ref_ch3
      path res_dir

   output:
      set val(name), path("${name}_sv_del") into sv_ch1

   script:
     """
	graphtyper genotype_sv ${ref_genome.get(0)} --sam=${bam.get(0)} --region=${region_a1} --output=${name}_sv_del ${res_dir}/sv_test.vcf.gz

     """

}

process call_sv_dup {
//   maxForks 5

   input:
      set val(name), file(bam) from data4
      set val(ref_name), file(ref_genome) from ref_ch4
      path res_dir

   output:
      set val(name), path("${name}_sv_dup") into sv_ch2

   script:
     """
	graphtyper genotype_sv ${ref_genome.get(0)} --sam=${bam.get(0)} --region=${region_a1} --output=${name}_sv_dup ${res_dir}/sv_test3.vcf.gz

     """
}


process get_depth {

   input:
      set val(name), file(bam) from data5
      path res_dir

   output:
      set val(name), file("${name}_2d6_ctrl.depth") into sv_ch3

   script:
     """   
	samtools bedcov ${res_dir}/test3.bed ${bam.get(0)} > ${name}_2d6_ctrl.depth      

     """

}


process format_snvs {

   input:
      set val(name), path("${name}_var_1") from var_ch1
      set val(name), path("${name}_var_2") from var_ch2

   output:
      set val(name), path("${name}_var") into (var_norm1, var_norm2)

   script:
    """
        bcftools isec -p ${name}_var -Oz ${name}_var_1/22/${region_a2}.vcf.gz ${name}_var_2/22/${region_a2}.vcf.gz
        bcftools concat -a -D -r ${region_b1} ${name}_var/0000.vcf.gz ${name}_var/0001.vcf.gz ${name}_var/0002.vcf.gz -Oz -o ${name}_var/${name}_${region_b2}.vcf.gz
        tabix ${name}_var/${name}_${region_b2}.vcf.gz
        bcftools norm -m - ${name}_var/${name}_${region_b2}.vcf.gz | bcftools view -e 'GT="1/0"' | bcftools view -e 'GT="0/0"' | bgzip -c > ${name}_var/${name}_all_norm.vcf.gz
        tabix ${name}_var/${name}_all_norm.vcf.gz
    """

}


process get_core_var {
//   maxForks 5
   
   errorStrategy 'ignore'
   tag "${name}"   

   input:
      set val(name), path("${name}_vars") from var_norm1
      path res_dir

   output:
      set val(name), path("${name}_int") into (core_vars1, core_vars2)

   script:
 
    """
      bcftools isec ${name}_vars/${name}_all_norm.vcf.gz ${res_dir}/allele_def_var.vcf.gz -p ${name}_int -Oz
      bcftools norm -m - ${name}_int/0002.vcf.gz | bcftools view -e 'GT="1/0"' | bcftools view -e 'GT="0/0"' | bgzip -c > ${name}_int/${name}_core.vcf.gz
      tabix ${name}_int/${name}_core.vcf.gz

     """

}



process analyse_1 {
//   maxForks 5

   errorStrategy 'ignore'
   tag "${name}"

   input:
      set val(name), path("${name}_gene_del") from sv_ch1

   output:
      set val(name), path("${name}_gene_del/${name}_gene_del_summary.txt") into del_ch

   script:
     """
      bcftools query -f'%ID\t%ALT\t[\t%GT\t%DP]\t%INFO/ABHet\t%INFO/ABHom\n' ${name}_gene_del/22/${region_a2}.vcf.gz > ${name}_gene_del/${name}_gene_del_summary.txt
     """

}


process analyse_2 {
//   maxForks 5

   errorStrategy 'ignore'
   tag "${name}"

   input:
      set val(name), path("${name}_gene_dup") from sv_ch2
      set val(name), path("${name}_int") from core_vars1

   output:
      set val(name), path("${name}_gene_dup/${name}_gene_dup_summary.txt") into dup_ch

   script:
     """
      bcftools query -f'%POS~%REF>%ALT\t[\t%GT\t%DP]\t%INFO/ABHet\t%INFO/ABHom\n' -i'GT="alt"' ${name}_gene_dup/22/${region_a2}.vcf.gz > ${name}_gene_dup/${name}_gene_dup_summary.txt
      bcftools query -f'%POS~%REF>%ALT\t[\t%GT\t%DP]\t%INFO/ABHet\t%INFO/ABHom\n' -i'GT="alt"' ${name}_int/${name}_core.vcf.gz >> ${name}_gene_dup/${name}_gene_dup_summary.txt

     """

}


process analyse_3 {
//   maxForks 5

   errorStrategy 'ignore'
   tag "${name}"

   input:
      set val(name), path("${name}_vars") from var_norm2
      set val(name), path("${name}_int") from core_vars2

   output:
      set val(name), path("${name}_vars/${name}_core_snvs.dip"), path("${name}_vars/${name}_full.dip"), path("${name}_vars/${name}_gt.dip") into prep_ch

   script:
     """
      bcftools query -f'[%POS~%REF>%ALT~%GT\n]' ${name}_int/${name}_core.vcf.gz  > ${name}_vars/${name}_core_snvs.dip
      bcftools query -f '%POS~%REF>%ALT\n' ${name}_vars/${name}_all_norm.vcf.gz > ${name}_vars/${name}_full.dip
      bcftools query -f'[%POS~%REF>%ALT~%GT\n]' ${name}_vars/${name}_all_norm.vcf.gz > ${name}_vars/${name}_gt.dip

     """

}


process call_stars {
//   maxForks 5

   publishDir params.out_dir, mode: 'copy', overwrite: 'true'

   errorStrategy 'ignore'
   tag "${name}"

   input:
      set val(name), path("${name}_vars/${name}_core_snvs.dip"), path("${name}_vars/${name}_full.dip"), path("${name}_vars/${name}_gt.dip") from prep_ch
      set val(name), path("${name}_gene_del/${name}_gene_del_summary.txt") from del_ch
      set val(name), path("${name}_gene_dup/${name}_gene_dup_summary.txt") from dup_ch
      set val(name), file("${name}_2d6_dp") from sv_ch3
      path db
      path caller_dir

   output:
      set val(name), file("${name}_2d6.alleles") into star_ch

   script:
   
    """
     python3 ${caller_dir}/cypgen_caller.py ${db}/diplo_db_debugged2.dbs ${name}_vars/${name}_core_snvs.dip ${name}_vars/${name}_full.dip ${name}_vars/${name}_gt.dip ${db}/genotypes4.dbs ${name}_gene_del/${name}_gene_del_summary.txt ${name}_gene_dup/${name}_gene_dup_summary.txt ${name}_2d6_dp ${db}/haps_var_new.dbs > ${name}_2d6.alleles  

    """

}

