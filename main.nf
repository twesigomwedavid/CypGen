#!/usr/bin/env nextflow


params.build='hg38'

if (params.build=='b37') {
    db = params.build37db
    res_dir = params.build37res
    caller_dir = params.build37call
    chrom = "22"
    region_a1 = "22:42522000-42542000"
    region_a2 = "042522000-042542000"
    region_b1 = "22:42522300-42528400"
    region_b2 = "042522300-042528400"
    debug38 = ""

} else {
    db = params.build38db
    res_dir = params.build38res
    caller_dir = params.build38call
    chrom = "chr22"
    region_a1 = "chr22:42126000-42137500"
    region_a2 = "042126000-042137500"
    region_b1 = "chr22:42126300-42132400"
    region_b2 = "042126300-042132400" 
    debug38 = "--minimum_extract_score_over_homref=0"
}

params.format='binary'

if (params.format=='compressed') {
    ext = "cram"
    ind = "crai"
    cram_options = "--force_use_input_ref_for_cram_reading" 

} else { 
    ext = "bam"
    ind = "bai"
    cram_options = ""

}


align_file = Channel.fromFilePairs(params.in_bam, type: 'file') {  file -> file.name.replaceAll(/.${ext}|.${ind}$/,'') }

//ref_ch = Channel.fromFilePairs(params.ref_genome, type: 'file') { file -> file.simpleName }

//ref_genome = file(params.ref_file, type: 'file', checkIfExists: true)

//align_file.subscribe {println "$it"}

// switch (params.ref_file) {
//     case [null]:
//         ref_file = "Unspecified!"
//         break
//     default:
//         ref_file = file(params.ref_file, type: 'file', checkIfExists: true)
// }


// switch (params.in_bam) {
//     case [null]:
//         align_file = "Unspecified!"
//         break
//     default:
//         align_file = Channel.fromFilePairs(params.in_bam, type: 'file', checkIfExists: true) { file -> file.simpleName }
// }



align_file.into { data1; data2; data3; data4; data5 }

//ref_genome.into {ref_ch1; ref_ch2; ref_ch3; ref_ch4 }

ref_dir_val = new File("${params.ref_file}").getParent()
ref_genome = new File("${params.ref_file}").getName()


//ref_dir = Channel . value ("${ref_dir_val}")


process call_snvs1 {
//   maxForks 5

   input:
      set val(name), file(bam) from data1
      path ref_dir from Channel.value("${ref_dir_val}")
      path res_dir

   output:	        
      set val(name), path("${name}_var_1") into var_ch1
      
   script:
    """
   	graphtyper genotype ${ref_dir}/${ref_genome} --sam=${name}.${ext} --region=${region_a1} --output=${name}_var_1 --prior_vcf=${res_dir}/common_plus_core_var.vcf.gz -a ${debug38} ${cram_options}
    """

}


process call_snvs2 {
//   maxForks 5

   input:
      set val(name), file(bam) from data2
      path ref_dir from Channel.value("${ref_dir_val}")

   output: 
      set val(name), path("${name}_var_2") into var_ch2

   script:
    """  
	graphtyper genotype ${ref_dir}/${ref_genome} --sam=${name}.${ext} --region=${region_a1} --output=${name}_var_2 -a ${debug38} ${cram_options}
    """

}


process call_sv_del {
//   maxForks 5

   input:
      set val(name), file(bam) from data3
      path ref_dir from Channel.value("${ref_dir_val}")
      path res_dir

   output:
      set val(name), path("${name}_sv_del") into sv_ch1

   script:
     """
	graphtyper genotype_sv ${ref_dir}/${ref_genome} --sam=${name}.${ext} --region=${region_a1} --output=${name}_sv_del ${res_dir}/sv_test.vcf.gz

     """

}

process call_sv_dup {
//   maxForks 5

   input:
      set val(name), file(bam) from data4
      path ref_dir from Channel.value("${ref_dir_val}")
      path res_dir

   output:
      set val(name), path("${name}_sv_dup") into sv_ch2

   script:
     """
	graphtyper genotype_sv ${ref_dir}/${ref_genome} --sam=${name}.${ext} --region=${region_a1} --output=${name}_sv_dup ${res_dir}/sv_test3.vcf.gz

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
	samtools bedcov ${res_dir}/test3.bed ${name}.${ext} > ${name}_2d6_ctrl.depth      

     """

}


var_ch1.join(var_ch2).set { var_ch_joined }


process format_snvs {

   input:
      set val(name), path("${name}_var_1"), path("${name}_var_2") from var_ch_joined

   output:
      set val(name), path("${name}_var") into (var_norm1, var_norm2)

   script:
    """
        bcftools isec -p ${name}_var -Oz ${name}_var_1/${chrom}/${region_a2}.vcf.gz ${name}_var_2/${chrom}/${region_a2}.vcf.gz
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
      bcftools query -f'%ID\t%ALT\t[\t%GT\t%DP]\t%INFO/ABHet\t%INFO/ABHom\n' ${name}_gene_del/${chrom}/${region_a2}.vcf.gz > ${name}_gene_del/${name}_gene_del_summary.txt
     """

}


sv_ch2.join(core_vars1).set {dup_int}

process analyse_2 {
//   maxForks 5

   errorStrategy 'ignore'
   tag "${name}"

   input:
      set val(name), path("${name}_gene_dup"), path("${name}_int") from dup_int
     // set val(name), path("${name}_int") from core_vars1

   output:
      set val(name), path("${name}_gene_dup/${name}_gene_dup_summary.txt") into dup_ch

   script:
     """
      bcftools query -f'%POS~%REF>%ALT\t[\t%GT\t%DP]\t%INFO/ABHet\t%INFO/ABHom\n' -i'GT="alt"' ${name}_gene_dup/${chrom}/${region_a2}.vcf.gz > ${name}_gene_dup/${name}_gene_dup_summary.txt
      bcftools query -f'%POS~%REF>%ALT\t[\t%GT\t%DP]\t%INFO/ABHet\t%INFO/ABHom\n' -i'GT="alt"' ${name}_int/${name}_core.vcf.gz >> ${name}_gene_dup/${name}_gene_dup_summary.txt

     """

}


var_norm2.join(core_vars2).set {dip_req}

process analyse_3 {
//   maxForks 5

   errorStrategy 'ignore'
   tag "${name}"

   input:
      set val(name), path("${name}_vars"), path("${name}_int") from dip_req
     // set val(name), path("${name}_int") from core_vars2

   output:
      set val(name), path("${name}_vars/${name}_core_snvs.dip"), path("${name}_vars/${name}_full.dip"), path("${name}_vars/${name}_gt.dip") into prep_ch

   script:
     """
      bcftools query -f'[%POS~%REF>%ALT~%GT\n]' ${name}_int/${name}_core.vcf.gz  > ${name}_vars/${name}_core_snvs.dip
      bcftools query -f '%POS~%REF>%ALT\n' ${name}_vars/${name}_all_norm.vcf.gz > ${name}_vars/${name}_full.dip
      bcftools query -f'[%POS~%REF>%ALT~%GT\n]' ${name}_vars/${name}_all_norm.vcf.gz > ${name}_vars/${name}_gt.dip

     """

}


prep_ch.join(del_ch).set {fin_files1}
fin_files1.join(dup_ch).set {fin_files2}
fin_files2.join(sv_ch3).set {fin_files}


process call_stars {
//   maxForks 5

   publishDir params.out_dir, mode: 'copy', overwrite: 'true'

   errorStrategy 'ignore'
   tag "${name}"

   input:
      set val(name), path("${name}_vars/${name}_core_snvs.dip"), path("${name}_vars/${name}_full.dip"), path("${name}_vars/${name}_gt.dip"), path("${name}_gene_del/${name}_gene_del_summary.txt"), path("${name}_gene_dup/${name}_gene_dup_summary.txt"), file("${name}_2d6_dp") from fin_files
     // set val(name), path("${name}_gene_del/${name}_gene_del_summary.txt") from del_ch
     // set val(name), path("${name}_gene_dup/${name}_gene_dup_summary.txt") from dup_ch
     // set val(name), file("${name}_2d6_dp") from sv_ch3
      path db
      path caller_dir

   output:
      set val(name), file("${name}_2d6.alleles") into star_ch

   script:
   
    """
     python3 ${caller_dir}/cypgen_caller.py ${db}/diplo_db_debugged2.dbs ${name}_vars/${name}_core_snvs.dip ${name}_vars/${name}_full.dip ${name}_vars/${name}_gt.dip ${db}/genotypes4.dbs ${name}_gene_del/${name}_gene_del_summary.txt ${name}_gene_dup/${name}_gene_dup_summary.txt ${name}_2d6_dp ${db}/haps_var_new.dbs > ${name}_2d6.alleles  

    """

}

