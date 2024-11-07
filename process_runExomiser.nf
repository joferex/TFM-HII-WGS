#!/usr/bin/env nextflow

// Define parameters
// Personalizar el análisis eligiendo el archivo .yml deseado
params.yaml_config = "/Volumes/EXTERNAL_USB/TFM/exomiser-cli-14.0.0/examples/HII_analysis_panel_nofilter.yml"
//////
params.exomiser_jar = "/Volumes/EXTERNAL_USB/TFM/exomiser-cli-14.0.0/exomiser-cli-14.0.0.jar"
params.spring_config = "file:/Volumes/EXTERNAL_USB/TFM/exomiser-cli-14.0.0/application.properties"
params.memory_opts = "-Xms2g -Xmx4g"
params.assembly = "hg38"
params.vcf_dir = "/Volumes/EXTERNAL_USB/TFM/vcf_samples"  // Directory containing VCF files cambiarlo a voluntad con --vcf_dir en la terminal
params.output_dir = "/Volumes/EXTERNAL_USB/TFM/results/HPO-nofilter-HII"  // Output directory

// Create a channel that picks up all VCF files in the specified directory
vcf_files = Channel.fromPath("${params.vcf_dir}/*.vcf.gz")

process runExomiser {
    
    input:
    path vcf_file

    script:
    // Extract the base name of the VCF file (without extension)
    def base_name = vcf_file.getBaseName().split('\\.')[0]
    
    // Create the specific output directory for this VCF file
    def specific_output_dir = "${params.output_dir}/sample-${base_name}"
    """
    mkdir -p ${specific_output_dir}
    
    java ${params.memory_opts} -jar ${params.exomiser_jar} \\
        --analysis ${params.yaml_config} \\
        --assembly ${params.assembly} \\
        --vcf ${vcf_file} \\
        --spring.config.location=${params.spring_config}
    
    # Move generated files to the specific output directory
    mv ${params.output_dir}/${base_name}* ${specific_output_dir}/

    echo "Analysis completed for: ${vcf_file}"
    echo "Files moved to: ${specific_output_dir}"
    """
}

workflow {
    runExomiser(vcf_files)
}