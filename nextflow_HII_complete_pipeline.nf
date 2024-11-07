

// Parámetros del flujo de trabajo
params.ref_dir = "/Volumes/EXTERNAL_USB/TFM/ref"
params.reads_dir = "/Volumes/EXTERNAL_USB/TFM/samples_fastq"
params.aligned_reads_dir = "/Volumes/EXTERNAL_USB/TFM/alineamientos"
params.results_dir = "/Volumes/EXTERNAL_USB/TFM/results"
params.data_dir = "/Volumes/EXTERNAL_USB/TFM/data"

// Parámetros de Exomiser
//Personalizar proceso priorizazion y anotación de variantes con el archivo yml 
params.yaml_config = 
"/Volumes/EXTERNAL_USB/TFM/exomiser-cli-14.0.0/examples/HII_analysis_panel_nofilter.yml"
params.exomiser_jar = 
"/Volumes/EXTERNAL_USB/TFM/exomiser-cli-14.0.0/exomiser-cli-14.0.0.jar"
params.spring_config = 
"file:/Volumes/EXTERNAL_USB/TFM/exomiser-cli-14.0.0/application.properties"
params.memory_opts = "-Xms2g -Xmx4g"
params.assembly = "hg38"
params.exomiser_output_dir = 
"/Volumes/EXTERNAL_USB/TFM/results/HPO-nofilter-HII"

// Proceso para descargar archivos necesarios
// Exomiser así como las bases de datos que utiliza deberá ser descargado de arcuerdo a las instrucciones de
// su documentación oficial https://exomiser.readthedocs.io/en/latest/installation.html
process DownloadFiles {
    output:
    path "${params.ref_dir}/hg38.fa"
    path "${params.ref_dir}/hg38.dict"
    path "${params.ref_dir}/hg38.fa.fai"
    path "${params.ref_dir}/Homo_sapiens_assembly38.dbsnp138.vcf"
    path "${params.ref_dir}/Homo_sapiens_assembly38.dbsnp138.vcf.idx"

    script:
    """
    mkdir -p ${params.ref_dir}

    # Descargar archivos de referencia solo si no existen
    if [ ! -f ${params.ref_dir}/hg38.fa ]; then
        wget -P ${params.ref_dir} 
https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
        gunzip ${params.ref_dir}/hg38.fa.gz
    fi

    # Indexar referencia solo si no existe el índice
    if [ ! -f ${params.ref_dir}/hg38.fa.fai ]; then
        samtools faidx ${params.ref_dir}/hg38.fa
    fi

    # Crear diccionario de referencia solo si no existe
    if [ ! -f ${params.ref_dir}/hg38.dict ]; then
        gatk CreateSequenceDictionary R=${params.ref_dir}/hg38.fa 
O=${params.ref_dir}/hg38.dict
    fi

    # Descargar archivos de sitios conocidos solo si no existen
    if [ ! -f ${params.ref_dir}/Homo_sapiens_assembly38.dbsnp138.vcf ]; 
then
        wget -P ${params.ref_dir} 
https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf
    fi

    if [ ! -f ${params.ref_dir}/Homo_sapiens_assembly38.dbsnp138.vcf.idx 
]; then
        wget -P ${params.ref_dir} 
https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx
    fi
    """
}

// Crear los directorios de salida
process SetupDirectories {
    output:
    path "${params.aligned_reads_dir}"
    path "${params.results_dir}"
    path "${params.data_dir}"
    path "${params.exomiser_output_dir}"

    script:
    """
    mkdir -p ${params.aligned_reads_dir} ${params.results_dir} 
${params.data_dir} ${params.exomiser_output_dir}
    """
}

// Define el proceso para ejecutar FastQC
process RunFastQC {
    input:
    path file1
    path file2

    output:
    path "${params.results_dir}/FASTQC_analysis/*.zip" into fastqc_reports

    script:
    """
    mkdir -p ${params.results_dir}/FASTQC_analysis/
    fastqc ${file1} -o ${params.results_dir}/FASTQC_analysis/
    fastqc ${file2} -o ${params.results_dir}/FASTQC_analysis/
    """
}

// Define el proceso para el mapeo con BWA-MEM
process BwaMem {
    input:
    path file1
    path file2
    val sample_name

    output:
    path "${params.aligned_reads_dir}/${sample_name}.paired.sam"

    script:
    """
    bwa index ${params.ref_dir}/hg38.fa
    bwa mem -t 4 -R 
"@RG\\tID:${sample_name}\\tPL:ILLUMINA\\tSM:${sample_name}" 
${params.ref_dir}/hg38.fa ${file1} ${file2} > 
${params.aligned_reads_dir}/${sample_name}.paired.sam
    """
}

// Define el proceso para marcar duplicados y ordenar
process MarkDuplicates {
    input:
    path "${params.aligned_reads_dir}/${sample_name}.paired.sam"
    val sample_name

    output:
    path 
"${params.aligned_reads_dir}/${sample_name}_sorted_dedup_reads.bam"

    script:
    """
    gatk MarkDuplicatesSpark -I 
${params.aligned_reads_dir}/${sample_name}.paired.sam -O 
${params.aligned_reads_dir}/${sample_name}_sorted_dedup_reads.bam
    """
}

// Define el proceso para la recalibración de calidad de base
process BaseRecalibrator {
    input:
    path 
"${params.aligned_reads_dir}/${sample_name}_sorted_dedup_reads.bam"
    val sample_name

    output:
    path "${params.data_dir}/${sample_name}_recal_data.table"

    script:
    """
    gatk BaseRecalibrator -I 
${params.aligned_reads_dir}/${sample_name}_sorted_dedup_reads.bam -R 
${params.ref_dir}/hg38.fa --known-sites 
${params.ref_dir}/Homo_sapiens_assembly38.dbsnp138.vcf -O 
${params.data_dir}/${sample_name}_recal_data.table
    """
}

// Define el proceso para aplicar la recalibración de calidad
process ApplyBQSR {
    input:
    path 
"${params.aligned_reads_dir}/${sample_name}_sorted_dedup_reads.bam"
    path "${params.data_dir}/${sample_name}_recal_data.table"
    val sample_name

    output:
    path 
"${params.aligned_reads_dir}/${sample_name}_sorted_dedup_bqsr_reads.bam"

    script:
    """
    gatk ApplyBQSR -I 
${params.aligned_reads_dir}/${sample_name}_sorted_dedup_reads.bam -R 
${params.ref_dir}/hg38.fa --bqsr-recal-file 
${params.data_dir}/${sample_name}_recal_data.table -O 
${params.aligned_reads_dir}/${sample_name}_sorted_dedup_bqsr_reads.bam
    """
}

// Define el proceso para el llamado de variantes
process HaplotypeCaller {
    input:
    path 
"${params.aligned_reads_dir}/${sample_name}_sorted_dedup_bqsr_reads.bam"
    val sample_name

    output:
    path "${params.results_dir}/${sample_name}.raw_variants.vcf.gz"

    script:
    """
    gatk HaplotypeCaller -R ${params.ref_dir}/hg38.fa -I 
${params.aligned_reads_dir}/${sample_name}_sorted_dedup_bqsr_reads.bam -O 
${params.results_dir}/${sample_name}.raw_variants.vcf
    bgzip -f ${params.results_dir}/${sample_name}.raw_variants.vcf
    """
}

// Proceso para ejecutar Exomiser en cada archivo VCF generado
process runExomiser {
    input:
    path vcf_file

    output:
    path "${params.exomiser_output_dir}"

    script:
    def base_name = vcf_file.getBaseName().split('\\.')[0]
    def specific_output_dir = 
"${params.exomiser_output_dir}/sample-${base_name}"

    """
    mkdir -p ${specific_output_dir}
    
    java ${params.memory_opts} -jar ${params.exomiser_jar} \\
        --analysis ${params.yaml_config} \\
        --assembly ${params.assembly} \\
        --vcf ${vcf_file} \\
        --spring.config.location=${params.spring_config}
    
    mv ${params.exomiser_output_dir}/${base_name}* ${specific_output_dir}/
    """
}



// Flujo de trabajo principal que ejecuta cada proceso en secuencia
workflow {
    // Primero, descargamos los archivos necesarios
    DownloadFiles()

    // Configuramos los directorios de salida
    setup_dirs = SetupDirectories()

    Channel
        .fromFilePairs("${params.reads_dir}/*_1.fastq.gz", size: 2)
        .map { pair -> [pair[0].baseName - '_1', pair[0], pair[1]] }
        .set { sample_channel }

    sample_channel
        | RunFastQC
        | BwaMem
        | MarkDuplicates
        | BaseRecalibrator
        | ApplyBQSR
        | HaplotypeCaller

    // Canal para recoger los archivos VCF generados
    Channel
        .fromPath("${params.results_dir}/*.raw_variants.vcf.gz")
        .set { vcf_files }

    // Ejecuta Exomiser en cada archivo VCF comprimido
    vcf_files | runExomiser
}

