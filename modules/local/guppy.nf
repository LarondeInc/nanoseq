include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process GUPPY {
    label 'process_medium'
    publishDir "${params.outdir}",
    	mode : params.publish_dir_mode,
    	saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process)) }

    if (params.guppy_gpu) {
        container = 'genomicpariscentre/guppy-gpu:6.4.6'
        clusterOptions = params.gpu_cluster_options
    } else {
        container = 'genomicpariscentre/guppy:6.4.6'
    }

    input:
    path(input_path), stageAs: 'input_path/*'
    val meta
    path guppy_config
    val barcode_kit

    output:
    path "fastq/*.fastq.gz"                    , emit: fastq
    tuple val(meta), path("basecalling/*.txt") , emit: summary
    path "basecalling/*"                       , emit: called
    path "versions.yml"                        , emit: versions

    script:
    def fast5_dir_path = workflow.profile.contains('test') ? "input_path" : "$input_path"
    def save_path = "./fastq"
    def barcode_kit = params.barcode_kit ? "--barcode_kits $params.barcode_kit" : ""
    def config   = ""
    if (params.guppy_config) config = file(params.guppy_config).exists() ? "--config ./$guppy_config" : "--config $params.guppy_config"
    """
    guppy_basecaller \\
        --input_path $fast5_dir \\
	--save_path ./basecalling \\
        --compress_fastq \\
        --recursive \\
	--do_read_splitting \\
	--min_score_read_splitting 58 \\
        $config
    """
}

