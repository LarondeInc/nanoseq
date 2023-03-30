process GUPPY {
    label 'process_medium'

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

    output:
    path "fastq/*.fastq.gz"                    , emit: fastq
    tuple val(meta), path("basecalling/*.txt") , emit: summary
    path "basecalling/*"                       , emit: called
    path "versions.yml"                        , emit: versions

    script:
    def fast5_dir_path = workflow.profile.contains('test') ? "input_path" : "$input_path"
    def save_path = "./fastq"
    def trim_barcodes = params.trim_barcodes ? "--trim_barcodes" : ""
    def barcode_ends = params.barcode_both_ends ? "--require_barcodes_both_ends" : ""
    def config   = ""
    if (params.guppy_config) config = file(params.guppy_config).exists() ? "--config ./$guppy_config" : "--config $params.guppy_config"
    """
    guppy_basecaller \\
        --input_path $fast5_dir_path \\
	    --save_path $save_path \\
        --compress_fastq \\
        --recursive \\
	    --do_read_splitting \\
	    --max_read_depth 2 \\
	    --min_score_read_splitting 58 \\
        $config
    """
}

