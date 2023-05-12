
params.options = [:]

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
    def barcode_kit = params.barcode_kit ? "--barcode_kits $params.barcode_kit" : ""
    def config   = ""
    if (params.guppy_config) config = file(params.guppy_config).exists() ? "--config ./$guppy_config" : "--config $params.guppy_config"
    """
    guppy_basecaller \\
        --input_path $fast5_dir_path \\
        --save_path ./basecalling \\
        --compress_fastq \\
        --recursive \\
	    --do_read_splitting \\
	    --min_score_read_splitting 58 \\
        $config
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        guppy: \$(echo \$(guppy_basecaller --version 2>&1) | sed -r 's/.{81}//')
    END_VERSIONS


    ## Concatenate fastq files
    mkdir fastq
    cd basecalling
    ls completed/*
    pwd
    if [ "\$(find . -type d -name "barcode*" )" != "" ]
    then
        for dir in barcode*/
        do
            dir=\${dir%*/}
            cat \$dir/*.fastq.gz > ../fastq/\$dir.fastq.gz
        done
    else
        cat *.fastq.gz > ../fastq/${meta}.fastq.gz
    fi

    """
}

