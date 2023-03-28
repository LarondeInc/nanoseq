/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

////////////////////////////////////////////////////
/* --          VALIDATE INPUTS                 -- */
////////////////////////////////////////////////////

// Check input path parameters to see if they exist
checkPathParamList = [ params.input, params.multiqc_config ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters (missing protocol or profile will exit the run.)
if (params.input) { 
    ch_input = file(params.input) 
} else {
    exit 1, 'Input samplesheet not specified!'
}

// Function to check if running offline
def isOffline() {
    try {
        return NXF_OFFLINE as Boolean
    }
    catch( Exception e ) {
        return false
    }
}

def ch_guppy_model  = Channel.empty()
def ch_guppy_config = Channel.empty()

if (params.protocol != 'DNA' && params.protocol != 'cDNA' && params.protocol != 'directRNA') {
    exit 1, "Invalid protocol option: ${params.protocol}. Valid options: 'DNA', 'cDNA', 'directRNA'"
}

if (!params.skip_basecalling) {
    if (!params.guppy_config) {
        if (!params.flowcell) { exit 1, "Please specify a valid flowcell identifier for basecalling!" }
        if (!params.kit)      { exit 1, "Please specify a valid kit identifier for basecalling!"      }
    } else if (file(params.guppy_config).exists()) {
        ch_guppy_config = Channel.fromPath(params.guppy_config)
    }

    if (params.guppy_model) {
        if (file(params.guppy_model).exists()) {
            ch_guppy_model = Channel.fromPath(params.guppy_model)
        }
    }
} else {
    if (!params.skip_demultiplexing) {
        if (!params.barcode_kit) {
            params.barcode_kit = 'Auto'
        }

        def qcatBarcodeKitList = ['Auto', 'RBK001', 'RBK004', 'NBD103/NBD104',
                                'NBD114', 'NBD104/NBD114', 'PBC001', 'PBC096',
                                'RPB004/RLB001', 'PBK004/LWB001', 'RAB204', 'VMK001', 'DUAL']

        if (params.barcode_kit && qcatBarcodeKitList.contains(params.barcode_kit)) {
            if (params.input_path) {
                ch_input_path = Channel.fromPath(params.input_path, checkIfExists: true)
            } else {
                exit 1, "Please specify a valid input fastq file to perform demultiplexing!"
            }
        } else {
            exit 1, "Please provide a barcode kit to demultiplex with qcat. Valid options: ${qcatBarcodeKitList}"
        }
    }
}

if (!params.skip_alignment) {
    if (params.aligner != 'minimap2' && params.aligner != 'graphmap2') {
        exit 1, "Invalid aligner option: ${params.aligner}. Valid options: 'minimap2', 'graphmap2'"
    }
    if (params.protocol != 'DNA' && params.protocol != 'cDNA' && params.protocol != 'directRNA') {
        exit 1, "Invalid protocol option: ${params.protocol}. Valid options: 'DNA', 'cDNA', 'directRNA'"
    }
}

if (params.call_variants) {
    if (params.protocol != 'DNA') {
        exit 1, "Invalid protocol option: ${params.protocol}. Valid options: 'DNA'"
    }
    if (!params.skip_vc && params.variant_caller != 'medaka' && params.variant_caller != 'deepvariant' && params.variant_caller != 'pepper_margin_deepvariant') {
        exit 1, "Invalid variant caller option: ${params.variant_caller}. Valid options: 'medaka', 'deepvariant' or 'pepper_margin_deepvariant'"
    }
    if (!params.skip_sv && params.structural_variant_caller != 'sniffles' && params.structural_variant_caller != 'cutesv') {
        exit 1, "Invalid structural variant caller option: ${params.structural_variant_caller}. Valid options: 'sniffles', 'cutesv"
    }
    if (!params.skip_vc && params.enable_conda && params.variant_caller != 'medaka') {
        exit 1, "Conda environments cannot be used when using the deepvariant or pepper_margin_deepvariant tools. Valid options: 'docker', 'singularity'"
    }
}

////////////////////////////////////////////////////
/* --          CONFIG FILES                    -- */
////////////////////////////////////////////////////

ch_multiqc_config        = file("$baseDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

////////////////////////////////////////////////////
/* --    IMPORT LOCAL MODULES/SUBWORKFLOWS     -- */
////////////////////////////////////////////////////

include { GET_TEST_DATA         } from '../modules/local/get_test_data'
include { GET_NANOLYSE_FASTA    } from '../modules/local/get_nanolyse_fasta'
include { GUPPY                 } from '../modules/local/guppy'
include { DEMUX_FAST5           } from '../modules/local/demux_fast5'
include { QCAT                  } from '../modules/local/qcat'
include { BAM_RENAME            } from '../modules/local/bam_rename'
include { BAMBU                 } from '../modules/local/bambu'
include { MULTIQC               } from '../modules/local/multiqc'

/*
 * SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
 */

include { INPUT_CHECK                      } from '../subworkflows/local/input_check'
include { PREPARE_GENOME                   } from '../subworkflows/local/prepare_genome'
include { QCBASECALL_PYCOQC_NANOPLOT       } from '../subworkflows/local/qcbasecall_pycoqc_nanoplot'
include { QCFASTQ_NANOPLOT_FASTQC          } from '../subworkflows/local/qcfastq_nanoplot_fastqc'
include { ALIGN_GRAPHMAP2                  } from '../subworkflows/local/align_graphmap2'
include { ALIGN_MINIMAP2                   } from '../subworkflows/local/align_minimap2'
include { BAM_SORT_INDEX_SAMTOOLS          } from '../subworkflows/local/bam_sort_index_samtools'
include { SHORT_VARIANT_CALLING            } from '../subworkflows/local/short_variant_calling'
include { STRUCTURAL_VARIANT_CALLING       } from '../subworkflows/local/structural_variant_calling'
include { BEDTOOLS_UCSC_BIGWIG             } from '../subworkflows/local/bedtools_ucsc_bigwig'
include { BEDTOOLS_UCSC_BIGBED             } from '../subworkflows/local/bedtools_ucsc_bigbed'
include { QUANTIFY_STRINGTIE_FEATURECOUNTS } from '../subworkflows/local/quantify_stringtie_featurecounts'
include { DIFFERENTIAL_DESEQ2_DEXSEQ       } from '../subworkflows/local/differential_deseq2_dexseq'
include { RNA_MODIFICATION_XPORE_M6ANET    } from '../subworkflows/local/rna_modifications_xpore_m6anet'
include { RNA_FUSIONS_JAFFAL               } from '../subworkflows/local/rna_fusions_jaffal'

////////////////////////////////////////////////////
/* --    IMPORT NF-CORE MODULES/SUBWORKFLOWS   -- */
////////////////////////////////////////////////////

/*
 * MODULE: Installed directly from nf-core/modules
 */
include { NANOLYSE                    } from '../modules/nf-core/modules/nanolyse/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'

/*
 * SUBWORKFLOW: Consisting entirely of nf-core/modules
 */

////////////////////////////////////////////////////
/* --           RUN MAIN WORKFLOW              -- */
////////////////////////////////////////////////////

// Info required for completion email and summary
def multiqc_report      = []

workflow NANOSEQ{

    // Pre-download test-dataset to get files for '--input_path' parameter
    // Nextflow is unable to recursively download directories via HTTPS
    if (workflow.profile.contains('test') && !workflow.profile.contains('vc')) {
        if (!params.skip_basecalling || !params.skip_modification_analysis) {
            if (!isOffline()) {
                GET_TEST_DATA ()
                if (params.skip_modification_analysis) {
                    GET_TEST_DATA.out.ch_input_fast5s_path
                        .set { ch_input_path }
                } else {
                    GET_TEST_DATA.out.ch_input_dir_path
                        .set { ch_input_path }
                }
            } else {
                exit 1, "NXF_OFFLINE=true or -offline has been set so cannot download and run any test dataset!"
            }
        } else {
            if (params.input_path) {
                ch_input_path = Channel.fromPath(params.input_path, checkIfExists: true)
            } else {
                ch_input_path = 'not_changed'
            }
        }
    } else {
        if (params.input_path) {
            ch_input_path = Channel.fromPath(params.input_path, checkIfExists: true)
        } else {
            ch_input_path = 'not_changed'
        }
    }

    /*
     * Create empty software versions channel to mix
     */
    ch_software_versions = Channel.empty()

    /*
     * SUBWORKFLOW: Read in samplesheet, validate and stage input files
     */
    INPUT_CHECK ( ch_input, ch_input_path )
        .set { ch_sample }

    if (!params.skip_basecalling) {
        ch_sample
            .first()
            .map { it[0] }
            .set { ch_sample_name }

        /*
         * MODULE: Basecalling and demultipexing using Guppy
         */
        GUPPY ( ch_input_path, ch_sample_name, ch_guppy_config.ifEmpty([]), ch_guppy_model.ifEmpty([]) )
        ch_guppy_summary = GUPPY.out.summary
        ch_software_versions = ch_software_versions.mix(GUPPY.out.versions.ifEmpty(null))

        if (params.skip_demultiplexing) {
            ch_sample
                .map { it -> [ it[0], it[0].id, it[2], it[3], it[4], it[5] ] }
                .set { ch_sample }
        }

        GUPPY.out.fastq
            .flatten()
            .map { it -> [ it, it.baseName.substring(0,it.baseName.lastIndexOf('.')) ] }
            .join(ch_sample, by: 1) // join on barcode
            .map { it -> [ it[2], it[1], it[3], it[4], it[5], it[6] ] }
            .set { ch_fastq }
        if (params.output_demultiplex_fast5) {

            /*
            * MODULE: Demultiplex fast5 files using ont_fast5_api/demux_fast5
            */
            DEMUX_FAST5 ( ch_input_path, ch_guppy_summary )
            ch_software_versions = ch_software_versions.mix(DEMUX_FAST5.out.versions.ifEmpty(null))
        }
    } else {
        ch_guppy_summary = Channel.empty()

        if (!params.skip_demultiplexing) {

            /*
             * MODULE: Demultipexing using qcat
             */
            QCAT ( ch_input_path )
            ch_fastq = Channel.empty()
            QCAT.out.fastq
                .flatten()
                .map { it -> [ it, it.baseName.substring(0,it.baseName.lastIndexOf('.'))] }
                .join(ch_sample, by: 1) // join on barcode
                .map { it -> [ it[2], it[1], it[3], it[4], it[5], it[6] ] }
                .set { ch_fastq }
            ch_software_versions = ch_software_versions.mix(QCAT.out.versions.ifEmpty(null))
        } else {
            if (!params.skip_alignment) {
                ch_sample
                    .map { it -> if (it[6].toString().endsWith('.gz')) [ it[0], it[6], it[2], it[1], it[4], it[5] ] }
                    .set { ch_fastq }
            } else {
                ch_fastq = Channel.empty()
            }
        }
    }

    ch_pycoqc_multiqc = Channel.empty()
    ch_fastqc_multiqc = Channel.empty()
    
    ch_samtools_multiqc = Channel.empty()
    
    ch_featurecounts_gene_multiqc       = Channel.empty()
    ch_featurecounts_transcript_multiqc = Channel.empty()
    if (!params.skip_multiqc) {
        workflow_summary    = WorkflowNanoseq.paramsSummaryMultiqc(workflow, summary_params)
        ch_workflow_summary = Channel.value(workflow_summary)

        /*
         * MODULE: MultiQC
         */
        MULTIQC (
        ch_multiqc_config,
        ch_multiqc_custom_config.collect().ifEmpty([]),
        ch_pycoqc_multiqc.collect().ifEmpty([]),
        ch_fastqc_multiqc.ifEmpty([]),
        ch_samtools_multiqc.collect().ifEmpty([]),
        ch_featurecounts_gene_multiqc.ifEmpty([]),
        ch_featurecounts_transcript_multiqc.ifEmpty([]),
        CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect(),
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml')
        )
    }
}

////////////////////////////////////////////////////
/* --              COMPLETION EMAIL            -- */
////////////////////////////////////////////////////

workflow.onComplete {
    if (params.email) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
        //Completion.email(workflow, params, params.summary_params, log, multiqc_report)
    }
//    Completion.summary(workflow, params, log)
    NfcoreTemplate.summary(workflow, params, log)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
