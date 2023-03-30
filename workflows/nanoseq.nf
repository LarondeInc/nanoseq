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
checkPathParamList = [ params.input ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters (missing protocol or profile will exit the run.)
if (params.input) { 
    ch_input = file(params.input) 
} else {
    exit 1, 'Input folder is not specified!'
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

def ch_guppy_config = Channel.empty()

if (!params.skip_basecalling) {
    if (!params.guppy_config) {
        exit 1, 'No guppy config file!'
    } else if (file(params.guppy_config).exists()) {
        ch_guppy_config = Channel.fromPath(params.guppy_config)
    }
}

////////////////////////////////////////////////////
/* --          CONFIG FILES                    -- */
////////////////////////////////////////////////////

////////////////////////////////////////////////////
/* --    IMPORT LOCAL MODULES/SUBWORKFLOWS     -- */
////////////////////////////////////////////////////

include { GUPPY                 } from '../modules/local/guppy'

/*
 * SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
 */


////////////////////////////////////////////////////
/* --    IMPORT NF-CORE MODULES/SUBWORKFLOWS   -- */
////////////////////////////////////////////////////

/*
 * MODULE: Installed directly from nf-core/modules
 */

/*
 * SUBWORKFLOW: Consisting entirely of nf-core/modules
 */

////////////////////////////////////////////////////
/* --           RUN MAIN WORKFLOW              -- */
////////////////////////////////////////////////////

// Info required for completion email and summary

workflow NANOSEQ{

    // Pre-download test-dataset to get files for '--input_path' parameter
    // Nextflow is unable to recursively download directories via HTTPS
    if (workflow.profile.contains('test') && !workflow.profile.contains('vc')) {
        if (!params.skip_basecalling) {
            if (!isOffline()) {
                GET_TEST_DATA ()
                {
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
    
    if (!params.skip_basecalling) {
        //ch_sample
        //    .first()
        //    .map { it[0] }
        //    .set { ch_sample_name }
        ch_sample_name = "test"
        /*
         * MODULE: Basecalling and demultipexing using Guppy
         */
        GUPPY ( ch_input_path, ch_sample_name, ch_guppy_config.ifEmpty([]))
        ch_guppy_summary = GUPPY.out.summary
        ch_software_versions = ch_software_versions.mix(GUPPY.out.versions.ifEmpty(null))

        //if (params.skip_demultiplexing) {
        //    ch_sample
        //        .map { it -> [ it[0], it[0].id, it[2], it[3], it[4], it[5] ] }
        //        .set { ch_sample }
        //}
        ch_sample = "test"

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
        }
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
