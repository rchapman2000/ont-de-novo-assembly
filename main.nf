#!/usr/bin/env/ nextflow

nextflow.enable.dsl = 2
def helpMessage() {
    log.info"""
    Nanopore Sequencing De Novo Assembly Pipeline

    Takes an input of raw reads produced by Oxford Nanopore sequencing, trims
    the reads using Porechop, and performs De Novo Assembly using the assembler
    Flye.

    USAGE: nextflow run main.nf [OPTIONS] --input INPUT_DIR --output OUTPUT_DIR --model MEDAKA_MODEL
    OPTIONS:

    --input INPUT_DIR - [Required] A directory containing paired-end fastq files

    --output OUTPUT_DIR - [Required] A directory to place output files (If not existing, pipeline will create)
    
    --model MEDAKA_MODEL - [Required] Medaka requires a 'model' which corresponds to the flow-cell type/basecaller parameters used to correct for errors typical to that technology. Use 'medaka tools list_models' to find this.

    OPTIONAL:
        --trimONTAdapters - Enables ONT Adapter/Barcode trimming using Porechop [Default = off]

        --minReadLen INT - If supplied, the pipeline will perform length filtering using Chopper excluding reads less than this size [Default = off]

        --maxReadLen INT - If supplied, the pipeline will perform legnth filtering using Chopper excluding reads greater than this size [Default = off] 

        --medakaBatchSize INT - Medaka uses a lot of GPU memory, and if you're assembly is large enough it may cause Medaka to crash. Reducing the batch size will help solve this issue. [Default = 100]

        --preGuppy5 - Flye handles data generated version of Guppy < 5.0 differently. Supply this parameter if your data was generated pre-Guppy 5.0

        --threads INT - The number of CPU threads that can be use to run pipeline tools in parallel
"""
}

// If the help parameter is supplied, link display the help message
// and quit the pipeline
params.help = false
if (params.help){
    helpMessage()
    exit 0
}

def createSummaryHeader (trim, minLenFilt, maxLenFilt) {
    FinalHeader = "Sample,Raw Reads,Average Raw Read Length,"

    if (trim != false) {
        FinalHeader = FinalHeader + "Read Post Porechop Trimming,Average Trimmed Read Length,"
    }

    if (minLenFilt != 0 || maxLenFilt != 0) {
        FinalHeader = FinalHeader + "Reads Post Length Filter,Average Filtered Read Length,"
    }

    FinalHeader = FinalHeader + "Draft Contigs,Average Draft Contig Length,Corrected Contigs,Average Corrected Contig Length"

    return FinalHeader
}

// Defines input parameters. Setting parameters
// to false by default allows us to check whether
// these parameters have been set by the user (useful
// for required parameters or binary parameters)
params.input = false
params.output = false
params.threads = 1
params.model = false
params.preGuppy5 = false
params.medakaBatchSize = 100
params.trimONTAdapters = false
params.minReadLen = 0
params.maxReadLen = 0

// Import pipeline modeuls from 'modules.nf'
include { Setup } from "./modules.nf"
include { Collect_Raw_Read_Stats } from "./modules.nf"
include { QC_Report } from "./modules.nf"
include { QC_Report_Filtered } from "./modules.nf"
include { Porechop_Trimming } from "./modules.nf"
include { Length_Filtering } from "./modules.nf"
include { Flye_Assembly } from "./modules.nf"
include { Medaka_Correct } from "./modules.nf"
include { Write_Summary } from "./modules.nf"

// Checks the input parameter
inDir = ''
if (params.input == false) {
    // If the parameter is not set, notify the user and exit.
    println "ERROR: No input directory provided. Pipeline requires an input directory."
    exit(1)
}
else if (!(file(params.input).isDirectory())) {
    // If the input directory is not set, notify the user and exit.
    println "ERROR: ${params.input} is not an existing directory."
    exit(1)
}
else {
    // If the parameter is set, convert the value provided to a file type
    // to get the absolute path, and then convert back to a string to be
    // used in the pipeline.
    inDir = file(params.input).toString()
    println "Input Directory: ${inDir}"
}

// Create a channel for hte input files.
inputFiles_ch = Channel
    // Pull from pairs of files (illumina fastq files denoted by having R1 or R2 in
    // the file name).
    .fromPath("${inDir}/*.fastq*")
    // The .fromPath() function spits out each file, and we can use
    // the build-in nextflow functionality to get the name of the file
    // and place it into a tuple along with its corresponding file. That way
    // we can use the file name in the pipeline to name any output files. The
    // getSimpleName() function is used as opposed to getName(), as the
    // getName() function only removes the last .extension and would
    // cause issues if the fastqs were gzipped.
    .map { it -> [it.getSimpleName(), it]}

// Checks the output parameter.
outDir = ''
if (params.output == false) {
    // If the parameter is not set, notify the user and exit.
    println "ERROR: No output directory provided. Pipeline requires an output directory."
    exit(1)
}
else {
    // If the parameter is set, ensure that the directory provided ends
    // in a trailing slash (to keep things consistent throughout) the
    // pipeline code.
    outDir = file(params.output).toString()
    println(outDir)
}

// Checks the model parameter.
model = ''
if (params.model == false) {
    // If the parameter is not set, notify the user and exit.
    //
    // This parameter requires a special input, so an example value is provided to the user
    // along with instructions to view the possible models using a medaka command.
    println "ERROR: no ONT model provided. Medaka requires a model in the format:"
    println ""
    println "{pore}_{device}_{caller variant}_{caller version}"
    println ""
    println "To see existing models enter: medaka tools list_models"
    exit(1)
}
else {
    // If the parameter is provided, set the model variable equal to the
    // provided value.
    model = params.model
}

// Flye handles reads produced by older versions of guppy (prior to v5.0)
// differently. The user denotes this using a command line option. Thus,
// we can handle this by containing the parameter in a string and then placing
// it into the flye command run by nextflow.
//
// By default, we assume that the reads were generated with guppy v5.0 or greater.
// Thus, we create a variable holding the flag "--nano-hq"
nanoporeReadType = "--nano-hq"
if (params.preGuppy5) {
    // If the user supplied the parameter denoting their reads
    // were generated prior to guppy v5.0, change the 
    // value of the nanoporeReadType variable to '--nano-raw.'
    nanoporeReadType = "--nano-raw"
}

// Handles the --trimONTAdapters parameter. Updates
// a variable to be written to the analysis parameters file
// if trimming was enabled.
trimmingEnabled = "DISABLED"
if (params.trimONTAdapters != false) {
    trimmingEnabled = "ENABLED"
}


// Handles the --minReadLen parameter.
minReadLenVal = "DISABLED"
minReadLenParam = ""
if (params.minReadLen != 0) {
    // Updates a variable to be written to the analysis parameters file
    // if a minimum read length threshold was supplied.
    minReadLenVal = params.minReadLen + " bp"

    // As well if a value was supplied, it creates a variable
    // containing teh string to be passed to chopper. By doing so,
    // this allows us to handle the case when a minimum length
    // was supplied but not a maximum length (as the maximum length
    // parameter string will be blank and not affect the command).
    minReadLenParam = "-l ${params.minReadLen}"
}


// Handles the --maxReadLen parameter.
maxReadLenVal = "DISABLED"
maxReadLenParam = ""
if (params.maxReadLen != 0) {

    // Updates a variable to be written to the analysis parameters file
    // if a maximum read length threshold was supplied.
    maxReadLenVal = params.maxReadLen + " bp"

    // As well if a value was supplied, it creates a variable
    // containing teh string to be passed to chopper. By doing so,
    // this allows us to handle the case when a minimum length
    // was supplied but not a maximum length (as the maximum length
    // parameter string will be blank and not affect the command).
    maxReadLenParam = "--maxlength ${params.maxReadLen}"
}

// Generates the summary header based on the parameters supplied. Allows for the header
// to be dynamic depending on the statistics collected.
summaryHeader = createSummaryHeader(params.trimONTAdapters, params.minReadLen, params.maxReadLen)

// To keep track of the summary, a string will be passed between each process that contains
// metrics generated by each. At the end of each process, the metrics collected will be
// appended to the end and returned as output from that process. At the end of the workflow,
// the metrics will be written to the summary file.
workflow {
    // Creates files to store analysis parameters and statistics.
    Setup( trimmingEnabled, minReadLenVal, maxReadLenVal, model, params.medakaBatchSize, summaryHeader, outDir )

    QC_Report( inputFiles_ch, outDir )

    Collect_Raw_Read_Stats( inputFiles_ch )

    if (params.trimONTAdapters != false && (params.minReadLen != 0 || params.maxReadLen != 0)) {
        Porechop_Trimming( Collect_Raw_Read_Stats.out[0], outDir, Collect_Raw_Read_Stats.out[1] )

        Length_Filtering( Porechop_Trimming.out[0], minReadLenParam, maxReadLenParam, outDir, Porechop_Trimming.out[2] )

        QC_Report_Filtered( Length_Filtering.out[0], outDir )

        // Performs De Novo Assembly of the filtered reads using Flye
        Flye_Assembly( Length_Filtering.out[0], nanoporeReadType, outDir, params.threads, Length_Filtering.out[1] )

    }
    else if (params.trimONTAdapters != false && params.minReadLen == 0 && params.maxReadLen == 0) {

        Porechop_Trimming( Collect_Raw_Read_Stats.out[0], outDir, Collect_Raw_Read_Stats.out[1] )

        QC_Report_Filtered( Porechop_Trimming.out[0], outDir )

       // Performs De Novo Assembly of the trimmed reads using Flye
        Flye_Assembly( Porechop_Trimming.out[0], nanoporeReadType, outDir, params.threads, Porechop_Trimming.out[2] )
    }
    else if (params.trimONTAdapters == false && (params.minReadLen != 0 || params.maxReadLen != 0)) {
        Length_Filtering( Collect_Raw_Read_Stats.out[0], minReadLenParam, maxReadLenParam, outDir, Collect_Raw_Read_Stats.out[1] )

        QC_Report_Filtered( Length_Filtering.out[0], outDir )

        // Performs De Novo Assembly of the filtered reads using Flye
        Flye_Assembly( Length_Filtering.out[0], nanoporeReadType, outDir, params.threads, Length_Filtering.out[1] )
    }
    else {
        // Performs De Novo Assembly of the reads using Flye
        Flye_Assembly( Collect_Raw_Read_Stats.out[0], nanoporeReadType, outDir, params.threads, Collect_Raw_Read_Stats.out[1] )
    }

    // Corrects the assembled contigs using Medaka
    Medaka_Correct( Flye_Assembly.out[0], outDir, params.threads, model, params.medakaBatchSize, Flye_Assembly.out[2] )

    // Writes statistics to an analysis summary file.
    Write_Summary( Medaka_Correct.out[1], outDir)
}