// Creates a parameters file and a summary file to 
// be added to later
process Setup {
    input:
        // The provided medaka model
        val model
        // The batch size option supplied to medaka
        val batchSize
        // The output directory to be used.
        val outDir
        
    output:
        // The parameters file created.
        file "analysis-parameters.txt"
        // The blank summary file to be added to.
        file "stats-summary.csv"

    publishDir "${outDir}", mode: 'copy'

    script:
    /*
    Creates a parameters file (in case the user needs to look back at how they ran the analysis)
    as well as a blank summary file, which will later be populated with data from each
    sample.

    The parameters file contains:
        1. The medaka model provided.
        2. The batch size parameter provided to medaka.

    The summary file will always contain:
        1. The sample
        2. Raw Reads
        3. Average Raw Read Length
        4. Trimmed Reads
        5. Average Trimmed Read Length
        6. Number of Draft Contigs
        7. Average Draft Contig Lengths
        8. Number of Corrected Contigs
        9. Average Corrected Contig Length
    */
    """
    #!/bin/bash

    touch analysis-parameters.txt

    echo "Medaka Model : ${model}" >> analysis-parameters.txt
    echo "Medaka Batch Size: ${batchSize}" >> analysis-parameters.txt

    touch stats-summary.csv

    echo "Sample,Raw Reads,Average Raw Read Length,Trimmed Reads,Average Trimmed Read Length,Draft Contigs,Average Draft Contig Length,Corrected Contigs,Average Corrected Contig Length" > stats-summary.csv
    """
}

// Detects and trims ONT adapter/barcode sequences from raw reads using Porechop
process Porechop_Trimming {
    input:
        // Tuple contains the sample base name and
        // the reads to be trimmed
        tuple val(base), file(reads)
        // The name of the output directory to write
        // output files to.
        val outDir

    output:
        // Tuple contains the sample base name and 
        // the trimmed reads.
        tuple val(base), file("${base}-trimmed.fastq")
        // A string containing summary metrics generated 
        // during this process.
        env summary

    publishDir "${outDir}", mode: 'copy'

    script:
    /*
    Uses Porechop to detect and trim ONT adapter/barcodes sequences found
    within the provided reads. 

    As well, the number of raw reads, average raw read length, number of 
    trimmed reads, and average trimmed read length are calculated
    and added to the summary data.
    */
    """
    #!/bin/bash

    raw_reads=\$((\$(zcat -f ${reads} | wc -l)/4))

    avg_raw_read_len=\$(bioawk -c fastx '{ totalBases += length(\$seq); totalReads++ } END{ print totalBases/totalReads }' ${reads})

    porechop -i ${reads} -o ${base}-trimmed.fastq --verbosity 3 > ${base}-porechop-report.txt

    trimmed_reads=\$((\$(gunzip -c ${base}-trimmed.fastq | wc -l)/4))

    avg_trimmed_read_len=\$(bioawk -c fastx '{ totalBases += length(\$seq); totalReads++ } END{ print totalBases/totalReads }' ${base}-trimmed.fastq)

    summary="${base},\$raw_reads,\$avg_raw_read_len,\$trimmed_reads,\$avg_trimmed_read_len"
    """
}


// Perfroms a De Novo assembly of the ONT reads using
// Flye
process Flye_Assembly {
    input:
        // Tuple contains the sample base name
        // and the reads to be assembled.
        tuple val(base), file(reads)
        // Flye requires that reads generated
        // from verions of Guppy prior to
        // v5.0 be handled separatedly noted via a command line option. 
        // This parameter specifies the type of reads present.
        val readParam
        // The name of the output directory to copy
        // generated output files.
        val outDir
        // The number of threads to be used.
        val threads
        // A string containing previously generated
        // statistics to be appended to.
        val existingSummary

    output:
        // Tuple contains the sample base name, the reads used to 
        // generate the assembly, and the draft assembly produced by
        // Flye
        tuple val(base), file(reads), file("${base}-draft-assembly.fasta")
        // The directory produced by Flye
        file("${base}-flye")
        // A string containing the previously generated assembly statistics
        // as well as statistics generated during this step.
        env summary
        

    publishDir "${outDir}", mode: 'copy'

    script:
    /*
    Uses flye to assemble the reads.

    If the reads can be assembled, there will be an 'assembly.fasta' file created
    in the flye output directory. If this file exists, then it will be moved to 
    a file called SAMPLE-draft-assembly.fasta. If not, an empty file will be
    created (to prevent an error from occurring).

    Finally, the number of contigs generated and average contig length are collected and
    appended to the existing summary statistics.
    */
    """
    #!/bin/bash

    flye ${readParam} ${reads} --out-dir ${base}-flye --thread ${threads}

    if [[ -f "${base}-flye/assembly.fasta" ]]; then
        mv ${base}-flye/assembly.fasta ./${base}-draft-assembly.fasta
    else
        touch ./${base}-draft-assembly.fasta
    fi
    
    num_contigs=\$(grep ">" ${base}-draft-assembly.fasta | wc -l)

    avg_contig_len=\$(bioawk -c fastx '{ totalBases += length(\$seq); totalReads++ } END{ print totalBases/totalReads }' ${base}-draft-assembly.fasta)

    summary="${existingSummary},\$num_contigs,\$avg_contig_len"
    """
}


// Uses Medaka to Polish the contigs.
process Medaka_Correct {
    input:
        // Tuple contains the sample base name,
        // the reads used to generate the draft assembly,
        // and the draft assembly.
        tuple val(base), file(reads), file(draft)
        // The name of the output directory to copy
        // generated output files.
        val outDir
        // The number of threads available to
        // compute using.
        val threads
        // The medaka model provided.
        val model
        // The batch size to be used
        // by medaka when polishing.
        val batchSize
        // A string containing previously generated
        // statistics to be appended to.
        val existingSummary

    output:
        // Tuple contains the sample base name and the
        // corrected assembly in fasta format.
        tuple val(base), file("${base}-corrected-assembly.fasta")
        // A string containing statistics including
        // those generated during this step.
        env summary

    publishDir "${outDir}", mode: "copy"

    script:
    /*
    Uses the medaka_consensus pipeline to correct the draft contigs using the raw reads.

    If the medaka pipeline worked correctly, a file named 'consensus.fasta' will be
    present in the medaka output directory. If this file exists, then it is moved into
    a file named "SAMPLE-corrected-assembly.fasta". If not, an empty file is created (to
    prevent file not found errors).

    Finally, statistics on the number of corrected contigs and average corrected contig length
    are appended to the string containing existing statistics.
    */
    """
    #!/bin/bash

    medaka_consensus -i ${reads} -d ${draft} -o ${base}-medaka -t ${threads} -m ${model} -b ${batchSize}

    if [[ -f "${base}-medaka/consensus.fasta" ]]; then
        mv ${base}-medaka/consensus.fasta ./${base}-corrected-assembly.fasta
    else
        touch ./${base}-corrected-assembly.fasta
    fi

    num_corrected_contigs=\$(grep ">" ${base}-corrected-assembly.fasta | wc -l)

    avg_corrected_contig_len=\$(bioawk -c fastx '{ totalBases += length(\$seq); totalReads++ } END{ print totalBases/totalReads }' ${base}-corrected-assembly.fasta)

    summary="${existingSummary},\$num_corrected_contigs,\$avg_corrected_contig_len"
    """
}

// Writes a line to the summary file for the sample.
process Write_Summary {
    input:
        // Tuple contains the sample basename and forward/reverse reads (the basename
        // is the only value important to this function).
        val summary
        // The output directory.
        val outDir

    script:
    /*
    The summary string containing the statistics collected as the pipeline
    was run are appended to the summary file.
    */
    """
    #!/bin/bash

    echo "${summary}" >> ${outDir}/stats-summary.csv
    """  

}