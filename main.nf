#! /usr/bin/env nextflow

/*************************************
 nextflow nanoQCtrim
 *************************************/

racon_container = 'quay.io/biocontainers/racon:1.4.13--he513fc3_0'
medaka_container = 'quay.io/biocontainers/medaka:1.0.3--py36hbecb4b7_1'


 def helpMessage() {
     log.info isuGIFHeader()
     log.info """
      Usage:
      The typical command for running the pipeline are as follows:

      nextflow run isugifNF/assemblyStats --genomes "*fasta" --outdir newStats3 --threads 16 -profile condo,singularity


      Mandatory arguments:

      --genome                      genome assembly fasta file to run stats on. (./data/*.fasta)
      -profile singularity           as of now, this workflow only works using singularity and requires this profile [be sure singularity is in your path]
      --reads                       Raw nanopore reads to use in polish
      Optional arguments:
      --outdir                       Output directory to place final output
      --threads                      Number of CPUs to use during the NanoPlot job [16]
      --queueSize                    Maximum number of jobs to be queued [18]

      --help                         This usage statement.
     """
}

     // Show help message
if (params.help) {
   helpMessage()
   exit 0
}

// create a channel for the genome
/*
  Channel
   .fromPath(params.genomes)
   .map { file -> tuple(file.simpleName, file) }
   .into { genome_runMinimap2; genome_runAssemblathonStats; genome_BUSCO }
*/
   Channel
    .fromPath(params.genomes)
    .map { file -> file.simpleName}
    .into { genomeLabel_runMinimap2; genomeLabel_runAssemblathonStats; genomeLabel_BUSCO }

    Channel
     .fromPath(params.genomes)
     .into { genome_runMinimap2; genome_runAssemblathonStats; genome_BUSCO }
/*
   process splitTuple {
     input:
     set val(label), file(genomeFile) from genome_runMinimap2

     output:
     val label into genomeLabel_ch
     path genomeFile into genomeFile_ch

     script:
     """
     echo "process requires a script"
     """
   }
*/
// chunk the fastq file and create a channel for the chunks
   Channel
       .fromPath(params.reads)
       .splitFastq(by: params.chunkSize, file:true)
       .set { read_chunks }


    process runMinimap2 {

      container = "$medaka_container"

      input:
      //set val(label), file(genomeFile) from genome_runMinimap2
      val label from genomeLabel_runMinimap2.val
      path genomeFile from into genome_runMinimap2.val
      path readsChunk from read_chunks

      output:
      file("${label}.sam") into alignment_output
      //publishDir "${params.outdir}/assemblyStats", mode: 'copy', pattern: '*.assemblyStat'

      script:
      """
      minimap2 -ax map-ont ${genomeFile} ${readsChunk} > ${label}.sam
      """
    }

    alignment_output
        .collectFile(name: 'aligned_combined.txt', storeDir: params.outdir)
        .subscribe {
            println "Entries are saved to file: $it"
        }






/*
    process runRacon {

      container = "$racon_container"

      input:
      set val(label), file(genomeFile) from genome_runAssemblathonStats

      output:
      file("${label}.assemblathonStats")
      publishDir "${params.outdir}/assemblathonStats", mode: 'copy', pattern: '*.assemblathonStats'

      script:
      """
      new_Assemblathon.pl  ${genomeFile} > ${label}.assemblathonStats
      """
    }


  process runMedaka {

    container = "$medaka_container"


    input:
    set val(label), file(genomeFile) from genome_BUSCO
    file(config) from config_ch.val

    output:
    file("${label}/short_summary.specific.*.txt")
    publishDir "${params.outdir}/BUSCOResults/${label}/", mode: 'copy', pattern: "${label}/short_summary.specific.*.txt"
    file("${label}/*")
    publishDir "${params.outdir}/BUSCO"

    script:
    """
    busco \
    -o ${label} \
    -i ${genomeFile} \
    ${params.options} \
    -m ${params.mode} \
    -c ${params.threads} \
    -f
    """

  }

*/

    def isuGIFHeader() {
        // Log colors ANSI codes
        c_reset = params.monochrome_logs ? '' : "\033[0m";
        c_dim = params.monochrome_logs ? '' : "\033[2m";
        c_black = params.monochrome_logs ? '' : "\033[1;90m";
        c_green = params.monochrome_logs ? '' : "\033[1;92m";
        c_yellow = params.monochrome_logs ? '' : "\033[1;93m";
        c_blue = params.monochrome_logs ? '' : "\033[1;94m";
        c_purple = params.monochrome_logs ? '' : "\033[1;95m";
        c_cyan = params.monochrome_logs ? '' : "\033[1;96m";
        c_white = params.monochrome_logs ? '' : "\033[1;97m";
        c_red = params.monochrome_logs ? '' :  "\033[1;91m";

        return """    -${c_dim}--------------------------------------------------${c_reset}-
        ${c_white}                                ${c_red   }\\\\------${c_yellow}---//       ${c_reset}
        ${c_white}  ___  ___        _   ___  ___  ${c_red   }  \\\\---${c_yellow}--//        ${c_reset}
        ${c_white}   |  (___  |  | / _   |   |_   ${c_red   }    \\-${c_yellow}//         ${c_reset}
        ${c_white}  _|_  ___) |__| \\_/  _|_  |    ${c_red  }    ${c_yellow}//${c_red  } \\        ${c_reset}
        ${c_white}                                ${c_red   }  ${c_yellow}//---${c_red  }--\\\\       ${c_reset}
        ${c_white}                                ${c_red   }${c_yellow}//------${c_red  }---\\\\       ${c_reset}
        ${c_cyan}  isugifNF/nanoQCtrim  v${workflow.manifest.version}       ${c_reset}
        -${c_dim}--------------------------------------------------${c_reset}-
        """.stripIndent()
    }
