#! /usr/bin/env nextflow

/*************************************
 nextflow nanoQCtrim
 *************************************/

racon_container = 'quay.io/biocontainers/racon:1.4.13--he513fc3_0'
medaka_container = 'quay.io/biocontainers/medaka:1.0.3--py36hbecb4b7_1'

//nextflow run isugifNF/nanoPolish --genomes tail.fasta --reads test.fastq -profile singularity,condo -resume


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

// Channels for the genome and its label
   Channel
    .fromPath(params.genomes)
    .map { file -> file.simpleName}
    .into { genomeLabel_runMinimap2; genomeLabel_runRacon; genomeLabel_BUSCO }

    Channel
     .fromPath(params.genomes)
     .into { genome_runMinimap2; genome_runRacon; genome_getScaffolds }

//Channels for reads and chunks of reads

  Channel
      .fromPath(params.reads)
      .into { read_file; read_file2 }

   Channel
       .fromPath(params.reads)
       .splitFastq(by: params.chunkSize, file:true)
       .set { read_chunks }


process runMinimap2 {

  container = "$medaka_container"

  input:
  //set val(label), file(genomeFile) from genome_runMinimap2
  val label from genomeLabel_runMinimap2.val
  path genomeFile from genome_runMinimap2.val
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
    .collectFile(name: 'aligned_combined.sam', storeDir: params.outdir)
    .set { overlaps_ch }



process runRacon {

  container = "$racon_container"

  input:
  path reads from read_file.val
  path genomeFile from genome_runRacon.val
  path overlaps from overlaps_ch
  //file overlaps from alignment_output.collectFile(name: 'aligned_combined.txt')
  //path overlaps from Channel.fromPath("${params.outdir}/aligned_combined.txt")
  val label from genomeLabel_runRacon.val

  output:
  file("${label}_racon.fasta") into raconGenome_ch
  publishDir "${params.outdir}", mode: 'copy', pattern: "${label}_racon.fasta"

  script:
  """
  racon -m 8 -x -6 -g -8 -w 500 -t ${params.threads} ${reads} ${overlaps} ${genomeFile} > ${label}_racon.fasta
  """
}

// https://nanoporetech.github.io/medaka/installation.html


//Run Medaka in 3 steps

// Step 1 Align reads to assembly
process medakaAlign {
container = "$medaka_container"

input:
path raconGenome from raconGenome_ch
path reads from read_file2.val

output:
file("calls_to_draft.bam") into medakaAlign_ch

script:
"""

mini_align -i ${reads} -r ${raconGenome} -m \
    -p calls_to_draft \
    -t ${params.threads}
"""

}

// Medaka Step 2 run consensus on each scaffold

////Get scaffold list

process getScaffolds {

input:
path genomeFile from genome_getScaffolds.val

output:
//publishDir "${params.outdir}"
//file("nfhead.txt")
stdout regions_ch

script:
"""
grep ">" ${genomeFile} | perl -pe 's/>//g'
"""

}


process medakaConsensus {
container = "$medaka_container"

input:
path inputAlign from medakaAlign_ch.val
val region from regions_ch.splitText()

output:
file("out.hdf") into medakaConsensus_ch
script:
"""
medaka consensus ${inputAlign} out.hdf \
    --model r941_min_high_g3210 --batch 200 --threads 8 \
    --region ${region.trim()}
"""

}
//--region contig1 contig2 contig3 contig4


// Medaka Step 3 collate results
process medakaStich {
container = "$medaka_container"

input:
path hdf from medakaConsensus_ch.collect()

output:

script:
"""
medaka stitch ${hdf} polished.assembly.fasta
"""

}


//r941_min_high_g3210


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
