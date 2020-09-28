#! /usr/bin/env nextflow

/*************************************
 nextflow nanoQCtrim
 *************************************/

racon_container = 'quay.io/biocontainers/racon:1.4.13--he513fc3_0'
medaka_container = 'quay.io/biocontainers/medaka:1.0.3--py36hbecb4b7_1'
samtools19_container = 'quay.io/biocontainers/samtools:1.9--h10a08f8_12'

//params
params.genome = "assembly_racon.fasta"
params.inputBam = "out_dir/calls_to_draft.bam"
params.inputBai = "out_dir/calls_to_draft.bam.bai"
params.outdir = "medakaOut"
params.model = "medaka/medaka/data/r941_min_high_g303_model.hdf5"
//process
process {
  executor = 'slurm'
  clusterOptions =  '-N 1 -n 16 -t 48:00:00'
}



executor {
  queueSize = params.queueSize
  submitRateLimit = '20 sec'
}

  docker {
    docker.enabled = true
  }

  singularity {
    singularity.enabled = true
    singularity.autoMounts = true
    // note: could probably go here... but I haven't tested this method
    // process {
    //   withLabel:'BUSCO' { container = 'list container name' } // <- if you use labels
    //   withName:runBUSCO { container = 'list container name' } // <- if you want to just use process name
    //}
    // note: otherwise I know putting above three lines into separate configs/singularity.config file works
    // includeConfig './configs/singularity.config'
  }
}

/* Not sure why this is necessary but nf-core/rnaseq had this line */
docker.runOptions = '-u \$(id -u):\$(id -g)'

// Capture exit codes from upstream processes when piping
process.shell = ['/bin/bash', '-euo', 'pipefail']

timeline {
  enabled = true
  file = "$params.outdir/timeline.html"
}

report {
  enabled = true
  file = "$params.outdir/report.html"
}





// Channels for the genome and its label
Channel
.fromPath(params.genome)
.map { file -> file.simpleName}
.into { genomeLabel_runMinimap2; genomeLabel_runRacon; genomeLabel_BUSCO }

Channel
.fromPath(params.genome)
.into { genome_runMinimap2; genome_runRacon; genome_getScaffolds }


//Channel for medaka consensus model





Channel
.fromPath(params.model)
.set { model_medaka }

Channel
.fromPath(params.inputBam)
.set { medakaAlign_ch }

Channel
.fromPath(params.inputBai)
.set { medakaAlignBai_ch }



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
   path inputBai from medakaAlignBai_ch.val
   val region from regions_ch.splitText()
   path modelIn from model_medaka.val

   output:
   path("out*.hdf") into medakaConsensus_ch

   script:
   """
   medaka consensus ${inputAlign} out_${region.trim()}.hdf \
       --model ${modelIn} --batch 200 --threads 8 \
       --region ${region.trim()}
   """

   }
   //--region contig1 contig2 contig3 contig4
   //r941_min_high_g303

   // Medaka Step 3 collate results
   process medakaStich {
   container = "$medaka_container"

   input:
   path hdf from medakaConsensus_ch.collect()

   output:
   path("medaka_polished.assembly.fasta")
   publishDir "${params.outdir}", mode: 'copy', pattern: "medaka_polished.assembly.fasta"

   script:
   """
   medaka stitch ${hdf} medaka_polished.assembly.fasta
   """

   }
