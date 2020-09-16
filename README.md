# isugifNF/nanoPolish

```
----------------------------------------------------
                                    \\---------//       
      ___  ___        _   ___  ___    \\-----//        
       |  (___  |  | / _   |   |_       \-//         
      _|_  ___) |__| \_/  _|_  |        // \        
                                      //-----\\       
                                    //---------\\       
      isugifNF/nanoPolish  v1.0.0       
    ----------------------------------------------------
```

[Genome Informatics Facility](https://gif.biotech.iastate.edu/) | [![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A519.10.0-brightgreen.svg)](https://www.nextflow.io/)

---

### Introduction

Polish a nanopore assembly using Racon and Medaka

**isugifNF/blast** is a [nextflow pipeline](https://www.nextflow.io/)


```
git clone git@github.com:nanoporetech/medaka.git
```

### Installation and running on Ceres HPCC

Nextflow is already installed on Ceres HPCC. Therefore, running **isugifNF/blast** involves (1) allocating a debug node `salloc -N 1 -p debug -t 01:00:00`, (2) loading nextflow `module load nextflow`, and (3) running the pipeline `nextflow run isugifNF/blast`. The `--help` flag prints out the usage statement.

```
salloc -N 1 -p debug -t 01:00:00
module load nextflow
nextflow run isugifNF/blast --help

```


```
nextflow run isugifNF/nanoPolish --genome tail.fasta --reads test.fastq --chunkSize 25000 --model "medaka/medaka/data/r941_min_high_g303_model.hdf5" -profile singularity,condo -resume
```


<details><summary>Usage Statement</summary>

<pre>

```
Usage:
The typical command for running the pipeline are as follows:

nextflow run isugifNF/nanoPolish --genome tail.fasta --reads test.fastq --chunkSize 25000 --model "medaka/medaka/data/r941_min_high_g303_model.hdf5" -profile singularity,condo


Mandatory arguments:

--genome                      genome assembly fasta file to run stats on. (./data/*.fasta)
-profile singularity           as of now, this workflow only works using singularity and requires this profile [be sure singularity is in your path]
--reads                       Raw nanopore reads to use in polish
Optional arguments:
--outdir                       Output directory to place final output
--threads                      Number of CPUs to use during the NanoPlot job [16]
--queueSize                    Maximum number of jobs to be queued [18]
--model                        Medaka hdf5 model used during polishing.
--chunkSize                    Number of fasta records to use when splitting the input fastq raw reads
--help                         This usage statement.

```

</pre>
</details>



## Common Errors

<details><summary>Unable to open file (file signature not found)</summary>

<pre>


The hd5 files that are used as models are stored in gits large file storage `lfs`, if you do a git clone those files will just be text pointers unless you have git lfs installed.  If you run into a `Unable to open file (file signature not found)` error it is because the hd5 file is not an actual model file but a pointer to git's lfs storage.

In the medaka repo execute the following commands.

```
git lfs install
git lfs pull
```

</pre>
</details>
