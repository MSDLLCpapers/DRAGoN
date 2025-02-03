
## Usage

```
Typical pipeline command:

  nextflow run DRAGoN -profile <docker/singularity/.../institute> --IO.metadata metadata.xlsx --IO.fastqdir s3://bucket/path/to/fastq/folder/ --IO.outdir <OUTDIR>

--help                                [boolean, string] Show the help message for all top level parameters. When a parameter is given to `--help`, the full help message of that parameter will be printed.
--helpFull                            [boolean]         Show the help message for all non-hidden parameters.
--showHidden                          [boolean]         Show all hidden parameters in the help message. This needs to be used in combination with `--help` or `--helpFull`.

Input/Output options
  --IO.fastqdir                       [string]  Path containing fastq files. Can be an empty or nonexistent directory, in which case these will be generated.
  --IO.bcldir                         [string]  Path to raw base calls off the sequencer. Required if --IO.fastqdir is an empty or nonexistent directory.
  --IO.sampleSheet                    [string]  Path to sample sheet defining sample names and I5/I7 indices for use with bcl2fastq. Required if --IO.fastqdir is an empty or nonexistent directory.
  --IO.reference                      [string]  Name of a reference genome for building STAR indices. Soft required if --IO.star is an empty or nonexistent directory. Hard required if --IO.star is missing.
  --IO.star                           [string]  Path to STAR genome indices for the target species. Can be an empty or nonexistent directory, in which case these will be generated. If omitted, --IO.reference is required.
  --IO.gtf                            [string]  Path to gene annotation GTF file for building STAR indices. Required if --IO.star is an empty or nonexistent directory and --IO.reference is not provided.
  --IO.fasta                          [string]  Path to genome sequence fasta file for building STAR indices. Required if --IO.star is an empty or nonexistent directory and --IO.reference is not provided.
  --IO.contaminant                    [string]  A STAR (versionGenome == 2.7.4a) index of a contaminant genome reference (rRNA, ncRNA, adapter sequences, etc.). If supplied, removes reads that align to the contaminant genome from further consideration. Such reads will be marked with flag 0x800 (SUPPLEMENTARY).
  --IO.metadata                       [string]  Path to Excel spreadsheet defining the experiment metadata
  --IO.barcodes                       [string]  Path to tab-separated text file defining barcode sequences and well names. This file should have no header row. The first column contains the barcode sequence, and the second column contains the name of the well. All other columns are ignored. Required if --IO.metadata is omitted.
  --IO.sampleName                     [string]  Prefix of fastq filenames within --IO.fastqdir to process. For example, `--IO.sampleName Sample1` will select all files matching `path/to/fastq/Sample1_*.fastq.gz`. Specify multiple prefixes using a quoted list. Required if --IO.metadata is omitted.
  --IO.outdir                         [string]  Path to store the analysis outputs. Will be created if it does not already exist. Previous analyses in this directory will be overwritten.
  --IO.keep_bam                       [string]  If set, a BAM file is included in the output.
  --IO.phred                          [string]  Phred quality score offset  (accepted: 33, 64) [default: 33]
  --IO.no_prestage_index              [boolean] When running locally with locally-mounted STAR index, disable copying it to the workdir

Options for QC and well demultiplexing
  --Demux.dmuxn                       [integer] Split input fastq into chunks of this size (0: no chunking) [default: 0]
  --Demux.split                       [integer] Subdivide the subsequent tasks into at most this many parallel jobs. Reads will be grouped by barcode, and barcodes will be assorted into downstream tasks in order to distribute the workload as evenly as possible. Depending on the plate layout, this can still result in unbalanced workloads. If 0, each well gets its own downstream task. [default: 0]
  --Demux.umilen                      [integer] Length of the UMI portion of the R1 read. For an experiment with well barcodes of length B, the UMI begins at position B+1. [default: 10]
  --Demux.mismatch                    [integer] Maximum number of mismatches (substitutions) allowed in the barcode sequence. If 0, only reads whose barcodes exactly match the predefined barcode sequence for a given well will be assorted to that well. Increasing this figure lets you capture more reads per well at the cost of some added noise. Setting this too high will cause an error. [default: 1]
  --Demux.bcNcountMax                 [integer] Maximum number of base calls in the barcode sequence that are allowed to be N. Defaults to --Demux.mismatch.
  --Demux.minQuality                  [integer] Base call quality threshold for up-front UMI filtering and cDNA base trimming. [default: 20]
  --Demux.minUmiQualBases             [integer] Minimum number of bases in the UMI portion of R1 whose quality scores must meet or exceed --Demux.minQuality. [default: 6]
  --Demux.maxATcontent                [number]  Maximum fraction of bases in the cDNA (R2) read allowed to be A or T. Reads that exceed this are discarded. [default: 0.9]
  --Demux.maxGCcontent                [number]  Maximum fraction of bases in the cDNA (R2) read allowed to be G or C. Reads that exceed this are discarded. [default: 0.9]
  --Demux.minLength                   [integer] Minimum allowed length of the cDNA (R2) read before or after trimming. Reads that do not meet this are discarded. [default: 20]
  --Demux.homoPolymer                 [integer] Define a homopolymer to be a run of at least this many of the same base in the cDNA (R2) read. If one is found, the read is trimmed from the start of the homopolymer to the 3' end of the read. [default: 10]
  --Demux.maxHomoMmatch               [integer] Maximum number of bases in a suspected homopolymer which are allowed to vary from that base. If set to 1 or higher, there can be a significant performance impact. [default: 0]
  --Demux.adapter                     [string]  PCR/sequencing adapter to trim from the 3' end of R2 reads [default: CTGTCTCTTATA]
  --Demux.data_quality_errors         [string]  Governs how the pipeline should handle wells with low quality.
raise (default): Aborts the pipeline on error
warn: Emits a warning to stderr
ignore: Don't run heuristics  (accepted: ignore, warn, raise) [default: raise]
  --Demux.downsample_to               [string]  Determine the threshold number of QC-pass reads above which a well will be downsampled. If a decimal between 0 and 1, the threshold is that fraction of the total QC-pass reads in the experiment. If no value is provided, the threshold will be determined adaptively as the lower of 15% of the total number of reads in the experiment, or the expression `floor(pct75+1.75*IQR)`, where pct75 is the 75th percentile of reads assigned to wells and IQR is the difference between the 25th and 75th percentiles. Default: no downsampling.
  --Demux.downsample_threshold        [string]  Downsample wells to a maximum of N reads. If a decimal between 0 and 1, will downsample to that fraction regardless of original size. If not provided, will choose the threshold given by --Demux.downsample_threshold. Default: no downsampling.
  --Demux.downsample_seed             [integer] Set a seed for the downsampling logic, for the sake of reproducibility. Default: random seed
  --Demux.plate_layout_xoff           [integer]
  --Demux.plate_layout_yoff           [integer]

Options for genome alignment
  --Align.outFilterMultimapNmax       [integer] Discard reads with more than this many valid alignments to the genome [default: 20]
  --Align.outFilterMatchNmin          [integer] Discard alignments with fewer than this many bases matching the reference [default: 20]
  --Align.outFilterMultimapScoreRange [integer] Discard secondary alignments whose scores are more than this far away from the best alignment for that read [default: 1]
  --Align.alignIntronMin              [integer] Minimum number of intron bases to map the read [default: 1]
  --Align.alignIntronMax              [integer] Maximum number of intron bases to map the read [default: 1]
  --Align.clip5pNbases                [integer] Clip this many bases from the 5' end of reads [default: 0]
  --Align.clip3pNbases                [integer] Clip this many bases from the 3' end of reads [default: 0]
  --Align.outFilterMismatchNoverLmax  [number]  Discard alignments with more than this fraction of mismatched bases [default: 0.1]
  --Align.twopassMode                 [string]  Whether to run STAR in single-pass (default) or Basic two-pass mode  (accepted: None, Basic) [default: None]
  --Align.extraOptions                [string]  A quoted string of extra flags to pass to STAR. Use this for options not enumerated above.
  --Align.MEM                         [string]  Amount of memory to reserve per CPU when running alignment processes. Interpreted as a number of bytes, or you can add a memory suffix e.g. 8GB. (def. 20GB; recommend increasing it for larger genomes) [default: 64GB]

Options for feature annotation
  --Count.fracOverlap                 [number]  Minimum fraction of overlapping bases in a read that is required for read assignment. Value should be within range [0,1]. 0 by default. Number of overlapping bases is counted from both reads if paired end. Both this option and '--minOverlap' option need to be satisfied for read assignment. [default: 0]
  --Count.strandness                  [integer] Denote whether reads are unstranded (0, def.), forward-strand (1) or reverse-strand (2) specific. [default: 0]
  --Count.extraOptions                [string]  A quoted string of extra flags to pass to featureCounts.
  --Count.adjust_effective_lengths    [boolean] If set, when normalizing counts, the gene or transcript length will be adjusted by the average read (fragment) length for each sample as part of the calculation.
  --Count.MEM                         [string]  Amount of memory to reserve per CPU when running alignment processes. Interpreted as a number of bytes, or you can add a memory suffix e.g. 8GB. (def. 8GB) [default: 8GB]

Options for duplicate detection
  --Dedup.umi_mismatch                [integer] Maximum edit distance between two UMIs to consider them duplicates. [default: 2]
  --Dedup.umiDist                     [integer] Search for duplicates within a sliding genomic window of this many bases. [default: 10]
  --Dedup.ambiguous                   [string]  Comma-separated list of strategies for handling multimapping reads. Accepted values are:
                                                  Unique: Only count reads assigned to a single gene.
                                                  Uniform: Distribute each multimapping read's count uniformly across all genes they map to.
                                                  PropUnique: Attempt to distribute multimapping reads proportional to the number of reads mapping uniquely to each gene. Falls back on Uniform.
                                                  EM: Computes the maximum likelihood distribution of multimapping reads using an expectation-maximization algorithm.
                                                  Rescue: Equivalent to running a single iteration of EM. [default: Unique]
  --Dedup.em_iter_max                 [integer] DRAGoN only: If --Dedup.ambiguous contains EM, the maximum number of iterations of EM [default: 100]
  --Dedup.em_max_diff                 [number]  DRAGoN only: If --Dedup.ambiguous contains EM, the tolerance for convergence, measured by the largest absolute change in counts for any gene in the current well [default: 0.01]
  --Dedup.em_small_thresh             [number]  DRAGoN only: If --Dedup.ambiguous contains EM, clamp counts smaller than this to 0 [default: 0.01]
  --Dedup.debug                       [boolean]

uncategorized_options
  --private.versionGenome             [string]  [default: 2.7.4a]
  --private.overhang                  [integer] [default: 89]
  --awsbatch_queue                    [string]  awsbatch queue to submit jobs to
  --awsbatch_jobrole                  [string]  IAM role to use for awsbatch worker instances
  --errorStrategy                     [string]  If a worker hits an error after two attempts, the workflow should either
    continue to "retry",
    allow all running tasks to "finish" before exiting (default), or
    "terminate" the workflow immediately.  (accepted: finish, terminate, retry) [default: finish]
  --maxRetries                        [integer] If --errorStrategy retry, the maximum number of times the task will be retried before the workflow terminates. [default: 3]
  --MINMEM                            [string]  Minimum amount of memory to reserve for any given job. Applies globally. [default: 2GB]
  --MAXMEM                            [string]  Maximum amount of memory to reserve for any given job. Applies globally. [default: 4GB]
  --MAXCPU                            [integer] Maximum number of CPUs to reserve for any given job. Applies globally. [default: 2]
  --PLATE_LAYOUT                      [string]  For MergeDemuxJSON, a CSV file representing the plate well barcodes in a grid layout. [default: assets/layout_384.csv]
  --notifyEmails                      [string]  Comma-separated list of emails to notify on workflow completion.
```

<br><br>
