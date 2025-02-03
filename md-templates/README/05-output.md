
<a name="output"></a>

## Output
All pipeline outputs will be placed in a new directory or S3 prefix defined by `--IO.outdir`. For each plate in the metadata excel, a separate subdirectory is created for each sample, and the name of this directory is taken from the bcl2fastq samplesheet if provided or from the names of the Excel tabs. If passing sample info via `--IO.barcodes` and `--IO.sampleName`, the latter argument is used for the name of the sole subdirectory of `--IO.outdir`.

<a name="starsolo"></a>

### STARsolo Output
The output of `STARsolo.nf` follows.
* `counts.tsv` - Tab-separated file containing counts of genes in wells. Rows are features (genes) and columns are barcodes (wells).
* `tpm.tsv`, `fpkm.tsv` - Tab-separated files containing normalized counts of genes in wells, with the same shape as `counts.tsv`.
* `demux.json` - JSON summary of demultiplexing stats with the following keys:
  * `read_count`: The total number of reads ingested.
  * `well_counts`: A list of per-barcode read counts. The order of this list is the same as the user-supplied barcode whitelist. An extra element at the end counts the number of reads discarded due to no matching barcode.
  * `written_counts`: A list of read counts per intermediate task. The length of this list is determined by `--Demux.split` plus an extra element for the noBC output.
  * `mismatch_counts`: A list of mappings from well name to the number of reads whose barcode was matched to that well with a given number of substitutions. The first element of this list counts the reads that were matched perfectly, the second those with one mismatch, and so on up to `nmismatch`.
  * `fail counts`: A mapping from QC fail reason (or PASS) to mappings from well name to the number of reads in that well passing or failing QC for the stated reason.
  * `unmatched_counts`: Observed frequencies of barcodes that failed to match one of the given well barcodes.
  * `plate_map`: 2d nested array representing the number of reads that are mapped to each well whether or not it has been assigned a sample ID. If the workflow complains about too many reads being discarded, you can use this as a diagnostic resource.
* `Aligned.sortedByCoord.out.bam` - The sorted BAM file, merged from all the STAR runs, sorted by barcode index (read group) and genomic coordinate.
* `Aligned.sortedByCoord.out.flagstat.json` - The result of running `samtools flagstat -O json Aligned.sortedByCoord.out.bam`. See the [`samtools` documentation](https://www.htslib.org/doc/samtools-flagstat.html) for details.
* `STARsoloReport.txt` - Click [here](#starsoloreport) for details on this file.
* `Solo.out` - A consolidation of the STARsolo reports from across the mapping tasks. Output is similar to that of a single run of STARsolo, except the three stats/summary tables (`Solo.out/Barcodes.stats`, `Solo.out/Gene/Features.stats`, and `Solo.out/Gene/Summary.csv`) are tables with each column representing a single STARsolo task. Column order is arbitrary. Additionally, contains subdirectory `Solo.out/Gene/raw` which includes three or more files:
  * `features.tsv`: A tab-separated file with 3 columns describing each annotated gene. Column 1 is the ENSEMBL gene ID (ENSG##########). Column 2 is the HGNC gene symbol (_PARP2L_, _ACTB_, etc.). Column 3 is the feature type (by default, this is "Gene Expression").
  * `barcodes.tsv`: A text file with one barcode sequence per line, copied from the first column of `--IO.barcodes`.
  * `matrix.mtx`: A MatrixMarket sparse matrix file where each entry is the number of uniquely-mapped, uniquely-assigned UMIs per gene per well.
  * `UniqueAndMult-*.mtx`: Depending on the argument to `--Dedup.ambiguous`, this file is a MatrixMarket sparse matrix where each entry is the number of UMIs per gene per well after correcting for multimapping sites using the user-selected strategy. The previously-mentioned `counts.tsv` will then be a dense version of this.
* `DRUGseq_Merged.bam`: This BAM file contains all the reads in the experiment, sorted by barcode. The tags `NH HI nM AS CR UR CY UY CB UB GX GN sS sQ sM NM MD jM jI` areq as described in the STAR manual, section 4.2.2. The following tags are custom to this workflow:
  * `BM`: Number of mismatches (substitutions) between `CR` and `CB`
  * `QC`: The reason why the read was discarded in the initial filtering step
  * `BI`: Barcode index or "noBC", used for properly splitting the job after the STAR run (DRAGoN workflow only)

Note: This workflow generates a large number of intermediate files. For a 1.7B read pair experiment with 384 barcodes and 1M read chunks for demultiplexing, a total of over 650k SAM files will be generated. You should periodically clean up temporary files from old nextflow runs using `nextflow clean -f`.

<br><br>

<a name="starsoloreport"></a>

#### STARsoloReport.txt
This is a TSV file with per-wellread counts summarizing the demultiplexing output. Column definition is as such:
* **TOTAL_IN** - Total number of reads assigned to this well.
* **PCT_LOSS_QC** - Percent of reads lost to initial QC failure.
* **PCT_LOSS_DOWNSAMPLING** - Percent of reads lost to downsampling. If no downsampling occurs, this will be all 0s.
* **PCT_LOSS_UNMAPPED** - Percent of reads surviving initial QC that fail map to the supplied reference genome.
* **PCT_LOSS_UMI** - Percent of reads in genes that are removed as duplicates.
* **FINAL_UMI** - Final number of UMIs, including multimappers.
###### Initial QC
* **QC_PASS** - Number of reads passing all QC filter criteria.
* **QC_FAIL_UMI_QUALITY** - Number of reads for which too few bases in the UMI sequence pass the configured quality threshold (default: &lt;6 bases with quality score &gt;= 20).
* **QC_FAIL_AT_CONTENT** - Number of reads for which the cDNA sequence had too many bases that were A or T (default: &gt;90% of base calls).
* **QC_FAIL_GC_CONTENT** - Number of reads for which the cDNA sequence had too many bases that were G or C (default: &gt;90% of base calls).
* **QC_FAIL_ADAPTER_TRIMMING** - Number of reads for which the cDNA sequence became too short (default:&lt; 20 bases) after trimming sequencing adapters (single sequence specified by --Demux.adapter).
* **QC_FAIL_QUALITY_TRIMMING** - Number of reads for which the cDNA sequence became too short (default: &lt;20 bases) after trimming homopolymer artifacts and low-quality bases.
* **QC_FAIL_LENGTH** - Number of reads for which the cDNA sequence was too short before trimming (default: &lt;20 bases).
###### Barcode matching
* **MATCH_EXACT** - Number of reads (both QC-PASS and QC-FAIL) whose barcode sequence exactly matches the given well barcode.
* **MATCH_1MM, MATCH_2MM, etc.** - Number of reads (both QC-PASS and QC-FAIL) whose barcode sequence matches the given well barcode with 1, 2, etc. base substitutions.
  * There will be as many such columns as the value of `--Demux.mismatch`.
###### Downsampling
* **DOWNSAMPLE_QC_PASS_READS** - Number of reads going into downsampling, should be the same as QC_PASS.
* **DOWNSAMPLE_LIMIT** - Threshold for downsampling, should be the same for the whole plate.
* **DOWNSAMPLE_FRACTION** - Downsample-to fraction. If no downsampling happens, this will be 1.0. If the well is discarded, it will be 0.
###### Mapping and annotation
* **UNMAP_COUNT** - Number of reads that were discarded because they failed to map to the genome reference.
* **SOLO_UNIQUE** - Number of uniquely-assigned UMIs.
* **SOLO_UNIFORM, SOLO_PROPUNIQUE, SOLO_EM, SOLO_RESCUE** - (optional) Number of unique and multimapped UMIs assigned according to the titular strategy.
<br><br>

<a name="dragon"></a>

### DRAGoN
The `DRAGoN.nf` output includes the same `counts.tsv` matrix, `flagstat.json` and `demux.json` as from STARsolo, plus additional tables reporting stats at each intermediate step.
* `DRAGoNreport.txt`: A tab-separated file reporting QC statistics from across the execution. Concatenates portions of the demultiplexing and deduplication reports.
* Directory `DRAGoN.out` mimics the sparse matrix structure of `Solo.out/Gene/raw`.
* `featureCounts.summary` is a tab-separated file with summary statistics for the assignment of reads to genes. The first column names the statistic, and the second column gives its value. This is aggregated across the entire experiment. See the featureCounts documentation for a description of this file.
* `dedup_stats.tsv` is a tab-separated file with summary statistics for post-hoc filtering and UMI counts. This is entirely redundant to [DRAGoNreport.txt](#dragon-report) and will be removed in a future update.
* `DRUGseq_Merged.bam`: Same as from STARsolo, but the following additional custom tags are set:
  * `MX`: Optional: The coordinates of all alternate alignments in `chr:start-end` format, separated by semicolons
  * `XS`: featureCounts assignment result ("Assigned" or "Unassigned_{reason}") for all alignments (primary first, then secondaries), separated by semicolons.
  * `XT`: Optional: featureCounts gene assignment for all alignments (primary first, then secondaries), separated by semicolons.
  * `XN`: Optional: Number of genes to which this read was assigned, taking into account all alignments. This may not equal the length of tag XT i.e. in the case where multiple alignments of the same read overlap the same gene. If missing, the read did not overlap any genes.
  * `DN`: Optional: If the read is marked as a duplicate, the name of the original UMI.

<br><br>

<a name="dragon-report"></a>

#### DRAGoNreport.txt
This is a tab-separated file reporting QC and quantification summary statistics over the DRAGoN execution. Each row is a sample well, with an extra row at the end describing reads not assigning to any well. Columns are the same as in [STARsoloReport.txt](#starsoloreport) plus the following:

###### Index
* **barcode** - Name of the well. "noBC" refers to reads not assigning to any well.
###### Filtering summary
* **TOTAL_IN** - Total number of reads assigned to this well.
* **PCT_LOSS_QC** - Percent of reads lost to initial QC failure.
* **PCT_LOSS_DOWNSAMPLING** - Percent of reads lost to downsampling. If no downsampling occurs, this will be all 0s.
* **PCT_LOSS_UNMAPPED** - Percent of reads surviving initial QC that fail map to the supplied reference genome.
* **PCT_LOSS_UNASSIGNED** - Percent of reads that map to the genome that do not overlap any annotated gene.
* **PCT_LOSS_UMI** - Percent of reads in genes that are removed as duplicates.
* **FINAL_UMI** - Final number of UMIs, including multimappers.
###### Initial QC
* **QC_PASS** - Number of reads passing all QC filter criteria.
* **QC_FAIL_UMI_QUALITY** - Number of reads for which too few bases in the UMI sequence pass the configured quality threshold (default: &lt;6 bases with quality score &gt;= 20).
* **QC_FAIL_AT_CONTENT** - Number of reads for which the cDNA sequence had too many bases that were A or T (default: &gt;90% of base calls).
* **QC_FAIL_GC_CONTENT** - Number of reads for which the cDNA sequence had too many bases that were G or C (default: &gt;90% of base calls).
* **QC_FAIL_ADAPTER_TRIMMING** - Number of reads for which the cDNA sequence became too short (default:&lt; 20 bases) after trimming sequencing adapters (single sequence specified by --Demux.adapter).
* **QC_FAIL_QUALITY_TRIMMING** - Number of reads for which the cDNA sequence became too short (default: &lt;20 bases) after trimming homopolymer artifacts and low-quality bases.
* **QC_FAIL_LENGTH** - Number of reads for which the cDNA sequence was too short before trimming (default: &lt;20 bases).
###### Barcode matching
* **MATCH_EXACT** - Number of reads (both QC-PASS and QC-FAIL) whose barcode sequence exactly matches the given well barcode.
* **MATCH_1MM, MATCH_2MM, etc.** - Number of reads (both QC-PASS and QC-FAIL) whose barcode sequence matches the given well barcode with 1, 2, etc. base substitutions.
  * There will be as many such columns as the value of `--Demux.mismatch`.
###### Downsampling
* **DOWNSAMPLE_QC_PASS_READS** - Number of reads going into downsampling, should be the same as QC_PASS.
* **DOWNSAMPLE_LIMIT** - Threshold for downsampling, should be the same for the whole plate.
* **DOWNSAMPLE_FRACTION** - Downsample-to fraction. If no downsampling happens, this will be 1.0. If the well is discarded, it will be 0.
###### Mapping and annotation
* **UNMAP_COUNT** - Number of reads that were discarded because they failed to map to the genome reference.
* **DEDUP_TOTAL** - Total number of reads assigning to this well.
* **DEDUP_DISCARD_QCFAIL** - Number of reads that were discarded due to a QC fail status.
* **DEDUP_DISCARD_UNMAPPED** - Number of reads that were discarded because they do not map to any genomic site.
* **DEDUP_DISCARD_NOFEATURE** - Number of reads that were discarded because they do not overlap any gene feature.
* **DEDUP_PASSED_QC** - Final number of reads passing all QC, alignment, and annotation filters, including multimappers.
* **DEDUP_N_MULTIMAPPED** - Number of reads mapping to multiple genomic sites.
###### Deduplication and disambiguation
* **DEDUP_UMI_COUNTS** - Final number of UMIs, including multimappers.
* **DEDUP_UNIQUE_COUNTS** - Final number of QC-pass reads mapping unambiguously to a single site/gene.
* **DEDUP_UNIQUE_UMIS** = Final number of uniquely-mapping UMIs.

<br><br>
