## Configuration profiles

This pipeline supports singularity, docker, and conda executors via the `-profile` argument. Additional executors can be enabled by supplying a custom config.

If no profile is provided:
- You must have a minimum of 80 GB RAM available.
- The `conda` executable must exist in the runtime environment.
- Nextflow will submit jobs to the default (`local`) executor. The number of jobs that can be executed in parallel will therefore be limited by the resources available on the current host.

New in 1.1: You can now select between the DRAGoN and STARsolo pipelines via the `main.nf` script. By default, DRAGoN will be selected. Use the `--pipeline` switch to select the STARsolo pipeline if so desired.

With that configured, the pipeline can be launched as such:
```bash
nextflow run -w /SFS/scratch/$USER/.nextflow -profile ctc_hpc_sge MSDLLCPapers/DRAGoN [--PARAM VALUE ...|-params-file PARAMS.yml]
```

<br><br>
