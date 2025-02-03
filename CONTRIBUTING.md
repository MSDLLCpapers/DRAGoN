This is a guide to contributing to this project.

## Compiling and testing

This workflow contains custom C++ sources which must be compiled. No extra actions are necessary if using the Docker container with the awsbatch executor (see below).

To compile, you should first create and activate the Conda environment defining the required packages, then use cmake to build the project.
```bash
cd /where/you/cloned/this/repository
conda env create --force --file conf/conda/drugseq-env.yml
conda activate drugseq-env
cmake . -B build -DCMAKE_MAKE_PROGRAM=qmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$(realpath .) -DCMAKE_PREFIX_PATH=${CONDA_PREFIX} -DCMAKE_MODULE_PATH=${CONDA_PREFIX}/unpacked_source/cmake
cmake --build build --target test install/strip --parallel
```

<br>

To test Python components, run the following inside this directory:
```bash
python -m unittest discover -s tests/python-unittest
```

To test the workflow, you will also need a compatible version of nextflow and a compute environment that supports the processes with high resource demand. You will need a minimum of 10 logical CPUs and 80 GB of RAM at your disposal. If you have less memory available, you can pass a smaller human-readable value via `--MAXMEM`, however this is potentially unstable.

```bash
nextflow run -profile [...],test --IO.outdir ./outs
```

## Making changes

This project is written primarily in three languages: Python, C++, and Nextflow (Groovy). The general organization principle for this project is as such:

- Workflow organization is done in Nextflow. `main.nf`, `workflows/DRAGoN.nf`, and `workflows/STARsolo.nf` are the main entrypoints. Subworkflows are imported from `components/`. See the Nextflow [documentation](https://www.nextflow.io/docs/latest/index.html) for language details.
- Heavy-lifting processes are performed using either third-party tools or custom C++ programs. For the former, the space of tools that are available to workers are defined in `conf/conda/drugseq-env.yml`. This also defines the Python packages and C/C++ libraries available to the C++ sources in `src/` and the Python scripts in `bin/`.
  - Every custom C++ program has an equivalent Python fallback script.
  - C++ binaries are also installed into `bin/` for local execution.
- Smaller processes such as merging intermediates are handled using Python scripts in `bin/`.

## CI/CD

This project uses GitHub for fine-grained version control. The main branch is `master`, and this has branch protections preventing direct push. Thus any change made to the pipeline should be made in a separate branch, and work logged in a JIRA ticket. The branch name convention is `feature/<ticket>` or `bugfix/<ticket>`.

All contributors should verify that their changes satisfy the following checks:

1. If C++ sources are modified, they should compile, and all Ctest tests should pass.
1. If python scripts are modified, all python tests should pass.
1. Running this workflow locally (using the `conda` profile and locally-built binaries) should pass.
1. All output files should be checked for proper formatting.
1. The Build Docker Image workflow on GitHub Actions should pass.
1. During the PR review process:
  1. The Run Unit Tests workflow should pass.
  1. version.txt should be updated and committed.
  1. nf-test snapshots should be updated (`nf-test test . --updateSnapshots`) and committed.
1. Running this workflow on the CTC cluster using the `singularity` profile should pass. You may need to manually update the container image into your cache directory (specified by `NXF_SINGULARITY_CACHEDIR`).
1. Running this workflow from a DRUG-seq head node should pass.
