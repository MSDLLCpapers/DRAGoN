
<a name="quick"></a>

## Quick start

### If your data are on S3...

You will need to explicitly tell nextflow which credentials to use. The recommendation is to use the following procedure:

1. Create `~/.aws/credentials` with mode 600 and define a profile there with your key pair. Example:

```
[my_aws_profile]
aws_access_key_id = YourAwsAccessKey
aws_secret_access_key = YourSecretAccessKey
```

2. Create `~/.nextflow/config` and add the following lines:

```
aws.profile = my_aws_profile
env.AWS_PROFILE = aws.profile
```

Replace `my_aws_profile` above with the actual name of your profile as defined in step 1.

3. Update `~/.aws/credentials` every time your access keys rotate.

### Generic Linux environment

First, install git and whichever engine you wish to use to run processes (anaconda, singularity, docker). Then install nextflow 23.10.1 (or later) from https://get.nextflow.io. Finally submit your job using an appropriate `nextflow run` command for your system. If you have your own local cluster, you may need to configure nextflow to use it yourself. You will also need to specify full paths to reference datasets at runtime.

If you want to run this workflow using the conda profile, you will need to compile the binary tools yourself. See [CONTRIBUTING.md](CONTRIBUTING.md) for instructions.
