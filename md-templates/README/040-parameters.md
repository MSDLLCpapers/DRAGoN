
<a name="parameters"></a>

## Parameters description

#### Parameters
Users have the option of using a YAML file to pass parameters to this workflow. [A template is provided](params.example.yml) for your convenience and for a brief description of available parameters. You may copy this file to a new path, add the missing reuqired fields, edit optional parameters as appropriate, and submit it to nextflow using the option `-params-file path/to/params.yml`. Options are arranged in a categorical hierarchy, with categories `IO` for file paths, `Demux` for demultiplexing and QC options, `Align` for STAR options, `Count` for featureCounts options, and `Dedup` for deduplication and final counting options. These same parameters can be passed from the commandline i.e. `--Demux.umilen 12 --IO.star path/to/genome-idxs` would be equivalent to setting in your YAML file
```yaml
Demux:
  umilen: 12
IO:
  star: path/to/genome-idxs
```

<br><br>

#### Accessing resources on AWS

Depending on your runtime environment and S3/batch policies, you may need to specify credentials. The pipeline can discover AWS credentials using any combination of environment variables, (AWS_ACCESS_KEY_ID, AWS_SECRET_ACCESS_KEY, AWS_PROFILE, AWS_DEFAULT_PROFILE, AWS_REGION, AWS_DEFAULT_REGION), configuration profiles (~/.aws/config, ~/.aws/credentials), or your EC2 headnode configuration. For more info, [click here](https://www.nextflow.io/docs/latest/aws.html#aws-page).

<br><br>

<a name="downsampling"></a>

#### Downsampling and wells masking

If the pipeline detects wells with excessive read counts, it will by default raise an error. You can reduce these to warnings and, optionally, invoke downsampling logic. By default, downsampling is disabled, however you can control the behavior using two flags and an optional metadata column.

<a name="downsampling-threshold-flag"></a>

##### `--Demux.downsample_threshold`

This flag controls the number of QC-passing reads a well must have in order for downsampling to apply. The behavior is as such:

- If not provided, downsampling is disabled.
- If provided without a value, the threshold is determined by taking the lower of 15% of the total QC-pass reads in the experiment or an outlier threshold determined by the formula `Q3 + 1.75 * (Q3 - Q1)` where `Q1` and `Q3` are the first and third quartiles of the per-well QC-pass read counts, respectively.
- If provided with a fraction between 0 and 1, the threshold is set at that fraction of the total QC-pass reads in the experiment.
- If provided with a whole number greater than 1, the threshold is set at that number of reads.

<a name="downsampling-to-flag"></a>

##### `--Demux.downsample_to`

This flag controls the target number of reads to downsample a well to if it is flagged for downsampling per the above logic. The behavior is as such:

- If not provided, downsampling is disabled.
- If provided without a value, the target number of reads is set to equal the threshold.
- If provided with a fraction between 0 and 1, the target number of reads is set per well to that percent of its QC-pass reads.
- If provided with a whole number greater than 1, the target number of reads is set to that number.

> [!NOTE]
> If either of these two options is set but not the other, the other option will be behave as though it was provided without a value.

<a name="downsampling-mask-wells"></a>

##### Metadata column "MaskWells"

This is an optional column in the metadata sheet. The intention is for this to be a boolean value. Set this to true to discard that well from downstream consideration. The well will still be processed through FastQC and initial QC, but will be discarded at the point of downsampling.

<br><br>
