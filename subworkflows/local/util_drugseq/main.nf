
// nf-core versions to yaml flow helpers

// Use is permitted under the condition that the pipeline is released with an MIT License.

//
// Generate workflow version string
//
def getGitHash() {
    if (workflow.commitId) {
        return workflow.commitId
    }
    def sout = new StringBuffer()
    def serr = new StringBuffer()
    def proc = "git -C ${workflow.projectDir} rev-parse --short=7 HEAD".execute()
    proc.waitForProcessOutput(sout, serr)
    return proc.exitValue() ? null : "${sout}"
}

def getWorkflowVersion() {
    def version_string = ""
    if (workflow.manifest.version) {
        def prefix_v = workflow.manifest.version[0] != 'v' ? 'v' : ''
        version_string += "${prefix_v}${workflow.manifest.version}"
    }

    def commitId = getGitHash()
    if (commitId) {
        def git_shortsha = commitId.substring(0, 7)
        version_string += "+g${git_shortsha}"
    }

    if (workflow.revision) {
        version_string += "@${workflow.revision}"
    }

    return version_string
}

//
// Get software versions for pipeline
//
def processVersionsFromYAML(yaml_file) {
    def yaml = new org.yaml.snakeyaml.Yaml()
    def versions = yaml.load(yaml_file).collectEntries { k, v -> [ k.tokenize(':')[-1], v ] }
    return yaml.dumpAsMap(versions).trim()
}

//
// Get workflow version for pipeline
//
def workflowVersionToYAML() {
    return """
    Workflow:
        $workflow.manifest.name: ${getWorkflowVersion()}
        Nextflow: $workflow.nextflow.version
    """.stripIndent().trim()
}

//
// Get channel of software versions used in pipeline in YAML format
//
def softwareVersionsToYAML(ch_versions) {
    return ch_versions
                .unique()
                .map { it -> processVersionsFromYAML(it) }
                .unique()
                .mix(Channel.of(workflowVersionToYAML()))
}

def usingProfile(profile) {
    workflow.profile.tokenize(',').contains(profile)
}
