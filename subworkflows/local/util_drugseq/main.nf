// Copyright © 2025 Merck & Co., Inc., Rahway, NJ, USA and its affiliates. All rights reserved.
// This file is part of DRAGoN.
//
// This source code is licensed under the MIT License found in the
// LICENSE file in the root directory of this source tree.

// nf-core versions to yaml flow helpers

//
// Get commit hash
//
def getGitHash() {
    if (workflow.commitId) {
        return workflow.commitId
    }
    def sout = new StringBuffer()
    def serr = new StringBuffer()
    def proc = "git rev-parse HEAD".execute(null, file(workflow.projectDir).toFile())
    proc.waitForProcessOutput(sout, serr)
    return proc.exitValue() ? null : "${sout}".trim()
}

//
// Get current branch
//
def getGitRevision() {
    if (workflow.revision) {
        return workflow.revision
    }
    def sout = new StringBuffer()
    def serr = new StringBuffer()
    def proc = "git rev-parse --abbrev-ref HEAD".execute(null, file(workflow.projectDir).toFile())
    proc.waitForProcessOutput(sout, serr)
    return proc.exitValue() ? null : "${sout}".trim()
}

//
// Generate version string
//
def getWorkflowVersion() {
    def version_string = "" as String
    if (workflow.manifest.version) {
        def prefix_v = workflow.manifest.version[0] != 'v' ? 'v' : ''
        version_string += "${prefix_v}${workflow.manifest.version}"
    }

    def commitId = getGitHash()
    if (commitId) {
        def git_shortsha = commitId.substring(0, 7)
        version_string += "-g${git_shortsha}"
    }

    def revision = getGitRevision()
    if (revision) {
        version_string += "@${revision}"
    }

    return version_string
}

//
// Get software versions for pipeline
//
def processVersionsFromYAML(yaml_file) {
    def yaml = new org.yaml.snakeyaml.Yaml()
    def versions = yaml.load(yaml_file).collectEntries { k, v -> [k.tokenize(':')[-1], v] }
    return yaml.dumpAsMap(versions).trim()
}

//
// Get workflow version for pipeline
//
def workflowVersionToYAML() {
    return """
    Workflow:
        ${workflow.manifest.name}: ${getWorkflowVersion()}
        Nextflow: ${workflow.nextflow.version}
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
