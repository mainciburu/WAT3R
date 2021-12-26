task wat3r {

    File bc_fastq
    File tcr_fastq

    Float memory
    Int disk_space
    Int num_threads

    command <<<
        wat3r -b ${bc_fastq} -t ${tcr_fastq} -p ${num_threads}
    >>>

    output {
        File umi_table = "wat3r/sample_igblast_db-pass.tsv"
        File metrics = "wat3r/wat3rMetrics.txt"
        File plot_filter_qscore = "wat3r/QC/QCplots_preFiltering.pdf"
        File plot_cluster_tcrs = "wat3r/QC/QCplot_clusters.pdf"
    }

    runtime {
        docker: "mainciburu/wat3r:1.0"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${num_threads}"
    }

    meta {
        author: "Peter van Galen"
        version: "1.0"
    }
}

workflow wat3r_workflow {
    call wat3r
}