workflow wat3r_workflow {

    String? sample_name
    Boolean annotation

    call wat3r

    if( annotation ) {
        call downstream_anno { input: umi_alignments=wat3r.umi_alignments, cluster_metrics=wat3r.cluster_metrics, filter_metrics=wat3r.filter_metrics, stats=wat3r.stats, sample_name=sample_name }
    }

    if( ! annotation ) {
        call downstream { input: umi_alignments=wat3r.umi_alignments, cluster_metrics=wat3r.cluster_metrics, filter_metrics=wat3r.filter_metrics, stats=wat3r.stats, sample_name=sample_name }
    }
}

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
        File umi_alignments = "wat3r/sample_igblast_db-pass.tsv"
        File cluster_metrics = "wat3r/QC/BC_UMI_cluster_metrics.txt"
        File filter_metrics = "wat3r/wat3rMetrics.txt"
        File stats = "wat3r/stats.log"
        File plot_filter_qscore = "wat3r/QC/QCplots_preFiltering.pdf"
        File plot_cluster_tcrs = "wat3r/QC/QCplot_clusters.pdf"
    }

    runtime {
        docker: "mainciburu/wat3r:1.1"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${num_threads}"
    }

    meta {
            author: "Peter van Galen"
        }
}

task downstream_anno {

    String? sample_name
    File? umi_alignments
    File? cluster_metrics
    File? filter_metrics
    File? stats
    File? cluster_annotation

    command <<<
        downstream -u ${umi_alignments} -c ${cluster_metrics} -f ${filter_metrics} -s ${stats} -n ${sample_name} -a ${cluster_annotation}
    >>>

    output {
        File barcode_table_filter = "downstream/${sample_name}_barcode_results.csv"
        File umi_alignments_filter = "downstream/${sample_name}_barcode_UMI_results.csv"
        File plot1 = "downstream/plots/db_histograms.pdf"
        File plot2 = "downstream/plots/ReadPercentage_FilteringSteps.pdf"
        File plot3 = "downstream/plots/scRNAseq_TCRrecovery_proportions.pdf"
        File plot4 = "downstream/plots/CDR3_UMIcount_distribution.pdf"
        File plot5 = "downstream/plots/valid_reads.pdf"
        File plot6 = "downstream/plots/CDR3_clones_heatmap.pdf"
        File plot7 = "downstream/plots/TRB_TRA_correspondence.pdf"
        File plot8 = "downstream/plots/TRA_TRB_clone_size.pdf"
        File plot9 = "downstream/plots/TRB_clone_size_celltype.pdf"
        File plot10 = "downstream/plots/trb_top_clones.pdf"
        File plot11 = "downstream/plots/trb_top_clones_norm.pdf"
        File plot12 = "downstream/plots/TRB_distance_heatmap.pdf"        
    }

    runtime {
        docker: "mainciburu/wat3r:1.1"
    }

    meta {
        author: "Peter van Galen"
    }
}

task downstream {

    String? sample_name
    File? umi_alignments
    File? cluster_metrics
    File? filter_metrics
    File? stats

    command <<<
        downstream -u ${umi_alignments} -c ${cluster_metrics} -f ${filter_metrics} -s ${stats} -n ${sample_name}
    >>>

    output {
        File barcode_table_filter = "downstream/${sample_name}_barcode_results.csv"
        File umi_alignments_filter = "downstream/${sample_name}_barcode_UMI_results.csv"
        File plot1 = "downstream/plots/db_histograms.pdf"
        File plot2 = "downstream/plots/ReadPercentage_FilteringSteps.pdf"
    }

    runtime {
        docker: "mainciburu/wat3r:1.1"
    }

    meta {
        author: "Peter van Galen"
    }
}
