rule pyDESeq2_interaction:
    input:
        in_adata="output/preprocessing/objects/QC/adata.h5ad"
    output:
        out_deg_interaction="output/DEA/objects/pyDESeq2_interaction/deg.tsv",
        out_deg_ifng_pos="output/DEA/objects/pyDESeq2_ifng_pos/deg.tsv"
    params:
        out_plots_interaction="output/DEA/plots/pyDESeq2_interaction",
        out_plots_ifng_pos="output/DEA/plots/pyDESeq2_ifng_pos"
    resources:
        mem_mb=32000
    conda:
        "../envs/pydeseq2_interaction.yaml"
    script:
        "../scripts/DEA/pyDESeq2_interaction.py"


rule calibrate_pvalues_interaction:
    input:
        in_deg="output/DEA/objects/pyDESeq2_interaction/deg.tsv"
    output:
        out_deg_cal="output/DEA/objects/calibrate_pvalues_interaction/deg_calibrated.tsv",
        out_summary="output/DEA/objects/calibrate_pvalues_interaction/calibration_summary.tsv",
        out_hist_raw_png="output/DEA/plots/calibrate_pvalues_interaction/P_hist_raw.png",
        out_hist_raw_pdf="output/DEA/plots/calibrate_pvalues_interaction/P_hist_raw.pdf",
        out_hist_cal_png="output/DEA/plots/calibrate_pvalues_interaction/P_hist_calibrated.png",
        out_hist_cal_pdf="output/DEA/plots/calibrate_pvalues_interaction/P_hist_calibrated.pdf"
    conda:
        "../envs/pydeseq2_interaction.yaml"
    script:
        "../scripts/DEA/calibrate_pvalues_interaction.R"


rule calibrate_pvalues_ifng_pos:
    input:
        in_deg="output/DEA/objects/pyDESeq2_ifng_pos/deg.tsv"
    output:
        out_deg_cal="output/DEA/objects/calibrate_pvalues_ifng_pos/deg_calibrated.tsv",
        out_summary="output/DEA/objects/calibrate_pvalues_ifng_pos/calibration_summary.tsv",
        out_hist_raw_png="output/DEA/plots/calibrate_pvalues_ifng_pos/P_hist_raw.png",
        out_hist_raw_pdf="output/DEA/plots/calibrate_pvalues_ifng_pos/P_hist_raw.pdf",
        out_hist_cal_png="output/DEA/plots/calibrate_pvalues_ifng_pos/P_hist_calibrated.png",
        out_hist_cal_pdf="output/DEA/plots/calibrate_pvalues_ifng_pos/P_hist_calibrated.pdf"
    conda:
        "../envs/pydeseq2_interaction.yaml"
    script:
        "../scripts/DEA/calibrate_pvalues_interaction.R"
