# path rel to root: input, output, singularity
# path rel to .smk file: script


# make plots summarizing the results from cellRanger QC 
rule plot_cellRangerQC:
    input: 
        in_cellRangerQC='data/full_dataset/cellRanger_summary.csv'
    output:
        out_plotsFlag='output/scRNAseq/A_preprocessing/plots/plot_cellRangerQC/flag.txt'
    params:
        out_plots = 'output/scRNAseq/A_preprocessing/plots/plot_cellRangerQC/', 
        custom_funs='workflow/scripts/funs/funs.R'
    resources:
        mem_mb=10000
    singularity:
        'workflow/envs/renv_0.2.sif'
    script:
        '../scripts/scRNAseq/A_preprocessing/plot_cellRangerQC.R'


# generate individual adata objects. One per batch. Add some metadata info.
rule perBatch_makeAdata:
    input:
        in_cellRanger_h5='data/full_dataset/{batch}/sample_filtered_feature_bc_matrix.h5',
        in_sampleMetadata='data/full_dataset/metadata.csv'
    output:
        out_adata=temp('output/scRNAseq/A_preprocessing/objects/{batch}/adata.h5ad')
    resources:
        mem_mb=10000
    singularity:
        'workflow/envs/scanpy_scirpy_0.4.sif'
    script:
        '../scripts/scRNAseq/A_preprocessing/perBatch_makeAdata.py'


# perform QC, sample filtering per batch and demultiplexing.
# includes figures with demultiplexing results.
rule perBatch_QCFilterDemultiplex:
    input:
        in_adata='output/scRNAseq/A_preprocessing/objects/{batch}/adata.h5ad'
    output:
        out_adata='output/scRNAseq/A_preprocessing/objects/{batch}/adata_filtered.h5ad',
        out_hashsolo='output/scRNAseq/A_preprocessing/objects/{batch}/hashsolo.tsv',
        out_QCres='output/scRNAseq/A_preprocessing/objects/{batch}/QCres.tsv',
        out_plotsFlag = 'output/scRNAseq/A_preprocessing/plots/perBatch_QCFilterDemultiplex/{batch}/flag.txt'
    params:
        out_plots = 'output/scRNAseq/A_preprocessing/plots/perBatch_QCFilterDemultiplex/{batch}/'
    resources:
        mem_mb=10000
    singularity:
        'workflow/envs/scanpy_scirpy_0.4.sif'
    script:
        '../scripts/scRNAseq/A_preprocessing/perBatch_QCFilterDemultiplex.py'
    
#make plots summarizing the results from scanpy QC 
rule plot_batchQC:
    input: 
        in_hashsolo=lambda w: expand('output/scRNAseq/A_preprocessing/objects/{batch}/hashsolo.tsv', batch=config['batches']),
        in_QCres=lambda w: expand('output/scRNAseq/A_preprocessing/objects/{batch}/QCres.tsv', batch=config['batches'])
    output:
        out_plotsFlag='output/scRNAseq/A_preprocessing/plots/plot_batchQC/flag.txt'
    params:
        out_plots = 'output/scRNAseq/A_preprocessing/plots/plot_batchQC/', 
        custom_funs='workflow/scripts/funs/funs.R'
    resources:
        mem_mb=10000
    singularity:
        'workflow/envs/renv_0.2.sif'
    script:
        '../scripts/scRNAseq/A_preprocessing/plot_batchQC.R'

# # create a single adata object from all the single-batch adatas
rule clean_metadata:
    input:
        in_OLINK_clinical='data/full_dataset/OLINK_clinical.csv',
        in_metadata='data/full_dataset/metadata.csv'
    output:
        out_OLINK_clinical='output/scRNAseq/A_preprocessing/objects/OLINK_clinical.tsv',
        out_metadata='output/scRNAseq/A_preprocessing/objects/metadata.tsv'
    resources:
        mem_mb=5000
    singularity:
        'workflow/envs/scanpy_scirpy_0.4.sif'
    script:
        '../scripts/scRNAseq/A_preprocessing/clean_metadata.py'

# # create a single adata object from all the single-batch adatas
rule clean_metadata_new:
    input:
        in_metadata='data/full_dataset/metadata.csv',
        # the eGFR_calc values for some patients are wrong 
        in_OLINK_clinical_old='data/full_dataset/OLINK_clinical.csv',
        # this matrix includes the correct eGFR_calc values
        in_clinical_simple='data/full_dataset/clinical_simple.csv'
    output:
        out_OLINK_clinical='output/scRNAseq/A_preprocessing/objects/OLINK_clinical_new.tsv'
    resources:
        mem_mb=5000
    singularity:
        'workflow/envs/scanpy_scirpy_0.4.sif'
    script:
        '../scripts/scRNAseq/A_preprocessing/clean_metadata_new.py'


# # create a single adata object from all the single-batch adatas
rule merge_batches:
    input:
        in_adata_lst=lambda w: expand('output/scRNAseq/A_preprocessing/objects/{batch}/adata_filtered.h5ad', batch=config['batches']),
        in_OLINK_clinical='output/scRNAseq/A_preprocessing/objects/OLINK_clinical.tsv'
    output:
        #out_adata=temp('output/scRNAseq/A_preprocessing/objects/merged/adata_all.h5ad')
        out_adata='output/scRNAseq/A_preprocessing/objects/merged/adata_all.h5ad'
    resources:
        mem_mb=20000
    singularity:
        'workflow/envs/scanpy_scirpy_0.4.sif'
    script:
        '../scripts/scRNAseq/A_preprocessing/merge_batches.py'

# # perform QC on merged adata. Integrate all batches into single projection.
rule merged_QCandHarmonize:
    input:
        in_adata='output/scRNAseq/A_preprocessing/objects/merged/adata_all.h5ad'
    output:
        out_adata='output/scRNAseq/A_preprocessing/objects/merged/adata_harmonized.h5ad',
        out_plotsFlag = 'output/scRNAseq/A_preprocessing/plots/merged_QCandHarmonize/flag.txt'
    params:
        out_plots = 'output/scRNAseq/A_preprocessing/plots/merged_QCandHarmonize/'
    resources:
        mem_mb=120000, 
        disk_mb=20000
    singularity:
        'workflow/envs/scanpy_scirpy_0.4.sif'
    script:
        '../scripts/scRNAseq/A_preprocessing/merged_QCandHarmonize.py'

# save subset of data for each cluster to facilitate manual cell annotation
## Todo: 
    ## - How to access the wildcard within the script? For now you don't even use the filenames in the out_adata_lst to save the oputput.
    ## as sub-optimal fix the cellCluster_filename in params is used.
rule save_adata_byCluster:
    input:
        in_adata='output/scRNAseq/A_preprocessing/objects/merged/adata_harmonized.h5ad'
    output:
        out_adata_lst=expand('output/scRNAseq/A_preprocessing/objects/merged/clust{cell_cluster}/adata.h5ad', cell_cluster=['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14'])
    params:
        cellCluster_filename={ cell_cluster:f"output/scRNAseq/A_preprocessing/objects/merged/clust{cell_cluster}/adata.h5ad" 
        for cell_cluster in ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14']},
    resources:
        mem_mb=120000, 
        disk_mb=20000
    singularity:
        'workflow/envs/scanpy_scirpy_0.4.sif'
    script:
        '../scripts/scRNAseq/A_preprocessing/save_adata_byCluster.py'
    
# cell annotation
rule merged_annotateCells_lvl1:
    input:
        in_adata='output/scRNAseq/A_preprocessing/objects/merged/adata_harmonized.h5ad',
        in_cellMarkers='data/markers_pbmc/cellMarkers_all_lvl1.txt',
        #in_cellMarkers_l2=lambda w: expand('data/markers_pbmc/DE_geneMarkers_Level2/{cell_type}_lvl2.txt', cell_type=config['cell_type'])
    output:
        ## subset of adata with 3 batches
        out_adata='output/scRNAseq/A_preprocessing/objects/merged/adata_lvl1_3b.h5ad',
        out_clusterMarkers='output/scRNAseq/A_preprocessing/objects/merged/clusterMarkers_25.tsv',
        out_gsea='output/scRNAseq/A_preprocessing/objects/merged/gsea_top3.tsv',
        out_ora='output/scRNAseq/A_preprocessing/objects/merged/ora_top3.tsv',
        out_plotsFlag = 'output/scRNAseq/A_preprocessing/plots/merged_annotateCells_lvl1/flag.txt'
    params:
        colname_clusters='leiden_0.4_new',
        out_plots = 'output/scRNAseq/A_preprocessing/plots/merged_annotateCells_lvl1/'
    resources:
        mem_mb=20000 
    singularity:
        'workflow/envs/scanpy_scirpy_0.4.sif'
    script:
        '../scripts/scRNAseq/A_preprocessing/merged_annotateCells_lvl1.py'

# cell annotation 2, with leiden_0.4
rule merged_annotateCells_lvl2:
    input:
        in_adata='output/scRNAseq/A_preprocessing/objects/merged/adata_harmonized.h5ad'
    output:
        # full anndata
        out_adata='output/scRNAseq/A_preprocessing/objects/merged/adata_lvl2.h5ad',
        # anndata for 3 batches only
        out_adata_3b='output/scRNAseq/A_preprocessing/objects/merged/adata_lvl2_3b.h5ad',
        out_clusterMarkers_all='output/scRNAseq/A_preprocessing/objects/merged/lvl2_clusterMarkers_all.tsv',
        out_clusterMarkers_25='output/scRNAseq/A_preprocessing/objects/merged/lvl2_clusterMarkers_25.tsv',
        out_ora='output/scRNAseq/A_preprocessing/objects/merged/lvl2_ora.tsv',
        out_plotsFlag = 'output/scRNAseq/A_preprocessing/plots/merged_annotateCells_lvl2/flag.txt'
    params:
        folder_custom_functions='workflow/scripts/funs/',
        out_plots = 'output/scRNAseq/A_preprocessing/plots/merged_annotateCells_lvl2/', 
        out_objects = 'output/scRNAseq/A_preprocessing/objects/merged/', 
        colname_clusters='leiden_0.4_new',
        # folder where the input data tables with the markers are stored (HH lab markers)
        folder_subclusters = 'data/markers_pbmc/DE_geneMarkers_Level2/',
        # HVGs are found within each batch. 
        # In what fraction of the batches should a gene be a HVG to be considered?
        fraction_batches_keep_HVG=0.7,
        # markers to examine in the subclusters of this cluster.
        subclusters = {
            '00': ['T_CD4_NonNaive'], # T_CD4_NonNaive
            '01': [''], # T_CD4_naive
            '02': ['T_CD8_NonNaive', 'UTC'], # T_CD8_NonNaive and some gd
            '03': [''], # T_CD8_naive
            '04': ['B', 'Progenitors', 'Plasma'], # naive and non-naive B
            '05': ['Plasma'], # Plasma
            '06': ['Mono'], # Mono Classical 
            '07': ['Mono'], # Mono non-classical
            '08': ['ILC'], # NKs
            '09': ['UTC'], # MAIT
            '10': [''], # DC
            '11': [''], # pDCs
            '12': [''], # ukn_01
            '13': [''], # ukn_02
            '14': ['']}, # ukn_03
        
        # clustering parameters for the subclusters
        clustering_params = {'00':{'do_subCluster': True, 'resol':0.18}, #
                        '01':{'do_subCluster': False, 'resol':0}, #
                        '02':{'do_subCluster': True, 'resol':0.2}, #
                        '03':{'do_subCluster': False, 'resol':0}, #
                        '04':{'do_subCluster': True, 'resol':0.05}, #
                        '05':{'do_subCluster': True, 'resol':0.01}, #
                        '06':{'do_subCluster': True, 'resol':0.2}, 
                        '07':{'do_subCluster': False, 'resol':0}, 
                        '08':{'do_subCluster': True, 'resol':0.1}, 
                        '09':{'do_subCluster': False, 'resol':0}, 
                        '10':{'do_subCluster': False, 'resol':0}, 
                        '11': {'do_subCluster': False, 'resol': 0}, 
                        '12': {'do_subCluster': False, 'resol': 0}, 
                        '13': {'do_subCluster': False, 'resol': 0}, 
                        '14': {'do_subCluster': False, 'resol': 0}}
    resources:
        mem_mb=20000 
    singularity:
        'workflow/envs/scanpy_scirpy_0.4.sif'
    script:
        '../scripts/scRNAseq/A_preprocessing/merged_annotateCells_lvl2.py'


# Defining cell-type annotation based on previous annotation scripts, in short:
## Define lvl2 cell-types (some sub-clusters are merged) and add annotations
## generate a file with the top 25 markers
## label low Q clusters
## generate plots: marker plots, UMAPS
rule merged_defineAnnot:
    input:
        # full anndata with lvl2 clusters
        in_adata='output/scRNAseq/A_preprocessing/objects/merged/adata_lvl2.h5ad'
    output:
        # anndata objects
        out_adata='output/scRNAseq/A_preprocessing/objects/merged/adata_defineAnnot.h5ad',
        out_adata_3b='output/scRNAseq/A_preprocessing/objects/merged/adata_defineAnnot_3b.h5ad', # anndata for 3 batches only
        # markers
        out_clusterMarkers_all='output/scRNAseq/A_preprocessing/objects/merged/defineAnnot_clusterMarkers_all.tsv',
        out_clusterMarkers_25='output/scRNAseq/A_preprocessing/objects/merged/defineAnnot_clusterMarkers_25.tsv',
        out_plotsFlag = 'output/scRNAseq/A_preprocessing/plots/merged_defineAnnot/flag.txt'
    params:
        folder_custom_functions='workflow/scripts/funs/',
        out_plots = 'output/scRNAseq/A_preprocessing/plots/merged_defineAnnot/', 
        out_objects = 'output/scRNAseq/A_preprocessing/objects/merged/', 
        colname_clusters='leiden_0.4_new',
        # HVGs are found within each batch. 
        # In what fraction of the batches should a gene be a HVG to be considered?
        fraction_batches_keep_HVG=0.7,
        # annotations 
        ## lvl1
        annotation_celltype_lvl1 = {
            '00': 'T_CD4_nonNaive',
            '01': 'T_CD4_naive',
            '02': 'T_CD8_nonNaive_UTC',
            '03': 'T_CD8_naive',
            '04': 'B',
            '05': 'Plasma',
            '06': 'Mono_c',
            '07': 'Mono_nonC',
            '08': 'NK',
            '09': 'MAIT',
            '10': 'DC',
            '11': 'pDC',
            '12': 'ukn_01',
            '13': 'ukn_02',
            '14': 'ukn_03'},
        ## lvl2
        annotation_celltype_lvl2 = {
            '00,0': 'T_CD4_CM', '00,1': 'T_CD4_eff', '00,2': 'T_CD4_EM', '00,3': 'T_CD4_Tregs',
            '01': 'T_CD4_naive',
            '02,0': 'T_CD8_CM', '02,1': 'T_CD8_EM', '02,2': 'T_gd', '02,3': 'T_CD8_EM', '02,4': 'T_CD8_activated',
            '03': 'T_CD8_naive',
            '04,0': 'B_transitional', '04,1': 'B_mem_switched', '04,2': 'B_lowQ', 
            '05,0': 'Plasma', '05,1': 'Plasma_lowQ', '05,2': 'Plasma_lowQ', '05,3': 'Plasma_lowQ', '05,4': 'Plasma_lowQ', '05,5': 'Plasma_lowQ', '05,6': 'Plasma_lowQ',
            '06,0': 'Mono_c', '06,1': 'Mono_c_AP', '06,2': 'Mono_c_IFN', '06,3': 'Mono_c_lowQ',
            '07': 'Mono_nonC',
            '08,0': 'NK_CD16', '08,1': 'NK_CD16', '08,2': 'NK_CD56',
            '09': 'MAIT',
            '10': 'DC',
            '11': 'pDC',
            '12': 'ukn_01',
            '13': 'ukn_02',
            '14': 'ukn_03'}
    resources:
        mem_mb=20000 
    singularity:
        'workflow/envs/scanpy_scirpy_0.4.sif'
    script:
        '../scripts/scRNAseq/A_preprocessing/merged_defineAnnot.py'

# UMAP density plots
rule merged_plotDensity:
    input:
        in_adata='output/scRNAseq/A_preprocessing/objects/merged/adata_defineAnnot.h5ad',
        in_clinical='data/full_dataset/clinical_simple_matchedGroups.txt'
    output:
        out_plotsFlag = 'output/scRNAseq/A_preprocessing/plots/merged_plotDensity/flag.txt'
    params:
        out_plots = 'output/scRNAseq/A_preprocessing/plots/merged_plotDensity/'
    resources:
        mem_mb=20000 
    singularity:
        'workflow/envs/scanpy_scirpy_0.4.sif'
    script:
        '../scripts/scRNAseq/A_preprocessing/merged_plotDensity.py'

# save adata.obs from the full data to debug rule below
rule merged_save_adata_obs:
    input:
        in_adata='output/scRNAseq/A_preprocessing/objects/merged/adata_defineAnnot.h5ad'
    output:
        out_adata_obs='output/scRNAseq/A_preprocessing/objects/merged/adata_obs.tsv'
    resources:
        mem_mb=20000 
    singularity:
        'workflow/envs/scanpy_scirpy_0.4.sif'
    script:
        '../scripts/scRNAseq/A_preprocessing/merged_save_adata_obs.py'

# prepare input for compositions
## the actual; compositional analysis is not run in snakemake but locally 
## because of the python package for clr transformation, it's not in the singularity image
## so a new one has to be created
rule merged_compositions_prepInput:
    input:
        in_adata='output/scRNAseq/A_preprocessing/objects/merged/adata_defineAnnot.h5ad',
        in_clinical='data/full_dataset/clinical_simple_matchedGroups.txt'
    output:
        out_tab = 'output/scRNAseq/A_preprocessing/objects/merged/compositions_prop.tsv'
    resources:
        mem_mb=20000 
    singularity:
        'workflow/envs/scanpy_scirpy_0.4.sif'
    script:
        '../scripts/scRNAseq/A_preprocessing/merged_compositions_prepInput.py'



# # # temporary rule to get the avr/median reads/cell across batches
# # rule tmp_seqDepth:
# #     input:
# #         in_adata='output/scRNAseq/A_preprocessing/objects/merged/adata_lvl2.h5ad'
# #     output:
# #         out_plotsFlag = 'output/scRNAseq/A_preprocessing/plots/tmp_seqDepth/flag.txt'
# #     params:
# #         out_plots = 'output/scRNAseq/A_preprocessing/plots/tmp_seqDepth/'
# #     resources:
# #         mem_mb=20000 
# #     singularity:
# #         'workflow/envs/scanpy_scirpy_0.4.sif'
# #     script:
# #         '../scripts/scRNAseq/A_preprocessing/tmp_seqDepth.py'




# # # optimize parameters at the clustering step:
# # ## HVGs: fraction of batches where a gene is a HVG to be kept.
# # ## leiden resolution.
# # rule merged_annotateCells_lvl2_optimClustParams:
# #     input:
# #         in_adata='output/scRNAseq/A_preprocessing/objects/merged/adata_harmonized.h5ad'
# #     output:
# #         out_plotsFlag = 'output/scRNAseq/A_preprocessing/plots/merged_annotateCells_lvl2_optimClustParams/flag.txt'
# #     params:
# #         folder_custom_functions='workflow/scripts/funs/',
# #         out_plots = 'output/scRNAseq/A_preprocessing/plots/merged_annotateCells_lvl2_optimClustParams/', 
# #         # clustering parameters for the subclusters
# #         clustering_params = {'0':{'fraction_batches_keep_HVG': [0.5, 0.6, 0.7], 'resol':[0.2, 0.17, 0.15, 0.1]}, # nHVGs: 443
# #                         '1':{'fraction_batches_keep_HVG': [0.5, 0.6, 0.7], 'resol':[0.17, 0.15, 0.1]}, # 987
# #                         '2':{'fraction_batches_keep_HVG': [0.5, 0.6, 0.7], 'resol':[0.1]}, # 359
# #                         '3':{'fraction_batches_keep_HVG': [0.5, 0.6, 0.7], 'resol':[0.1, 0.05]}, # 433
# #                         '4':{'fraction_batches_keep_HVG': [0.5, 0.6, 0.7], 'resol':[0.1, 0.05]}, # 322
# #                         '5':{'fraction_batches_keep_HVG': [0.5, 0.6, 0.7], 'resol':[0.1, 0.06]}, # 168
# #                         '10':{'fraction_batches_keep_HVG': [0.5, 0.6, 0.7], 'resol':[0.01, 0.008, 0.005]}} # 246
# #     resources:
# #         mem_mb=20000 
# #     singularity:
# #         'workflow/envs/scanpy_scirpy_0.4.sif'
# #     script:
# #         '../scripts/scRNAseq/A_preprocessing/merged_annotateCells_lvl2_optimClustParams.py'

# # cell-typist on lvl2 clusters (high Q cells only)
# # rule merged_cellTypist_lvl2:
# #     input:
# #         in_adata='output/scRNAseq/A_preprocessing/objects/merged/adata_harmonized.h5ad'
# #     output:
# #         # full anndata
# #         out_adata='output/scRNAseq/A_preprocessing/objects/merged/adata_lvl2.h5ad',
# #         # anndata for 3 batches only
# #         out_adata_3b='output/scRNAseq/A_preprocessing/objects/merged/adata_lvl2_3b.h5ad',
# #         out_clusterMarkers_all='output/scRNAseq/A_preprocessing/objects/merged/lvl2_clusterMarkers_all.tsv',
# #         out_clusterMarkers_25='output/scRNAseq/A_preprocessing/objects/merged/lvl2_clusterMarkers_25.tsv',
# #         out_ora='output/scRNAseq/A_preprocessing/objects/merged/lvl2_ora.tsv',
# #         out_plotsFlag = 'output/scRNAseq/A_preprocessing/plots/merged_annotateCells_lvl2/flag.txt'
# #     params:
# #         folder_custom_functions='workflow/scripts/funs/',
# #         out_plots = 'output/scRNAseq/A_preprocessing/plots/merged_annotateCells_lvl2/', 
# #         out_objects = 'output/scRNAseq/A_preprocessing/objects/merged/', 
# #         # folder where the input data tables with the markers are stored (HH lab markers)
# #         folder_subclusters = 'data/markers_pbmc/DE_geneMarkers_Level2/',
# #         # HVGs are found within each batch. 
# #         # In what fraction of the batches should a gene be a HVG to be considered?
# #         fraction_batches_keep_HVG=0.7,
# #         # markers to examine in the subclusters of this cluster.
# #         subclusters = {
# #             '0': ['T_CD4_NonNaive'],
# #             '1': ['Mono', 'DC'], # very likely just Mono
# #             '2': ['T_CD8_NonNaive', 'UTC', 'ILC'], # likely different T_CD8_NonNaive
# #             '3': ['ILC', 'UTC'], # very likely NKs, explore UTC just in case
# #             '4': ['B', 'DC', 'Progenitors', 'Plasma'], # very likely B, maybe a small subcluster of DC in the outskirts?
# #             '5': ['Progenitors'], # how to distinguish CD4-CD8 naive?
# #             '6': ['Mono', 'DC'], # very likely just Mono
# #             '7': [''], # (subcluster and make QC plots, a lot of MT genes)
# #             '8': ['DC', 'Mono'], # very likely DC, explore Mono just in case
# #             '9': [''], # (subcluster and make QC plots, a lot of MT genes)
# #             '10': ['DC'] },
# #         # clustering parameters for the subclusters
# #         clustering_params = {'0':{'do_subCluster': True, 'resol':0.2}, #
# #                         '1':{'do_subCluster': True, 'resol':0.17}, #
# #                         '2':{'do_subCluster': True, 'resol':0.18}, #
# #                         '3':{'do_subCluster': True, 'resol':0.1}, #
# #                         '4':{'do_subCluster': True, 'resol':0.05}, #
# #                         '5':{'do_subCluster': True, 'resol':0.06}, #
# #                         '6':{'do_subCluster': False, 'resol':0.1}, #
# #                         '7':{'do_subCluster': False, 'resol':0.2}, 
# #                         '8':{'do_subCluster': False, 'resol':0.2}, 
# #                         '9':{'do_subCluster': False, 'resol':0.2}, 
# #                         '10':{'do_subCluster': True, 'resol':0.015}}
# #     resources:
# #         mem_mb=20000 
# #     singularity:
# #         'workflow/envs/scanpy_scirpy_0.4.sif'
# #     script:
# #         '../scripts/scRNAseq/A_preprocessing/merged_annotateCells_lvl2.py'

# # add cell-types/sub-types to adata (after manually checking enriched markers)


# # rule make_adata:
# #     input:
# #         in_cellRanger='data/scRNAseq/',
# #         in_cellMetadata='data/scRNAseq/cell_metadata.tsv'
# #     output:
# #         out_adata='output/scRNAseq/A_preprocessing/objects/adata.h5ad'
# #     resources:
# #         mem_mb=10000
# #     singularity:
# #         'workflow/envs/scanpy_scirpy_0.2.sif'
# #     script:
# #         '../scripts/scRNAseq/A_preprocessing/01_make_adata.py'

# # rule QCandFilter_adata:
# #     input:
# #         in_adata='output/scRNAseq/A_preprocessing/objects/adata.h5ad'
# #     output:
# #         out_adata='output/scRNAseq/A_preprocessing/objects/01_adata.h5ad',
# #         out_plotsFlag = 'output/scRNAseq/A_preprocessing/plots/02_QC&filter_adata/flag.txt'
# #     params:
# #         out_plots = 'output/scRNAseq/A_preprocessing/plots/02_QC&filter_adata/'
# #     resources:
# #         mem_mb=10000
# #     singularity:
# #         'workflow/envs/scanpy_scirpy_0.2.sif'
# #     script:
# #         '../scripts/scRNAseq/A_preprocessing/02_QC&filter_adata.py'

# # rule harmony:
# #     input:
# #         in_adata='output/scRNAseq/A_preprocessing/objects/01_adata.h5ad'
# #     output:
# #         out_adata='output/scRNAseq/A_preprocessing/objects/03_adata.h5ad',
# #         out_plotsFlag = 'output/scRNAseq/A_preprocessing/plots/03_harmony/flag.txt'
# #     params:
# #         out_plots = 'output/scRNAseq/A_preprocessing/plots/03_harmony/'
# #     resources:
# #         mem_mb=10000
# #     singularity:
# #         'workflow/envs/scanpy_scirpy_0.2.sif'
# #     script:
# #         '../scripts/scRNAseq/A_preprocessing/03_harmony.py'

# # rule format_metadata:
# #     input:
# #         in_metadata='data/IgAN_OLINK_4x4.csv'
# #     output:
# #         out_metadata = 'output/scRNAseq/A_preprocessing/objects/metadata_OLINK.tsv'
# #     resources:
# #         mem_mb=1000
# #     singularity:
# #         'workflow/envs/scanpy_scirpy_0.2.sif'
# #     script:
# #         '../scripts/scRNAseq/A_preprocessing/04_format_metadata.py'
