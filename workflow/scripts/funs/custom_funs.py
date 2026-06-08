# Custom functions

# Required packages for: Comparing ranked vectors
from __future__ import division

import math
from bisect import bisect_left
from collections import namedtuple


############################
##### Cell annotation ######
############################

# These functions streamline cell annotation and make the code more readable.

# function to pre-process adata including: 
## identify HVGs, PCA, harmonization, neighborhood graph, umap embedding and leiden sub-clustering.
## returns: adata object with updated slots, saves top HVGs by .
def perform_cluster_subclustering(adata_subset, cluster_id='0',
HVGs_batchKey=None, fraction_batches_keep_HVG=0.8, HVG_save_tsv=False, HVG_tsv_path=None, 
harmonize=True, harmony_key=[], 
neighbors_nNeighbors=15, neighbors_nPCs=10, 
leiden_resol=0.1):
    """ 
    Function to perform data preprocessing in the context of cluster (cell-type) subclustering 
    Input: sliced adata object for the cell cluster of interest.
    Performs: 1) identify HVGs (saves tsv with top 20 HVGs by norm_dispersion), 2) PCA, 3) batch harmonization,
              4) neighborhood graph, 5) umap embedding, 6) leiden sub-clustering.
    Returns: adata with updated slots.
    """
    import pandas as pd
    import numpy as np
    import re
    import scanpy as sc
    import scanpy.external as sce
    # identify highly variable genes
    print('\n### Identify HVGs ###')
    sc.pp.highly_variable_genes(adata_subset, batch_key=HVGs_batchKey, min_mean=0.0125, max_mean=3, min_disp=0.5)
    ## assign HVG if it is highly variable in at least fraction_batches_keep_HVG
    nbatches_keep_HVG=round(max(adata_subset.var['highly_variable_nbatches'].unique())*fraction_batches_keep_HVG)
    adata_subset.var['highly_variable']=adata_subset.var['highly_variable_nbatches']>=nbatches_keep_HVG
    ## remove TCR/BCR from HVGs
    adata_subset.var['is_TCR_BCR']=[re.search(pattern='^(IG[KHL])|(TR[ABGD])', string=x)!=None for x in adata_subset.var.index]
    adata_subset.var['highly_variable']=(adata_subset.var['highly_variable'])&(adata_subset.var['is_TCR_BCR']==False)
    print(f"- HVGs in equalmore than {nbatches_keep_HVG} batches were kept.")
    print(f"- n HVGs to be used downstream (PCA, etc): {sum(adata_subset.var['highly_variable'])}")
    ## save tsv
    if HVG_save_tsv and HVG_tsv_path:
        print(f"- Saving tsv with top 20 HVGs by descending norm_dispersion")
        df=adata_subset.var.copy()
        df=df[df['highly_variable']==True].reset_index(drop=False).sort_values(by='dispersions_norm', ascending=False).head(20)
        df.to_csv(HVG_tsv_path, sep='\t', index=False)

    # run PCA
    print('\n### Run PCA ###\n')
    sc.tl.pca(adata_subset, svd_solver='arpack', use_highly_variable=True)

    # harmonize (adjusts PCs)
    if harmonize and len(harmony_key)>=1:
        print('\n### Run harmony ###')
        print(f"- Bathes to harmonize: {harmony_key}")
        for h_key in harmony_key:
            sce.pp.harmony_integrate(adata_subset, key=h_key, max_iter_harmony=50)
            adata_subset.obsm['X_pca'] = adata_subset.obsm['X_pca_harmony']

    # compute the neighborhood graph
    print('\n### Compute the neighborhood graph ###\n')
    sc.pp.neighbors(adata_subset, n_neighbors=neighbors_nNeighbors, n_pcs=neighbors_nPCs, use_rep='X_pca_harmony') 
    # embedding the neighborhood graph
    print('\n### Embed the neighborhood graph as UMAP###\n')
    sc.tl.umap(adata_subset)
    # get finer clusters
    print('\n### Get finer clusters with Leiden ###')
    print(f"- resolution to use: {leiden_resol}")
    sc.tl.leiden(adata_subset, key_added='subclust', resolution=leiden_resol)
    ## name them accordingly
    adata_subset.obs['subclust']=[cluster_id+','+x for x in adata_subset.obs['subclust']]
    ## how many clusters?
    print(f"- Total subclusters: {adata_subset.obs['subclust'].nunique()}")
    print(f"- Subclusters with > 10 cells: {sum(adata_subset.obs['subclust'].value_counts() > 10)}")
    print(f"- n cells per subcluster:\n{adata_subset.obs['subclust'].value_counts().head(10)}")
    
    return(adata_subset)



# function to perform cell-type annotation 
def annotate_cells(adata_subset, cell_cluster='0',groups=[], groupby='subclust', 
                   method='t-test_overestim_var', markers=None, topn_markers=20, topn_DEGs=40):
    """ 
    Function to perform cell type annotation bases on overrepresentation 
    analysis.
    Input: anndata object
    Performs: First it gets DEGs for each group in 'groups'. 
              Then, it performs ORA for each group using the DEGs and the reference genes in markers.
    Returns: (adata_subset, DEGs, ORA)
    Arguments:
        - adata_subset: adata object
        - cell_cluster: str. id of the cell cluster on which the subclustering was performed
        - groups: list of strings. List of ids for which the DEGs will be generated. Must be in .obs column specified in 'groupby'.
        - groupby: column in .obs where the group ids are to perform the DEA.
        - method: DEA method in sc.tl.rank_genes_groups().
        - markers: reference gene markers for cell annotation. This is the query set for ORA.
        - topn_markers: how many top genes to use in the query set? (markers). Rank is performed on 'scores' column of markers.
        - topn_DEGs: how many top genes to use in the foreground? (DEGs). Rank is performed on 'scores' column.
    """
    import pandas as pd
    import numpy as np
    import scanpy as sc
    import decoupler as dc

    # get markers per subcluster
    sc.tl.rank_genes_groups(adata_subset, groupby=groupby, groups=groups, method=method)
    res_lst=[]
    for cell_subcluster in groups:
        #### get result for each subcluster
        df=sc.get.rank_genes_groups_df(adata_subset, group=cell_subcluster).sort_values(by='scores', ascending=False)
        #### label subcluster
        df['cell_subcluster']=cell_subcluster
        #### add to list
        res_lst.append(df)
    degs=pd.concat(res_lst).rename(columns={'names': 'gene_name'})
    degs['cell_cluster']=cell_cluster

    ## perform ORA enrichment
    print(f"- Perform ORA enrichment")
    res_ora = pd.DataFrame()
    ### iterate over each cell subcluster
    for cell_subcluster in degs['cell_subcluster'].unique():
        #### filtered df only including degs for current subcluster
        df=degs.copy()
        df=df[df['cell_subcluster']==cell_subcluster].set_index('gene_name').head(topn_DEGs)
        #### n background genes (rest of the genes)
        n_background=degs[degs['cell_subcluster']==cell_subcluster].shape[0] - topn_DEGs
        #### markers
        markers=markers.sort_values(by=['cell_type', 'scores'],ascending=[True, False]).groupby('cell_type').head(topn_markers).reset_index(drop=True)
        #### perform enrichment
        res = dc.get_ora_df(df=df, 
        net = markers, source='cell_type', target='symbol', n_background=n_background,
        verbose=True)
        #### format and store result
        res['cell_subcluster']=cell_subcluster
        res_ora = pd.concat([res_ora, res])
    ### save ORA table for plots 
    res_ora['cell_cluster']=cell_cluster
    res_ora=res_ora.sort_values(by=['cell_subcluster', 'Combined score'], ascending=[True, False])

    return((adata_subset, degs, res_ora))


# function to make diagnostic plots that help: 
## to choose pre-processing parameters,
## to visualize sub-clusters and determine if overclustering happened
## saves: plots
## returns: nothing
def make_diagnostic_plots(adata_subset, path_folder:'', plot_suffix=''):
    """ 
    Function to generate and save diagnostic plots. These are intended to help:
    - choose pre-processing parameters,
    - visualize sub-clusters,
    - determine if overclustering happened.
    perform_cluster_subclustering() and annotate_cells() should have been run first.
    """
    import pandas as pd
    import numpy as np
    import scanpy as sc
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import seaborn as sns
    
    print(f"\n### Make diagnostic plots ###\n")
    # PCA
    ## PCA variance ratio
    _=plt.figure()
    _=sc.pl.pca_variance_ratio(adata_subset, log=False, n_pcs=50, show=False)
    plt.savefig(path_folder+plot_suffix+"_pca_01.pdf", bbox_inches="tight")
    ## PCAs
    _=plt.figure()
    _=sc.pl.pca(adata_subset, color=['subclust', 'batch'], annotate_var_explained=True, components = ['1,2', '3,4', '5,6', '7,8', '9,10'], show=False, ncols=2)
    plt.savefig(path_folder+plot_suffix+"_pca_02.pdf", bbox_inches="tight")
    # _=plt.figure()
    # _=sc.pl.pca(adata_subset, color=['subclust'], annotate_var_explained=True, components = ['1,2', '3,4', '5,6'], show=False )
    # plt.savefig(path_folder+plot_suffix+"_pca_02.pdf", bbox_inches="tight")
    # _=plt.figure()
    # _=sc.pl.pca(adata_subset, color=['subclust'], annotate_var_explained=True, components = ['7,8', '9,10'], show=False )
    # plt.savefig(path_folder+plot_suffix+"_pca_03.pdf", bbox_inches="tight")

    # UMAP
    keys_lst=['subclust', 'batch', 'n_genes_by_counts', 'log1p_n_genes_by_counts', 
    'total_counts', 'log1p_total_counts', 'pct_counts_in_top_20_genes', 'log10nGenesPerTranscripts',
    'pct_counts_mt', 'total_counts_mt', 'pct_counts_ribo', 'total_counts_ribo',
    'pct_counts_hb', 'total_counts_hb',
    'doublet_score']
    keys_lst=[x for x in keys_lst if x in adata_subset.obs.columns]
    _=plt.figure()
    _=sc.pl.umap(adata_subset, color=keys_lst, frameon=False, show=False, ncols=2)
    plt.savefig(path_folder+plot_suffix+"_UMAP_01.pdf", bbox_inches="tight")

    # top 5 markers per subcluster
    _=plt.figure()
    _=sc.pl.rank_genes_groups_dotplot(adata_subset, groupby='subclust', standard_scale="var", n_genes=8, return_fig=True)
    _.savefig(path_folder+plot_suffix+"_markers.pdf")
    
    # n cells per patient-cluster
    plotTab=pd.DataFrame(adata_subset.obs[['subclust', 'patient_number']].value_counts()).reset_index()
    plotTab.columns.values[2]='count'
    print(f"\n- plotTab.head():\n{plotTab.head()}\n")
    _=plt.figure()
    _=sns.barplot(plotTab, x='subclust', y='count', hue='patient_number', palette='viridis', legend=False)
    plt.savefig(path_folder+plot_suffix+"_nCells_01.pdf", bbox_inches="tight")

    # n cells per condition-cluster
    plotTab=pd.DataFrame(adata_subset.obs[['subclust', 'condition']].value_counts()).reset_index()
    plotTab.columns.values[2]='count'
    _=plt.figure()
    _=sns.barplot(plotTab, x='subclust', y='count', hue='condition', palette='viridis', legend=True)
    plt.savefig(path_folder+plot_suffix+"_nCells_02.pdf", bbox_inches="tight")

    # proportion of cells
    ## prepare data
    df=adata_subset.obs.copy()
    ### n cells per cluster
    df['n_cells_clust']=df.groupby('subclust')['patient_number'].transform('count')
    ### prop cells per cluster-batch
    df['n_cells_clust-batch']=df.groupby(['subclust', 'batch'])['patient_number'].transform('count')
    df['prop_cells_clust-batch']=df['n_cells_clust-batch']/df['n_cells_clust']
    ### prop cells cluster-pat
    df['n_cells_clust-pat']=df.groupby(['subclust', 'patient_number'])['batch'].transform('count')
    df['prop_cells_clust-pat']=df['n_cells_clust-pat']/df['n_cells_clust']
    ## Plot: prop per cluster-patient
    plotTab = df[['subclust', 'patient_number','prop_cells_clust-pat']].copy().drop_duplicates()
    _=plt.figure()
    _=sns.barplot(plotTab, x='subclust', y='prop_cells_clust-pat', hue='patient_number', palette='viridis', legend=False)
    plt.savefig(path_folder+plot_suffix+"_propCells_01.pdf", bbox_inches="tight")
    ## Plot: prop per cluster-batch
    plotTab = df[['subclust', 'batch','prop_cells_clust-batch']].copy().drop_duplicates()
    _=plt.figure()
    _=sns.barplot(plotTab, x='subclust', y='prop_cells_clust-batch', hue='batch')
    plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
    plt.savefig(path_folder+plot_suffix+"_propCells_02.pdf", bbox_inches="tight")

    # QC figures in subclusters
    keys_lst=[('n_genes_by_counts', 'log1p_n_genes_by_counts'), 
    ('total_counts', 'log1p_total_counts'), ('pct_counts_in_top_20_genes', 'log10nGenesPerTranscripts'),
    ('pct_counts_mt', 'total_counts_mt'), ('pct_counts_ribo', 'total_counts_ribo'),
    ('pct_counts_hb', 'total_counts_hb'),'doublet_score']
    keys_lst=[x for x in keys_lst if all(pd.Series(x).isin(adata_subset.obs.columns))]
    for i,var in enumerate(keys_lst):
        if type(var) is tuple:
            var=list(var)
        with mpl.rc_context({'figure.figsize': (6, 6)}):
            _=plt.figure()
            _=sc.pl.violin(adata_subset, keys=var, groupby='subclust', stripplot=True, show=False)
            plt.savefig(path_folder+plot_suffix+f"_violin_0{i}.pdf", bbox_inches="tight")



# function to make diagnostic plots that help: 
## to choose pre-processing parameters,
## to visualize sub-clusters and determine if overclustering happened
## saves: plots
## returns: nothing
## No need to run annotate_cells() beforehand.
def make_diagnostic_plots2(adata_subset, path_folder:'', plot_suffix=''):
    """ 
    Function to generate and save diagnostic plots. These are intended to help:
    - choose pre-processing parameters,
    - visualize sub-clusters,
    - determine if overclustering happened.
    perform_cluster_subclustering() should have been run first.
    """
    import pandas as pd
    import numpy as np
    import scanpy as sc
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import seaborn as sns
    
    print(f"\n### Make diagnostic plots ###\n")
    # PCA
    ## PCA variance ratio
    _=plt.figure()
    _=sc.pl.pca_variance_ratio(adata_subset, log=False, n_pcs=50, show=False)
    plt.savefig(path_folder+plot_suffix+"_pca_01.pdf", bbox_inches="tight")
    ## PCAs
    _=plt.figure()
    _=sc.pl.pca(adata_subset, color=['subclust', 'batch'], annotate_var_explained=True, components = ['1,2', '3,4', '5,6', '7,8', '9,10'], show=False, ncols=2)
    plt.savefig(path_folder+plot_suffix+"_pca_02.pdf", bbox_inches="tight")

    # UMAP
    _=plt.figure()
    _=sc.pl.umap(adata_subset, color=['subclust', 'batch'], frameon=False, show=False)
    plt.savefig(path_folder+plot_suffix+"_UMAP_01.pdf", bbox_inches="tight")
    
    # n cells per patient-cluster
    plotTab=pd.DataFrame(adata_subset.obs[['subclust', 'patient_number']].value_counts()).reset_index()
    plotTab.columns.values[2]='count'
    _=plt.figure()
    _=sns.barplot(plotTab, x='subclust', y='count', hue='patient_number', palette='viridis', legend=False)
    plt.savefig(path_folder+plot_suffix+"_nCells_01.pdf", bbox_inches="tight")

    # n cells per condition-cluster
    plotTab=pd.DataFrame(adata_subset.obs[['subclust', 'condition']].value_counts()).reset_index()
    plotTab.columns.values[2]='count'
    _=plt.figure()
    _=sns.barplot(plotTab, x='subclust', y='count', hue='condition', palette='viridis', legend=True)
    plt.savefig(path_folder+plot_suffix+"_nCells_02.pdf", bbox_inches="tight")

    # proportion of cells
    ## prepare data
    df=adata_subset.obs.copy()
    ### n cells per cluster
    df['n_cells_clust']=df.groupby('subclust')['patient_number'].transform('count')
    ### prop cells per cluster-batch
    df['n_cells_clust-batch']=df.groupby(['subclust', 'batch'])['patient_number'].transform('count')
    df['prop_cells_clust-batch']=df['n_cells_clust-batch']/df['n_cells_clust']
    ### prop cells cluster-pat
    df['n_cells_clust-pat']=df.groupby(['subclust', 'patient_number'])['batch'].transform('count')
    df['prop_cells_clust-pat']=df['n_cells_clust-pat']/df['n_cells_clust']
    ## Plot: prop per cluster-patient
    plotTab = df[['subclust', 'patient_number','prop_cells_clust-pat']].copy().drop_duplicates()
    _=plt.figure()
    _=sns.barplot(plotTab, x='subclust', y='prop_cells_clust-pat', hue='patient_number', palette='viridis', legend=False)
    plt.savefig(path_folder+plot_suffix+"_propCells_01.pdf", bbox_inches="tight")
    ## Plot: prop per cluster-batch
    plotTab = df[['subclust', 'batch','prop_cells_clust-batch']].copy().drop_duplicates()
    _=plt.figure()
    _=sns.barplot(plotTab, x='subclust', y='prop_cells_clust-batch', hue='batch')
    plt.savefig(path_folder+plot_suffix+"_propCells_02.pdf", bbox_inches="tight")

    # QC figures in subclusters
    for i,var in enumerate(['n_genes_by_counts', 'total_counts', 'pct_counts_mt']):
        with mpl.rc_context({'figure.figsize': (6, 6)}):
            _=plt.figure()
            _=sc.pl.violin(adata_subset, keys=var, groupby='subclust', stripplot=True, show=False)
            plt.savefig(path_folder+plot_suffix+f"_violin_0{i}.pdf", bbox_inches="tight")



############################
####### Pseudobulks ########
############################

def which_genes_pass_expr(pdata, celltype_colname='celltype', celltype='B',
grp_colname='condition', min_ct_grp=10, min_total_ct=15, min_prop_grp=0.5):
    """ 
    Determine which genes have sufficient counts to be kept in a statistical analysis.
    Returns a boolean array essentially answering: 
        which gene has counts >= min_ct_grp in at least min_prop_grp of at least one group and total counts >= min_total_ct?
    Inspired by decoupler's filter_by_expr
    Parameters: 
        - pdata: AnnData. pseudobulk data in AnnData format.
        - celltype_colname: str. name of .obs column where celltype categories are stored.
        - celltype: str. celltype category to restrict the filtering.
        - grp_colname: str. name of .obs column where the group categories are stored.
        - min_ct_grp: int. minimum number of counts a gene should have for at least some samples (min_prop_grp) to be retained.
        - min_prop_grp: float. minimum proportion of samples per group that expresses the gene.
    Returns:
        genes: ndarray. Genes to be kept True/False.
    """
    import pandas as pd
    import scanpy as sc

    # dictionary with obs_names for each group (e.g., patients for each condition)
    ## get the name for each group
    grps=pdata.obs[grp_colname].unique()
    ## populate dictionary
    obsnames_by_group={grp:[] for grp in grps}
    for grp in grps:
        mask=(pdata.obs[celltype_colname]==celltype)&(pdata.obs[grp_colname]==grp)
        obs=list(pdata.obs.loc[mask,:].index)
        obsnames_by_group[grp]=obs
    
    # df with raw gene counts for each observation (gene x obs)
    counts_df = pd.DataFrame(pdata.X.copy(), index=list(pdata.obs_names), columns=list(pdata.var_names)).T

    # genes passing expr thresholds
    ## list of matrices, one per each group
    mat_lst=[counts_df.loc[:,counts_df.columns.isin(obsnames_by_group[grp])] for grp in grps]
    ## get the genes that pass the filter for EACH group (boolean)
    ## (which gene has counts >= min_ct_grp in at least min_prop_grp)
    gene_lst=[(mat>=min_ct_grp).sum(axis=1)>=round(mat.shape[1]*min_prop_grp) for mat in mat_lst]
    ## get the genes that pass the filter for AT LEAST ONE group 
    genes=pd.DataFrame(gene_lst).T.sum(axis=1) >= 1
    ## and also have total counts >= min_total_ct
    genes=genes.mul(counts_df.sum(axis=1)>=min_total_ct)

    return pd.DataFrame(genes, columns=['pass_expr'])


############################
########### DEA ############
############################

def run_mixedlm(df, var_y_lst=[], var_x_lst=[], formula_rightSide = None, grps=''):
    """ 
    Runs a linear mixed model for each dependent variable. 
    var_x_lst is ignored if a formula expressing the independent variables is provided (formula_rightSide is not None).
    This is useful when you want to analyze interaction effects or more complex designs.
    Returns a list with two dfs: results and residuals.
    If an independent variable is categorical, it assumes it is already type:<category> and ordered 
    such that the resulting coeficient is relative to the baseline.

    Parameters:     df: DataFrame
                        Wide-format dataframe with the variables to be included in the ols as columns. 
                    var_y_lst: list
                        list with dependent variables. One model will be fitted for each element of the list.
                    var_x_lst: list
                        list with independent variables used to explain each variable in var_y_lst. 
                    formula_rightSide: str
                        string with the formula for the independent variables only.
                    grps: str
                        indicates the random effect
    Returns:        list with two elements:
                        - DataFrame with results.
                        - DataFrame with residuals.
    """
    import pandas as pd
    import numpy as np
    import re
    import statsmodels.formula.api as smf
    ## silence warnings
    import warnings
    from statsmodels.tools.sm_exceptions import ConvergenceWarning
    warnings.simplefilter('ignore', ConvergenceWarning)

    if formula_rightSide:
        formula_right=formula_rightSide
    else:
        formula_right=re.sub(pattern='[\'\"\[\]]', repl='', string=str(var_x_lst)).replace(', ', ' + ')
    ## lst of dfs with results
    results_lst = [] 
    residuals_lst = []
    
    for var in var_y_lst:
        ## define formula 
        formula=f"{var} ~ {formula_right}"
        try:
            ## fit model
            model = smf.mixedlm(formula, groups = grps, data = df).fit()
            ## retrieve results
            ### regression 
            res=model.summary().tables[1].reset_index()
            res=res.iloc[1:-1,:res.shape[1]-2] # get rid of the intercept and random effect
            res.columns = ['regressor', 'coef', 'stdErr', 'z', 'p']
            res['p']= list(model.pvalues)[1:-1] # these p values have better precision, the others in the summary table have few decimal positions.
            ### residuals
            df_residuals=pd.DataFrame(model.resid, columns=['resid'])
            df_residuals['dependent']=var
            df_residuals['formula']=formula; df_residuals['randomEff'] = grps
        except:
            ## empty results
            res=pd.DataFrame({colname:[np.nan] for colname in ['regressor', 'coef', 'stdErr', 'z', 'p']})
            df_residuals=pd.DataFrame()
        
        res['formula'] = formula; res['randomEff'] = grps; res['dependent']=var
    
        results_lst.append(res)
        residuals_lst.append(df_residuals)
    
    return([pd.concat(results_lst), pd.concat(residuals_lst)])


def run_lm(df, var_y_lst=[], var_x_lst=[], formula_rightSide = None, robust=False):
    """ 
    Runs ols regression for each dependent variable. If robust=True, it runs robust linear regression.
    var_x_lst is ignored if a formula expressing the independent variables is provided (formula_rightSide is not None).
    This is useful when you want to analyze interaction effects or more complex designs.
    Returns a list with two dfs: ols results and residuals.
    If an independent variable is categorical, it assumes it is already type:<category> and ordered 
    such that the resulting coeficient is relative to the baseline.

    Parameters:     df: DataFrame
                        Wide-format dataframe with the variables to be included in the ols as columns. 
                    var_y_lst: list
                        list with dependent variables. One ols model will be fitted for each element of the list.
                    var_x_lst: list
                        list with independent variables used to explain each variable in var_y_lst. 
                    formula_rightSide: str
                        string with the formula for the independent variables only.
    Returns:        list with two elements:
                        - DataFrame with ols results.
                        - DataFrame with residuals.
    """
    import pandas as pd
    import re
    import statsmodels.formula.api as smf
    ## silence warnings
    import warnings
    from statsmodels.tools.sm_exceptions import ConvergenceWarning
    warnings.simplefilter('ignore', ConvergenceWarning)

    if formula_rightSide:
        formula_right=formula_rightSide
    else:
        formula_right=re.sub(pattern='[\'\"\[\]]', repl='', string=str(var_x_lst)).replace(', ', ' + ')
    ## lst of dfs with results
    results_lst = [] 
    residuals_lst = []
    
    for var in var_y_lst:
        ## define formula and fit model
        formula=f"{var} ~ {formula_right}"
        if robust:
            model = smf.rlm(formula, data = df).fit()
        else:
            model = smf.ols(formula, data = df).fit()
        ## retrieve results
        ### regression 
        res=model.summary2().tables[1].reset_index()
        res=res.iloc[1:,:res.shape[1]-2] # get rid of the intercept
        if robust:
            res.columns = ['regressor', 'coef', 'stdErr', 'z', 'p']
        else: 
            res.columns = ['regressor', 'coef', 'stdErr', 't', 'p']
            res['R2adj']=model.rsquared_adj
        res['p']=list(model.pvalues[model.pvalues.index!='Intercept']) # these p values have better precision, the others in the summary table have few decimal positions.
        res['formula'] = formula; res['dependent']=var
        ### residuals
        df_residuals=pd.DataFrame(model.resid, columns=['resid'])
        df_residuals['dependent']=var
        df_residuals['formula']=formula
    
        results_lst.append(res)
        residuals_lst.append(df_residuals)
    
    return([pd.concat(results_lst), pd.concat(residuals_lst)])

def run_lrt(df, var_y_lst=[], var_x_lst=[]):
    """ 
    Runs likelihood ratio test (LRT) between two ols models for each dependent variable. 
    Returns a df with LRT stat, pvalue, and design of the complex and simpler models.

    Parameters:     df: DataFrame
                        Wide-format dataframe with the variables to be included in the ols model as columns. 
                    var_y_lst: list
                        list with dependent variables. One ols model will be fitted for each element of the list.
                    var_x_lst: list
                        list with independent variables used to explain each variable in var_y_lst. 
                        The LRT will be between a model with all variables vs a model without the last variable
    Returns:        DataFrame:
                        Lrt results: lr stat, pvalue, complex and simpler models.
    """
    import pandas as pd
    import re
    import statsmodels.formula.api as smf
    ## silence warnings
    import warnings
    from statsmodels.tools.sm_exceptions import ConvergenceWarning
    warnings.simplefilter('ignore', ConvergenceWarning)

    # check input
    assert (type(var_x_lst)==type([]))&((type(var_x_lst)==type([]))), 'var_x_lst and var_x_lst should be lists'  
    assert len(var_y_lst)>=1, 'var_y_lst should contain at least one element'
    assert len(var_x_lst)>=2, 'var_x_lst should contain at least two elements'

    formula_right=re.sub(pattern='[\'\"\[\]]', repl='', string=str(var_x_lst)).replace(', ', ' + ')
    formula_right_simple=re.sub(pattern='[\'\"\[\]]', repl='', string=str(var_x_lst[:-1])).replace(', ', ' + ')
    ## lst of dfs with results
    results_lst = [] 
    
    for var in var_y_lst:
        ## define formula and fit model
        formula=f"{var} ~ {formula_right}"
        formula_simple=f"{var} ~ {formula_right_simple}"
        model = smf.ols(formula, data = df).fit()
        model_simple = smf.ols(formula_simple, data = df).fit()
        ## retrieve results
        stat, p, _ = model.compare_lr_test(restricted=model_simple)
        res=pd.DataFrame({'stat': [stat], 'p': [p]})
        res['dependent']=var; res['formula']=formula; res['formula_simple']=formula_simple

        results_lst.append(res)
    
    return(pd.concat(results_lst))


##########################################
######## Comparing ranked vectors ########
##########################################

# determine the weight parameter (p) of the Rank biased overlap (RBO)
# Modified from: https://towardsdatascience.com/rbo-v-s-kendall-tau-to-compare-ranked-lists-of-items-8776c5182899
# credit: Preetam Joshi
## The res

def wRBO(p, depth):
    """ 
    This function helps one to set the p parameter of the Rank Biased Overlap (RBO) measure.
    The RBO measure is used to compare the similarity of two ranked lists. In ranked lists, the 
    first elements are more important to the last ones (e.g., DEGs in DEA). The p parameter controls
    how important the first n elements are in the final score.
    The function requires two arguments: p and depth. The resulting number indicates the fraction of 
    the weight at a given depth (e.g., the fraction of the weight that will be given to the first n elements of the list).
    For example, if p=0.9937, then 20 % of the weight will be given to the first 10 elements 
    (by running: wRBO(p=0.9937, depth=10)), 54 % of the weight to the top 50 elements, 
    73 % to the first 100 and 99% to the first 500.
    """
    import math
    import numpy as np
    # tail recursive helper function
    def helper(ret, p, depth, i):
        term = math.pow(p, i)/i
        if depth == i:
            return ret + term
        return helper(ret + term, p, depth, i+1)
    def sum_series(p, depth):
        return helper(0, p, depth, 1)

    wrbo1_d = 1 - math.pow(p, depth-1) + (((1-p)/p) * depth *(np.log(1/(1-p)) - sum_series(p, depth-1)))
    return wrbo1_d


# Modified from: https://github.com/dlukes/rbo/blob/master/rbo.py
# credit: 2016 ÚČNK/David Lukeš

"""Rank-biased overlap, a ragged sorted list similarity measure.

See http://doi.acm.org/10.1145/1852102.1852106 for details. All functions
directly corresponding to concepts from the paper are named so that they can be
clearly cross-identified.

The definition of overlap has been modified to account for ties. Without this,
results for lists with tied items were being inflated. The modification itself
is not mentioned in the paper but seems to be reasonable, see function
``overlap()``. Places in the code which diverge from the spec in the paper
because of this are highlighted with comments.

The two main functions for performing an RBO analysis are ``rbo()`` and
``rbo_dict()``; see their respective docstrings for how to use them.

"""

RBO = namedtuple("RBO", "min res ext")
RBO.__doc__ += ": Result of full RBO analysis"
RBO.min.__doc__ = "Lower bound estimate"
RBO.res.__doc__ = "Residual corresponding to min; min + res is an upper bound estimate"
RBO.ext.__doc__ = "Extrapolated point estimate"


def _round(obj):
    if isinstance(obj, RBO):
        return RBO(_round(obj.min), _round(obj.res), _round(obj.ext))
    else:
        return round(obj, 3)


def set_at_depth(lst, depth):
    ans = set()
    for v in lst[:depth]:
        if isinstance(v, set):
            ans.update(v)
        else:
            ans.add(v)
    return ans


def raw_overlap(list1, list2, depth):
    """Overlap as defined in the article.

    """
    set1, set2 = set_at_depth(list1, depth), set_at_depth(list2, depth)
    return len(set1.intersection(set2)), len(set1), len(set2)


def overlap(list1, list2, depth):
    """Overlap which accounts for possible ties.

    This isn't mentioned in the paper but should be used in the ``rbo*()``
    functions below, otherwise overlap at a given depth might be > depth which
    inflates the result.

    There are no guidelines in the paper as to what's a good way to calculate
    this, but a good guess is agreement scaled by the minimum between the
    requested depth and the lengths of the considered lists (overlap shouldn't
    be larger than the number of ranks in the shorter list, otherwise results
    are conspicuously wrong when the lists are of unequal lengths -- rbo_ext is
    not between rbo_min and rbo_min + rbo_res.

    >>> overlap("abcd", "abcd", 3)
    3.0

    >>> overlap("abcd", "abcd", 5)
    4.0

    >>> overlap(["a", {"b", "c"}, "d"], ["a", {"b", "c"}, "d"], 2)
    2.0

    >>> overlap(["a", {"b", "c"}, "d"], ["a", {"b", "c"}, "d"], 3)
    3.0

    """
    return agreement(list1, list2, depth) * min(depth, len(list1), len(list2))
    # NOTE: comment the preceding and uncomment the following line if you want
    # to stick to the algorithm as defined by the paper
    # return raw_overlap(list1, list2, depth)[0]


def agreement(list1, list2, depth):
    """Proportion of shared values between two sorted lists at given depth.

    >>> _round(agreement("abcde", "abdcf", 1))
    1.0
    >>> _round(agreement("abcde", "abdcf", 3))
    0.667
    >>> _round(agreement("abcde", "abdcf", 4))
    1.0
    >>> _round(agreement("abcde", "abdcf", 5))
    0.8
    >>> _round(agreement([{1, 2}, 3], [1, {2, 3}], 1))
    0.667
    >>> _round(agreement([{1, 2}, 3], [1, {2, 3}], 2))
    1.0

    """
    len_intersection, len_set1, len_set2 = raw_overlap(list1, list2, depth)
    return 2 * len_intersection / (len_set1 + len_set2)


def cumulative_agreement(list1, list2, depth):
    return (agreement(list1, list2, d) for d in range(1, depth + 1))


def average_overlap(list1, list2, depth=None):
    """Calculate average overlap between ``list1`` and ``list2``.

    >>> _round(average_overlap("abcdefg", "zcavwxy", 1))
    0.0
    >>> _round(average_overlap("abcdefg", "zcavwxy", 2))
    0.0
    >>> _round(average_overlap("abcdefg", "zcavwxy", 3))
    0.222
    >>> _round(average_overlap("abcdefg", "zcavwxy", 4))
    0.292
    >>> _round(average_overlap("abcdefg", "zcavwxy", 5))
    0.313
    >>> _round(average_overlap("abcdefg", "zcavwxy", 6))
    0.317
    >>> _round(average_overlap("abcdefg", "zcavwxy", 7))
    0.312

    """
    depth = min(len(list1), len(list2)) if depth is None else depth
    return sum(cumulative_agreement(list1, list2, depth)) / depth


def rbo_at_k(list1, list2, p, depth=None):
    # ``p**d`` here instead of ``p**(d - 1)`` because enumerate starts at
    # 0
    depth = min(len(list1), len(list2)) if depth is None else depth
    d_a = enumerate(cumulative_agreement(list1, list2, depth))
    return (1 - p) * sum(p ** d * a for (d, a) in d_a)


def rbo_min(list1, list2, p, depth=None):
    """Tight lower bound on RBO.

    See equation (11) in paper.

    >>> _round(rbo_min("abcdefg", "abcdefg", .9))
    0.767
    >>> _round(rbo_min("abcdefgh", "abcdefg", .9))
    0.767

    """
    depth = min(len(list1), len(list2)) if depth is None else depth
    x_k = overlap(list1, list2, depth)
    log_term = x_k * math.log(1 - p)
    sum_term = sum(
        p ** d / d * (overlap(list1, list2, d) - x_k) for d in range(1, depth + 1)
    )
    return (1 - p) / p * (sum_term - log_term)


def rbo_res(list1, list2, p):
    """Upper bound on residual overlap beyond evaluated depth.

    See equation (30) in paper.

    NOTE: The doctests weren't verified against manual computations but seem
    plausible. In particular, for identical lists, ``rbo_min()`` and
    ``rbo_res()`` should add up to 1, which is the case.

    >>> _round(rbo_res("abcdefg", "abcdefg", .9))
    0.233
    >>> _round(rbo_res("abcdefg", "abcdefghijklmnopqrstuvwxyz", .9))
    0.239

    """
    S, L = sorted((list1, list2), key=len)
    s, l = len(S), len(L)
    x_l = overlap(list1, list2, l)
    # since overlap(...) can be fractional in the general case of ties and f
    # must be an integer --> math.ceil()
    f = int(math.ceil(l + s - x_l))
    # upper bound of range() is non-inclusive, therefore + 1 is needed
    term1 = s * sum(p ** d / d for d in range(s + 1, f + 1))
    term2 = l * sum(p ** d / d for d in range(l + 1, f + 1))
    term3 = x_l * (math.log(1 / (1 - p)) - sum(p ** d / d for d in range(1, f + 1)))
    return p ** s + p ** l - p ** f - (1 - p) / p * (term1 + term2 + term3)


def rbo_ext(list1, list2, p):
    """RBO point estimate based on extrapolating observed overlap.

    See equation (32) in paper.

    NOTE: The doctests weren't verified against manual computations but seem
    plausible.

    >>> _round(rbo_ext("abcdefg", "abcdefg", .9))
    1.0
    >>> _round(rbo_ext("abcdefg", "bacdefg", .9))
    0.9

    """
    S, L = sorted((list1, list2), key=len)
    s, l = len(S), len(L)
    x_l = overlap(list1, list2, l)
    x_s = overlap(list1, list2, s)
    # the paper says overlap(..., d) / d, but it should be replaced by
    # agreement(..., d) defined as per equation (28) so that ties are handled
    # properly (otherwise values > 1 will be returned)
    # sum1 = sum(p**d * overlap(list1, list2, d)[0] / d for d in range(1, l + 1))
    sum1 = sum(p ** d * agreement(list1, list2, d) for d in range(1, l + 1))
    sum2 = sum(p ** d * x_s * (d - s) / s / d for d in range(s + 1, l + 1))
    term1 = (1 - p) / p * (sum1 + sum2)
    term2 = p ** l * ((x_l - x_s) / l + x_s / s)
    return term1 + term2


def rbo(list1, list2, p):
    """Complete RBO analysis (lower bound, residual, point estimate).

    ``list`` arguments should be already correctly sorted iterables and each
    item should either be an atomic value or a set of values tied for that
    rank. ``p`` is the probability of looking for overlap at rank k + 1 after
    having examined rank k.

    >>> lst1 = [{"c", "a"}, "b", "d"]
    >>> lst2 = ["a", {"c", "b"}, "d"]
    >>> _round(rbo(lst1, lst2, p=.9))
    RBO(min=0.489, res=0.477, ext=0.967)

    """
    if not 0 <= p <= 1:
        raise ValueError("The ``p`` parameter must be between 0 and 1.")
    args = (list1, list2, p)
    return RBO(rbo_min(*args), rbo_res(*args), rbo_ext(*args))


def sort_dict(dct, *, ascending=False):
    """Sort keys in ``dct`` according to their corresponding values.

    Sorts in descending order by default, because the values are
    typically scores, i.e. the higher the better. Specify
    ``ascending=True`` if the values are ranks, or some sort of score
    where lower values are better.

    Ties are handled by creating sets of tied keys at the given position
    in the sorted list.

    >>> dct = dict(a=1, b=2, c=1, d=3)
    >>> list(sort_dict(dct)) == ['d', 'b', {'a', 'c'}]
    True
    >>> list(sort_dict(dct, ascending=True)) == [{'a', 'c'}, 'b', 'd']
    True

    """
    scores = []
    items = []
    # items should be unique, scores don't have to
    for item, score in dct.items():
        if not ascending:
            score *= -1
        i = bisect_left(scores, score)
        if i == len(scores):
            scores.append(score)
            items.append(item)
        elif scores[i] == score:
            existing_item = items[i]
            if isinstance(existing_item, set):
                existing_item.add(item)
            else:
                items[i] = {existing_item, item}
        else:
            scores.insert(i, score)
            items.insert(i, item)
    return items


def rbo_dict(dict1, dict2, p, *, sort_ascending=False):
    """Wrapper around ``rbo()`` for dict input.

    Each dict maps items to be sorted to the score according to which
    they should be sorted. The RBO analysis is then performed on the
    resulting sorted lists.

    The sort is descending by default, because scores are typically the
    higher the better, but this can be overridden by specifying
    ``sort_ascending=True``.

    >>> dct1 = dict(a=1, b=2, c=1, d=3)
    >>> dct2 = dict(a=1, b=2, c=2, d=3)
    >>> _round(rbo_dict(dct1, dct2, p=.9, sort_ascending=True))
    RBO(min=0.489, res=0.477, ext=0.967)

    """
    list1, list2 = (
        sort_dict(dict1, ascending=sort_ascending),
        sort_dict(dict2, ascending=sort_ascending),
    )
    return rbo(list1, list2, p)


