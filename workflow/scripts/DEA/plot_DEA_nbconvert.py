import numpy as np
import pandas as pd
import scanpy as sc
import decoupler as dc
import sys
import math
from pathlib import Path
from scipy import sparse

# When this notebook is exported and run as a script, add the project root to sys.path
# so local workflow.* imports behave like they do in the notebook kernel.
if "__file__" in globals():
    PROJECT_ROOT = Path(__file__).resolve().parents[3]
    if str(PROJECT_ROOT) not in sys.path:
        sys.path.insert(0, str(PROJECT_ROOT))

## helper functions
from workflow.scripts.funs.helper_funs import get_project_root

## plotting
# import marsilea as ma
# import marsilea.plotter as mp
# from matplotlib.colors import NoNorm # for marsilea dotplot
import plotnine as p9
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl

from pathlib import Path

if "snakemake" in globals():
    # in_metadata = Path(snakemake.input["in_metadata"])
    # out_metadata = Path(snakemake.output["out_metadata"])
    print('Hello World')
else:
    ## get paths
    ROOT = get_project_root(globals())
    ## input
    in_adata = ROOT / "output/preprocessing/objects/QC/adata.h5ad"
    in_deg = ROOT / "output/DEA/objects/calibrate_pvalues/deg_calibrated.tsv"
    ## output
    ### plots
    out_plots = ROOT / "output/DEA/plots/plot_DEA"
    out_plots.mkdir(parents=True, exist_ok=True)


# adata

adata = sc.read_h5ad(in_adata)
## remove group B
adata = adata[adata.obs['QC_group']=='A'].copy()
## make sure counts are norm
adata.X = adata.layers['cts.deseq2_norm'].copy()
## col with ENSEMBL ID
adata.var['ENS']=adata.var.index
## rename var_names
adata.var_names = [symb if (not pd.isna(symb)) else ENS for ENS,symb in zip(adata.var.index, adata.var['genesymbol'])]
## because some genesymbol won't be unique
adata.var_names_make_unique()
## conv dictionary to re-name degs
conv_dict = {ENS:new_symb for ENS,new_symb in zip(adata.var['ENS'],adata.var.index)}

adata

# add individualID

col_dict = {
    'SD1': 'individual_01','SD2': 'individual_01','SD3': 'individual_02','SD4': 'individual_02',
    'SD5': 'individual_03','SD6': 'individual_03','SD7': 'individual_04','SD8': 'individual_04',
    'SD9': 'individual_05','SD10': 'individual_05','SD11': 'individual_06','SD12': 'individual_06',
    'SD13': 'individual_07','SD14': 'individual_07','SD15': 'individual_08','SD16': 'individual_08',
    'SD17': 'individual_09','SD18': 'individual_09','SD19': 'individual_10','SD20': 'individual_10',
    'SD21': 'individual_11','SD22': 'individual_11','SD23': 'individual_12','SD24': 'individual_12',
    'SD25': 'individual_13','SD26': 'individual_13','SD27': 'individual_14','SD28': 'individual_14',
    'SD29': 'individual_15','SD30': 'individual_15','SD31': 'individual_16','SD32': 'individual_16',
    'SD33': 'individual_17','SD34': 'individual_17','SD35': 'individual_18','SD36': 'individual_18',
    'SD37': 'individual_19','SD38': 'individual_19','SD39': 'individual_20','SD40': 'individual_20',
    'SD41': 'individual_21','SD42': 'individual_21','SD43': 'individual_22','SD44': 'individual_22',
    'SD45': 'individual_23','SD46': 'individual_23','SD47': 'individual_24','SD48': 'individual_24',
    'SD49': 'individual_25','SD50': 'individual_25','SD51': 'individual_26','SD52': 'individual_26',
    'SD53': 'individual_27','SD54': 'individual_27','SD55': 'individual_28','SD56': 'individual_28',
    'SD57': 'individual_29','SD58': 'individual_29','SD59': 'individual_30','SD60': 'individual_30',
    'SD61': 'individual_31','SD62': 'individual_31','SD63': 'individual_32','SD64': 'individual_32',
    'SD65': 'individual_33','SD66': 'individual_33',
}
adata.obs['individualID']=adata.obs.index.map(col_dict)
adata.obs['individualID']=adata.obs['individualID'].astype('category')

# DEA results

deg = pd.read_csv(in_deg,sep='\t',index_col=0).set_index('gene_id')
## rename genesymbol as in adata.var
deg['genesymbol'] = deg.index.map(conv_dict)

deg.head()

# QC plot: volcano plot with decoupler

volcano_tab = deg.dropna(subset=["log2FoldChange", "padj"]).copy()
volcano_tab = volcano_tab.set_index("genesymbol", drop=False)

fig = dc.pl.volcano(
    data=volcano_tab,
    x="log2FoldChange",
    y="padj",
    top=20,
    thr_stat=0.5,
    thr_sign=0.05,
    figsize=(7, 6),
    dpi=120,
    return_fig=True,
)
fig.axes[0].set_title("HR_case__IFNG_pos vs HR_control__IFNG_pos")
fig.savefig(out_plots / "volcano.png", bbox_inches="tight")
fig.savefig(out_plots / "volcano.pdf", bbox_inches="tight")

plt.close(fig)


topn=20
## subset adata to samples in DEA
select_grp=['HR_control__IFNG_pos','HR_case__IFNG_pos']
tmp=adata[adata.obs['grp'].isin(select_grp)].copy()
tmp.obs['grp']=tmp.obs['grp'].astype('category').cat.reorder_categories(select_grp,ordered=True)
## palette
pal = [c for c in sns.color_palette('Oranges',n_colors=2).as_hex()]
tmp.uns['grp_colors']=pal
## get topn genes in comparison
keys=deg.sort_values('padj',ascending=True)
keys=keys[keys['padj']<0.05]
keys=keys.head(topn)['genesymbol'].tolist()

## get table with norm counts and metadata in long format
ptab = pd.DataFrame(tmp.layers['cts.deseq2_norm'], index=tmp.obs.index,columns=tmp.var.index)
ptab = (
    ptab[keys].reset_index().
    melt(id_vars=['sample_id'],var_name='genesymbol',value_name='expr'))
### add grp from adata
ptab['grp']=ptab['sample_id'].map( tmp.obs['grp'].to_dict() )
ptab['grp']=['HR case, IFNG+' if grp=='HR_case__IFNG_pos' else 'HR ctrl, IFNG+' for grp in ptab['grp']]
ptab['grp']=ptab['grp'].astype('category').cat.reorder_categories(['HR ctrl, IFNG+','HR case, IFNG+'],ordered=True)
#ptab['grp']=ptab['grp'].astype('category').cat.reorder_categories(tmp.obs['grp'].cat.categories,ordered=True)
## add stat and padj
cols=['stat','padj']
for col in cols:
    col_dict={x:y for x,y in zip(deg['genesymbol'],deg[col])}
    ptab[col]=ptab['genesymbol'].map(col_dict)
### plot title (ordered by significance)
#ptab=ptab.sort_values(by='padj',ascending=True)
ptab['padj']=['{:.2e}'.format(x) for x in ptab['padj']]
ptab['p_title']=[f"{a}\nStat: {round(float(b),2)}\nPadj: {c}" for a,b,c in zip(ptab['genesymbol'],ptab['stat'],ptab['padj'])]
ptab['p_title']=ptab['p_title'].astype('category').cat.reorder_categories(ptab['p_title'].unique(), ordered=True)

ptab.head()

# as boxplot


## palette
pal = {grp:color for grp,color in zip(ptab['grp'].cat.categories,sns.color_palette('Oranges',n_colors=2).as_hex())}

plots_per_col=8
ncols=math.ceil(len(keys)/plots_per_col)
w=2.3*plots_per_col if plots_per_col>=4 else 6
h=4*ncols

contrast='grp'
p=(
    p9.ggplot(ptab, p9.aes(x=contrast, y='expr', group=contrast)) +
    p9.geom_boxplot(outlier_alpha=0,size=0.8,width=0.3,show_legend=False) +
    p9.geom_jitter(p9.aes(fill=contrast), color='black', size=3, width=0.1,alpha=1.0,show_legend=True) + 
    p9.scale_fill_manual(pal) +
    p9.guides(fill=p9.guide_legend(override_aes={'size':6})) +
    p9.theme_classic() +
    p9.labs(title="Expression of selected genes") +
    p9.theme(
        plot_title=p9.element_text(hjust=0.5, size=14), 
        strip_text=p9.element_text(size=12), 
        axis_title=p9.element_text(size=14), 
        axis_text_x=p9.element_blank(), 
        legend_text=p9.element_text(size=12),
        figure_size=(w, h),
    ) + 
    p9.facet_wrap('p_title', scales = 'free', ncol=plots_per_col) 
)
p.save(out_plots / "DEA_boxplot_01.png", width=w, height=h, limitsize=False)
p.save(out_plots / "DEA_boxplot_01.pdf", width=w, height=h, limitsize=False)
plt.close()

## subset adata to samples in DEA
tmp=adata.copy() # no change since we plot all samples
## get topn genes in comparison
keys=[
    'APOL1', 'APOL6', 'IRF1', 'CCL2', 'CXCL9', 'CXCL10', 
    'ICAM1', 'IDO1', 'PDCD1LG2', 'GBP1', 'PSMB9', 'CCDC68', 
    'CT5S', 'TMEM140', 'TRIM21', 'STAT1', 'STAT2', 'STAT3', 
    'IRF2', 'IRF8', 'IRF9', 'CXCL11']
keys=[k for k in keys if k in tmp.var_names]

## get table with norm counts and metadata in long format
ptab = pd.DataFrame(tmp.layers['cts.deseq2_norm'], index=tmp.obs.index,columns=tmp.var.index)
ptab = (
    ptab[keys].reset_index().
    melt(id_vars=['sample_id'],var_name='genesymbol',value_name='expr'))
### add grp from adata
ptab['grp']=ptab['sample_id'].map( tmp.obs['grp'].to_dict() )
col_dict={
    'HR_case__IFNG_neg':'HR case, IFNG-', 'HR_case__IFNG_pos':'HR case, IFNG+', 
    'HR_control__IFNG_neg':'HR control, IFNG-','HR_control__IFNG_pos':'HR control, IFNG+'}
ptab['grp']=ptab['grp'].map(col_dict)
ptab['grp']=ptab['grp'].astype('category').cat.reorder_categories(['HR case, IFNG-','HR case, IFNG+','HR control, IFNG-','HR control, IFNG+'],ordered=True)
# ## add stat and padj
# cols=['stat','padj']
# for col in cols:
#     col_dict={x:y for x,y in zip(deg['genesymbol'],deg[col])}
#     ptab[col]=ptab['genesymbol'].map(col_dict)
# ### plot title (ordered by significance)
# ptab['padj']=['{:.2e}'.format(x) for x in ptab['padj']]
# ptab['p_title']=[f"{a}\nStat: {round(float(b),2)}\nPadj: {c}" for a,b,c in zip(ptab['genesymbol'],ptab['stat'],ptab['padj'])]
# ptab['p_title']=ptab['p_title'].astype('category').cat.reorder_categories(ptab['p_title'].unique(), ordered=True)
ptab['p_title']=ptab['genesymbol'].astype('category').cat.reorder_categories(ptab['genesymbol'].unique(), ordered=True)

ptab.head()

# as boxplot


## palette
oranges=[c for c in sns.color_palette('Oranges',n_colors=2).as_hex()] # HR case
blues=[c for c in sns.color_palette('Blues',n_colors=2).as_hex()] # HR control
pal = {grp:color for grp,color in zip(ptab['grp'].cat.categories,oranges+blues)}

plots_per_col=8
ncols=math.ceil(len(keys)/plots_per_col)
w=2.3*plots_per_col if plots_per_col>=4 else 6
h=4*ncols

contrast='grp'
p=(
    p9.ggplot(ptab, p9.aes(x=contrast, y='expr', group=contrast)) +
    p9.geom_boxplot(outlier_alpha=0,size=0.8,width=0.3,show_legend=False) +
    p9.geom_jitter(p9.aes(fill=contrast), color='black', size=3, width=0.1,alpha=1.0,show_legend=True) + 
    p9.scale_fill_manual(pal) +
    p9.guides(fill=p9.guide_legend(override_aes={'size':6})) +
    p9.theme_classic() +
    p9.labs(title="Expression of selected genes") +
    p9.theme(
        plot_title=p9.element_text(hjust=0.5, size=14), 
        strip_text=p9.element_text(size=12), 
        axis_title=p9.element_text(size=14), 
        axis_text_x=p9.element_blank(), 
        legend_text=p9.element_text(size=12),
        figure_size=(w, h),
    ) + 
    p9.facet_wrap('p_title', scales = 'free', ncol=plots_per_col) 
)
p.save(out_plots / "Fig5g_boxplot_01.png", width=w, height=h, limitsize=False)
p.save(out_plots / "Fig5g_boxplot_01.pdf", width=w, height=h, limitsize=False)
plt.close()

## subset adata to samples in DEA
tmp=adata.copy() # no change since we plot all samples
## get topn genes in comparison
keys=[
    'APOL1', 'APOL6', 'IRF1', 'CCL2', 'CXCL9', 'CXCL10', 
    'ICAM1', 'IDO1', 'PDCD1LG2', 'GBP1', 'PSMB9', 'CCDC68', 
    'CT5S', 'TMEM140', 'TRIM21', 'STAT1', 'STAT2', 'STAT3', 
    'IRF2', 'IRF8', 'IRF9', 'CXCL11']
keys=[k for k in keys if k in tmp.var_names]

# get fold change
## by individualID, IFNG+ / IFNG-

## get expression
ptab = pd.DataFrame(tmp.layers['cts.deseq2_norm'], index=tmp.obs.index,columns=tmp.var.index)
ptab = (
    ptab[keys].reset_index().
    melt(id_vars=['sample_id'],var_name='genesymbol',value_name='expr'))
## add individualID
col_dict=tmp.obs['individualID'].to_dict()
ptab['individualID']=ptab['sample_id'].map(col_dict)
## add other columns
cols=['condition_group']
for col in cols:
    col_dict=tmp.obs[col].to_dict()
    ptab[col]=ptab['sample_id'].map(col_dict)
## get fold change
ptab = (
    ptab
    .groupby(["individualID", "genesymbol", "condition_group"], observed=True)["expr"]
    .mean()
    .reset_index()
    .pivot(index=["individualID", "genesymbol"], columns="condition_group", values="expr")
    .reset_index()
)
eps = 1e-6
ptab["fold_change"] = (ptab["IFNG_pos"] + eps) / (ptab["IFNG_neg"] + eps)
ptab["log2_fold_change"] = np.log2(ptab["fold_change"])

ptab.head()

# prepare table for plotting

### add disease_group
col_dict=tmp.obs.set_index('individualID')['disease_group'].to_dict()
ptab['grp']=ptab['individualID'].map( col_dict )
col_dict={'HR_case':'HR case, IFNG+ vs IFNG-', 'HR_control':'HR control, IFNG+ vs IFNG-'}
ptab['grp']=ptab['grp'].map(col_dict)
ptab['grp']=ptab['grp'].astype('category').cat.reorder_categories(['HR case, IFNG+ vs IFNG-', 'HR control, IFNG+ vs IFNG-'],ordered=True)
## add stat and padj
cols=['stat','padj']
for col in cols:
    col_dict={x:y for x,y in zip(deg['genesymbol'],deg[col])}
    ptab[col]=ptab['genesymbol'].map(col_dict)
# ### plot title (ordered by significance)
# ptab['padj']=['{:.2e}'.format(x) for x in ptab['padj']]
# ptab['p_title']=[f"{a}\nStat: {round(float(b),2)}\nPadj: {c}" for a,b,c in zip(ptab['genesymbol'],ptab['stat'],ptab['padj'])]
# ptab['p_title']=ptab['p_title'].astype('category').cat.reorder_categories(ptab['p_title'].unique(), ordered=True)
ptab['p_title']=ptab['genesymbol'].astype('category').cat.reorder_categories(ptab['genesymbol'].unique(), ordered=True)

ptab.head()

# as boxplot


## palette
pal = [c for c in sns.color_palette('tab10',n_colors=2).as_hex()] # blue, orange
pal = pal[::-1]
pal = {cat:color for cat,color in zip(ptab['grp'].cat.categories,pal)}

plots_per_col=8
ncols=math.ceil(len(keys)/plots_per_col)
w=2.3*plots_per_col if plots_per_col>=4 else 6
h=4*ncols

contrast='grp'
p=(
    p9.ggplot(ptab, p9.aes(x=contrast, y='log2_fold_change', group=contrast)) +
    p9.geom_boxplot(outlier_alpha=0,size=0.8,width=0.3,show_legend=False) +
    p9.geom_jitter(p9.aes(fill=contrast), color='black', size=3, width=0.1,alpha=1.0,show_legend=True) + 
    p9.scale_fill_manual(pal) +
    p9.guides(fill=p9.guide_legend(override_aes={'size':6})) +
    p9.theme_classic() +
    p9.labs(title="Log2 fold change of selected genes",y='Log2FC',x='disease group') +
    p9.theme(
        plot_title=p9.element_text(hjust=0.5, size=14), 
        strip_text=p9.element_text(size=12), 
        axis_title=p9.element_text(size=14), 
        axis_text_x=p9.element_blank(), 
        legend_text=p9.element_text(size=12),
        figure_size=(w, h),
    ) + 
    p9.facet_wrap('p_title', scales = 'free', ncol=plots_per_col) 
)
p.save(out_plots / "Fig5g_boxplot_02.png", width=w, height=h, limitsize=False)
p.save(out_plots / "Fig5g_boxplot_02.pdf", width=w, height=h, limitsize=False)
plt.close()
