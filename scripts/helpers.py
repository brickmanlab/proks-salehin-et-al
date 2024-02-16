from typing import Dict, List, Optional, Tuple

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from scvi.model import SCANVI

###################################################################################################
## Core helpers


def get_ensembl_genes(organism: str, save: bool = True):
    """
    Get gene ids and symbols from Ensembl databse.

    organism
        examples: hsapiens, mmusculus
    """
    from pybiomart import Server

    server = Server(host="http://www.ensembl.org")
    print(server.marts["ENSEMBL_MART_ENSEMBL"])

    dataset = server.marts["ENSEMBL_MART_ENSEMBL"].datasets[f"{organism}_gene_ensembl"]
    genes = dataset.query(
        attributes=["ensembl_gene_id", "external_gene_name"]
    ).set_index("Gene stable ID")
    genes.loc[genes["Gene name"].isna(), "Gene name"] = genes[
        genes["Gene name"].isna()
    ].index
    genes = genes.drop_duplicates()

    if save:
        genes.to_csv(f"../data/external/{organism}_gene_ensembl.csv")

    return genes


def get_gene_length(genes: pd.DataFrame, filename: str):
    """
    Create dataframe with gene symbols as index and average length per gene.

    genes:
    filename
        output from gtftools
    """
    gtf = pd.read_table(filename, index_col=0)
    common_genes = genes.index.intersection(gtf.index)
    gtf = gtf.loc[common_genes].dropna()
    return pd.DataFrame(
        gtf.loc[common_genes, "mean"].values,
        index=genes.loc[common_genes, "Gene name"],
        columns=["length"],
    )


def convert_to_gene_symbols(
    adata: sc.AnnData, conversion: pd.DataFrame, upper: bool = False
) -> None:
    conversion.index = conversion.index.str.lower()
    adata.var_names = adata.var_names.str.lower()

    common_genes = adata.var_names.intersection(conversion.index)
    adata.var["symbol"] = adata.var_names
    adata.var.loc[common_genes, "symbol"] = conversion.loc[common_genes, "Gene name"]
    adata.var["symbol"] = adata.var["symbol"].astype(str)
    adata.var.loc[adata.var.symbol == "nan", "symbol"] = adata.var[
        adata.var.symbol == "nan"
    ].index
    adata.var = adata.var.set_index("symbol")
    adata.var_names_make_unique()

    if upper:
        adata.var_names = adata.var_names.str.upper()


def filter_markers(df: pd.DataFrame, n_genes: int = 5, upper: bool = False):
    # significant only
    df = df[
        (df["is_de_fdr_0.05"])
        & (df["bayes_factor"] > 3)
        & (df["non_zeros_proportion1"] > 0.1)
        & (df["lfc_median"] > 0)
    ]
    comparisons = df.comparison.unique()

    deg_df = {}
    for comparison in comparisons:
        cluster = comparison.split(" ")[0]
        markers = (
            df.query("comparison == @comparison")
            .sort_values(by="lfc_median", ascending=False)
            .head(n_genes)
        )
        deg_df[cluster] = (
            markers.index.str.upper().tolist() if upper else markers.index.tolist()
        )

    return deg_df


###################################################################################################
## Plotting helpers


def get_summary(adata) -> None:
    from IPython.display import Markdown, display

    adata.var["mt"] = adata.var_names.str.startswith("mt-")
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
    )
    # sc.pp.calculate_qc_metrics(adata, inplace=True)

    print(f"Number of MITO genes: {np.sum(adata.var['mt'])}")
    print(f"Number of ERCC genes: {np.sum(adata.var_names.str.startswith('ercc-'))}")
    print(
        f"Number of RIBO genes: {np.sum(adata.var_names.str.startswith(('rps','rpl')))}"
    )

    display(Markdown("### Number of cells per experiment"))
    fig, ax = plt.subplots(1, 2, figsize=[20, 5])
    sns.barplot(
        x="experiment",
        y="index",
        data=adata.obs.experiment.value_counts().reset_index(),
        ax=ax[0],
    )
    sc.pl.highest_expr_genes(adata, n_top=20, ax=ax[1])
    plt.show()

    display(Markdown("### Number of cells cell type"))
    plt.figure(figsize=[25, 4])
    counts_stats = (
        # adata.obs.groupby(["experiment", "ct"]).count()["batch"].unstack(fill_value=0)
        adata.obs.groupby(["experiment", "ct"]).size().unstack(fill_value=0)

    )
    sns.heatmap(counts_stats, annot=True, cmap="viridis", fmt="g")
    plt.show()

    display(Markdown("### QC plots"))
    fig, ax = plt.subplots(1, 2, figsize=[25, 4])
    sns.histplot(x="n_genes_by_counts", hue="experiment", data=adata.obs, ax=ax[0])
    sns.histplot(x="total_counts", hue="experiment", data=adata.obs, ax=ax[1])
    ax[1].set_yscale("log")
    plt.show()

    fig, ax = plt.subplots(1, 2, figsize=(25, 4))
    sns.kdeplot(x="total_counts", hue="experiment", data=adata.obs, ax=ax[0])
    ax[0].set_xscale("log")
    plt.show()

    sns.histplot(x="total_counts", hue="experiment", data=adata.obs, ax=ax[1])
    # ax[1].set_xscale('log')
    plt.show()

    fig, ax = plt.subplots(1, 2, figsize=(25, 4))
    sns.violinplot(y=adata.obs["n_genes_by_counts"], orient="v", ax=ax[0])
    sns.violinplot(y=adata.obs["total_counts"], orient="v", ax=ax[1])
    plt.show()

    fig, ax = plt.subplots(1, 2, figsize=(25, 4))
    sns.violinplot(
        x="experiment", y="n_genes_by_counts", orient="v", data=adata.obs, ax=ax[0]
    ).tick_params(axis="x", rotation=90)
    sns.violinplot(
        x="experiment", y="total_counts", orient="v", data=adata.obs, ax=ax[1]
    ).tick_params(axis="x", rotation=90)
    plt.show()

    fig, ax = plt.subplots(figsize=(20, 5), sharey=True)
    sns.scatterplot(
        x="total_counts", y="n_genes_by_counts", ax=ax, data=adata.obs, hue="ct"
    )
    plt.show()


def plt_boxplot(
    adata: sc.AnnData,
    groupby: str,
    gene: str,
    layer: str,
    order: Optional[List[str]] = None,
):
    if gene not in adata.var_names:
        raise ValueError(f"Gene `{gene}` not present")

    tmp = sc.AnnData(adata.layers[layer], obs=adata.obs, var=adata.var)

    groups = tmp.obs.groupby(groupby).groups
    order = order if order else tmp.obs[groupby].cat.categories
    expr_df = {key: np.ravel(tmp[groups[key].values, gene].X) for key in order}
    expr_df = pd.DataFrame({key: pd.Series(v) for key, v in expr_df.items()})
    expr_df = pd.melt(expr_df)

    sns.boxplot(x="variable", y="value", data=expr_df, orient="v")
    plt.title(gene)


def plt_pseudotime(
    adata: sc.AnnData,
    genes: List[str],
    groupby: str,
    layer: str = "latent_normalized",
    key="dpt_pseudotime",
):
    tmp = adata.copy()
    tmp.X = adata.layers[layer]
    tmp = tmp[tmp.obs.sort_values(by=key).index, genes]

    tmp_df = tmp.to_df()
    tmp_df.index = tmp.obs[key].values

    lut = dict(zip(tmp.obs[groupby].cat.categories, tmp.uns[f"{groupby}_colors"]))
    norm = matplotlib.colors.Normalize(vmin=tmp.obs[key].min(), vmax=tmp.obs[key].max())
    sm = plt.cm.ScalarMappable(cmap="viridis", norm=norm)
    print(sm)

    g = sns.clustermap(
        tmp_df.T,
        row_cluster=False,
        col_cluster=False,
        col_colors=[
            tmp.obs[groupby].map(lut).tolist(),
            tmp.obs[groupby].map(lut).tolist(),
        ],
        xticklabels=False,
        dendrogram_ratio=0.01,
        standard_scale="row",
        figsize=(20, 5),
    )
    g.cax.set_visible(False)
    plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0)


def plt_trend(
    adata: sc.AnnData,
    genes: List[str],
    groupby: str,
    layer: str = "latent_normalized",
    key="dpt_pseudotime",
):
    tmp = adata.copy()
    tmp.X = tmp.layers[layer]
    tmp_df = tmp[tmp.obs.sort_values(by=key).index, genes].to_df()

    n = 2
    rm = tmp_df.rolling(n).mean()[n:]
    rm.index = np.arange(rm.shape[0])

    sns.lineplot(data=rm)


def plt_cat_heatmaps(
    df: pd.DataFrame, cell_order: List[str], prefix: str
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    df_heat = pd.pivot_table(df, index=["source", "dest"], values="dist mean").unstack(
        fill_value=0
    )
    df_heat.index = df_heat.index.str.replace(prefix, "")
    df_heat.columns = df_heat.columns.droplevel().str.replace(prefix, "")

    df_heat2 = pd.crosstab(df.source, df.dest)
    df_heat2.index = df_heat2.index.str.replace(prefix, "")
    df_heat2.columns = df_heat2.columns.str.replace(prefix, "")

    common_indexes = [x for x in cell_order if x in df_heat.index]
    common_columns = [x for x in cell_order if x in df_heat.columns]

    return (
        df_heat.loc[common_indexes, common_columns],
        df_heat2.loc[common_indexes, common_columns],
    )


###################################################################################################
## Integration helpers


def compute_scib_metrics(adata, emb_key, label_key, batch_key, model_name):
    from scib.metrics.lisi import lisi_graph
    from scib.metrics.silhouette import silhouette, silhouette_batch

    emb_key_ = "X_emb"
    adata.obsm[emb_key_] = adata.obsm[emb_key]
    sc.pp.neighbors(adata, use_rep=emb_key_)
    df = pd.DataFrame(index=[model_name])
    df["ilisi"], df["clisi"] = lisi_graph(adata, batch_key, label_key, type_="embed")
    df["sil_batch"] = silhouette_batch(adata, batch_key, label_key, emb_key_)
    df["sil_labels"] = silhouette(adata, label_key, emb_key_)

    return df


###################################################################################################
## Classification helpers
import xgboost as xgb
from sklearn.model_selection import train_test_split, StratifiedKFold, cross_val_score


def train_xgboost(
    df: pd.DataFrame,
    X_train: np.ndarray,
    y_train: np.ndarray,
    X_test: np.ndarray,
    y_test: np.ndarray,
) -> xgb.XGBClassifier:
    # Ref: https://machinelearningmastery.com/evaluate-gradient-boosting-models-xgboost-python/
    # Doc: https://xgboost.readthedocs.io/en/latest/parameter.html#learning-task-parameters
    # objective='multi:softmax'
    # objective='multi:softprob'
    # by default it uses objective='multi:softprob'
    
    # Setup
    xgb_clf = xgb.XGBClassifier(
        objective="multi:softmax",
        num_class=df.target.cat.categories.size,
        tree_method="gpu_hist",
        predictor="gpu_predictor",
        eval_metric=["merror", "mlogloss"],
        random_state=42,
    )

    # Fit
    xgb_clf.fit(
        X_train, y_train, eval_set=[(X_train, y_train), (X_test, y_test)], 
        early_stopping_rounds=10, verbose=True
    )

    # Plots
    results = xgb_clf.evals_result()
    epochs = len(results["validation_0"]["mlogloss"])
    x_axis = range(0, epochs)

    fig, ax = plt.subplots(1, 2, figsize=(20, 5))
    for idx, metric in enumerate(["merror", "mlogloss"]):
        ax[idx].plot(x_axis, results["validation_0"][metric], label="Train")
        ax[idx].plot(x_axis, results["validation_1"][metric], label="Test")
        ax[idx].legend()
        ax[idx].set_ylabel(metric)
    plt.show()

    return xgb_clf


def split_group(group):
    train_data, test_data = train_test_split(group, test_size=0.2)
    return train_data, test_data


def train_test_split_by_group(
    df: pd.DataFrame,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    train_test_data = df.groupby("target").apply(split_group)

    X_train = pd.concat([item[0] for item in train_test_data])
    y_train = X_train["target"].cat.codes
    X_train = X_train.values[:, :-1]

    X_test = pd.concat([item[1] for item in train_test_data])
    y_test = X_test["target"].cat.codes
    X_test = X_test.values[:, :-1]
    # X_test.groupby('target').apply(len) / df.groupby('target').apply(len)
    # X_train.groupby('target').apply(len) / df.groupby('target').apply(len)

    return X_train, y_train, X_test, y_test


# def train_test_group_split(adata: anndata.AnnData, groupby: str):
def train_test_split_by_group_torch(adata: sc.AnnData, groupby: str):
    """
    Function to split anndata object 80/20 per group in format
    required for SCANVIDeep explainer.
    """
    import torch
    groups = adata.obs.groupby(groupby)
    train, test = [], []
    for _, cells in groups.groups.items():
        train_test = train_test_split(cells.values, test_size=0.1)
        
        train.append(train_test[0])
        test.append(train_test[1])

    train, test = np.concatenate(train), np.concatenate(test)
    
    X_train = {
        'X': torch.from_numpy(adata[train].layers['counts'].A).type(torch.DoubleTensor),
        'batch': torch.from_numpy(adata[train].obs.batch.cat.codes.values[:, np.newaxis]),
        'labels': torch.from_numpy(adata[train].obs.ct.cat.codes.values[:, np.newaxis])
    }

    X_test = {
        'X': torch.from_numpy(adata[test].layers['counts'].A).type(torch.DoubleTensor),
        'batch': torch.from_numpy(adata[test].obs.batch.cat.codes.values[:, np.newaxis]),
        'labels': torch.from_numpy(adata[test].obs.ct.cat.codes.values[:, np.newaxis])
    }
    
    return train, X_train, test, X_test


def predict(model: SCANVI, threshold: int = 0.85):
    predictions = model.predict(soft=True)
    df = pd.DataFrame(
        zip(predictions.idxmax(axis=1), predictions.max(axis=1)),
        columns=["pred", "score"],
    )
    df["pred_filt"] = "Unknown"
    df.loc[df.score >= threshold, "pred_filt"] = df.loc[df.score >= threshold, "pred"]

    return df.pred_filt.values


def naive_classification(adata: sc.AnnData, markers: Dict[str, list[str]]):
    suffix: str = "_score"
    ctrl_size = min(list(map(len, markers.values())))
    for ct, genes in markers.items():
        sc.tl.score_genes(
            adata, genes, score_name=f"{ct}{suffix}", ctrl_size=ctrl_size, use_raw=False
        )

    adata.obs["naive"] = (
        adata.obs[[f"{x}{suffix}" for x in markers.keys()]]
        .idxmax(axis=1)
        .str.replace(suffix, "")
    )

    return adata


###################################################################################################
## SHAP helpers


def feature_plot(df: pd.DataFrame, shap_values: np.ndarray, classes: pd.Index):
    fig, ax = plt.subplots(7, 2, sharex=True, figsize=[20, 40])

    for idx, ct in enumerate(classes):
        avg = (
            pd.DataFrame(shap_values[idx], index=df.index, columns=df.columns)
            .abs()
            .mean(axis=0)
            .sort_values(ascending=False)
            .reset_index()
            .rename(columns={"index": "feature", 0: "weight"})
            .query("weight > 0")
            .head(10)
        )

        sns.barplot(x="weight", y="feature", data=avg, ax=ax[idx // 2, idx % 2])
        ax[idx // 2, idx % 2].set_title(
            f"Mean(|SHAP value|) average importance for: {ct}"
        )
