import numpy as np
import pandas as pd


def convert_10x_h5_to_h5ad(path_10x_h5, output, batch_name=None):
    import scanpy
    adata_all = scanpy.read_10x_h5(path_10x_h5, gex_only=False)
    adata_gex_only = scanpy.read_10x_h5(path_10x_h5)
    for feature_type in [
            x for x in np.unique(adata_all.var['feature_types'])
            if x != 'Gene Expression'
    ]:
        mask = adata_all.var['feature_types'] == feature_type
        # TODO: evaluate whether better to save data/metadata together as a single dataframe?
        adata_gex_only.obsm[feature_type] = adata_all.X[:, mask]
        adata_gex_only.uns[feature_type + '_metadata'] = adata_all.var[mask]
    if batch_name is None:
        batch_name = output
    adata_gex_only.obs.rename(lambda x: x + '_' + batch_name, inplace=True)
    adata_gex_only.var['gene_name'] = adata_gex_only.var.index
    adata_gex_only.var.index = adata_gex_only.var['gene_ids']
    adata_gex_only.write_h5ad(output)


def predict_guide_status_from_grna_counts(counts,
                                          nguide_threshold_for_cell=0,
                                          ncell_threshold_for_guide=10):
    from pomegranate import GeneralMixtureModel, PoissonDistribution
    np.seterr(divide='ignore')
    logcounts = np.log(counts) / np.log(2)
    if (~np.isfinite(logcounts)).sum() > 0:
        cellmask = np.isfinite(logcounts)
        if cellmask.sum(
        ) >= ncell_threshold_for_guide:  # have minimum cells for guide
            try:
                model = GeneralMixtureModel.from_samples(
                    [PoissonDistribution, PoissonDistribution],
                    n_components=2,
                    X=logcounts[cellmask].reshape(-1, 1))
                predictions = []
                for val in logcounts:
                    if not np.isfinite(val):
                        predictions.append(np.nan)
                    else:
                        predictions.append(
                            model.predict(np.array(val).reshape(-1, 1))[0])
            except:
                pass
                predictions = [0] * len(counts)

        else:
            predictions = [0] * len(counts)
    else:
        predictions = [0] * len(counts)
    predictions = np.array(predictions)
    if predictions.sum() > 0:
        if counts[predictions == 0].mean() > counts[predictions == 1].mean():
            predictions = 1 - predictions
    predictions = np.nan_to_num(predictions, nan=0.0)
    threshold_for_cell_mask = counts >= nguide_threshold_for_cell
    predictions *= threshold_for_cell_mask
    return predictions


def generate_crispr_grna_predictions(adata,
                                     write_to_obsm=True,
                                     write_to_obs=False,
                                     obsm_key='CRISPR Guide Capture'):
    from panopticrispr.utilities import predict_guide_status_from_grna_counts
    from tqdm import tqdm
    if (not write_to_obs) and (not write_to_obsm):
        raise Exception("One of write_to_obs or write_to_obsm must be True")
    dfs = []
    if type(adata.obsm[obsm_key]) == pd.core.frame.DataFrame:
        for col, series in tqdm(adata.obsm[obsm_key].items(),
                                total=adata.obsm[obsm_key].shape[1]):
            predictions = predict_guide_status_from_grna_counts(series.values)

            if predictions.sum() > 0:
                dfs.append(
                    pd.DataFrame(predictions,
                                 columns=[col],
                                 index=adata.obs.index))
                if write_to_obs:
                    adata.obs[guide + '_prediction'] = predictions
    else:
        for iguide, guide in enumerate(
                tqdm(adata.uns[obsm_key + '_metadata'].index.values)):
            predictions = predict_guide_status_from_grna_counts(
                np.array(adata.obsm[obsm_key][:, iguide].todense())[:, 0])

            if predictions.sum() > 0:
                dfs.append(
                    pd.DataFrame(predictions,
                                 columns=[guide],
                                 index=adata.obs.index))
                if write_to_obs:
                    adata.obs[guide + '_prediction'] = predictions
    if write_to_obsm:
        adata.obsm[obsm_key + '_predictions'] = pd.concat(dfs, axis=1)


def get_masked_differential_expression(adata,
                                       mask1=None,
                                       mask2=None,
                                       verbose=False,
                                       mask1_downsample_size=None,
                                       mask2_downsample_size=None,
                                       min_cluster_size=0,
                                       gene_alternate_name=None):

    from scipy.stats import mannwhitneyu
    from statsmodels.stats.multitest import fdrcorrection
    from time import time
    from tqdm import tqdm

    if mask1_downsample_size:
        p = np.min([mask1_downsample_size, np.sum(mask1)]) / np.sum(mask1)
        mask1 *= np.random.choice([True, False],
                                  p=[p, 1 - p],
                                  size=mask1.shape[0])
    if mask2_downsample_size:
        p = np.min([mask2_downsample_size, np.sum(mask2)]) / np.sum(mask2)
        mask2 *= np.random.choice([True, False],
                                  p=[p, 1 - p],
                                  size=mask2.shape[0])
    if verbose:
        print('Group 1 size: ', np.sum(mask1), ', group 2 size: ',
              np.sum(mask2))
    pvalues = []
    uvalues = []
    genes = []
    meanexpr1 = []
    meanexpr2 = []
    meanexpexpr1 = []
    meanexpexpr2 = []
    fracexpr1 = []
    fracexpr2 = []
    start = time()
    data1 = np.array(adata.X[mask1, :].todense().T)
    if verbose:
        print('First matrix extracted in', time() - start, 'seconds')
    start = time()
    data2 = np.array(adata.X[mask2, :].todense().T)
    print(data1.shape, data2.shape)

    if verbose:
        print('Second matrix extracted', time() - start, 'seconds')
    for igene, gene in enumerate(
            tqdm(adata.var['gene_name'],
                 desc='Computing Mann-Whitney p-values')):
        genes.append(gene)
        if np.std(data1[igene, :]) + np.std(data2[igene, :]) < 1e-14:
            pvalues.append(1)
            uvalues.append(np.nan)
        else:
            mw = mannwhitneyu(data1[igene, :],
                              data2[igene, :],
                              alternative='two-sided')
            pvalues.append(mw.pvalue)
            uvalues.append(mw.statistic)
        meanexpr1.append(data1[igene, :].mean())
        meanexpr2.append(data2[igene, :].mean())
        meanexpexpr1.append(np.mean(2**data1[igene, :]))
        meanexpexpr2.append(np.mean(2**data2[igene, :]))
        fracexpr1.append((data1[igene, :] > 0).mean())
        fracexpr2.append((data2[igene, :] > 0).mean())
    output = pd.DataFrame(genes)
    output.columns = ['gene']
    output['pvalue'] = pvalues
    output['CommonLanguageEffectSize'] = np.array(uvalues) / (data1.shape[1] *
                                                              data2.shape[1])
    output['MeanExpr1'] = meanexpr1
    output['MeanExpr2'] = meanexpr2
    output['MeanExpExpr1'] = meanexpexpr1
    output['MeanExpExpr2'] = meanexpexpr2
    output['FracExpr1'] = fracexpr1
    output['FracExpr2'] = fracexpr2
    if gene_alternate_name is not None:
        gene2altname = {
            gene: altname
            for gene, altname in zip(loom.ra['gene'],
                                     loom.ra[gene_alternate_name])
        }
        altnames = [gene2altname[x] for x in genes]
        output['GeneAlternateName'] = altnames

    output = output.sort_values('CommonLanguageEffectSize', ascending=False)
    output['BenjaminiHochbergQ'] = fdrcorrection(output['pvalue'],
                                                 is_sorted=False)[1]
    return output


def generate_grna_umi_merge(adata, umi_str='_umi', grna_umi_index='index'):
    if grna_umi_index == 'index':
        adata.uns['CRISPR Guide Capture_metadata']['target'] = [
            x.split(umi_str)[0]
            for x in adata.uns['CRISPR Guide Capture_metadata'].index.values
        ]
    else:
        adata.uns['CRISPR Guide Capture_metadata']['target'] = adata.uns[
            'CRISPR Guide Capture_metadata'][grna_umi_index].apply(
                lambda x: x.split(umi_str)[0])
    grna_umi_collapsed_counts = []
    grnas = []
    for grna in adata.uns['CRISPR Guide Capture_metadata']['target'].unique():
        grnas.append(grna)
        mask = adata.uns['CRISPR Guide Capture_metadata']['target'] == grna
        grna_umi_collapsed_counts.append(
            np.array(adata.obsm['CRISPR Guide Capture'][:,
                                                        mask].sum(axis=1)[:,
                                                                          0]))
    df = pd.DataFrame(np.hstack(grna_umi_collapsed_counts),
                      columns=grnas,
                      index=adata.obs.index)
    adata.obsm['CRISPR Guide Capture_UMI collapsed'] = df


def get_grna_umi_prediction_summary(adata, target_suffix='_crispr'):
    df = adata.obsm['CRISPR Guide Capture_predictions'].copy()
    df['nGuide'] = df.sum(axis=1)
    df = df[df['nGuide'] == 1].T
    df['nCell'] = df.sum(axis=1).astype(int)
    df['target'] = [x.split(target_suffix)[0] for x in df.index.values]
    return df


def generate_grna_umi_labeling(adata, target_suffix='_crispr'):
    df = get_grna_umi_prediction_summary(adata, target_suffix=target_suffix)
    df = df[[x for x in df.columns if x not in ['nCell', 'target']]].iloc[0:-1]
    hits = np.where(df == 1)
    labeldict = {
        cell: grna_umi
        for cell, grna_umi in zip(df.columns[hits[1]], df.index[hits[0]])
    }
    adata.obs['gRNA_UMI'] = [
        labeldict[x] if x in labeldict.keys() else np.nan
        for x in adata.obs.index
    ]

def get_differential_expression_over_continuum(adata,
                                               mask,
                                               covariate,
                                               method='spearman',
                                               gene_alternate_name=None):

    if method not in ['spearman', 'kendall']:
        raise Exception(
            "Requested method not implemented.  Only `kendall` or `spearman` currently available."
        )
    if np.sum(mask) != len(covariate):
        raise Exception(
            "Length of covariate vector does not match mask length.")
    from tqdm import tqdm
                                     
    if method == 'spearman':
        from scipy.stats import spearmanr
        method_module = spearmanr
    elif method == 'kendall':
        from scipy.stats import kendalltau
        method_module = kendalltau
    corrs = []
    pvals = []
    genes = []
    X = np.array(adata.X[mask,:].todense())
    for i, gene in enumerate(tqdm(adata.var['gene_name'], desc='looping over genes')):
        if np.std(X[:,i]) < 1e-14:
            pvals.append(1)
            corrs.append(0)
        else:
            result = method_module(X[:,i], covariate, nan_policy='omit')
            pvals.append(result.pvalue)
            corrs.append(result.correlation)
        genes.append(gene)
    df = pd.DataFrame(genes)
    df.columns = ['gene']
    df['pval'] = pvals
    df['corr'] = corrs
    return df
