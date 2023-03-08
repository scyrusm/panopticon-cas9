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
        adata.obsm[obsm_key+'_predictions'] = pd.concat(dfs, axis=1)
