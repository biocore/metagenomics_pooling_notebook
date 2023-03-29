import pandas as pd
import numpy as np


def _to_column_percentage(df):
    return df / df.sum() * 100


def _to_row_percentage(df):
    return df.div(df.sum(axis=1), axis=0) * 100


def _to_log10_cpm(df):
    aligned_reads = df.sum()
    df_log10_cpm = np.log10(df / aligned_reads * 1000000)
    return df_log10_cpm


def _apply_model(df_log10_cpm, linear_models):
    cols = []
    failed_colnames = []
    for c in df_log10_cpm:
        r = linear_models[linear_models["sample_name"] == c]
        if r.empty:
            # TODO FIXME HACK:  What do we do if this sample can't be rescaled
            #  because no linear model was created?
            failed_colnames.append(c)
            continue
        converted = r["a_intercept"].iloc[0] + df_log10_cpm[c] * \
            r["b_slope"].iloc[0]
        converted = np.power(10, -converted)
        cols.append(converted)

    converted_reads = pd.DataFrame(data=cols).T
    return converted_reads, failed_colnames


def convert_read_count_to_cell_count_per_g_input(df, metadata_features, prep_info_samples):
    lengths = metadata_features.loc[df.index, "total_length"]
    sample_weights = prep_info_samples.loc[df.index, "calc_mass_sample_aliquot_input_g"]
    # Careful with types here.  If you use ints,
    # the length * 650 * 10**9 can overflow integers with very long genomes
    mult_row = (6.022 * (10.0 ** 23)) / (lengths * (650 * 10.0 ** 9) * sample_weights)
    df = df.mul(mult_row, axis=0)
    return df

def to_absolute_abundance_read_count(df, linear_models):
     df = _to_log10_cpm(df)
     df, failed_cols = _apply_model(df, linear_models)
     df.index.name = "sample_name"
     return df, failed_cols

def to_absolute_abundance_cell_count(df, linear_models, metadata_features):
    df = _to_log10_cpm(df)
    df, failed_cols = _apply_model(df, linear_models)
    df = convert_read_count_to_cell_count_per_g_input(df, metadata_features, prep_info_samples)
    df.index.name = "sample_name"
    return df, failed_cols
