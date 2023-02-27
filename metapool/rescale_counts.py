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
        converted = r["a_intercept"].iloc[0] + df_log10_cpm[c] * r["b_slope"].iloc[0]
        converted = np.power(10, -converted)
        cols.append(converted)

    converted_reads = pd.DataFrame(data=cols).T
    return converted_reads, failed_colnames


def convert_read_count_to_cell_count(df, metadata_features):
    lengths = metadata_features.loc[df.index, "total_length"]
    # Careful with types here.  If you use ints,
    # the length * 650 * 10**9 can overflow integers with very long genomes
    mult_row = (6.022 * (10.0**23)) / (lengths * (650 * 10.0**9))
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
    df = convert_read_count_to_cell_count(df, metadata_features)
    df.index.name = "sample_name"
    return df, failed_cols


# Port of R code's coverage filtering that we aren't planning to use at the moment:
# Assumes 150 bp per read, if this changes, need to update to modify coverage calcuation.
# BP_PER_READ = 150

# feature_counts_raw = table_community.melt(id_vars="OTUID", var_name="SampleID", value_name="count_raw")
# feature_counts_raw = feature_counts_raw[feature_counts_raw["count_raw"] > 0]
# feature_counts_raw = feature_counts_raw.merge(aligned_reads, left_on="SampleID", right_index=True)
# feature_counts_raw.columns = ["OGU", "SampleID", "count_raw", "AlignedReads"]
# feature_counts_raw = feature_counts_raw[["SampleID", "AlignedReads", "OGU", "count_raw"]]

# R's coverage filter states: Keep only those gotus with an expected sequencing depth of 1x
# based on the number of reads, the assumed 150bp read length, and the genome length.
# We are skipping this filter entirely as it can be done externally.
# genome_length_average = []
# total_reads = []

# TODO FIXME HACK: R code is kind of unclear to me here: it seems that it just wants the length of each OGU
#  but it calculates a mean... mean of what?  All OTUIDs only have a single length...
# feature_counts_raw = feature_counts_raw.merge(metadata_features[["OTUID", "total_length"]], left_on="OGU", right_on="OTUID")
# feature_counts_raw = feature_counts_raw.merge(metadata_pools[["read_count_total", "sample_name"]], left_on="SampleID", right_on="sample_name")
# feature_counts_raw = feature_counts_raw.drop(["OTUID", "sample_name"], axis=1)
# feature_counts_raw.columns = ["SampleID", "AlignedReads", "OGU", "count_raw", "genome_length", "total_reads"]

# Coverage is the total number of basepairs covered by reads (assuming 150bp per read) divided by the genome length
# feature_counts_raw["coverage"] = feature_counts_raw["count_raw"] * BP_PER_READ / feature_counts_raw["genome_length"]
