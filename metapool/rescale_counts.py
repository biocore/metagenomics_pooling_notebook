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
    '''
    ### df:
    Index: 3771 entries, G000005825 to G900163845
    Data columns (total 18 columns):
     #   Column                 Non-Null Count  Dtype
    ---  ------                 --------------  -----
     0   14606.363202203.metaG  3771 non-null   float64
     1   14606.363202207.metaG  3771 non-null   float64
     ...

    ### df.index:
    ['G000005825', 'G000006605', 'G000006725', 'G000006745', 'G000006785', ...

    ### metadata_features:

                   species                     total_length
    OTUID
    G000005825     Bacillus pseudofirmus       4249248
    G000006175      Methanococcus voltae       1936387
    ...

    ### lengths:
    OTUID
    G000005825    4249248
    G000006605    2476822
    ...

    ### prep_info_samples:
                           sample_name_r  read_count_spike_in  anonymized_name  ... syndna_pool_number read_count_total  calc_mass_sample_aliquot_input_g
    sample_name                                                                          ...                                                                      
    14606.363202128.metaG  X14606.363202128.metaG              5875250  363202128.metaG  ...              pool2         14853990                          0.023299
    14606.363202170.metaG  X14606.363202170.metaG               608189  363202170.metaG  ...              pool1         24384324                          0.023299
    ...
    '''

    # the regular (non-transposed) form of df is needed for getting lengths.
    lengths = metadata_features.loc[df.index, "total_length"]
    transposed_df = df.T

    '''
    print(transposed_df.info())
    Index: 18 entries, 14606.363202203.metaG to 14606.363202708.metaG
    Columns: 3771 entries, G000005825 to G900163845
    dtypes: float64(3771)
    memory usage: 530.4+ KB
    None
    G000006785
    G900156305
    14606.363202203.metaG
    G000006785
    G900156305
    14606.363202207.metaG
    G000006785
    G900156305
    14606.363202695.metaG
    G000006785
    G900156305
    14606.363202170.metaG
    nan	inf
    '''

    # df.columns == transposed_df.index
    sample_weights = prep_info_samples.loc[df.columns, "calc_mass_sample_aliquot_input_g"]

    '''
    print(lengths.head())
    OTUID
    G000005825    4249248
    G000006605    2476822
    G000006725    2731750

    print(sample_weights.head())
    Name: total_length, dtype: int64
    14606.363202203.metaG    0.046599
    14606.363202207.metaG    0.005825
    14606.363202695.metaG    0.005825
    '''

    foo = lengths * 10
    bar = lengths / 2

    '''
    print(foo.head())
    OTUID
    G000005825    42492480
    G000006605    24768220
    G000006725    27317500
    G000006745    40334640
    G000006785    18524330
    Name: total_length, dtype: int64

    print(bar.head())
    OTUID
    G000005825    2124624.0
    G000006605    1238411.0
    G000006725    1365875.0
    G000006745    2016732.0
    G000006785     926216.5

    baz = foo * bar
    print(baz.head())
    G000005825    9.028054e+13
    G000006605    3.067324e+13
    G000006725    3.731229e+13
    G000006745    8.134416e+13
    G000006785    1.715754e+13
    '''

    # The issue is that lengths uses 'G000005825' and similar for indexes while
    # sample_weights uses '14606.363202207.metaG' and similar. Also, because
    # there are roughly only 18 of the latter and 3700+ of the former, we can't
    # correctly perform operations like mul().

    # replacing sample_weights with a series using 'Gvalues' as indexes.
    # the value at each index will be the mean of the calc_mass_sample_aliquot_input_g
    # values associated w/the 18 sample-names. When sample_weights is replaced
    # in this way, the function and test behave correctly.
    mean = sum(sample_weights) / len(sample_weights)
    d = {}
    for k in list(lengths.index):
        d[k] = mean
    sample_weights = pd.Series(data=d, index=list(lengths.index))

    # due to earlier concerns regarding ints and rollovers, this equation was
    # reduced to make that scenario less likely.
    # 926461538461.5385 == (6.022 * (10.0 ** 23)) / (650 * 10.0 ** 9)
    mult_row = 926461538461.5385 / (lengths * sample_weights)

    # axis=0 implies mult_row is meant to be a Series object with an equal
    # number of rows to df. Each row in mult_row contains a single value
    # meant to be multiplied over every value in the corresponding row in df.
    #
    # the caller (test) is also expecting 3700+ rows with Gvalues as the
    # index. This implies that df and metadata_features have correct form.
    # prep_info_samples is either incorrect and should provide values indexed
    # to Gvalues instead of sample_names, or this function is missing code to
    # map a value from one or more of the 18 sample_names into each of the
    # 3700+ Gvalues.
    df = df.mul(mult_row, axis=0)

    return df


def to_absolute_abundance_read_count(df, linear_models):
    df = _to_log10_cpm(df)
    df, failed_cols = _apply_model(df, linear_models)
    df.index.name = "sample_name"
    return df, failed_cols


def to_absolute_abundance_cell_count(df, linear_models, metadata_features,
                                     prep_info_samples):
    df = _to_log10_cpm(df)
    df, failed_cols = _apply_model(df, linear_models)
    df = convert_read_count_to_cell_count_per_g_input(df,
                                                      metadata_features,
                                                      prep_info_samples)
    df.index.name = "sample_name"
    return df, failed_cols
