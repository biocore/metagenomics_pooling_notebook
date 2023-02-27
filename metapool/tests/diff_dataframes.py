import pandas as pd
import numpy as np


def dataframes_are_equal(df1, df2, verbose=False):
    df1 = df1[sorted(df1.columns)]
    df2 = df2[sorted(df2.columns)]

    df1 = df1.loc[sorted(df1.index), :]
    df2 = df2.loc[sorted(df2.index), :]

    indices_match = (df1.index == df2.index).all()
    columns_match = (df1.columns == df2.columns).all()
    if verbose:
        print("Indices Match:", indices_match)
        print("Columns Match:", columns_match)

    fuzzy_match_sum = 0
    for gotu in df1.columns:
        if df1.dtypes[gotu] == np.float64:
            for sample in df1.index:
                l = df1.loc[sample, gotu]
                r = df2.loc[sample, gotu]
                if -0.0001 < l - r < 0.0001:
                    fuzzy_match_sum += 1
                else:
                    if verbose:
                        print("UNMATCHED:", gotu, sample, l, r)
        else:
            for sample in df1.index:
                if df1.loc[sample, gotu] == df2.loc[sample, gotu]:
                    fuzzy_match_sum += 1
                else:
                    if verbose:
                        print("UNMATCHED:", gotu, sample, df1.loc[sample, gotu], df2.loc[sample, gotu])


    fuzzy_match_exp = df2.shape[0] * df2.shape[1]
    if verbose:
        print("Matched cells", fuzzy_match_sum)
        print("Unmatched cells", fuzzy_match_exp - fuzzy_match_sum)
        if fuzzy_match_sum != fuzzy_match_exp:
            print("Error: SOME ARE UNMATCHED!!")
    return (fuzzy_match_exp == fuzzy_match_sum) and indices_match and columns_match


if __name__ == "__main__":
    import sys
    df1 = pd.read_csv(sys.argv[1], sep="\t", index_col="sample_name")
    df2 = pd.read_csv(sys.argv[2], sep="\t", index_col="sample_name")
    dataframes_are_equal(df1, df2)

