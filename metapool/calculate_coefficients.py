import scipy
import seaborn as sns
import pandas as pd
import numpy as np


def plot_with_fit(xcol, ycol, groups, xmax=4, ymax=4,
                  textpos=(.6, .2)):
    from matplotlib import pyplot as plt
    import matplotlib
    matplotlib.use("TkAgg")
    x = 0
    y = 0
    fig, axes = plt.subplots(xmax, ymax, sharex=True, sharey=True)
    for name, group in groups:
        ax = axes[y, x]
        # ax.scatter(group["CPMlog"], group["Dilution2"])

        # TODO FIXME HACK: Unclear what R^2 value is being printed by ggplot
        #  but it definitely isn't a pearson computed per group.  *shrug*
        no_nan = group.dropna()
        r, p = scipy.stats.pearsonr(no_nan[xcol], no_nan[ycol])
        sns.regplot(x=xcol, y=ycol, data=group, ax=ax)
        ax.text(textpos[0], textpos[1], 'r={:.2f}\np={:.2g}'.format(r, p),
                transform=ax.transAxes)
        if x > 0:
            ax.set_ylabel("")
        if y < ymax - 1:
            ax.set_xlabel("")

        ax.set_title("\n".join(name))
        x += 1
        if x == xmax:
            x = 0
            y += 1
        if y == ymax:
            plt.show()
            x = 0
            y = 0

    if x != 0 or y != 0:
        plt.show()


def calculate_coefficients(table_synthetic_hits, metadata_pools, dilutions,
                           plot_fit=False):
    # Strip leading X in sample names which are automatically removed by R
    # table reading.
    metadata_pools["sample_name_r"] = metadata_pools["sample_name_r"].str[1:]

    # Parsing the files
    # Calculate the total number of reads aligned
    # to the plasmid sequences for each sample
    table_synthetic_hit_totals = table_synthetic_hits.sum()
    table_synthetic_hit_totals = pd.DataFrame(
        data={"sample_name_r": table_synthetic_hit_totals.index.to_list(),
              "table_synthetic_hit_totals":
                  table_synthetic_hit_totals.values.tolist()},
        index=None)

    # Convert feature-table to long-form and merge with total number of hits to
    # plasmid sequences
    table_synthetic_hits_long_with_totals_beta = pd.melt(
        table_synthetic_hits.reset_index(), id_vars="OTUID")
    table_synthetic_hits_long_with_totals_beta.columns = ["plasmid_id",
                                                          "sample_name_r",
                                                          "count"]
    table_synthetic_hits_long_with_totals = \
        table_synthetic_hits_long_with_totals_beta.merge(
        table_synthetic_hit_totals, on="sample_name_r")

    # Merge updated feature-table with sample-pool information
    table_synthetic_hits_long_with_totals_and_pools = \
        table_synthetic_hits_long_with_totals.merge(
        metadata_pools, on="sample_name_r")

    # Aggregate across lanes (not needed)
    table_synthetic_hits_long_with_totals_and_pools_all_lanes = \
        table_synthetic_hits_long_with_totals_and_pools

    # Calculate counts per million
    table_synthetic_hits_long_with_totals_and_pools_all_lanes[
        "counts_per_million"] = (
            table_synthetic_hits_long_with_totals_and_pools_all_lanes[
                "count"] /
            table_synthetic_hits_long_with_totals_and_pools_all_lanes[
                "read_count_total"]
        ) * 1000000

    # Calculate percentages
    table_synthetic_hits_long_with_totals_and_pools_all_lanes["percentage"] = (
       table_synthetic_hits_long_with_totals_and_pools_all_lanes[
          "count"] /
       table_synthetic_hits_long_with_totals_and_pools_all_lanes[
          "read_count_total"]
    ) / 100

    # Calculate log10 of counts per million
    table_synthetic_hits_long_with_totals_and_pools_all_lanes[
        "counts_per_million_log10"] = np.log10(
        table_synthetic_hits_long_with_totals_and_pools_all_lanes[
            "counts_per_million"])

    # Merge table with plasmid dilution information
    table_synthetic_hits_with_dilutions = \
        table_synthetic_hits_long_with_totals_and_pools_all_lanes.merge(
            dilutions, on=["plasmid_id", "pool"])

    final_table_1 = table_synthetic_hits_with_dilutions

    dilution_log = final_table_1["dilution_id"].copy()
    # TODO FIXME HACK:  Check that the dilutions we aim for are 5%, .5%, .05%
    #  ... or this calculation is wrong!
    dilution_log[dilution_log == 1] = 1 / 20
    dilution_log[dilution_log == 2] = 0.1 / 20
    dilution_log[dilution_log == 3] = 0.01 / 20
    dilution_log[dilution_log == 4] = 0.001 / 20
    dilution_log[dilution_log == 5] = 0.0001 / 20
    final_table_1["dilution_log"] = np.log10(dilution_log)

    if plot_fit:
        plot_with_fit(
            "counts_per_million_log10",
            "dilution_log",
            final_table_1.groupby(["pool", "sample_name"]),
            xmax=4,
            ymax=4
        )

    result = {
        "pool": [],
        "sample_name_pool": [],
        "a_intercept": [],
        "b_slope": []
    }

    for name, group in final_table_1.groupby(["pool", "sample_name_pool"]):
        m, b, r, p, std_err = scipy.stats.linregress(
            group["counts_per_million_log10"], group["dilution_id"])
        result["pool"].append(name[0])
        result["sample_name_pool"].append(name[1])
        result["a_intercept"].append(b)
        result["b_slope"].append(m)

    coefall = pd.DataFrame(result)
    coefall = coefall.sort_values(by=["pool", "sample_name_pool"])
    return coefall
