import scipy
import seaborn as sns
import pandas as pd
import numpy as np


DILUTIONS = pd.DataFrame({
    'dilution': {
        0: 1.0, 1: 0.1, 2: 0.01, 3: 0.001, 4: 0.0001, 5: 1.0, 6: 0.1, 7: 0.01,
        8: 0.001, 9: 0.0001, 10: 0.25, 11: 0.025, 12: 0.0025, 13: 0.00025,
        14: 0.000025, 15: 0.25, 16: 0.025, 17: 0.0025, 18: 0.00025,
        19: 0.000025},
    'dilution_id': {
        0: 1, 1: 2, 2: 3, 3: 4, 4: 5, 5: 1, 6: 2, 7: 3, 8: 4, 9: 5, 10: 1,
        11: 2, 12: 3, 13: 4, 14: 5, 15: 1, 16: 2, 17: 3, 18: 4, 19: 5},
    'plasmid_id': {
        0: 'p126', 1: 'p136', 2: 'p146', 3: 'p156', 4: 'p166', 5: 'p266',
        6: 'p256', 7: 'p246', 8: 'p236', 9: 'p226', 10: 'p126', 11: 'p136',
        12: 'p146', 13: 'p156', 14: 'p166', 15: 'p266', 16: 'p256', 17: 'p246',
        18: 'p236', 19: 'p226'},
    'syndna_pool_number': {
        0: 'pool1', 1: 'pool1', 2: 'pool1', 3: 'pool1', 4: 'pool1', 5: 'pool1',
        6: 'pool1', 7: 'pool1', 8: 'pool1', 9: 'pool1', 10: 'pool2',
        11: 'pool2', 12: 'pool2', 13: 'pool2', 14: 'pool2', 15: 'pool2',
        16: 'pool2', 17: 'pool2', 18: 'pool2', 19: 'pool2'},
    'plasmid_sequence_id': {
        0: 'seq_1_gc=0.26_PmeI_NcoI', 1: 'seq_1_gc=0.36_PmeI_XhoI',
        2: 'seq_1_gc=0.46_PmeI_NcoI', 3: 'seq_1_gc=0.56_PmeI_XmaI',
        4: 'seq_1_gc=0.66_PmeI_XhoI', 5: 'seq_2_gc=0.66_XhoI_PacI',
        6: 'seq_2_gc=0.56_XmaI_PacI', 7: 'seq_2_gc=0.46_NcoI_PacI',
        8: 'seq_2_gc=0.36_XhoI_XmaI', 9: 'seq_2_gc=0.26_NcoI_PcaI',
        10: 'seq_1_gc=0.26_PmeI_NcoI', 11: 'seq_1_gc=0.36_PmeI_XhoI',
        12: 'seq_1_gc=0.46_PmeI_NcoI', 13: 'seq_1_gc=0.56_PmeI_XmaI',
        14: 'seq_1_gc=0.66_PmeI_XhoI', 15: 'seq_2_gc=0.66_XhoI_PacI',
        16: 'seq_2_gc=0.56_XmaI_PacI', 17: 'seq_2_gc=0.46_NcoI_PacI',
        18: 'seq_2_gc=0.36_XhoI_XmaI', 19: 'seq_2_gc=0.26_NcoI_PcaI'}
    })


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


def calculate_coefficients(table_synthetic_hits, metadata_pools,
                           plot_fit=False):
    # Calculate the total number of reads aligned
    # to the plasmid sequences for each sample
    table_synthetic_hit_totals = table_synthetic_hits.sum().to_frame()
    table_synthetic_hit_totals.index.name = 'sample_name'
    # table_synthetic_hit_totals.columns = ['table_synthetic_hit_totals']

    # table_synthetic_hits.loc['table_synthetic_hit_totals'] = table_synthetic_hits.sum()

    # Convert feature-table to long-form and merge with total number of hits to
    # plasmid sequences
    table_synthetic_hits_long_with_totals_beta = pd.melt(
        table_synthetic_hits.reset_index(), id_vars="OTUID")
    # table_synthetic_hits_long_with_totals_beta.columns = ["plasmid_id",
    #                                                       "sample_name",
    #                                                       "count"]

    table_synthetic_hits_long_with_totals = \
        table_synthetic_hits_long_with_totals_beta.merge(
            table_synthetic_hit_totals, on="sample_name")
    print (table_synthetic_hits_long_with_totals)

    # Merge updated feature-table with sample-pool information
    table_synthetic_hits_long_with_totals_and_pools = \
        table_synthetic_hits_long_with_totals.merge(
            metadata_pools, on="sample_name")

    # Aggregate across lanes (not needed)
    table_synthetic_hits_long_with_totals_and_pools_all_lanes = \
        table_synthetic_hits_long_with_totals_and_pools

    print (table_synthetic_hits_long_with_totals_and_pools_all_lanes)
    table_synthetic_hits_long_with_totals_and_pools_all_lanes = table_synthetic_hits_long_with_totals_and_pools_all_lanes[['plasmid_id', 'sample_name', 'count', 'syndna_pool_number', 'read_count_total']]

    # Calculate counts per million
    table_synthetic_hits_long_with_totals_and_pools_all_lanes[
        "counts_per_million"] = (
            table_synthetic_hits_long_with_totals_and_pools_all_lanes[
                "count"] /
            table_synthetic_hits_long_with_totals_and_pools_all_lanes[
                "read_count_total"]
        ) * 1000000

    # Calculate log10 of counts per million
    table_synthetic_hits_long_with_totals_and_pools_all_lanes[
        "counts_per_million_log10"] = np.log10(
        table_synthetic_hits_long_with_totals_and_pools_all_lanes[
            "counts_per_million"])

    # Merge table with plasmid dilution information
    table_synthetic_hits_with_dilutions = \
        table_synthetic_hits_long_with_totals_and_pools_all_lanes.merge(
            DILUTIONS, on=["plasmid_id", "syndna_pool_number"])

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
            final_table_1.groupby(["syndna_pool_number", "sample_name"]),
            xmax=4,
            ymax=4
        )

    result = {
        "syndna_pool_number": [],
        "sample_name": [],
        "a_intercept": [],
        "b_slope": []
    }

    for name, group in final_table_1.groupby(["syndna_pool_number",
                                              "sample_name"]):
        m, b, r, p, std_err = scipy.stats.linregress(
            group["counts_per_million_log10"], group["dilution_id"])
        result["syndna_pool_number"].append(name[0])
        result["sample_name"].append(name[1])
        result["a_intercept"].append(b)
        result["b_slope"].append(m)

    coefall = pd.DataFrame(result)
    coefall = coefall.sort_values(by=["syndna_pool_number",
                                      "sample_name"])

    return coefall
