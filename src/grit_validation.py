#!/usr/bin/env python
import os.path as osp

import numpy as np
import pandas as pd
import geopandas as gpd


def compare_didnt(*input, output=None, bad_azimuth=135, good_azimuth=45):
    ids = [osp.splitext(osp.basename(f))[0] for f in input]
    dat = [gpd.read_file(f, index_col=0) for f in input]
    dat = pd.concat(dat, keys=ids, names=["domain_id", dat[0].index.name])
    dat["invalid"] = (dat.snap_dist > dat.avgwidth) | (dat.avgwidth < 30) | dat.azimuth_diff.between(good_azimuth, bad_azimuth)
    dat["good"] = (~ dat["invalid"]) & (dat.azimuth_diff < good_azimuth)
    dat["bad"] = (~ dat["invalid"]) & (dat.azimuth_diff > bad_azimuth)
    # sum if true and convert to km length and to km2 river area
    sumifkm = lambda x: (x * dat.loc[x.index, "length_split"]).sum() * 1e-3
    sumifkm2 = lambda x: (x * dat.loc[x.index, "length_split"] * dat.loc[x.index, "avgwidth"]).sum() * 1e-6
    outtbl = []
    for f in (sumifkm, sumifkm2):
        deltas = dat.groupby("delta_name")["invalid good bad".split()].agg(f)
        # add total row at end
        deltas.loc["TOTAL"] = deltas.sum()
        # invalid as percentage of total
        deltas["invalid"] = deltas.invalid * 100 / deltas.sum(axis=1)
        # absolute valid
        deltas.insert(0, "valid", deltas[["good", "bad"]].sum(axis=1))
        # good & bad as percentage of valid
        deltas["good"] = deltas.good * 100 / deltas["valid"]
        deltas["bad"] = deltas.bad * 100 / deltas["valid"]
        outtbl.append(deltas)
    outtbl = pd.concat(outtbl, axis=1, keys=["length", "area"])
    if output:
        outtbl.round(2).to_csv(output)
    else:
        return outtbl


def select_random_bifurcations(segments, n=1000, output=None, include_lines=True):
    """
    Plot hist:
    ax=bif_sel.width.hist(bins=list(range(30, 5000, 30))+[bif_sel.width.max()],
            cumulative=True, density=True, label="width-weighted sample", histtype="step")
    ax=bif.width.hist(bins=list(range(30, 5000, 30))+[bif.width.max()],
            cumulative=True, density=True, label="all", histtype="step", ax=ax)
    """
    nodes = gpd.read_file(segments, layer="nodes").set_index("global_id")
    lines = gpd.read_file(segments, layer="lines").set_index("global_id")
    # only take bifurcations that are not on lakes
    bif = nodes[(nodes["node_type"] == "bifurcation") & (nodes["grwl_value"] != 180)]
    bif_sel = bif.sample(n=n, weights="width")
    if include_lines:
        bif_lines = lines[(lines.upstream_node_id.isin(bif_sel.index)) |
                          (lines.downstream_node_id.isin(bif_sel.index))]
    if output:
        bif_sel.to_file(output, layer="nodes")
        if include_lines:
            bif_lines.to_file(output, layer="lines")
    else:
        return (bif_sel, bif_lines) if include_lines else bif_sel


def plot_manual_uncertainty_scores(input, output=None, ax=None):
    import proplot as pplt

    valid = pd.read_csv(input, index_col="fid")

    valid["validation"] = valid['Louise_Yinxue_Laurence_Boen_Validation'].copy()

    check = valid['doublecheck_validation'].apply(lambda x: str(int(x)) if np.isfinite(x) else "")
    valid["validation_error"] = check != ""
    valid.loc[valid.validation_error, "validation"] = check[valid.validation_error]
    valid.validation = valid.validation.replace({"?": 2, np.nan: 2}).astype(int)

    # expand integer meanings
    # lump canals into rivers
    valid.replace({
        "grwl_value": {255: "inland", 86: "inland", 126: "coastal"},
        "validation": {0: "incorrect", 1: "correct", 2: "uncertain"},
    }, inplace=True)

    # aggregate
    counts = valid.groupby(["grwl_value", "validation"]).global_id.count()
    grwl_counts = counts.groupby(level="grwl_value").sum()
    shares = counts.div(grwl_counts) * 100


    if ax is None:
        fig, (ax,) = pplt.subplots(nrows=1, ncols=1)
    icshares = shares.unstack() #.drop("correct", axis=1)
    icshares.index = [f"{c} (N={grwl_counts[c]})" for c in icshares.index]
    bars = ax.bar(icshares, cycle="ggplot")
    ax.legend(title=None, ncols=1, location="ur")
    ax.format(ylabel="%", xlabel="", title="GRIT bifurcations")
    
    if output:
        fig.save(output, dpi=300)
    else:
        return fig, ax


def compare_hydrography_flow_accumulation(input_grit, input_hydro, output=None):
    buffer = 2.5/60  # in degrees, ie roughly 5x5 kernel

    grit = gpd.read_file(input_grit).to_crs(epsg=4326)
    hs = gpd.read_file(input_hydro).to_crs(epsg=4326)

    gb = grit.copy()
    gb.geometry = gb.buffer(buffer)

    jnd = gb.sjoin(hs).reset_index()
    jnd["diff"] = (jnd.upa_left - jnd.upa_right) * 100 / jnd.upa_left
    jnd["abs_diff"] = jnd["diff"].abs()
    idx = jnd.groupby("index")["abs_diff"].idxmin()
    sel = jnd.loc[idx].set_index("index")
    grit["diff"] = sel["diff"]
    grit["upa_other"] = sel["upa_right"]

    if output:
        grit.to_file(output, layer="drainage_area_difference")
    else:
        return grit
    return


def plot_hydrography_flow_accumulation(hydrosheds, merit, output):
    import proplot as pplt
    import matplotlib.cm as cm
    import matplotlib.colors as mcolors
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature

    cmap = pplt.Colormap('Balance')

    fig, axs = pplt.subplots(nrows=2, ncols=1, figsize=(6, 7), projection=ccrs.EqualEarth())

    points = [gpd.read_file(hs) for hs in (hydrosheds, merit)]
    titles = ["GRIT - HydroSHEDS", "GRIT - MERIT-Hydro"]
    for ax, pts, ti in zip(axs, points, titles):
        ax.add_feature(cfeature.LAND, facecolor='lightgrey', edgecolor='black', linewidth=0.1)
        pts = pts.sort_values("upa")
        pts["siz"] = (pts["upa"].values**0.32)*6e-4
        sc = ax.scatter(
            pts.geometry.x.values,
            pts.geometry.y.values,
            c = pts["diff"].values,
            s = pts.siz.values,
            cmap=cmap, vmin=-100, vmax=100, extend="both", edgecolor="none", alpha=0.5, absolute_size=True,
        )
        ax.set_title(ti)
    axs.format(
        lonlim=(-170, 170),  # Longitude limits
        latlim=(-60, 80),  # Latitude limits
    )
    # Create the colorbar separately
    cbar = fig.colorbar(sc, loc='b', length=0.8, label='Deviation in drainage area (%)')


    if output:
        fig.save(output, dpi=300)
    else:
        return fig, ax


if __name__ == "__main__":
    from commandline import interface as _cli
    _cli()
