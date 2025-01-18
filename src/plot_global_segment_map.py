#!/usr/bin/env python
import numpy as np
import pandas as pd
import geopandas as gpd
import proplot as pplt
from tqdm import tqdm
from shapely.geometry import box
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from matplotlib.lines import Line2D
from matplotlib.path import Path
from matplotlib.patches import PathPatch
from matplotlib_scalebar.scalebar import ScaleBar
import matplotlib as mpl
import matplotlib.patheffects as path_effects


mpl.rcParams['agg.path.chunksize'] = 10000

def plot(input_segments=None, input_sinks=None, lines=None, nodes=None, sinks=None,
         output=None, dpi=600, nolines=False):
    """Plot global map of segments with insets of deltas."""

    lines = lines if lines is not None else gpd.read_file(input_segments, layer='lines')
    #nodes = nodes if nodes is not None else gpd.read_file(input_segments, layer='nodes')
    sinks = gpd.read_file(input_sinks) if input_sinks is not None else sinks

    fig, (ax,) = pplt.subplots(
        projection="eqearth", refwidth="210mm", dpi=dpi, reso="xx-hi",
    )
    lw_func = lambda da: da**0.43 * 0.0018
    lines["line_width"] = lw_func(lines.drainage_area_out)
    #lines.loc[lines["line_width"] < 0.07, "line_width"] = 0.07
    #lines.loc[lines["line_width"] > 4, "line_width"] = 4
    print(lines.line_width.describe())

    # suble sink blobs
    if sinks is not None:
        sinkhdl = sinks.plot(ax=ax, edgecolor="none", facecolor="0.5")

    # color
    colorm = pplt.Colormap(["red", "blue", "green"])
    lines["category"] = pd.cut(lines["line_width"], 25)
    line_kw = dict(path_effects=[path_effects.Stroke(capstyle="round", joinstyle="round")], zorder=1,
                   vmax=100, vmin=100)  #cmap=colorm, column="bifurcation_balance",
    for bin, lix in tqdm(lines.groupby("category").groups.items()):
        if len(lix):
            if not nolines:
                lines.loc[lix].plot(ax=ax, lw=(bin.left+bin.right)/2, **line_kw)
            pass

    legvls = [5000, 50000, 500000, 5000000]
    leghdls, leglab = darea_legend_handles(legvls, map(lw_func, legvls))
    ax.legend(leghdls, leglab, loc="lower center", ncols=4, align="c", frame=False,
              markerscale=6, title="Drainage area [10$^3$km$^2$]", order="C", columnspacing=0.5)

    mapkw = dict(coast=True, borders=False,  lakes=True, lakescolor="#2B72AD", lakesalpha=0.5,
                 borderscolor='gray5', coastcolor='gray6')
    ax.format(
        suptitle="GRIT segments (>30m wide)",
        lonlim=(-179.9, 179.9), latlim=(-60, 75), **mapkw,
    )
    # insets
    keys = ["title", "axes_location", "scale_lines", "lonlim", "latlim", "extent_loc", "scale_loc"]
    values = [
        ("Amazon", [0.07, 0.073, 0.2, 0.45], 2, (-52.8, -49.3), (-2, 0.9), (1, 4), 2),
        ("Pearl", [0.85, 0.49, 0.12, 0.25], 2, (112.55, 113.7), (22.1, 23.45), (2, 3), 3),
        ("Padma-Brahmap.", [0.65, 0.29, 0.1, 0.2], 2, (89.1, 90.05), (23.54, 24.15), (1, 2), 3),
        ("Mekong", [0.64, 0.00, 0.17, 0.27], 2, (104.88, 106.7), (9.4, 11.65), (1, 4), 3),
        ("Rhine-Meuse", [0.32, 0.62, 0.13, 0.2], 4,(3.5, 6.1), (51.4, 52.1), (1, 4), 4),
        ("Congo", [0.41, 0.16, 0.11, 0.26], 2., (15.25, 15.6), (-4.02, -4.35), (1, 4), 2),
        ("Fraser", [0.058, 0.65, 0.12, 0.2], 3., (-123.15, -122.6), (49.04, 49.27), (1, 4), 4),
    ]
    insets = {v[0]: dict(zip(keys, v)) for v in values}
    # extent boxes in right crs
    boxes = gpd.GeoSeries({i: box(a["lonlim"][0], a["latlim"][0], a["lonlim"][1], a["latlim"][1])
                           for i, a in insets.items()}, crs="EPSG:4326").to_crs(lines.crs)
    # plot insets
    inset_mapkw = dict(longrid=False, latgrid=False, ocean=False, coast=False, title_kw=dict(size='small'), titlepad=0.1)
    for i, a in insets.items():
        a.update(inset_mapkw)
        zoom_inset(ax, lines.loc[lines.sindex.query(boxes[i])], **a)

    if output:
        fig.save(output, dpi=dpi)
    return ax


def zoom_inset(ax, lines, axes_location, scale_lines=1., **fmtkw):
    inset = ax.inset(axes_location, proj="eqearth", reso="xx-hi")
    amaz = lines.groupby("category").groups.items()
    for bin, lix in tqdm(amaz):
        if len(lix):
            lw = (bin.left+bin.right) / 2 * scale_lines
            lines.loc[lix].plot(ax=inset, lw=lw, zorder=1,
                                path_effects=[path_effects.Stroke(capstyle="round", joinstyle="round")])
            pass
    scloc = fmtkw.pop("scale_loc")
    mark_inset(ax, inset, alpha=0.3, **dict(zip(["loc1", "loc2"], fmtkw.pop("extent_loc"))))
    inset.format(**fmtkw)
    x, y, w, h = axes_location
    sb = ScaleBar(
        dx=1,
        font_properties=dict(size=6.5),
        frameon=False,
        location=scloc,
        width_fraction=0.0025/h,
        length_fraction=0.3,
        sep=1,
        label="(DA %ix)" % scale_lines,
        label_loc="right"
    )
    inset.add_artist(sb)
    return inset


def darea_legend_handles(values, width):
    path_data = [
        (Path.MOVETO, (-0.5, -0.2)),
        (Path.CURVE3, (-0.3, 0.5)),
        (Path.CURVE3, (0, 0)),
        (Path.CURVE3, (0.3, -0.5)),
        (Path.CURVE3, (0.5, 0.2)),
        ]
    codes, verts = zip(*path_data)
    path = Path(verts, codes)

    leghdls = [Line2D([0], [0], lw=0, markeredgewidth=d, marker=path, markerfacecolor="none", markersize=3) for d in width]
    leglab = [f"{int(l/1000)}" for l in values]
    return leghdls, leglab


if __name__ == "__main__":
    import fire
    fire.Fire(plot)
