#!/usr/bin/env python
import gc
import os.path as osp
import warnings

import numpy as np
import grass.script as grass
import pandas as pd
from grass_attribute_table import GrassAttributeTable
from tqdm import tqdm

from river_network_utils import (
    network_raster_segments,
    read_raster_xyz,
)

# enable pandas progress reporting
tqdm.pandas()


def route_value(
    value_raster,
    network_vector,
    segments_raster,
    output_raster=None,
    output=None,
    output_segments_csv=None,
    output_nodes_csv=None,
    output_column_prefix=None,
    centreline_subbasins=None,
    method="sum",
    ncells_raster=None,
    width_column="width",
    order_column="strahler_order",
    centreline_value_raster=None,
    network_cache=None,
):
    """Route raster values through a directed network.

    Routing methods
    ---------------
    sum :
        Cumulative sum of centreline values, partitioned at bifurcations.
    mean :
        Average catchment conditions for each centerline cell after the `sum`
        method. Requires the `ncells_raster` to devide the cumulative sum by.
    mainstem :
        Same as `sum` but without partitioning at bifurcations.

    Arguments
    ---------
    value_raster : str
        Name of raster with values to route by method. NULL values are assumed 0.
    network_vector : str
        Name of directed river network vector with the needed line attributes given by
        `[angle, width, order]_column`. Vertices not in `segment_raster` are ignored.
    segments_raster : str
        Name of the network raster, i.e. segments, the values have no significance.
    output_raster, output (deprecated): str
        Name of output raster.
    output_segments_csv : str
        Path to write segment input and output values to.
    output_nodes_csv : str
        Path to write node input and output values to.
    output_column_prefix : str
        Upload segment and node in/out values.
    centreline_subbasins : None | str
        Subbasins to aggregate values to. If not given, `value_raster` is assumed to be
        already aggregated to the centreline subbasins.
    method : str, "sum" (default) | "mean" | "mainstem"
        Aggregation method along the network, see above for explanation. If 'mean',
        the `ncells_raster` is required to created the weighted means.
    ncells_raster : str | None
        Upstream drainage area in cells of the centreline cells to be used with `method='mean'`
        to create the weighted means.
    centreline_value_raster : None | str
        The name of the centreline aggregated values raster if centreline_subbasins is given.
        If None, `value_raster`__centreline__`method` is used.
    width_column, order_column : str
        Columns of the network lines table needed for the routing and partitioning at
        bifurcations.
    network_cache : path | None
        File to read or write (if non existent path) dense network to.
    """
    if output and (output_raster is None):
        warnings.warn("Use output_raster instead of output.", DeprecationWarning)
        output_raster = output

    if method == "mean":
        assert ncells_raster, "If method=mean, `ncells_raster` must be given."
    # import ipdb;ipdb.set_trace()
    # sum values to centerline subbasins
    if centreline_subbasins:
        cl_values = (
            centreline_value_raster
            or f"{value_raster.split('@')[0]}__centreline__{method}"
        )
        grass.run_command(
            "r.stats.zonal",
            cover=value_raster,
            base=centreline_subbasins,
            output=cl_values,
            method="sum",
        )
    else:
        cl_values = value_raster

    # read in centerline values
    rasters = [segments_raster, cl_values] + ([ncells_raster] if method == "mean" else [])
    print("Reading raster values...", rasters)
    raster_values = read_raster_xyz(rasters)

    # get topology and raster cell coordinates without node cells
    segments, nodes, links = network_raster_segments(
        network_vector,
        check_missing=raster_values.index,
        cache_segments=network_cache,
        trim_nodes=True,
    )
    print("Reading segments is done!")
    default_value = -9999
    # import ipdb;ipdb.set_trace()
    # get values for each centerline cell
    # ISSUE: some segments are not presented in the raster_values
    segment_values = segments.progress_apply(
        lambda s: raster_values.loc[s, cl_values].values if all(val in raster_values.index for val in s) else default_value
    )
    # get values for nodes separately
    node_values = raster_values.loc[nodes["coor"], cl_values].set_axis(nodes.index)

    # prepare stores
    routed = segment_values.copy()
    node_in = pd.Series(0.0, index=nodes.index)
    node_out = pd.Series(
        {i: np.zeros(len(s)) for i, s in nodes["downstream_segments"].items()}
    )
    print(
        f"Routing {len(segments)} segments with {len(nodes)} nodes and "
        f"{nodes.bifurcation.sum()} bifurcations."
    )
    # link ids ordered by strahler order
    lordergrps = sorted(links.groupby(order_column).groups.items())
    for o, lix in tqdm(lordergrps, "Routing orders"):
        # get upstream nodes of all links in this order
        onodes = nodes.loc[links.loc[lix, "start_node"].unique()]
        for n, ds in onodes.downstream_segments.items():
            # add subbasin value of node
            node_in[n] += node_values[n]
            # no bifurcation
            if len(ds) == 1:
                node_out[n] = [node_in[n]]
            else:  # bifurcation partitioning
                wid = links.loc[ds, width_column]
                if method == "mainstem":
                    # set partition to 1 for mainstem and 0 for others
                    partition = np.zeros(len(wid))
                    partition[np.argmax(wid.values)] = 1
                else:
                    # default width partioning
                    # width - d. area relationship inversion (Frasson et al., 2019)
                    # parts = (wid / 9.6) ** (1/0.32)
                    partition = (wid/wid.sum()).values
                # partitioning
                node_out[n] = node_in[n] * partition
                # print("Bifurcation at node %i, partitioning [%%]:" % n)
                # print(parts.mul(100).round(1).to_string())
        # route
        for li in lix:
            usn, dsn = links.loc[li, ["start_node", "end_node"]]
            # Q going out of upstream node into line
            usq = node_out[usn][nodes.downstream_segments[usn].index(li)]
            # if line has pi
            if li in routed.index:
                routed[li] = routed[li].cumsum() + usq
                isnan = np.isnan(routed[li])
                if isnan.sum():
                    unroutednan = np.isnan(segment_values[li][isnan])
                    miss = isnan.sum()*100/len(routed[li])
                    warnings.warn(f"Found {miss}% nans in routed segment {li}. Loop? Using unrouted values instead.")
                    if unroutednan.sum():
                        warnings.warn(f"Found {unroutednan.sum()} nans in unrouted segment {li}. Gaps in predictor raster?")
                    routed[li][isnan] = segment_values[li][isnan]
                usq = routed[li][-1]
            node_in[dsn] += usq

    # divide accumulated sums by ncells if mean
    if method == "mean":
        # get values for each centerline cell
        segment_ncells = segments.progress_apply(
            lambda s: raster_values.loc[s, ncells_raster].values
        )
        # get values for nodes separately
        node_ncells = raster_values.loc[nodes["coor"], ncells_raster].set_axis(nodes.index)
        assert (node_ncells == 0).sum() == 0
        node_in = node_in / node_ncells
        # sum the node_out at bifurcation -YL 2023/08/08
        # NEED TO CHECK WHAT node_out is for
        node_out = node_out.apply(lambda x: sum(x) if len(x) > 0 else 0)
        node_out = node_out / node_ncells
        routed = routed / segment_ncells

    # write node input (single) / output (csv list)
    if output_nodes_csv or output_column_prefix:
        out_nod = node_in.to_frame("in")
        out_nod["out"] = node_out.progress_apply(lambda x: ", ".join(map(str, x)))
        if output_nodes_csv:
            out_nod.to_csv(output_nodes_csv, index_label="cat")
    # write out segment in (from node output) / out (last element in routed)
    if output_segments_csv or output_column_prefix:
        out_segm = pd.DataFrame(columns=["in", "out"], index=links.index)
        for l, li in tqdm(links.iterrows(), total=len(links)):
            idx = nodes.loc[li.start_node, "downstream_segments"].index(l)
            out_segm.loc[l, "in"] = node_out[li.start_node][idx]
            # zero-len segments have out=in
            out_segm.loc[l, "out"] = routed[l][-1] if l in routed.index else out_segm.loc[l, "in"]
        if output_segments_csv:
            out_segm.to_csv(output_segments_csv, index_label="cat")
    # upload to vector table
    if output_column_prefix:
        segm = GrassAttributeTable(network_vector, layer=1)
        nod = GrassAttributeTable(network_vector, layer=2)
        segm[[f"{output_column_prefix}_{c}" for c in out_segm.columns]] = out_segm
        nod[[f"{output_column_prefix}_{c}" for c in out_nod.columns]] = out_nod
        segm.write()
        nod.write()
    # write full centerline raster out
    if output_raster:
        print("Writing values to raster...")
        tmpf = grass.tempfile()
        with open(tmpf, "w") as f:
            node_in.set_axis(pd.MultiIndex.from_tuples(nodes.coor)).to_csv(f, header=False)
            for li, l in routed.iteritems():
                df = pd.DataFrame(segments[li])
                df.loc[:, 2] = l
                df.to_csv(f, header=False, index=False)
        # clear memory before reading in data
        del segments, raster_values, routed, segment_values
        gc.collect()
        grass.run_command(
            "r.in.xyz",
            input=tmpf,
            output=output_raster,
            separator="comma",
            percent=5,
        )
    return


def create_dense_segments_cache(network_vector, segments_raster, output):
    """Just create the network cache needed for the route_value/network_raster_segments.
    """
    assert not osp.exists(output)
    raster_values = read_raster_xyz([segments_raster])
    segments, nodes, links = network_raster_segments(
        network_vector, check_missing=raster_values.index, cache_segments=output
    )
    return


if __name__ == "__main__":
    from commandline import interface as _cli

    _cli()
