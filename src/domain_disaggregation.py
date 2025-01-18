#!/usr/bin/env python
from collections import defaultdict

import geopandas as gpd
from shapely.ops import unary_union
import pylab as pl
from tqdm import tqdm
import pandas as pd
import numpy as np

pd.options.mode.chained_assignment = None


def add_in_envelope(feature, otherfeatures, return_added_ix=True):
    missed = otherfeatures.index[
        otherfeatures.envelope.apply(lambda g: feature.envelope.contains(g))
    ]
    idgeom = unary_union([feature] + list(otherfeatures.geometry[missed]))
    return (idgeom, missed) if return_added_ix else idgeom


def rectangular_domains_grow_largest_within_distance(
    polygons, within_column=None, max_distance=None
):
    """Deprecated version in favour of rectangular_domains_grow_largest."""

    unassigned = polygons.copy()
    unassigned["area"] = unassigned.area

    # pick largest
    largest_area = unassigned.envelope.area.max()
    id = 0
    domains = {"largest_feature": {}, "geometry": {}}
    pbar = tqdm(total=len(polygons))
    while len(unassigned) > 0:
        nextix = unassigned.area.idxmax()
        nextpoly = unassigned.loc[nextix]
        unassigned.drop(nextix, inplace=True)
        idgeom = nextpoly.geometry
        if within_column:
            candidates = unassigned[
                unassigned[within_column] == nextpoly[within_column]
            ]
        else:
            candidates = unassigned.copy()
        candidates["next_distance"] = candidates.convex_hull.distance(
            idgeom.convex_hull
        )
        if max_distance:
            candidates = candidates[candidates["next_distance"] <= max_distance]
        candidates.sort_values(
            ["next_distance", "area"], ascending=[True, False], inplace=True
        )
        add = []
        for a, feat in candidates.iterrows():
            ngeom = unary_union([idgeom, feat.geometry])
            if ngeom.envelope.area <= largest_area:
                idgeom = ngeom
                add.append(a)
        # add those enclosed in new domain but missed by distance
        idgeom, missed = add_in_envelope(idgeom, unassigned)
        mergeix = set(add + list(missed))
        unassigned.drop(mergeix, inplace=True)
        # store domains and largest feature
        domains["geometry"][id] = idgeom
        domains["largest_feature"][id] = nextix
        # increase id
        id += 1
        pbar.update(len(mergeix) + 1)
    pbar.close()

    return gpd.GeoDataFrame(domains)


def envelope_compactness(feature):
    return feature.area / feature.envelope.area


def rectangular_domains_grow_largest(
    polygons, within_column="region", min_compactness=None, max_distance=None
):
    """Deprecated version in favour of rectangular_domains_grow_most_compact."""

    unassigned = polygons.copy()
    unassigned["area"] = unassigned.area

    # pick largest
    largest_area = unassigned.envelope.area.max()
    id = 0
    domains = {"largest_feature": {}, "geometry": {}}
    pbar = tqdm(total=len(polygons))
    while len(unassigned) > 0:
        nextix = unassigned.area.idxmax()
        nextpoly = unassigned.loc[nextix]
        unassigned.drop(nextix, inplace=True)
        idgeom = nextpoly.geometry
        # within
        candidates = (
            unassigned.loc[unassigned[within_column] == nextpoly[within_column]]
            if within_column
            else unassigned.copy()
        )
        candidates["next_distance"] = candidates.convex_hull.distance(
            idgeom.convex_hull
        )
        # search max_distance
        if max_distance:
            candidates = candidates.loc[candidates["next_distance"] <= max_distance]
        # intersecting area
        candidates["intersect_area"] = candidates.geometry.apply(
            lambda g: idgeom.intersection(g).area
        )
        # compactness
        candidates["compactness"] = candidates.geometry.apply(
            lambda g: envelope_compactness(unary_union([idgeom, g]))
        )
        candidates.sort_values(
            ["next_distance", "intersect_area", "compactness"],
            ascending=[True, True, False],
            inplace=True,
        )
        add = []
        for a, feat in candidates.iterrows():
            ngeom = unary_union([idgeom, feat.geometry])
            if ngeom.envelope.area <= largest_area:
                if min_compactness and envelope_compactness(ngeom) < min_compactness:
                    continue
                idgeom = ngeom
                add.append(a)
        # add those enclosed in new domain but missed by distance
        idgeom, missed = add_in_envelope(idgeom, unassigned)
        mergeix = set(add + list(missed))
        unassigned.drop(mergeix, inplace=True)
        # store domains and largest feature
        domains["geometry"][id] = idgeom
        domains["largest_feature"][id] = nextix
        # increase id
        id += 1
        pbar.update(len(mergeix) + 1)
    pbar.close()

    return gpd.GeoDataFrame(domains)


def rectangular_domains_grow_most_compact(
    polygons, region_column="region", min_compactness=None, debug=True
):
    """Deprecated version in favour of rectangular_domains_global_compactness."""

    unassigned = polygons.copy()
    unassigned["area"] = unassigned.area
    unassigned["compactness"] = unassigned.geometry.apply(envelope_compactness)

    # pick largest
    largest_area = unassigned.envelope.area.max()
    id = 0
    domains = {"largest_feature": {}, "geometry": {}}
    pbar = tqdm(total=len(polygons))
    while len(unassigned) > 0:
        nextix = unassigned.compactness.idxmin()
        nextpoly = unassigned.loc[nextix]
        unassigned.drop(nextix, inplace=True)
        idgeom = nextpoly.geometry
        # within region
        candidates = unassigned.loc[
            unassigned[region_column] == nextpoly[region_column]
        ]
        intersects = candidates[candidates.intersects(idgeom.envelope)]
        nadded = 1
        all_intersects_checked = False
        while len(intersects) > 0 and not all_intersects_checked:
            print("N intersecting: %s" % len(intersects))

            idgeom, added = add_in_envelope(idgeom, intersects)
            if len(added) > 0:
                intersects.drop(added, inplace=True)
                unassigned.drop(added, inplace=True)
                candidates.drop(added, inplace=True)
                nadded += len(added)
                if debug:
                    print(f"added {len(added)} polygons already within {id}")

            # compactness and potential new areas
            intersects["combined_geometry"] = intersects.geometry.apply(
                lambda g: unary_union([idgeom, g])
            )
            intersects["combined_area"] = intersects.combined_geometry.apply(
                lambda g: g.envelope.area
            )
            intersects["combined_compactness"] = intersects.combined_geometry.apply(
                lambda g: envelope_compactness(g)
            )

            intersects = intersects[intersects.combined_area <= largest_area]

            if min_compactness and len(intersects) > 0:
                intersects = intersects[
                    intersects.combined_compactness >= min_compactness
                ]

            if len(intersects) == 0:
                all_intersects_checked = True
                break

            intersects["next_distance"] = intersects.convex_hull.distance(
                idgeom.convex_hull
            )
            intersects.sort_values(
                ["next_distance", "combined_compactness"], ascending=False, inplace=True
            )

            addf = intersects.iloc[0]
            idgeom = addf.combined_geometry
            unassigned.drop(addf.name, inplace=True)
            candidates.drop(addf.name, inplace=True)
            nadded += 1
            if debug:
                print(f"added {addf.name} to {id}")

            # update intersects
            intersects = candidates[candidates.intersects(idgeom.envelope)]

        # store domains and largest feature
        domains["geometry"][id] = idgeom
        domains["largest_feature"][id] = nextix
        # increase id
        id += 1
        if not debug:
            pbar.update(nadded)
    pbar.close()

    return gpd.GeoDataFrame(domains)


def split_polygons_by_lines(polygons, lines):
    """Split polygon dataframe by line geometries."""
    from shapely.ops import split

    intersects = [polygons[polygons.intersects(l)] for l in lines]
    newpolys = []
    for ints, sl in zip(intersects, lines):
        splitgeoms = ints.geometry.apply(lambda g: split(g, sl)).to_frame()
        splitgeoms = splitgeoms.explode(index_parts=False)
        splitgeoms = splitgeoms[splitgeoms.index.duplicated(keep=False)]
        newpolys.append(ints.loc[splitgeoms.index].set_geometry(splitgeoms.geometry))
    newpolys = pd.concat(newpolys)
    # swap polygons
    newpolygons = polygons.drop(newpolys.index.unique())
    newpolys.reset_index(drop=True, inplace=True)
    newpolys.index = newpolys.index + 1 + newpolygons.index.max()
    return newpolygons.append(newpolys)


def neighbour_areas(poly, others, buffer=0):
    pl = (
        poly.geometry.convex_hull.buffer(buffer)
        if buffer
        else poly.geometry.convex_hull
    )
    nb = others[~others.convex_hull.disjoint(pl)]
    combgeom = nb.geometry.apply(lambda g: unary_union([poly.geometry, g])).to_frame()
    combgeom["area"] = combgeom.area
    combgeom["envelope_area"] = combgeom.envelope.area
    combgeom["compactness"] = combgeom["area"] / combgeom["envelope_area"]
    return combgeom


def neighbour_envelope_overlap(polygons, neighbours, dissolve=True):
    """Return areas where neighbours have overlaps with other polygons."""
    neiix = neighbours.envelope.to_frame("geometry").reset_index()
    polix = polygons.envelope.to_frame("geometry").reset_index()
    ol = neiix.overlay(polix, how="intersection")
    oll = ol[(ol["level_0"] != ol["index"]) & (ol["level_1"] != ol["index"])]
    if dissolve:
        oll = oll.dissolve(["level_0", "level_1"])
    return oll


def non_overlapping_ratio(polygons, neighbours):
    nor = pd.Series(1, index=neighbours.index)
    olld = neighbour_envelope_overlap(polygons, neighbours)
    nor[olld.index] = 1 - olld.area/neighbours.loc[olld.index].envelope.area
    return nor


def rectangular_domains_global_compactness(
    polygons,
    region_column=None,
    neighbour_buffer=0,
    largest_size=None,
    ignore_size=0,
    objective="compactness_change",
    return_diagnostics=False,
):
    """Aggregate polygons to most compact groups up to largest_size.

    Polygons are successively aggregated by merging neighbouring polygons that offer the
    largest increase in envelope compactness until largest_size is reached.

    Arguments
    =========
    polygons : geopandas.GeoDataFrame
        Polygons to be aggregated.
    region_column : str
        Column in polygons to preserve regions. If None, all polygons may be aggregated.
    largest_size : value in map units
        Largest domain size to create. Default is the domain size of the largest polygon.
    ignore_size : value in map units
        Area of polygons to filter out.
    objective : str | list
        The objective to chose next neighbours to merge, either a single or a
        combination, in which case their product is used. All values are normalised to 0-1:
        Options:
        - non_overlapping_ratio : ratio of overlap free area
        - compactness : area / envelope area
        - compactness_change : change in envelope compactness (maybe negative)
        - area_ratio : envelope area / largest_size
        - area_deficit : 1 - area_ratio
    return_diagnostics : bool
        Also return a DataFrame with stats about the aggregation iterations.
    """
    if type(objective) == str:
        objective = [objective]
    polygons = polygons[polygons.geometry.area > ignore_size]
    newpolygons = polygons.copy(deep=True)
    largesta = largest_size or newpolygons.envelope.area.max()
    regions = (
        polygons[region_column] if region_column else pd.Series(0, index=polygons.index)
    )
    regids = polygons.groupby(regions).groups

    # find neighbours
    alreadynb = defaultdict(lambda: [])
    neighbours = []
    for id, poly in tqdm(polygons.iterrows(), total=len(polygons)):
        reg = polygons.loc[regids[regions[id]]].drop([id] + alreadynb[id])
        nbrs = neighbour_areas(poly, reg, buffer=neighbour_buffer)
        neighbours.append(nbrs)
        [alreadynb[i].append(id) for i in nbrs.index]
    neighbours = pd.concat(neighbours, keys=polygons.index)
    n_neighbours = pd.Series(0, index=polygons.index)
    nn = neighbours.groupby(level=0).geometry.count()
    n_neighbours[nn.index] = nn
    # calculate other objectives
    ar, ea = polygons.area, polygons.envelope.area
    l0, l1 = neighbours.index.get_level_values(0), neighbours.index.get_level_values(1)
    cptbf = (ar[l0].values + ar[l1].values) / (ea[l0].values + ea[l1].values)
    neighbours["compactness_change"] = neighbours["compactness"] - cptbf
    neighbours["non_overlapping_ratio"] = non_overlapping_ratio(polygons, neighbours)
    neighbours["area_ratio"] = neighbours.envelope.area / largesta
    neighbours["area_deficit"] = 1 - neighbours["area_ratio"]

    # aggregate by domain compactness
    nlast = len(neighbours)
    pbar = tqdm(total=nlast)
    diagnostics = [[len(polygons), polygons.area.sum() / polygons.envelope.area.sum()]]
    while True:
        neighbours = neighbours[neighbours.envelope_area <= largesta]
        neighbours = neighbours[neighbours.compactness_change >= 0]
        if len(neighbours) == 0:
            break
        newids = neighbours[objective].product(axis=1).idxmax()
        # remove smaller, replace bigger
        sm, big = newpolygons.loc[list(newids)].area.sort_values().index
        oldcpt = newpolygons.area[[sm, big]].sum() / newpolygons.envelope.area[[sm, big]].sum()
        newpoly = newpolygons.loc[big]
        newpoly.geometry = neighbours.loc[newids].geometry
        newpolygons.loc[big] = newpoly
        newpolygons.drop(sm, inplace=True)
        # remove all neighbours with big and small
        neighbours.drop([big, sm], level=0, inplace=True, errors="ignore")
        neighbours.drop([big, sm], level=1, inplace=True, errors="ignore")
        # recalculate neighbours
        reg = newpolygons[regions[newpolygons.index] == regions[big]].drop(big)
        newnb = neighbour_areas(newpoly, reg, buffer=neighbour_buffer)
        if len(newnb):  # may be 0
            newnb.index = pd.MultiIndex.from_tuples([(big, i) for i in newnb.index])
            newnb["compactness_change"] = newnb["compactness"] - oldcpt
            newnb["non_overlapping_ratio"] = non_overlapping_ratio(newpolygons, newnb)
            newnb["area_ratio"] = newnb.area / largesta
            newnb["area_deficit"] = 1 - newnb["area_ratio"]
            neighbours = neighbours.append(newnb)
        pbar.update(nlast - len(neighbours))
        nlast = len(neighbours)
        if return_diagnostics:
            totalcpt = newpolygons.area.sum() / newpolygons.envelope.area.sum()
            diagnostics.append([len(newpolygons), totalcpt])

    if return_diagnostics:
        diagnostics = pd.DataFrame(diagnostics, columns=["n", "mean_compactness"])

    # aggregate all domains that are fully enclosed in other envelopes
    for id, poly in newpolygons.iterrows():
        contained = newpolygons[newpolygons.envelope.contains(poly.geometry)].drop(id)
        if len(contained):
            nid = contained.distance(poly.geometry).idxmin()
            newpoly = newpolygons.loc[nid]
            newpoly.geometry = unary_union([newpoly.geometry, poly.geometry])
            newpolygons.loc[nid] = newpoly
            newpolygons.drop(id, inplace=True)

    # add compactness to final polygons
    newpolygons["envelope_compactness"] = newpolygons.geometry.apply(
        envelope_compactness
    )
    return (newpolygons, diagnostics) if return_diagnostics else newpolygons


def plot_domains_map(basins, domains, output=None, ax=None):
    fig, ax = (ax.get_figure(), ax) if ax else pl.subplots()
    map = domains.plot(column="name", cmap="Paired", ax=ax)
    basins.plot(color="none", edgecolor="k", alpha=0.5, ax=map)
    domains.envelope.plot(color="none", edgecolor="k", alpha=0.6, ax=map)

    if output:
        fig.savefig(output, dpi=300)
    return map


def plot_domain_area_histograms(basins, domains, output=None, ax=None):
    fig, ax = (ax.get_figure(), ax) if ax else pl.subplots()
    # area histogram
    basins.envelope.area.hist(
        bins=np.logspace(8, 13), density=True, color="k", histtype="step", ax=ax
    )
    domains.envelope.area.hist(
        bins=np.logspace(8, 13), density=True, color="b", histtype="step", ax=ax
    )
    ax.set_xscale('log')
    ax.set_yscale('log')
    if output:
        fig.savefig(output, dpi=300)
    return ax


def domains_with_analysis(basins_path, output_dir, domains_output=None, **kw):
    from shapely.geometry import LineString
    from pathlib import Path

    outdir = Path(output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    basins = gpd.read_file(basins_path).to_crs("EPSG:8857")

    args = dict(
        region_column="region",
        neighbour_buffer=100e3,
        ignore_size=1000e6,
        return_diagnostics=True,
    )
    args.update(kw)
    domains, diagnostics = rectangular_domains_global_compactness(basins, **args)

    domains.to_crs("EPSG:4326").to_file(domains_output or str(outdir / "domains.gpkg"))
    diagnostics.to_csv(str(outdir / "diagnostics.csv"))
    envelopes = domains.copy(deep=True)
    envelopes.set_geometry(domains.envelope.to_crs("EPSG:4326"), inplace=True)
    envelopes.to_file(str(outdir / "envelopes.gpkg"))

    plot_domain_area_histograms(basins, domains,
        output=str(outdir / "plot_domain_histograms.png"))
    plot_domains_map(basins, domains, output=str(outdir / "plot_domain_map.png"))

    # stats
    st = ["basins", "continents", "domains"]
    stv = [basins, domains.dissolve("region"), domains]
    stats = pd.concat([s.envelope.area.mul(1e-6).describe() for s in stv], keys=st).unstack()
    for s, v in zip(st, stv):
        stats.loc[s, "total_area"] = v.envelope.area.mul(1e-6).sum()
        stats.loc[s, "envelope_compactness"] = v.area.mul(1e-6).sum()/stats.loc[s, "total_area"]
    stats.to_csv(str(outdir / "stats.csv"))
    return domains, stats


def tweak_basins(basins, exclude_areas=None, split_lines=None, output=None):
    """Exclude and/or split basins using other geometries.
    """
    basins = gpd.read_file(basins)

    if exclude_areas:
        deselect = gpd.read_file(exclude_areas)
        basins = basins[basins.geometry.apply(lambda g: ~deselect.contains(g).any())]
    if split_lines:
        splitlines = gpd.read_file(split_lines)
        basins = split_polygons_by_lines(basins, splitlines.geometry)

    if output:
        basins.to_file(output)
    return basins


def plot_tweaks(basins, exclude_areas=None, split_lines=None, output=None):
    """Plot split and exclusion geometries.
    """

    basins = gpd.read_file(basins)

    fig, ax = pl.subplots()
    basins.plot(color='none', edgecolor="0.5", ax=ax)

    if exclude_areas:
        deselect = gpd.read_file(exclude_areas)
        deselect.plot(ax=ax, color='r', alpha=0.5)
    if split_lines:
        splitlines = gpd.read_file(split_lines)
        splitlines.centroid.plot(ax=ax, color="r")

    if output:
        fig.savefig(output, dpi=300, bbox_inches='tight')
    return ax


def test(**kw):
    print(kw)
    for k, v in kw.items():
        print(k, v, type(v))


if __name__ == "__main__":
    from commandline import interface as _cli

    _cli()
