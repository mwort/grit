#!/usr/bin/env python
import os

import pandas as pd
import geopandas as gpd

DEBUG = True


def log(*args):
    if DEBUG:
        print(*args)
    return


def read_domains(domain_polygons_file, domains_selected_file):
    domain_polygons = gpd.read_file(domain_polygons_file)

    domains_selected = pd.read_csv(domains_selected_file, index_col=0)

    # filter by selected and merge multiple polygons to single cat
    domains = domain_polygons[domain_polygons.cat.isin(domains_selected.cat)].dissolve("cat")

    # set 4 letter ids and retain cat
    cat_id_map = pd.Series(domains_selected.index, index=domains_selected.cat)
    domains.set_index(cat_id_map[domains.index], inplace=True)
    domains["cat"] = domains_selected.loc[domains.index, "cat"]

    # envelope overlapping with others
    envs = domains.envelope
    for id, env in envs.items():
        domains.loc[id, "ids_in_envelope"] = ",".join(envs[envs.overlaps(env)].index)
    return domains

def read_nodes(domains, node_file_pattern):
    nodes_all = {}
    print("Reading nodes...")
    for id in domains.index:
        fl = node_file_pattern.format(domain=id)
        pts = gpd.read_file(fl).set_index("cat")[["comp", "geometry"]].dropna()
        pts["comp"] = pts.comp.astype(int)
        nodes_all[id] = pts
        print(id)
    return nodes_all

def allocate_components(nodes_all, domains):
    components = dict()
    ignore = dict()
    print("Allocating components (N / all):")
    for id, dom in domains.iterrows():
        nodes = nodes_all[id]
        # if no neighbours, take all components
        if not dom["ids_in_envelope"]:
            components[id] = nodes["comp"].unique()
            continue
        others = dom["ids_in_envelope"].split(",")
        log(f"{id}: Neighbours {others}")
        # ignore already allocated components
        if id in ignore:
            log(f"{id}: not using {ignore[id]}")
            nodes = nodes[nodes.comp.isin(ignore[id])]
        # check nodes overlaps
        other_nodes = pd.concat([nodes_all[i] for i in others], keys=others, names=["domain_id", "cat"])
        # chuck out already allocated ones
        for i in other_nodes.index.get_level_values(0).unique():
            if i in ignore:
                log(f"{id}: not using {ignore[i]} of {i}")
                icats = other_nodes.loc[i].index[other_nodes.loc[i, "comp"].isin(ignore[i])]
                other_nodes.drop(list(zip([i] * len(icats), icats)), inplace=True)
        # resulting df cols: domain_id  cat_1  comp_1  cat_2  comp_2  geometry
        ol = other_nodes.reset_index().overlay(nodes.reset_index())
        # get domain_id and (possibly multple) comps of neighbouring domains, indexed by component of this domain
        comps = ol.groupby(["domain_id", "comp_2"])["comp_1"].agg(lambda x: list(x.unique())).reset_index(0)
        # keep all that dont have overlaps
        ccomps = nodes.comp.unique()
        components[id] = list(set(ccomps) - set(comps.index))
        log(f"{id}: Adding {len(components[id])} without overlaps.")
        # loop over components and allocate
        for c, (did, ocomps) in comps.iterrows():
            ignore.setdefault(did, [])
            # if neighbour has more than 1 component, then this domain will have the complete component
            if len(ocomps) > 1:
                components[id].append(c)
                #ignore[did].extend(ocomps)
                log(f"{id}: Adding {c}, ignoring {ocomps} of {did}")
                continue
            # check which comp is bigger
            oc = ocomps[0]
            cnodes = nodes[nodes.comp == c]
            ocnodes = other_nodes.loc[did][other_nodes.loc[did, "comp"] == oc]
            # if current component has more or the same nodes, keep it
            # Since the neighbour comp is ignored, it's not allocated twice
            if len(cnodes) >= len(ocnodes):
                components[id].append(c)
                ignore[did].append(oc)
                log(f"{id}: Adding {c}, ignoring {oc} of {did}")
        print(f"{id}: {len(components[id])}/{len(ccomps)}")
    return components, ignore


def write_out(dictoflists, outfilepattern):
    print("Writing out ignored components...")
    for id, cps in dictoflists.items():
        if cps:
            with open(outfilepattern.format(domain=id), "w") as f:
                f.write("\n".join(map(str, cps)))
    return


def plot(nodes_all, ignore, output=None):
    import pylab as pl

    fig, ax = pl.subplots()
    for id, nds in nodes_all.items():
        if id in ignore:
            nds.drop(nds.index[nds.comp.isin(ignore[id])], inplace=True)
        nds.plot(ax=ax, markersize=1)
    if output:
        fig.savefig(output, dpi=300)
    return fig


def run(domain_polygons_file, domains_selected_file, node_file_dir, output_dir):
    domains = read_domains(domain_polygons_file, domains_selected_file)
    nodes_all = read_nodes(domains, os.path.join(node_file_dir, "{domain}.gpkg"))
    components, ignore = allocate_components(nodes_all, domains)
    write_out(ignore, os.path.join(output_dir, "{domain}.txt"))
    plot({k: nodes_all[k] for k in components}, ignore, "tmp.png")


if __name__ == "__main__":
    from commandline import interface as _cli
    _cli()
