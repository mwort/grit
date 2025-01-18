#!/usr/bin/env python
import os.path as osp
import os
import shutil
import warnings

import numpy as np
import pandas as pd
from tqdm import tqdm

from rivgraph.networks import RiverNetwork


def prune_single_flow_network(nodes, links):
    # prune single flow reaches and check the direction by accumulation
    sfd_in = [None]
    flips = []
    order = 1
    while True:
        # inlet nodes and lines with single flow
        sfd_inl = nodes[(nodes.n_lines == 1) & nodes.single_flow & ~nodes.outlet]
        sfd_lns = links.loc[sfd_inl.conn.apply(lambda x: x[0])]

        # receiving node also needs to be a single flow node, does it though?
        # No, because if the inlet is single flow it can only go one way into the next node, even if that is not single flow
        # Yes, because it prevents outlet only components after pruning in coastal areas
        nxt_nodes = nodes.loc[
            np.where(sfd_lns.end_node == sfd_inl.index, sfd_lns["start_node"], sfd_lns["end_node"])
        ]
        sfd_inl = sfd_inl[(nxt_nodes.single_flow).values]
        sfd_lns = sfd_lns[(nxt_nodes.single_flow).values]

        # make sure inlet - outlet links have the right direction (includes outlet - inlet lines)
        flip = (sfd_lns.end_node == sfd_inl.index)
        flips.extend(sfd_lns.index[flip].tolist())
        if len(sfd_inl) == 0:
            break
        # recalculate next nodes
        nxt_nodes = nodes.loc[
            np.where(sfd_lns.end_node == sfd_inl.index, sfd_lns["start_node"], sfd_lns["end_node"])
        ]
        # recalculate node connections
        for n, l in zip(nxt_nodes.index, sfd_lns.index):
            if l in nodes.loc[n, 'conn']:
                nodes.loc[n, 'conn'].remove(l)
            else:
                print(nodes.loc[n])
                print(links.loc[l])
                raise RuntimeError(f"line {l} is not connected to node {n}")
        nodes.loc[nxt_nodes.index, 'n_lines'] = nodes.loc[nxt_nodes.index, 'conn'].apply(lambda c: len(c))
        # drop from nodes and lines
        nodes.drop(sfd_inl.index, inplace=True)
        nlines0 = nodes.loc[nxt_nodes.index, 'n_lines'] == 0
        nodes.drop(nxt_nodes.index[nxt_nodes.outlet & nlines0], inplace=True)
        links.drop(sfd_lns.index, inplace=True)
        print(f'Pruning {len(sfd_lns)} single flow lines of order {order}')
        order += 1

    # set new inlets
    nodes.loc[(nodes.n_lines == 1) & ~nodes.outlet, "inlet"] = True
    # remove lines that have undefined start/end nodes
    undefined = (links[["start_node", "end_node"]] <= 0).any(axis=1)
    if undefined.sum():
        warnings.warn(f"Pruning {undefined.sum()} lines with undefined start/end nodes.")
        print(links[undefined])
        lns = undefined.index[undefined]
        nds = links.loc[undefined, ["start_node", "end_node"]].values.flatten()
        for n in set(nds[nds > 0]):
            nodes.at[n, "conn"] = list(set(nodes.loc[n, "conn"]) - set(lns))
            if len(nodes.loc[n, "conn"]) == 0:
                nodes.drop(n, inplace=True)
        links.drop(lns, inplace=True)
    return nodes, links, flips


def multi_network_rivgraph(links, nodes, idx_shape, results_folder='.', res=30):
    """Run RivgraphNetwork on multiple networks.
    """
    os.makedirs(results_folder, exist_ok=True)

    print(f"Reading {nodes}")
    nodes = pd.read_csv(nodes, index_col=0) if type(nodes) == str else nodes
    if nodes.isna().sum().max() > 0:
        warnings.warn("Ignoring some nodes with NAs: N\n %s" % (nodes.isna().sum()))
        nodes.dropna(inplace=True)
    nodes['conn'] = nodes.conn.apply(eval)
    print(f"Reading {links}")
    links_dir = osp.dirname(links)
    links = pd.read_csv(links, index_col=0) if type(links) == str else links
    if links[["start_node", "end_node", "idx"]].isna().sum().max() > 0:
        warnings.warn("Ignoring some links with NAs: N\n %s" % 
            (links[["start_node", "end_node", "idx"]].isna().sum()))
        links.dropna(inplace=True, subset=["start_node", "end_node", "idx"])
    links['idx'] = links.idx.apply(eval)
    links['wid_pix'] = links.wid_pix.apply(eval)


    # get rid of single flow path segments and calculate their directions based on accumulation
    nodes, links, flips = prune_single_flow_network(nodes, links)

    # in the rare case that all links were set by single flow
    if len(nodes) == 0 and len(links) == 0:
        write_flipped_link_ids(flips, osp.join(results_folder, 'flipped_links.csv'))
        print("All links set by single flow. Only flips but no link or node attributes written out.")
        return

    # sanity check input
    #  segments with less than 2 pixels?
    links_length = links.idx.apply(lambda s: len(s))
    assert links_length.min() >= 2
    #  no inlet also assigned as outlet
    assert len(set(nodes.index[nodes.inlet]) & set(nodes.index[nodes.outlet])) == 0
    # invalid start/end nodes
    assert (links[["start_node", "end_node"]] < 1).sum().sum() == 0, "Invalid start/end nodes."
    #  all node idx are at start/end of line idx
    assert (nodes.loc[links.start_node.values, "idx"].values == links['idx'].apply(lambda l: l[0])).all()
    assert (nodes.loc[links.end_node.values, "idx"].values == links['idx'].apply(lambda l: l[-1])).all()
    # no loops
    hasloops = links.start_node == links.end_node
    if (hasloops).sum():
        warnings.warn(f"Disregarding {hasloops.sum()} loops.")
        loopnodes = links[hasloops]["start_node"]
        for ln, nd in loopnodes.items():
            nodes.loc[nd, "conn"] = list(set(nodes.loc[nd, "conn"]) - set([ln]))
        links = links[~hasloops]
    # all components should have at least one inlet and one outlet
    n_inlets = nodes.groupby("component").inlet.sum()
    assert (n_inlets > 0).all(), f"These components dont have inlets: {n_inlets.index[n_inlets == 0]}"
    assert (nodes.groupby("component").outlet.sum() > 0).all(), "Some components dont have an outlet"

    # loop over all components
    nix, lix = nodes.groupby("component").groups, links.groupby("component").groups
    assert len(nix) == len(lix)
    print(f"Found {len(nix)} networks.")
    links_info, nodes_info = [], []
    # iterate over components starting with the smallest
    for n, _ in sorted(nix.items(), key=lambda item: len(item[1])):
        name = "component_%04i" % n
        fixedlinksfile = name+"_fixlinks.csv"
        # put link fixes where rivgraph will pick them up
        if osp.exists(osp.join(links_dir, fixedlinksfile)):
            shutil.copy(osp.join(links_dir, fixedlinksfile), osp.join(results_folder, fixedlinksfile))
        rgn = RiverNetwork(links.loc[lix[n]], nodes.loc[nix[n]], idx_shape,
                           results_folder=results_folder, res=res, name=name)
        rgn.run(output=False)
        links_info.append(rgn.links_direction_info.copy())
        nodes_info.append(rgn.nodes_direction_info.copy())
        flips.extend(rgn.flipped_links)
    linfo = pd.concat(links_info)
    linfo.loc[:, "component"] = links["component"]
    ninfo = pd.concat(nodes_info)
    # report cycles
    gb = linfo[linfo.cycle>0].groupby(["component", "cycle"]).certain_alg
    cygb = gb.agg(lambda x: (x == 29).sum())
    badln = linfo.set_index(["component", "cycle"]).loc[cygb.index]
    print("Cycle lines without randomly assigned directions in cycle:")
    print(badln["certain_alg"].to_string())
    print(f"N cycles: {len(cygb)}, N links: {len(badln)}, N discontinuity nodes: {(ninfo.continuity_violated == 1).sum()}")
    # write link/node info
    linfo.to_csv(osp.join(results_folder, 'links_direction_info.csv'))
    ninfo.to_csv(osp.join(results_folder, 'nodes_direction_info.csv'))
    write_flipped_link_ids(flips, osp.join(results_folder, 'flipped_links.csv'))
    return


def write_flipped_link_ids(flips, output):
    with open(output, 'w') as f:
        f.writelines(['%s\n' % i for i in flips])
    return


if __name__ == "__main__":
    from commandline import interface as _cli
    _cli()
    # import pstats
    # import cProfile

    # cProfile.run("_cli()", f"my_func_stats_{os.getpid()}")

    # p = pstats.Stats(f"my_func_stats_{os.getpid()}")
    # p.sort_stats("cumulative").print_stats("evoflood/src")

