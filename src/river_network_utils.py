#!/usr/bin/env python
import os.path as osp
import sys
import gc
from collections import defaultdict
import warnings

import numpy as np
import pandas as pd
from tqdm import tqdm
import grass.script as grass
from grass_attribute_table import GrassAttributeTable
import networkx as nx

# pandas apply progress reporting
tqdm.pandas()


def read_line_segments(vector):
    """Read a lines vector coordinates (np.array) into a pd.Series."""

    tmpf = grass.tempfile(False)
    grass.run_command("v.out.ascii", input=vector, format="standard", output=tmpf)
    cats, coords = [], []
    with open(tmpf) as tf:
        [tf.readline() for _ in range(10)]
        while True:
            nxt = tf.readline()
            if not nxt:
                break
            typ, n, nl = nxt.strip().split()
            if typ != "L":
                [tf.readline() for _ in range(int(n)+int(nl))]
                continue
            line = [tf.readline().strip().split() for _ in range(int(n))]
            il, cat = tf.readline().strip().split()
            cats.append(int(cat))
            try:
                coords.append(np.array(line, dtype=int))
            except ValueError:
                coords.append(np.array(line, dtype=float))
    return pd.Series(coords, index=cats)


def clean_network(network_vector, output=None, node_min_length_column='river_halfwidth',
                  node_lake_column=None, line_length_column="length", factor=1,
                  lake_factor=1, min_accum=None, node_accum_column="darea",
                  node_outlet_column=None, cleaned_cats_file=None,
                  preserve_longest=True, ncycles=1):

    lines = network_lines(network_vector, attr=True)
    nodes = network_nodes(network_vector, attr=True)
    nodes['min_len'] = nodes[node_min_length_column]
    if node_lake_column:
        nodes['is_lake'] = nodes[node_lake_column].astype(bool)
    if node_accum_column and min_accum:
        nodes["accum"] = nodes[node_accum_column]
    if node_outlet_column:
        nodes['is_outlet'] = nodes[node_outlet_column].astype(bool)
    drop_lines = []
    # remove loops
    loops = lines[lines.start_node == lines.end_node]
    nodes.loc[loops.start_node, "n_lines"] -= 1
    drop_lines.extend(list(loops.index))
    length = lines[line_length_column]
    print(f'Removing {len(loops)} loops.')

    cycle = 0
    while True:
        # get lines who have only one line (itself) attached at either nodes
        danglenodes = nodes[nodes.n_lines == 1]
        # only continue if unvisited dangle nodes remain
        if len(danglenodes) == 0 or cycle >= ncycles:
            break
        # get line ids of dangles,
        dangles_lines = danglenodes.line_cats.apply(lambda l: l[0])
        # exclude dangles that have dangle nodes on both sides
        #dangles = set(danglenodes.line_cats.apply(lambda l: l[0]))
        dangle_n_nodes = dangles_lines.groupby(dangles_lines).count()
        dangles = set(dangle_n_nodes[dangle_n_nodes == 1].index)
        # dangle line ids at each node
        danglelineids = nodes.line_cats.apply(lambda l: list(dangles.intersection(set(l))))
        danglelineids = [np.array([[i]*len(l), l]).T for i, l in danglelineids.items() if len(l) > 0]
        nl = pd.DataFrame(np.concatenate(danglelineids), columns=['node', 'line'])
        nl['min_len'] = nodes.loc[nl.node, 'min_len'].values
        maxlen = nl.loc[nl.groupby('line').min_len.idxmax().values].set_index("line")
        minl = maxlen["min_len"]
        llen = length[minl.index]
        tooshort = llen < minl*factor
        if node_lake_column:
            nl['is_lake'] = nodes.loc[nl.node, 'is_lake'].values
            islake = nl.groupby("line").is_lake.all()
            tooshort = tooshort | (islake & (llen < minl*lake_factor))
        # Preserve shortest outlet (also if just single)
        maxlen["length"] = length.loc[maxlen.index]
        if node_outlet_column:
            nl['is_outlet'] = nodes.loc[nl.node, 'is_outlet'].values
            isoutl = nl.groupby("line").is_outlet.any()
            ngb = maxlen[tooshort & isoutl].groupby("node")
            shortest_lines = ngb.length.idxmin().values
            tooshort[shortest_lines] = False
        # if there are more than one dangle that is too short from the same node, leave the longest
        if preserve_longest:
            ngb = maxlen[tooshort].groupby("node")
            longest_lines = ngb.length.idxmax()[ngb.length.count() > 1].values
            tooshort[longest_lines] = False
        print(f"Removing {tooshort.sum()} too short dangles.")
        drops = set(llen.index[tooshort].tolist())
        # Preserve dangles that have enough accumulation at start
        if min_accum and node_accum_column:
            minacn = danglenodes[danglenodes["accum"] >= min_accum]
            # but dont preserve dangles if they are only in lakes
            if node_lake_column:
                minacn = minacn[~minacn["is_lake"]]
            minacl = set(minacn.line_cats.apply(lambda n: n[0]).tolist())
            print(f"Preserving {len(minacl & drops)} dangles that have >= {min_accum} accumulation.")
            drops = drops - minacl
        if len(drops) == 0:
            break
        drop_lines.extend(drops)
        # remove all visited dangle nodes and update lines
        nodes.drop(index=danglenodes.index, inplace=True)
        nodes['line_cats'] = nodes.line_cats.apply(lambda ls: list(set(ls) - drops))
        nodes['n_lines'] = nodes.line_cats.apply(len)
        cycle += 1
    if output:
        print(f'Will remove {len(drop_lines)} ({int(len(drop_lines)*100/len(lines))}%) lines.')
        drop_cats(network_vector, output, drop_lines, id_file=cleaned_cats_file)
    return drop_lines


def clean_short_dangles(network_vector, output, min_length, node_outlet_column="is_outlet",
                        line_length_column="length"):
    """
    Remove short dangles that are not outlets and preserve longest at nodes with multiple dangles.
    """
    nodes = network_nodes(network_vector, attr=True)
    lines = network_lines(network_vector, attr=True)
    # pick dangles nodes but not outlets
    dangle_nodes = nodes[(nodes.n_lines == 1) & (~nodes[node_outlet_column].astype(bool))]
    dangle_lines = lines.loc[dangle_nodes.line_cats.apply(lambda c: c[0])]
    # make sure start_node is dangle_node
    stnd = dangle_lines.start_node.copy()
    flip = stnd != dangle_nodes.index
    dangle_lines.loc[flip, "start_node"] = dangle_lines.loc[flip, "end_node"]
    dangle_lines.loc[flip, "end_node"] = stnd[flip]

    short_dangles = dangle_lines[dangle_lines[line_length_column] < min_length]
    # make sure longest dangle remains where two short dangles exist
    sdgb = short_dangles.groupby("end_node")
    threeway = (sdgb.start_node.count() > 1) & (nodes.loc[sdgb.count().index, "n_lines"] == 3)
    keep = sdgb[line_length_column].idxmax()[threeway].tolist()
    nkeep = len(keep)
    short_dangles.drop(keep, inplace=True)
    # prevent loops to be created downstream
    dslns = nodes.loc[short_dangles.end_node.values, "line_cats"]
    drop = []
    for en, lns in dslns.items():
        nxtnds = set(lines.loc[lns, ["start_node", "end_node"]].values.flatten())
        nxtnds.remove(en)
        drop.append(len(nxtnds) == len(lns))
        if len(nxtnds) != len(lns):
            nkeep += 1
    short_dangles = short_dangles[np.array(drop)]
    print(f"Removing {len(short_dangles)} short dangles, preserving {nkeep}.")
    drop_cats(network_vector, output, short_dangles.index)
    return


def clean_small_loops(network_vector, output, factor=2, node_width_column="river_width",
                      line_length_column="length"):
    """Remove the longest line of direct line loops if both lines are shorter
    than node_width * factor.
    """
    nodes = network_nodes(network_vector, attr=True)
    lines = network_lines(network_vector, attr=True)
    # make sure start/end nodes are sorted
    sen = lines["start_node"].copy()
    flip = sen > lines.end_node
    lines.loc[flip, "start_node"] = lines.loc[flip, "end_node"]
    lines.loc[flip, "end_node"] = sen.loc[flip]
    # identify all loops
    dups = (lines[lines[["start_node", "end_node"]]
            .duplicated(keep=False)].copy()
            .sort_values(["start_node", "end_node"]))
    dups["mean_node_width"] = np.mean(
        (nodes.loc[dups.start_node, node_width_column],
         nodes.loc[dups.end_node, node_width_column]),
        axis=0)
    bad = dups[dups[line_length_column] < dups.mean_node_width * factor]
    badgb = bad.groupby(["start_node", "end_node"])[line_length_column]
    # only take the longest if both lines are short
    longestbad = badgb.idxmax()[badgb.count() > 1].values
    print(f"Removing {len(longestbad)} lines from small loops.")
    drop_cats(network_vector, output, longestbad)
    return


def read_raster_xyz(rasters):
    tmpf = grass.tempfile(False)
    grass.run_command("r.out.xyz", input=rasters, separator="comma", output=tmpf)
    kw = dict(names=["x", "y"]+(rasters.split(",") if type(rasters)==str else list(rasters)))
    return pd.read_csv(tmpf, header=None, index_col=(0, 1), **kw)


def dense_lines(line_vector, res=30, check_missing=None):
    """Creates lines from line_vector that are densified by the raster resolution.
    lines_vector : str | pd.Series with coordinates
    Returns pd.Series
    """
    undense = read_line_segments(line_vector) if type(line_vector) == str else line_vector
    dense = undense.copy()
    for i, line in tqdm(undense.items(), "Densify lines", total=len(undense)):
        gaps = np.linalg.norm(line[1:] - line[:-1], axis=1)
        # .round needed to avoid floating point problem with % in large gaps
        diag = (gaps / np.sqrt(res**2 * 2)).round(14)
        bdiag = np.all([diag > 1, diag % 1 == 0], axis=0)
        straight = gaps / res
        bstr = np.all([straight > 1, straight % 1 == 0], axis=0)
        # line arcs not complying with grid
        bodd = (straight % 1 != 0) & (diag % 1 != 0)
        ncells = np.ones_like(diag)
        ncells[bdiag] = diag[bdiag]
        ncells[bstr] = straight[bstr]
        newline = [np.array([line[0]])]
        stx, sty = line[0]
        for xy, nc in zip(line[1:], ncells.round(0).astype(int)):
            # raster conform arc
            x, y = np.linspace(stx, xy[0], int(nc)+1)[1:], np.linspace(sty, xy[1], int(nc)+1)[1:]
            newline.append(np.column_stack([x, y]))
            stx, sty = xy
        dense[i] = list(map(tuple, np.concatenate(newline, dtype=line.dtype, casting="unsafe")))
    del undense
    gc.collect()
    # check if coordinates are in predefined index
    if check_missing is not None:
        nunchckd = dense.apply(len).sum()
        # make sure all dense vertices are in rast index and missing ends/starts are consistent with nodes
        rix = set(check_missing)
        missing = defaultdict(lambda: [])
        for i, line in dense.items():
            mis = set(line) - rix
            if len(mis) > 0:
                for v in mis:
                    if line.index(v) in [0, len(line)-1]:  # at start or end
                        missing[i].append(line.index(v))
                    line.remove(v)
                    print(f"Missing vertex at {v} in line {i}")

        nnotin = nunchckd - dense.apply(len).sum()
        print('%1.3f%% (%i/%i) of vertices not in input rasters.' % (nnotin*100/nunchckd, nnotin, nunchckd))
        return dense, missing
    return dense


def network_nodes(network_vector, attr=False):
    conn = grass.read_command("v.net", input=network_vector, operation="nreport")
    lines = [i.split() for i in conn.split("\n")[:-1]]
    nds = pd.Series({int(i): list(map(int, c.split(","))) for i, c in lines}).to_frame("line_cats")
    if attr:
        nds = nds.join(GrassAttributeTable(network_vector, layer=2))
    if not (attr and "n_lines" in nds.columns):
        nds["n_lines"] = nds.line_cats.apply(len)
    return nds


def network_lines(network_vector, attr=False):
    conn = grass.read_command("v.net", input=network_vector, operation="report")
    lines = pd.DataFrame([i.split() for i in conn.split("\n")[:-1]], columns=["id", "start_node", "end_node"])
    lines = lines.astype(int).set_index("id")
    if attr:
        lines = lines.join(GrassAttributeTable(network_vector, layer=1))
    return lines


def river_networkx(network_vector, directed=True, edge_attributes=["id"]):
    import networkx as nx

    lines = network_lines(network_vector, attr=True)
    lines['id'] = lines.index
    G = nx.from_pandas_edgelist(lines, 'start_node', 'end_node',
                                edge_attr=edge_attributes,
                                create_using=nx.MultiDiGraph if directed else nx.MultiGraph)
    return G


def upload_line_counts(network_vector, column='n_lines'):
    nodes = GrassAttributeTable(network_vector, layer=2)
    nodes[column] = network_nodes(network_vector).n_lines
    nodes.write()
    return


def snap_network_ends(network_vector, output, snap_min_max=(60, 120), cleaning_factor=2,
                      node_min_length_column="river_halfwidth", line_length_column="length",
                      min_loop_length=1000):
    """Snap ends within threshold (only in rivers) if the network path traverses more than 1 node.
    Creates several intermediate files with *__* in name.
    """
    import grass.script as gs
    from grass.pygrass.vector import VectorTopo
    from grass.pygrass.vector.geometry import Line

    gs.run_command("g.copy", vector=[network_vector, output])

    upload_line_counts(output)
    gs.run_command("v.extract", input=output, layer=2,
                   where="n_lines=1 AND grwl_value in (86, 126, 255)", output="line__ends")
    gs.run_command("v.db.addcolumn", map=output, layer=2,
                   column="dangle2dangle_snap_dist double,dangle2dangle_snap_cat int")
    gs.run_command("v.distance", from_=output, to="line__ends", from_layer=2, to_layer=2,
                   dmin=snap_min_max[0], dmax=snap_min_max[1], upload="cat,dist",
                   col="dangle2dangle_snap_cat,dangle2dangle_snap_dist")
    gs.run_command("v.to.db", map=output, option="coor", columns="x,y", layer=2)

    nodes = GrassAttributeTable(output, layer=2)
    en = nodes[(nodes.n_lines == 1) & (~nodes.dangle2dangle_snap_cat.isna())]
    fromto = {min(f, t): max(f, t) for f, t in en.dangle2dangle_snap_cat.astype(int).items()}
    nodes.loc[:, 'dangle_path_cost'] = 1
    nodes.loc[fromto.keys(), 'dangle_path_cost'] = 0
    nodes.loc[fromto.values(), 'dangle_path_cost'] = 0
    nodes.write()
    lines = GrassAttributeTable(output, layer=1)
    lines.loc[:, 'dangle_path_cost'] = 0
    lines.write()

    # get paths
    tmpf = gs.tempfile(False)
    with open(tmpf, 'w') as f:
        f.write('\n'.join([f'{f} {f} {t}' for f, t in fromto.items()]))
    # cost = length along network
    gs.run_command('v.net.path', input=output, output='links__paths', file=tmpf)
    # cost = number of nodes along network
    gs.run_command('v.net.path', input=output, output='links__paths__nnodes', file=tmpf,
                   node_column='dangle_path_cost', arc_column='dangle_path_cost')
    paths = GrassAttributeTable('links__paths').rename({"cost": "length"}, axis=1)
    paths['nnodes'] = GrassAttributeTable('links__paths__nnodes')["cost"]
    paths["loop_length"] = paths['length'] + nodes.loc[paths.fcat, "dangle2dangle_snap_dist"].values
    # check if any of the paths start/end lines would remain after cleaning
    non_viable = clean_network(output, factor=cleaning_factor, preserve_longest=False,
                               node_min_length_column=node_min_length_column, line_length_column=line_length_column)
    node_line_cats = network_nodes(output)["line_cats"].apply(lambda l: l[0])
    paths["has_viable_line"] = [
        any([f not in non_viable, t not in non_viable])
        for f, t in zip(node_line_cats[paths.fcat], node_line_cats[paths.tcat])
    ]
    # snap conditions: no network connection
    nocon = (paths.sp == 1)
    print(f"Snapping {nocon.sum()} ends without network path.")
    # nodes traversed > 1 and min length
    minna = ((paths.sp == 0) & (paths.nnodes > 1) & (paths["loop_length"] > min_loop_length))
    print(f"Snapping {minna.sum()} ends with loops >{min_loop_length}m and more than 1 node.")
    paths_valid = paths[nocon | (minna & paths.has_viable_line)]
    gs.message(f"Found {len(paths_valid)} snapping points.")
    # write out connections to patch
    gs.run_command("v.db.droptable", map=output, layer=1, flags="f")
    gs.run_command("v.db.addtable", map=output, layer=1, columns="is_valid int")
    line_vect = VectorTopo(output)
    line_vect.open('rw', tab_cols=[(u'cat', 'INTEGER PRIMARY KEY'), (u'is_valid', 'int')])
    cats = lines.index.max()+paths_valid.index
    for c, (fc, tc) in zip(cats, paths_valid[['fcat', 'tcat']].values):
        fxy, txy = nodes.loc[fc, ['x', 'y']], nodes.loc[tc, ['x', 'y']]
        dist, dird = (fxy - txy).abs().astype(float), (fxy - txy).astype(float)
        if (dist.min() == 0) or (dist.min() == dist.max()):
            # streight or 45 diagonal line
            lnvrt = nodes.loc[[fc, tc], ['x', 'y']].values
        else:
            # non-45dg diagonal line
            print(dist)
            long, short = dist.idxmax(), dist.idxmin()
            interm_vrt = fxy - dist[short] * dist/dird
            lnvrt = [tuple(fxy), tuple(interm_vrt), tuple(txy)]
        ln = Line(lnvrt)
        line_vect.write(ln, cat=c, attrs=(1,))
    line_vect.table.conn.commit()
    line_vect.close()
    return


def prune_triangles(short_line_network, network_vector, output, node_darea_column="darea",
                    line_length_column="length", node_component_column="comp"):
    """Prune links of small triangles that are the result of merging the
    skeletonised mask and r.watershed streams.

    Removes the link connecting the nodes with the two lower drainage areas.
    """
    lines = GrassAttributeTable(short_line_network, layer=1)
    nlines = network_lines(short_line_network)
    nlines['length'] = lines[line_length_column]
    nodes = GrassAttributeTable(short_line_network, layer=2)
    ncounts = nodes.groupby(node_component_column).count()[node_darea_column]
    nodes.loc[:, "n_in_comp"] = ncounts[nodes[node_component_column]].values
    nodes.write()
    nnodes = network_nodes(short_line_network)

    # remove the dent in 2-node triangles
    lines_remove = []
    for c, ft in nodes[nodes.n_in_comp == 2].groupby(node_component_column).groups.items():
        lns = set(nnodes.line_cats[ft].sum())
        if len(lns) == 2:
            lines_remove.append(lines.loc[lns, line_length_column].idxmax())
        else:
            nodes.loc[ft, 'n_in_comp'] = 4
    # triangles with confluence, remove link btw nodes with least darea
    for c, nds in nodes[nodes.n_in_comp == 3].groupby(node_component_column).groups.items():
        lns = set(nnodes.line_cats[nds].sum())
        if len(lns) == 3:
            # 2 nodes with lowest darea
            nds_darea = nodes.loc[nds, node_darea_column].sort_values().index[:2]
            lns = nnodes.loc[nds_darea, "line_cats"].tolist()
            lines_remove.append(list(set(lns[0]) & set(lns[1]))[0])
        else:
            nodes.loc[nds, 'n_in_comp'] = 4
    # network approach on all other clusters, i.e. only preserve lines that belong to shortest path to outlet
    for c, nds in nodes[nodes.n_in_comp > 3].groupby(node_component_column).groups.items():
        lns = nlines.loc[set(nnodes.line_cats[nds].sum())]
        ndis = nodes.loc[nds]
        ton = [nodes.loc[nds, node_darea_column].idxmax()]
        frn = set(nds) - set(ton)
        # get all shortest lines to largest darea node
        paths = shortest_paths(lns, frn, ton, return_lines=True)
        obsolete = set(lns.index) - set(l for p in paths.values() for l in p)
        lines_remove.extend(obsolete)

    print(f"Removing {len(lines_remove)} out of {len(lines)} of triangle lines.")
    drop_cats(network_vector, output, lines_remove)
    return


def drop_cats(input, output, cats, id_file=None):
    # remove lines by inverse extraction
    tmpf = id_file or grass.tempfile(False)
    with open(tmpf, "w") as f:
        f.write("\n".join(map(str, cats)))
    grass.run_command("v.extract", flags="r", input=input, output=output, file=tmpf)
    return


def shortest_paths(lines, from_nodes, to_nodes, weight_column='length', return_lines=False):
    import networkx as nx
    from itertools import product
    lines['id'] = lines.index
    G = nx.from_pandas_edgelist(lines, 'start_node', 'end_node', edge_attr=[weight_column, 'id'])
    paths = {}
    for f, t in product(from_nodes, to_nodes):
        paths[(f, t)] = nx.dijkstra_path(G, f, t, weight=weight_column)
        if return_lines:
            paths[(f, t)] = [G[u][v]['id'] for u, v in zip(paths[(f, t)][0:-1], paths[(f, t)][1:])]
    return paths


def remove_cycles(network_vector, output, node_darea_column='darea',
                  line_length_column='length', drop_lines_file=None):
    import networkx as nx

    lines = network_lines(network_vector)
    lines['id'] = lines.index
    lines['length'] = GrassAttributeTable(network_vector, layer=1)[line_length_column]
    G = nx.from_pandas_edgelist(lines, 'start_node', 'end_node', edge_attr=['id'])
    nodes = network_nodes(network_vector)
    nodes['darea'] = GrassAttributeTable(network_vector, layer=2)[node_darea_column]
    lines_remove = []
    # find all cycles in undirected graph
    cycle_graph = nx.Graph(set(G.edges)-set(nx.bridges(G)))
    cycle_comps = [G.subgraph(c).copy() for c in nx.connected_components(cycle_graph)]
    # find all paths to cycles largest node
    for g in tqdm(cycle_comps, 'Breaking cycles'):
        nds = set(g.nodes)
        outl = set([nodes.loc[nds, 'darea'].idxmax()])
        lns = lines.loc[[G[u][v]['id'] for u, v in g.edges]]
        frn = set(nds) - outl
        # get all shortest lines to largest darea node
        paths = shortest_paths(lns, frn, outl, return_lines=True)
        obsolete = set(lns.index) - set(l for p in paths.values() for l in p)
        lines_remove.extend(obsolete)
    # add parallel edges that are not captured by the above method
    lines[['start_node', 'end_node']] = np.sort(lines[['start_node', 'end_node']])
    # sort by length to keep the shortes duplicated (first)
    lines.sort_values("length", inplace=True)
    dups = lines[lines[['start_node', 'end_node']].duplicated()]
    lines_remove.extend(dups.id.to_list())
    print(f'Removing {len(lines_remove)} lines of {len(cycle_comps)} cycles'
          f' and {len(dups)} parallels.')
    drop_cats(network_vector, output, lines_remove, id_file=drop_lines_file)
    return


def valid_network_components(segments_vector, valid_components=None,
                             invalid_components=None, node_type_column="type",
                             component_column="comp", grwl_column="grwl_value",
                             max_coastal_node_proportion=0.95):
    # connections
    nodes = network_nodes(segments_vector)
    nodestbl = GrassAttributeTable(segments_vector, layer=2)
    nodes["outlet"] = nodestbl[node_type_column] == "outlet"
    # outlets must be end points
    #nodes.loc[(nodes.n_lines > 1) & nodes.outlet, "outlet"] = False
    # all other nodes with single line are inlets
    nodes["inlet"] = False
    nodes.loc[(nodes.n_lines == 1) & (~nodes.outlet), "inlet"] = True
    print(f"Found {nodes.outlet.sum()} outlets and {nodes.inlet.sum()} inlets.")
    # coastal only components
    coprop = nodestbl.groupby(component_column)[grwl_column].agg(lambda x: (x == 126).sum()/len(x))
    # identify network components that dont have outlets
    nodes['comp'] = nodestbl[component_column]
    comps = nodes.groupby('comp')[['inlet', 'outlet']].any()
    comps["not_coastal"] = coprop <= max_coastal_node_proportion
    validcomps = comps.all(axis=1)
    print(f'Found {(~comps.inlet).sum()} network components without inlets.')
    print(f'Found {(~comps.outlet).sum()} network components without outlets.')
    print(f'Found {(~comps.not_coastal).sum()} network components with mostly coastal nodes.')
    print(f'Found {validcomps.sum()} network components with both.')

    for outc, cmps in [(invalid_components, validcomps), (valid_components, ~validcomps)]:
        if outc:
            grass.run_command('g.copy', vector=[segments_vector, outc])
            prune_network_components(outc, component_column, ids=validcomps.index[cmps])
    return


def network_components_in_domain(network_vector, domain_vector, component_column="comp",
    in_domain_column="in_domain", type_column="type", output=None):
    """Select components from network if majority of nodes is inside domain.

    Nodes of type='outlet' are ignored to avoid 2-node coastal components to be
    deselected.
    """
    # check which nodes are in domain
    grass.run_command("v.db.addcolumn", map=network_vector,
        column=f"{in_domain_column} int", layer=2)
    grass.run_command("v.db.update", map=network_vector,
        column=in_domain_column, layer=2, value=0)
    grass.run_command("v.what.vect", map=network_vector, layer=2,
        column=in_domain_column, query_map=domain_vector, query_col="cat")

    nodes = GrassAttributeTable(network_vector, layer=2)
    lines = GrassAttributeTable(network_vector, layer=1)

    # make sure any domain cat is 1
    nodes[in_domain_column] = nodes[in_domain_column].astype(bool).astype(int)

    # avoid two-node components at coast to lead to 50/50
    nsel = nodes[nodes[type_column] != 'outlet'].dropna(subset=[component_column])

    comp = nsel.groupby(component_column)
    nincomp = comp[in_domain_column].count()
    ratio = comp[in_domain_column].sum()/nincomp

    if len(ratio[ratio==0.5]) > 0:
        warnings.warn("Some components are half in, half outside of domain:\n %s" % (nincomp[ratio==0.5]))

    comp_sel = ratio[ratio >= 0.5].index
    grass.message(f"Found {len(comp_sel)}/{len(ratio)} components.")

    nodes[in_domain_column] = nodes[component_column].isin(comp_sel).astype(int)
    lines[in_domain_column] = lines[component_column].isin(comp_sel).astype(int)
    nodes.write()
    lines.write()

    if output:
        grass.run_command('g.copy', vector=[network_vector, output])
        comps_prune = list(set(ratio.index) - set(comp_sel))
        if comps_prune:
            prune_network_components(output, component_column, ids=comps_prune)
    return


def prune_network_components(network_vector, column, ids=None, ids_file=None):
    """Removes parts of network via attribute selection in both the nodes and lines table.
    
    `column` needs to exist in both tables.
    """
    assert ids is not None or ids_file
    if ids is None:
        with open(ids_file) as f:
            ids = list(map(int, f.read().split()))

    for i in (1, 2):
        tbl = GrassAttributeTable(network_vector, layer=i)
        tbl["remove__tmp"] = tbl[column].isin(ids).astype(int)
        tbl.write()
        grass.run_command('v.edit', map=network_vector, tool='delete', layer=i, where="remove__tmp=1")
        grass.run_command('db.execute', sql=f"DELETE FROM {tbl.table} WHERE remove__tmp=1;")
        grass.run_command("v.db.dropcolumn", map=network_vector, layer=i, columns="remove__tmp")
    return


def rivgraph_nodes_input(segments_vector, idx_raster, output=None,
                         accumulation_column="darea", type_column="type",
                         grwl_column='grwl_value', components_column='comp',
                         synthetic_elevation_column="outlet_distance",
                         mainchannel_darea_column="mainchannel_darea"):

    idx = grass.read_command("v.what.rast", map=segments_vector, layer=2,
                             raster=idx_raster, flags="p")
    nodes = (pd.DataFrame([l.split("|") for l in idx.split()], columns=["id", "idx"])
             .astype(int)
             .set_index("id")
             .sort_index())
    nodestbl = GrassAttributeTable(segments_vector, layer=2)
    # connections
    nodeslines = network_nodes(segments_vector)
    nodes['conn'] = nodeslines.line_cats
    nodes['n_lines'] = nodeslines.n_lines
    nodes["outlet"] = nodestbl[type_column] == "outlet"
    nodes["accumulation"] = nodestbl[accumulation_column]
    nodes['single_flow'] = nodestbl[grwl_column].isna()
    nodes["coastal"] = nodestbl[grwl_column] == 126
    nodes['component'] = nodestbl[components_column]
    nodes['synthetic_elevation'] = nodestbl[synthetic_elevation_column]
    nodes['mainchannel_darea'] = nodestbl[mainchannel_darea_column]
    # outlets must be end points
    #nodes.loc[(nodes.n_lines > 1), "outlet"] = False
    # in single flow only components, only the outlet with the largest accumulation can be an outlet
    sglfcomps = nodes.groupby("component").single_flow.all()
    sglfcomps = sglfcomps.index[sglfcomps]
    sglfoutl = nodes[nodes.component.isin(sglfcomps) & nodes.outlet]
    nodes.loc[sglfoutl.index, "outlet"] = False
    nodes.loc[sglfoutl.groupby("component").accumulation.idxmax(), "outlet"] = True
    # all other nodes with single line are inlets
    nodes["inlet"] = False
    nodes.loc[(nodes.n_lines == 1) & (~nodes.outlet), "inlet"] = True
    # interbasin nodes should also be inlets
    nodes.loc[nodestbl[type_column] == "interbasin", "inlet"] = True

    if output:
        nodes.to_csv(output)
    return nodes


def rivgraph_links_input(segments_vector, segments_raster, idx_raster, halfwidth_raster,
                         dem_raster, accumulation_raster, nodes, output=None, default_width=30):
    """Produce the links input file for rivgraph.

    The start and end nodes are checked in the lines and added to the lines if they
    are not part of it. This happens if nodes are not located on raster cell centers. default_width
    is used if no vertices of the line has a width.
    """
    import networkx as nx

    lines = network_lines(segments_vector)
    rast = read_raster_xyz([segments_raster, idx_raster, halfwidth_raster,
                            dem_raster, accumulation_raster])
    segments, missing = dense_lines(segments_vector, res=grass.region()["nsres"], check_missing=rast.index)

    lines["wid_pix"] = segments.apply(lambda s: rast.loc[map(tuple, s), halfwidth_raster].tolist())
    lines["idx"] = segments.apply(lambda s: rast.loc[map(tuple, s), idx_raster].tolist())
    # fix missing start/ends
    nodesinfo = pd.read_csv(nodes, index_col=0)
    for lix, ix in missing.items():
        for i in ix:
            mw = np.nanmean(lines.loc[lix, "wid_pix"])
            lines.loc[lix, "wid_pix"].insert(i, default_width if np.isnan(mw) else mw)
            stn, enn = lines.loc[lix, ["start_node", "end_node"]]
            lines.loc[lix, "idx"].insert(i, nodesinfo.loc[(stn if i == 0 else enn), "idx"])
    lines["length_cells"] = lines.idx.apply(len)

    # elevation based attributes
    print("Deriving accumulation gradient...")
    accum_lines = segments.apply(lambda s: rast.loc[map(tuple, s), accumulation_raster].values)
    # positive if correct direction, negative otherwise
    lines["darea_start"] = accum_lines.apply(lambda l: np.median(l[:max(1, len(l)//2)]))
    lines["darea_end"] = accum_lines.apply(lambda l: np.median(l[ - max(1, len(l)//2):]))
    lines["darea_gradient"] = lines["darea_end"] - lines["darea_start"]
    del accum_lines
    print("Deriving elevation gradient...")
    dem_lines = segments.apply(lambda s: rast.loc[map(tuple, s), dem_raster].values)
    lines["elevation_start"] = dem_lines.apply(lambda l: np.median(l[:max(1, len(l)//2)]))
    lines["elevation_end"] = dem_lines.apply(lambda l: np.median(l[ - max(1, len(l)//2):]))
    lines["elevation_gradient"] = lines["elevation_end"] - lines["elevation_start"]
    del dem_lines, rast, segments
    gc.collect()
    # add network attributes
    linestbl = GrassAttributeTable(segments_vector, layer=1)
    lines['component'] = linestbl.comp
    G = nx.from_pandas_edgelist(lines.rename_axis('id').reset_index(), 'start_node', 'end_node', edge_attr=['id'])
    lines['bridge'] = False
    # bridges excluding dangles
    lines.loc[[G[v][u]['id'] for v, u in nx.bridges(G) if min(len(G.adj[v]), len(G.adj[u])) > 1], 'bridge'] = True
    if output:
        lines.to_csv(output)
    return lines


def flip_segments(input_segments, flipped_links, output, node_attr=None, link_attr=None):
    """Flip segments and write to new output vector. Optionally add node and link attributes.

    flipped_links, node_attr, link_attr can be both an iterable/DataFrame or a file path.
    """
    if type(flipped_links) == str:
        with open(flipped_links) as f:
            flipped_links = map(int, f.read().split())
    grass.run_command("g.copy", vector=[input_segments, output])
    lines = GrassAttributeTable(output)
    lines.loc[:, "rivgraph_flipped"] = 0
    lines.loc[flipped_links, "rivgraph_flipped"] = 1
    if link_attr:
        lattr = pd.read_csv(link_attr, index_col=0) if type(link_attr) == str else link_attr
        lines.loc[lattr.index, lattr.columns] = lattr
    lines.write()
    grass.run_command("v.edit", map=output, tool="flip", where="rivgraph_flipped=1")

    nodes = GrassAttributeTable(output, layer=2)
    if node_attr:
        nattr = pd.read_csv(node_attr, index_col=0) if type(node_attr) == str else node_attr
        nodes.loc[nattr.index, nattr.columns] = nattr
        nodes["cycle"] = nodes["cycle"].fillna(0).astype(int)
        nodes.write()

    return


def remove_directed_cycles(segments, output, cycles_vect=None):
    """
    Deal with cycles after setting directions with RivGraph.

    """
    import networkx as nx


    lines = GrassAttributeTable(segments)
    linesnet = network_lines(segments)
    nodes = GrassAttributeTable(segments, layer=2)
    nlines = network_lines(segments)


    # check lines with invalid topology (missing nodes)
    stennode = linesnet[["start_node", "end_node"]].min(axis=1)
    stendbad = linesnet.index[stennode < 0]
    if len(stendbad):
        warnings.warn(f'Found {len(stendbad)} lines that are missing either start/end nodes. Will remove them!')
        lines.loc[stendbad, "cycle"] = 1

    # find cycles independently of RivGraph (as it excludes single flow and self loops)
    # and make cycle counting for region, rather than component
    G = nx.from_pandas_edgelist(nlines, "start_node", "end_node", create_using=nx.MultiDiGraph)
    # cycle clusters
    cycles_simple = nx.simple_cycles(G)
    Gcy = G.subgraph(set([n for c in cycles_simple for n in c]))
    cycles = list(nx.connected_components(Gcy.to_undirected()))
    lines["cycle"] = 0
    nodes["cycle"] = 0
    lines["cycle_io"] = ""
    nodes["cycle_io"] = ""
    for i, c in enumerate(cycles):
        # check if cycles are sources or sinks
        Gedg = G.edge_subgraph([e for n in c for e in
                                list(G.in_edges(n, keys=True)) +
                                list(G.out_edges(n, keys=True))])
        outflows = [n for n in Gedg if len(Gedg.out_edges(n)) == 0 and len(Gedg.in_edges(n)) == 1]
        inflows = [n for n in Gedg if len(Gedg.out_edges(n)) == 1 and len(Gedg.in_edges(n)) == 0]
        # cycle counters and source/sink info
        nodes.loc[c, "cycle"] = i + 1
        lnids = nlines.start_node.isin(c) & nlines.end_node.isin(c)
        lines.loc[lnids, "cycle"] = i + 1
        if outflows and not inflows:
            io = "source"
            if len(outflows):
                io = "single-"+io
        elif inflows and not outflows:
            io = "sink"
            if len(inflows):
                io = "single-"+io
        else:
            io = "in-out"
        nodes.loc[c, "cycle_io"] = io
        lines.loc[lnids, "cycle_io"] = io
    # check where continuity violated
    nodes['continuity_violated_rivgraph'] = (
        nodes['continuity_violated'] if "continuity_violated" in nodes else 0
    )
    nodes['continuity_violated'] = 0
    nodes['discontinuity_io'] = ''
    single = nodes[nodes["cycle_io"].str.startswith("single")].index
    for n in G:
        inedg, outedg = G.in_edges(n), G.out_edges(n)
        nin, nout = len(inedg), len(outedg)
        io = 'single-' if n in single else ''
        if (nin == 0 and nout > 1 and nodes.loc[n, 'type'] != "outlet"):
            nodes.loc[n, 'continuity_violated'] = 1
            nodes.loc[n, 'discontinuity_io'] = io+'source'
        elif (nin > 1 and nout == 0 and nodes.loc[n, 'type'] != "outlet"):
            nodes.loc[n, 'continuity_violated'] = 1
            nodes.loc[n, 'discontinuity_io'] = io+'sink'
    # report
    nsingle = len([list(cy)[0] for cy in cycles if len(cy) == 1])
    cyio = nodes[nodes.cycle > 0].groupby("cycle").cycle_io.first()
    nsource, nsink = (cyio == "source").sum(), (cyio == "sink").sum()
    print(f"Found {len(cyio)} cycles ({nsingle} single-node, {nsource} sources, {nsink} sinks)")
    disnd = nodes.groupby('discontinuity_io').continuity_violated.sum().drop('')
    print(f"Found discontinuity in node types: {disnd.to_dict()}")

    lines.write()
    nodes.write()

    # temporarily remove cycles
    grass.run_command("g.copy", vector=[segments, output])
    cynodes, cylines = nodes[nodes.cycle > 0].index, lines[lines.cycle != 0].index
    if len(cylines) > 0:
        store_cycnod = {}
        nodes = GrassAttributeTable(output, layer=2)
        nodesnet = network_nodes(output)
        cyclset = set(cylines)
        for n, lns in nodesnet.loc[cynodes].line_cats.items():
            if len(set(lns) - cyclset) > 0:
                store_cycnod.setdefault(n, nodes.loc[n, "cycle"])
                nodes.loc[n, "cycle"] = 0
        nodes.write()
        prune_network_components(output, "cycle", ids=list(range(1, int(lines.cycle.max())+1)))
        # reset cycle ids
        nodes = GrassAttributeTable(output, layer=2)
        store_cycnod = pd.Series(store_cycnod)
        nodes.loc[store_cycnod.index, "cycle"] = store_cycnod
        nodes.loc[store_cycnod.index, "jtype"] = "cycle"
        nodes.write()

    return


def compute_strahler_order(start_nodes, end_nodes, loop_raise=False):
    """Compute the Strahler order of lines segments defined by fromto nodes starting at inlet nodes.
    """
    import networkx as nx

    enuniq = np.unique(end_nodes)
    orderX = start_nodes[~start_nodes.isin(enuniq)]
    # registry of links upstream of each node
    endix = pd.Series({e: set(v) for e, v in end_nodes.groupby(end_nodes).groups.items()})
    # loops need a separate registry
    loops = pd.Series(0, index=start_nodes.index)
    df = start_nodes.rename_axis('id').reset_index()
    df['end_node'] = end_nodes.values
    G = nx.from_pandas_edgelist(df, 'start_node', 'end_node', edge_attr=['id'],
                                create_using=nx.DiGraph)
    udg = G.to_undirected()
    for i, c in enumerate([zip(n,(n[1:]+n[:1])) for n in nx.simple_cycles(G)]):
        loops[[udg[u][v]['id'] for u, v in c]] = i + 1
    loopends = pd.Series(loops[loops > 0].values, index=end_nodes[loops[loops > 0].index])
    loopix = pd.Series({i: set(end_nodes.index[end_nodes.isin(nds)])
                        for i, nds in loopends.groupby(loopends).groups.items()},
                       dtype=object)
    # take loop nodes out of general registry and loop lines out of the loop registry
    for l, lid in loops[loops > 0].items():
        endix[end_nodes[l]] -= {l}
        loopix[lid] -= {l}

    i = 1
    strahler_order = pd.Series(-1, index=start_nodes.index)
    traversed = []
    # as long as there is one routed
    while len(orderX) >= 1: # and any([n not in self.nodes["outlets"] for n in orderX[:, 1]]):
        # loop over each node and correct order in order dict
        sys.stdout.write('\rCalculating stream order %s' % i)
        for sb in orderX.index:
            strahler_order[sb] = max([strahler_order[sb], i])
        # get downstream node ids
        ds = end_nodes[orderX.index]
        # remove line ids from node inventory
        for lix, nix in ds.items():
            if lix in endix[nix]:
                endix[nix] -= {lix}
        # get next lines with all upstream lines traversed
        nnodes = endix[ds.unique()].apply(len)
        nxtnodes = set(nnodes[nnodes == 0].index) - set(traversed)
        # any loops?
        nxtloopnds = loopends.index.intersection(nxtnodes)
        # only go into loop nodes if all upstream have been done
        if len(nxtloopnds):
            olends = loopends[nxtloopnds]
            uplines = pd.Series(ds.index, index=ds)
            # remove lines going into the loop nodes from the registry
            for lid, enids in olends.groupby(olends).groups.items():
                loopix[lid] -= set(uplines[enids])
            # dont include nodes in next unless all upstream lines been done
            nloopl = loopix[olends.unique()].apply(len)
            loopcont = olends.index[olends.isin(nloopl.index[nloopl > 0])]
            nxtnodes -= set(loopcont)
        # store nodes we have already gone through
        traversed.extend(nxtnodes)
        orderX = start_nodes[start_nodes.isin(nxtnodes)]
        # increase order
        i += 1
        # raise or warn
        if i > len(start_nodes) + 1:
            msg = "Order is exceeding the number of links, which indicates loops. Remaining nodes: %s" % orderX
            if loop_raise:
                raise RuntimeError(msg)
            else:
                warnings.warn(msg)
                break

    sys.stdout.write('\n')
    return strahler_order


def tidy_segments_columns(segments_vect, global_id_offset=0, write=True):
    """Remove all columns and add them back selectively and tidy.
    
    - adds topological info to table
    - makes ids global and int
    """
    links = GrassAttributeTable(segments_vect, layer=1)
    nodes = GrassAttributeTable(segments_vect, layer=2)
    lnks, nds = links.copy(), nodes.copy()
    links.drop(links.columns, axis=1, inplace=True)
    nodes.drop(nodes.columns, axis=1, inplace=True)
    linkstop, nodestop = add_global_topology_ids(segments_vect, global_id_offset=global_id_offset, write=False)

    # node type in single column
    ntype = nds.type
    ntype.loc[nds.jtype == 'b'] = 'bifurcation'
    ntype.loc[nds.jtype == 'c'] = 'confluence'
    ntype.loc[~nds.coast.isna() & (nds.n_segments == 1)] = 'coastal_outlet'
    ntype.loc[~nds.sink_id.isna() & (nds.n_segments == 1)] = 'sink_outlet'
    ntype.loc[ntype.isna() & (nds.n_segments == 1)] = 'inlet'
    ntype.loc[ntype.isna()] = 'confluence'

    # tidy nodes
    nodes['global_id'] = nodes.index + global_id_offset
    nodes['catchment_id'] = nds.comp + global_id_offset
    nodes[nodestop.columns] = nodestop
    nodes['node_type'] = ntype
    nodes['grwl_value'] = nds.grwl_value.fillna(0).astype(int)
    nodes['width'] = nds.river_width
    nodes['grwl_transition'] = nds.grwl_change 
    nodes['cycle'] = nds.cycle.fillna(0).astype(int)
    nodes['continuity_violated'] = nds.continuity_violated.fillna(0).astype(int)

    # tidy lines
    links['global_id'] = links.index + global_id_offset
    links['catchment_id'] = lnks.comp + global_id_offset
    links[linkstop.columns] = linkstop
    links['direction_algorithm'] = lnks.certain_alg.fillna(1)
    links['width_adjusted'] = lnks.wid_adj
    links['length_adjusted'] = lnks.len_adj.fillna(lnks.length)
    links['is_mainstem'] = lnks.is_mainstem.fillna(1).astype(int)
    links['strahler_order'] = lnks.strahler_order.fillna(-1).astype(int)
    links['cycle'] = lnks.cycle.fillna(0).astype(int)
    links[['length', 'azimuth']] = lnks[['length', 'azimuth']]
    links['sinuousity'] = lnks.sinuous

    if write:
        nodes.write()
        links.write()
    else:
        return links, nodes


def add_global_topology_ids(network_vect, global_id_offset=0, write=True):
    """Create upstream/downstream topology columns with global ids.
    Optionally, write to table.
    """
    linksnet = network_lines(network_vect)
    nodesnet = network_nodes(network_vect)

    # node up/downstream line ids
    nusids = pd.Series(index=nodesnet.index, dtype=str)
    ndsids = pd.Series(index=nodesnet.index, dtype=str)
    for n, lids in nodesnet.line_cats.items():
        lns = linksnet.loc[lids]
        usi = lns.index[lns.end_node == n] + global_id_offset
        dsi = lns.index[lns.start_node == n] + global_id_offset
        nusids[n] = ','.join(usi.astype(str))
        ndsids[n] = ','.join(dsi.astype(str))

    if write:
        links = GrassAttributeTable(network_vect, layer=1)
        nodes = GrassAttributeTable(network_vect, layer=2)
    else:
        links, nodes = pd.DataFrame(index=linksnet.index), pd.DataFrame(index=nodesnet.index)

    nodes['upstream_line_ids'] = nusids
    nodes['downstream_line_ids'] = ndsids
    links[['upstream_node_id', 'downstream_node_id']] = linksnet[['start_node', 'end_node']] + global_id_offset
    links['upstream_line_ids'] = nusids[linksnet.start_node].values
    links['downstream_line_ids'] = ndsids[linksnet.end_node].values

    if write:
        nodes.write()
        links.write()
    else:
        return links, nodes


def widest_main_channel(self):
    lix = self.links['id']
    self.links["main_channel_order"] = np.zeros_like(lix)
    for inl in self.nodes["inlets"]:
        cn = inl
        ord = 1
        # travers away from inlet node via widest 
        while cn not in self.nodes["outlets"]:
            ncon = self.nodes["conn"][self.nodes["id"].index(cn)]
            width = np.array([self.links["wid_adj"][lix.index(i)] for i in ncon])
            widest = ncon[np.argmax(width)]
            wid = lix.index(widest)
            # dont continue if order already assigned
            # inefficiency here because channels will be traversed multiple
            # times until highest order
            if self.links["main_channel_order"][wid] >= ord:
                break
            self.links["main_channel_order"][wid] = ord
            lncon = self.links["conn"][wid]
            lncon.remove(cn)
            cn = lncon[0]
    return

def upload_strahler_order(directed_network, column="strahler_order"):
    nlines = network_lines(directed_network)
    linestbl = GrassAttributeTable(directed_network)
    order = compute_strahler_order(nlines["start_node"], nlines["end_node"])
    linestbl.loc[nlines.index, column] = order
    linestbl.write()
    return


def split_interbasin_network(network_vector, node_type_column="type"):
    """Remove from interbasin nodes outward facing lines to split the network.

    Outward facing is the part the part of the network with fewer nodes.
    """
    import networkx as nx
    import grass.script as gs

    lines = network_lines(network_vector, attr=True)
    lines['id'] = lines.index
    nodes = network_nodes(network_vector, attr=True)

    # ignore interbasin nodes connected to dangles as they will leave single node components behind
    dngls = lines.loc[nodes[nodes.n_lines == 1].line_cats.apply(lambda l: l[0]).unique()]
    nodes.drop(dngls[["start_node", "end_node"]].values.flatten(), inplace=True)

    lns_drop = []
    for id, ibn in nodes[(nodes.type == "interbasin") & (nodes.n_lines > 1)].iterrows():
        nxtlns = lines.loc[ibn.line_cats, ["start_node", "end_node"]]
        nxtnds = pd.Series(-1, index=nxtlns.index, dtype=int)
        nxtnds[nxtlns.start_node != id] = nxtlns[nxtlns.start_node != id].start_node
        nxtnds[nxtlns.end_node != id] = nxtlns[nxtlns.end_node != id].end_node
        G = nx.from_pandas_edgelist(lines.drop(ibn.line_cats), 'start_node', 'end_node', edge_attr=['id'])
        # find components that include nxtnds and count their nodes
        nnds = pd.Series(-1, index=nxtlns.index)
        for l, n in nxtnds.iteritems():
            for c in nx.connected_components(G):
                if n in c:
                    nnds[l] = len(c)
                    break
        # drop line with fewest nodes in component
        lns_drop.append(nnds.idxmin())
    #
    print(f"Dropping {len(lns_drop)} lines...")
    lnsdrop = ','.join(map(str, lns_drop))
    grass.run_command('v.edit', map=network_vector, tool='delete', layer=1, where="cat in (%s)" % lnsdrop)
    tbl = gs.vector_db(network_vector)[1]["table"]
    grass.run_command('db.execute', sql=f"DELETE FROM {tbl} WHERE cat in (%s);" % lnsdrop)
    return


def zonal_stats(raster, zones, method="mean", upload_vect=None, upload_column=None):
    """Return series of statistics over base raster categories.

    method may be any combination (comma separated) of r.univar -e output.
    upload_vect/upload_column allows upload to a vector table with the same categories as base.
    """
    method = [method] if type(method) == str else list(method)
    if upload_column is not None:
        upload_column = [upload_column] if type(upload_column) == str else list(upload_column)
        assert len(method) == len(upload_column)

    tmpfile = grass.tempfile(False)
    grass.run_command("r.univar", map=raster, zones=zones, output=tmpfile, separator="comma", flags="et")
    stats = pd.read_csv(tmpfile, index_col=0)[method]
    if upload_vect and upload_column:
        print(f"Uploading {upload_column} to table of {upload_vect}...")
        tbl = GrassAttributeTable(vector=upload_vect)
        tbl[upload_column] = stats
        tbl.write()
    return stats



def network_raster_segments(
    network_vector, check_missing=None, cache_segments=None, trim_nodes=True,
):
    # import ipdb;ipdb.set_trace()
    links = network_lines(network_vector, attr=True)

    nodes = network_nodes(network_vector)
    nodes["upstream_segments"] = [
        links.index[links.end_node == i].tolist() for i in nodes.index
    ]
    nodes["downstream_segments"] = [
        links.index[links.start_node == i].tolist() for i in nodes.index
    ]
    # ignore nodes that have neither up/down segments
    noseg = nodes.index[(nodes.upstream_segments.apply(len) == 0) &
                        (nodes.downstream_segments.apply(len) == 0)]
    nodes.drop(noseg, inplace=True)
    nodes["bifurcation"] = nodes["downstream_segments"].apply(len) > 1

    # segments with each raster cell either produce new or read from cached file
    if cache_segments is None or not osp.exists(cache_segments):
        segments_with_nodes = dense_lines(
            network_vector, res=grass.region()["nsres"], check_missing=check_missing
        )
        segments_with_nodes = segments_with_nodes if check_missing is None else segments_with_nodes[0]
        if cache_segments:
            segments_with_nodes.to_pickle(cache_segments)
    else:
        warnings.warn(f"Using cached dense lines from {cache_segments}!")
        segments_with_nodes = pd.read_pickle(cache_segments)

    nodes["coor"] = [
        tuple(
            segments_with_nodes[us[0]][-1] if len(us) else segments_with_nodes[ds[0]][0]
        )
        for _, (us, ds) in nodes[
            ["upstream_segments", "downstream_segments"]
        ].iterrows()
    ]
    if trim_nodes:
        # first and last cell removed
        segments = segments_with_nodes.progress_apply(lambda l: list(map(tuple, l[1:-1])))
        # remove any segments that only consist of nodes
        segments.drop([i for i, l in segments.items() if len(l) == 0], inplace=True)
    else:
        segments = segments_with_nodes
    # drop any line that doesnt connect to nodes on both ends
    discon = links[(links.start_node <= 0) | (links.end_node <= 0)]
    if len(discon):
        warnings.warn(f"Ignoring {len(discon)} lines that are not connected to nodes.")
        links.drop(discon.index, inplace=True)

    return segments, nodes, links


def drainage_from_network(network_vector, output, network_cache=None, resolution=30):
    """Create a D8 drainage raster (grass encoded) from a directed river network.
    """
    segments, nodes, links = network_raster_segments(
            network_vector, cache_segments=network_cache, trim_nodes=False,
        )
    re = resolution
    grass_drainage_map = {(-re,  re): 3, (  0,  re): 2, ( re,  re): 1,
                          (-re,   0): 4, (  0,   0): 0, ( re,   0): 8,
                          (-re, -re): 5, (  0, -re): 6, ( re, -re): 7}

    for ds in nodes.downstream_segments:
        if len(ds) > 1:
            ds.remove(links.loc[ds, "wid_adj"].idxmax())
            segments[ds[0]] = segments[ds[0]][1:]
    # segments with less than 2 cells will not work
    segments = segments[segments.apply(len) >= 2]

    def drain(segment):
        seg = np.array(segment)
        direction = seg[1:, :] - seg[:-1, :]
        direction[direction<-resolution] = -resolution
        direction[direction>resolution] = resolution
        return np.apply_along_axis(lambda x: grass_drainage_map[tuple(x)], 1, direction)

    # the actual translation
    drainage = segments.progress_apply(drain)

    print("Writing values to raster...")
    tmpf = grass.tempfile()
    outlets = nodes.loc[nodes.downstream_segments.apply(len) == 0, "coor"]
    with open(tmpf, "w") as f:
        for li, l in drainage.iteritems():
            # ignore end node coordinate
            xyz = pd.DataFrame(segments[li][: -1])
            xyz.loc[:, 2] = l
            xyz.to_csv(f, header=False, index=False)
        # set outlets to 0 as sinks to avoid loops
        pd.DataFrame([l + (0,) for l in outlets]).to_csv(f, header=False, index=False)
    grass.run_command(
        "r.in.xyz",
        input=tmpf,
        output=output,
        separator="comma",
        percent=5,
        type="CELL",
        method="min",  # ensures that 0 for outlets are always preserved
    )
    return


def upstream_bifurcations(segments, output_column=None):
    """For each line and node, report the number of unique bifurcations upstream."""
    G = river_networkx(segments)
    nodes, links = GrassAttributeTable(segments, layer=2), GrassAttributeTable(segments, layer=1)
    # to look upstream, we need to reverse
    Grev = G.reverse()
    bifurcations = {n for n in G if len(G.out_edges(n)) > 1}
    node_count = pd.Series(0, index=nodes.index)
    for n in tqdm(G):
        node_count[n] = len(nx.descendants(Grev, n) & bifurcations)
    # add the count of the upstream node to links plus 1 if the node is a bf
    link_count = pd.Series({
        l: node_count[s] + (1 if s in bifurcations else 0)
        for s, _, l in G.edges(data="id")
    })
    if output_column:
        nodes[output_column] = node_count
        links[output_column] = link_count
        nodes.write()
        links.write()
    return node_count, link_count


def snap_osm_river_names(segments, osm_lines, split_length=200, max_dist=1000, output=None):
    from scipy.stats import pearsonr

    # split into shorter segments and upload segment id and length
    grass.run_command("v.split", input=segments, output="segments__split", length=split_length)
    grass.run_command("v.category", input="segments__split", option="add", layer=3, output="segments__split__reach_ids", type="line")
    grass.run_command("v.db.addtable", map="segments__split__reach_ids", layer=3,
                      columns="segment_id int,dist double,osm_cat int")
    grass.run_command("v.to.db", map="segments__split__reach_ids", layer=3, option="query", query_layer=1,
                      query_column="cat", columns="segment_id")
    grass.run_command("v.to.db", map="segments__split__reach_ids", layer=3, option="length", units="meters", column="length")
    # get nearest osm attributes to split segments
    grass.run_command("v.distance", from_="segments__split__reach_ids", from_layer=3, to=osm_lines,
                      upload="dist,cat", column="dist,osm_cat")

    seg = GrassAttributeTable('segments__split__reach_ids', layer=3)
    osm = GrassAttributeTable(osm_lines)
    names = pd.DataFrame(columns=["local", "en"], dtype=str)
    names.index.name = "cat"
    for s in seg.segment_id.unique():
        se=seg[seg.segment_id == s]
        # linear correlation?
        if len(se) > 3:
            r, p = pearsonr(se.length.cumsum(),se.dist)
            if abs(r) > 0.5 and p < 0.005:
                continue
        # too far away?
        if ((se.dist > max_dist).sum() / len(se)) > 0.5:
            continue
        # if we have made it this far, take the most common id
        oml = osm.loc[se.osm_cat.mode()[0]]
        names.loc[s] = None
        names.loc[s, 'local'] = oml['name']
        if oml.other_tags:
            try:
                # [1:-1] removes first and last ""
                ot = dict(map(lambda s: s.split('"=>"'), oml.other_tags[1:-1].split('","')))
            except Exception:
                warnings.warn(f"Cant parse other tags: {oml.other_tags}")
            # checks preferred tags in order
            for k in ["name:en", "int_name", "name:fr"]:
                if k in ot:
                    names.loc[s, 'en'] = ot[k]
                    break
    # write out
    if output:
        names.to_csv(output)
        return
    return names


def fillna_shortest_path(network_vect, column, comp_column="catchment_id", data=None, priority=None,
                         weight="not implemented", max_length=20, write=False):
    from itertools import combinations, chain
    seg = GrassAttributeTable(network_vect) if data is None else data
    if False:
        return seg[column]
    lines = network_lines(network_vect)
    # we dont care if they are connected up/downstream
    gr = river_networkx(network_vect)
    sortnds = list(nx.topological_sort(gr))
    gr = gr.to_undirected()
    # set of line ids that can be filled
    nanlines = set(seg.index[seg[column].isna()])
    nnans = len(nanlines)

    # decide order
    gb = seg.groupby([comp_column, column])
    prior = gb[priority].sum() if priority else gb.count().iloc[:, 0]
    order = prior.sort_values(ascending=False).index
    # updating progress bar
    pbar = tqdm(order)
    for comp, name in pbar:
        idxs = gb.groups[(comp, name)]
        pbar.set_description(f"{name:20s}")
        # all up/down nodes of group
        named_nodes = set(lines.loc[idxs].values.flatten())
        # make subgraph to avoid calculating shortest paths that are already connected
        named_sg = gr.subgraph(named_nodes)
        named_sparse = set([n for n in named_sg if len(named_sg.edges(n)) != len(gr.edges(n))])
        # split into nan components
        nannodsincomp = lines.loc[nanlines & set(seg[seg[comp_column]==comp].index)].values.flatten()
        nangr = gr.subgraph(nannodsincomp).to_undirected().copy()
        # use order to approximate distance to limit number of path combinations
        order = {n: sortnds.index(n) for n in named_sparse if n in sortnds}
        # get all shortest paths
        all_paths_nodes = set()
        comps = list(nx.connected_components(nangr))
        for comp in tqdm(comps):
            nmdnds = comp & named_sparse
            if len(nmdnds) < 2:
                continue
            # combinations will only fire if more than 1 node in comp
            combos = set(combinations(nmdnds, 2))
            for s, e in tqdm(combos):
                if abs(order[s] - order[e]) <= max_length and nx.has_path(nangr, s, e):
                    pths = nx.all_simple_paths(nangr, s, e, cutoff=max_length)
                    all_paths_nodes |= set(chain.from_iterable(pths))
                    nangr.remove_edges_from([(p[i], p[i+1]) for p in pths for i in range(0, len(p)-1)])
        # unpack lines
        minimum_sg = gr.subgraph(set(all_paths_nodes))
        line_ids = set(i for _, _, i in minimum_sg.edges(data="id"))
        fillids = line_ids & nanlines
        seg.loc[fillids, column] = name
        # remove from nans
        nanlines -= fillids

    print(f"Filled {(nnans-len(nanlines))} of NaNs in {column} column.")
    if write:
        seg.write()
    else:
        return seg[column]


def rivgraph_fixlinks_from_lines(fixlinks_file, segments, outputdir,
                                 buffer=1, min_length_share=0.9, azimuth_tolarance=10):
    """Takes a selection of correctly directed segment lines in gpkg, matches them to undirected
    segments (without using ids, making it transferable across versions) and writes out Rivgraph
    fixlinks files for each component.
    """
    grass.run_command("v.import", input=fixlinks_file, output="fixlinks")
    # make sure we have a clean table
    grass.run_command("v.db.droptable", map="fixlinks", flags="f")
    grass.run_command("v.db.addtable", map="fixlinks")
    # length and azimuth of fixlinks
    grass.run_command("v.to.db", map="fixlinks", option="azimuth", columns="azimuth")
    grass.run_command("v.to.db", map="fixlinks", option="length", columns="length")
    # match via overlapping length with buffered fixlinks
    grass.run_command("v.buffer", flags="t", input="fixlinks", output="fixlinks__polygons", distance=buffer)
    grass.run_command("v.overlay", ainput=segments, binput="fixlinks__polygons", output="segments__overlay", operator="and")
    # length and azimuth of filtered segments
    grass.run_command("v.to.db", map="segments__overlay", option="azimuth", column="azimuth")
    grass.run_command("v.to.db", map="segments__overlay", option="length", column="length")

    seg = GrassAttributeTable("segments__overlay")
    # filter by min overlap
    seg = seg.loc[seg["length"]/seg["b_length"] > min_length_share]
    print(f"Found {len(seg)} matching fix lines.")

    seg.loc[:, "azimuth_diff"] = ((seg.b_azimuth - seg.azimuth + 180) % 360 - 180).abs()

    # assign upstream node based on direction
    seg["upstream_node"] = np.nan
    links = network_lines(segments)
    # same direction, take start node as upstream node
    samedir = seg.azimuth_diff < azimuth_tolarance
    seg.loc[samedir, "upstream_node"] = links.loc[seg.loc[samedir, "a_cat"], "start_node"].values
    # wrong direction, take end node as upstream node
    wrongdir = seg.azimuth_diff > (180 - azimuth_tolarance)
    seg.loc[wrongdir, "upstream_node"] = links.loc[seg.loc[wrongdir, "a_cat"], "end_node"].values

    nomatching = seg.loc[seg.upstream_node.isna()]
    if len(nomatching):
        nom = nomatching[["a_cat", "b_azimuth", "azimuth", "azimuth_diff"]].to_string()
        raise RuntimeError(f"The following fix links dont match by angle:\n{nom}")

    # write out rivgraph fixlink files by component
    for comp, ix in seg.groupby("a_comp").groups.items():
        pth = osp.join(outputdir, "component_%04i_fixlinks.csv" % comp)
        seg.loc[ix, ["a_cat", "upstream_node"]].to_csv(pth, index=False, header=["link_id", "usnode"])
        print(f"Wrote {len(ix)} fix links to {pth}")


def check_domain_drainage_share(segments, component_catchments, domain):
    """Return the proportion of the domain that drains into outlets of the domain including
    coastal catchments drainage area."""
    # drainage area of nodes needs to be attached from darea_node_info table
    nodes=GrassAttributeTable(segments, layer=2)
    catchments=GrassAttributeTable(component_catchments)
    domain=GrassAttributeTable(domain)

    outlet_da = nodes.loc[nodes["node_type"].isin(["coastal_outlet", "sink_outlet"]), "drainage_area"].sum()
    coastal_da = catchments.loc[catchments["is_coastal"] == 1, "area"].sum()
    da_share = (outlet_da + coastal_da)*100/domain.area.sum()

    return da_share


if __name__ == "__main__":
    from commandline import interface as _cli
    _cli()
