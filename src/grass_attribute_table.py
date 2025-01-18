#!/usr/bin/env python

import os.path as osp
import os
import sqlite3

import pandas as pd
from pandas import DataFrame


class GrassAttributeTable(DataFrame):
    """A plugin to represent a grass vector attribute table.
    This plugin can either be used with a project instance as the first
    argument to get the GRASS connection parameters or by specifically stating
    the database path. E.g.::
        dbpath = 'grassloc/PERMANENT/sqlite/sqlite.db'
        table = GrassAttributeTable(vector='test', database=dbpath)
    Specify `database` and `table`, if you dont want to rely on grass to get
    the table connection parameters.
    """
    vector = None
    #: Optional
    layer = 1
    #: Specify a different key index, default is the first column
    key = None
    #: list, optional subset of columns to read, writing reads+writes full tbl
    subset_columns = None
    #: optional if it shouldnt call grass to find out
    database = None
    #: optional if it shouldnt call grass to find out
    table = None
    #: needed to stop exposing all pd.DataFrame methods
    plugin = []

    def __init__(self, project=None, **override):
        super(GrassAttributeTable, self).__init__()
        if type(project) == str:
            project, override['vector'] = None, project
        self.__dict__.update(override)
        self.project = project
        em = 'vector (str) needed (others are optional).'
        assert type(self.vector) == str, em
        nms = self.vector.split('@')
        if project and not self.database:
            self.mapset = nms[1] if len(nms) > 1 else project.grass_mapset
            self.database = osp.join(project.grass_db, project.grass_location,
                                     self.mapset, 'sqlite', 'sqlite.db')
        # retrive info in a grass session active
        elif os.environ.get('GISRC', None) and not self.database:
            import grass.script as grass
            con = grass.vector_db(self.vector).get(self.layer, None)
            msg = '%s has no sqlite database connection.' % self.vector
            assert con and con['driver'] == 'sqlite', msg
            self.__dict__.update(
                {k: con[k] for k in ['table', 'database', 'key']})
        else:
            assert self.database, 'Without project, database must be given.'
        self.table = self.table or nms[0]
        self.key = self.key or 0
        # fill dataframe
        self.read()
        return

    def read(self):
        """Read table from db."""
        sc = self.subset_columns
        cols = ','.join(sc) if sc else '*'
        with self.dbconnection as con:
            tbl = pd.read_sql('select %s from %s;' % (cols, self.table), con)
        self.key = tbl.columns[self.key] if type(self.key) == int else self.key
        tbl.set_index(self.key, inplace=True, verify_integrity=True)
        # fill DataFrame
        super(GrassAttributeTable, self).__init__(tbl)
        return

    @property
    def dbconnection(self):
        return sqlite3.connect(self.database, timeout=20)

    def write(self, table_name=None):
        """Save table back to GRASS sqlite3 database.
        """
        cleantbl = self
        with self.dbconnection as con:
            if self.subset_columns:  # read other columns
                tbl = pd.read_sql('select * from %s;' % self.table, con)
                tbl.set_index(self.key, inplace=True, verify_integrity=True)
                tbl[cleantbl.columns] = cleantbl
                cleantbl = tbl
            cleantbl.to_sql(table_name or self.table, con, if_exists='replace')
        return


def join(from_vect, to_vect, from_ix="cat", to_ix="cat", columns=None, **kw):
    """Copy columns between arbitrary vectors (assuming cat match).

    Arguments
    ---------
    from_vect, to_vect : str
        Any vector string spec accepted by GrassAttributeTable.
    from_ix, to_ix : str
        Index columns to use for join.
    columns : list | None
        Subset of columns to copy or None for all (default).
    **kw :
        Additional keywords to parse to both from/to tables.
    """
    fromtbl = GrassAttributeTable(
        vector=from_vect, subset_columns=columns, key=from_ix, **kw
    )
    totbl = GrassAttributeTable(vector=to_vect, key=to_ix, **kw)

    for n, c in fromtbl.items():
        totbl.loc[c.index, n] = c

    totbl.write()
    return


if __name__ == "__main__":
    from commandline import interface
    interface()
