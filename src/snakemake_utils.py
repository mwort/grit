import os
import os.path as osp
from pathlib import Path as P

import pandas as pd
try:
    import pysftp
    from snakemake.remote.SFTP import RemoteProvider, RemoteObject
except ImportError:
    RemoteProvider = object


PROJECT_ROOT = P(os.environ['PROJECT_ROOT'])
DATA_DIR = PROJECT_ROOT/'data'


class CondaEnvs:
    def __init__(self, path_pattern):
        import glob
        envs = glob.glob(path_pattern)
        envsdict = {osp.basename(osp.splitext(f)[0]).split('_')[2]: f for f in envs}
        self.__dict__.update(envsdict)


conda_envs = CondaEnvs(os.environ['SNAKEMAKE_CONDA_ENVFILES'])


def file_or_expand(file, expand_pattern, column):
    """Return file if it doesnt exist or if it does, expand pattern using the csv column in file."""
    if osp.exists(file):
        ids = pd.read_csv(file, dtype=str)
        file = [expand_pattern.format(**{column: i}) for i in ids[column]]
    return file


class ExploreMichelWortmannRemote(RemoteProvider):
    """This allows files to be uploaded to Michel's server via Snakemake.
    It requires password-less ssh keys to be setup for the remote host.
    server = ExploreMichelWortmannRemote()
    rule:
        output: server.remote('path')
    """
    #host = 'h2895971.stratoserver.net'
    host = "michelwortmann.com"
    prefix = '/var/www/vhosts/michelwortmann.com/explore.michelwortmann.com/'

    def __init__(self, username="michel", **kwargs):
        super().__init__(username=username, **kwargs)

    def remote(self, path, **kwargs):
        lpth = osp.join(self.prefix, path)
        if not lpth.startswith("/"):
            lpth = "/" + lpth
        server_path = super().remote(self.host + lpth, **kwargs)
        return server_path


def group_permissions_hook(logfile=None, output=[".snakemake"]):
    """Intended to be used in onsuccess/onerror snakemake hooks to ensure.
    ouput has the right group permissions.
    """
    import re
    from snakemake.shell import shell
    if logfile:
        parsed_out = re.findall(r"\s*output: (.*)\n", open(logfile).read())
        output += [o.replace(",", "") for o in parsed_out]
    if output:
        print(f"Postprocessing output: {' '.join(output)}")
        shell("chmod --quiet -R g+w {output} || true")
    return output
