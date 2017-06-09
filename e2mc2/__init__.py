###############################################################################
# Copyright 2017 Adam Jackson
###############################################################################
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
###############################################################################

# Input files that need to be managed
# - lat.in
# - clusters.out
# - eci.out
# - gs_str.out

import os
from os.path import isdir, join, curdir
import io
import tarfile
import subprocess
import json
from contextlib import contextmanager
import glob
import re

import pandas as pd

@contextmanager
def run_directory(run_path):
    """
    Temporarily work in another directory, creating it if necessary.

    Inspired by https://pythonadventures.wordpress.com/2013/12/15/chdir-a-context-manager-for-switching-working-directories

    """
    home = os.getcwd()
    if not os.path.isdir(run_path):
        os.makedirs(run_path)
    os.chdir(run_path)
    yield
    os.chdir(home)


class MonteCarloCalc:
    default_params = {
        "T0": 100,
        "T1": 1500,
        "dT": 100,
        "db": None,
        "k": 8.617e-5,
        "abs": False,
        "mu0": None,
        "mu1": None,
        "dmu": None,
        "gs": 0,  # Initial state (-1 for disordered)
        "phi0": None,
        "er": 12,
        "innerT": False,
        "eq": 1e4,
        "n": 1e4,
        "x": None,  # Concentration
        "dx": None,  # Convergence threshold
        "aq": None,  # Convergence param: 0:energy 1:conc
        "tstat": 0,  # Discontinuity test (recommended value 3)
        "sigdig": 6,
        "o": "mc.out",  # Output filename
        "oss": "mcsnapshot.out",
        "opss": "psnapshot.out",  # periodic snapshots
        "is": None,  # Initial state file
        "ks": None,
        "cm": True,  # Canonical (fixed composition)
        "g2c": True,  # Report canonical values instead of GC
        "q": False,  # Quiet
        "sd": None,  # Random number seed
        "dl": False,  # Drop last point after phase transition
    }

    _bool_params = {"innerT", "cm", "q", "abs", "dl", "g2c"}
    _int_params = {"eq", "n", "gs"}

    def __init__(self, cluster_expansion, **kwargs):
        """
        Monte Carlo from cluster expansion using ATAT/EMC2

        Args:
            cluster_expansion (e2mc2.ClusterExpansion or str): 
                Object containing CE parameters or path to input files or path
                to serialised data (.tar or .json file)

            **kwargs: Optional arguments corresponding to EMC2 command line 
                parameters. To change defaults, you can modify 
                MonteCarloCalc.default_params before instantiating object.

        Properties:
            params (dict): Dictionary of EMC2 command line parameters
            cluster_expansion (e2mc2.ClusterExpansion): CE data for calculation

        Methods:
            set: Set EMC2 parameters
            set_ce: Update cluster expansion data
        """

        if isinstance(cluster_expansion, ClusterExpansion):
            self.cluster_expansion = cluster_expansion
        elif type(cluster_expansion) == str:
            self.cluster_expansion = ClusterExpansion(cluster_expansion)
        else:
            raise ValueError("Not appropriate CE data")

        self.params = MonteCarloCalc.default_params

        self.set(**kwargs)

    def set(self, **kwargs):
        """Update EMC2 run parameters with keyword args"""
        self.params.update(kwargs)

    def calc_write_files(self, calc_directory="emc2_run"):
        """Prepare directory for EMC2 calculation, with log of params"""
        if isdir(calc_directory):
            raise IOError("Calculation directory {0} already "
                          "exists".format(calc_directory))

        self.cluster_expansion.write_dir(calc_directory)

        with open(join(calc_directory, "emc2_params.json"), 'w') as f:
                json.dump(self.params, f)

    def run(self, calc_directory="emc2_run"):
        """Call EMC2 in a new directory with current parameter settings"""
        self.calc_write_files(calc_directory)

        command = ['emc2']

        for k, v in self.params.items():
            if v is None or v is False:
                pass
            elif k in MonteCarloCalc._int_params:
                command.append('-' + k)
                command.append(str(int(v)))
            elif k in MonteCarloCalc._bool_params:
                command.append('-' + k)
            else:
                command.append('-' + k)
                command.append(str(v))

        with run_directory(calc_directory):
            subprocess.call(command)

    def read_output(self, output="emc2_run"):
        """Read output files from EMC2 run and attach to object"""
        if isdir(output):

            n_clusters = len(self.ce.eci.split())
            data = pd.read_table(join(output, 'mc.out'), skipinitialspace=True,
                                 usecols=range(17 + n_clusters), header=None)
            clusters_headers = ['C' + str(i+1) for i in range(n_clusters)]
            data.columns = ['T', 'mu', 
                            'E', 'x', 'phi', 'E2', 'x2', 
                            'E_lte', 'x_lte', 'phi_lte', 
                            'E_mf', 'x_mf', 'phi_mf', 
                            'E_hte', 'x_hte', 'phi_hte', 
                            'lro'] + clusters_headers

            self.mc_data = data

            # emc2 is a bit inconsiderate in how it handles snapshots, and the
            # counter eats into the end of the filename. We need to use globs
            # and regex to find and sort these. We assume that the first four
            # characters are distinctive and intact.
            if "opss" in self.ce.params and self.ce.params["opss"] is not None:
                basename = self.ce.params["opss"]
                snapshots = glob.glob(basename[:4] + '*.out')

                num_parts = re.compile(r'\d+').findall
                def sortkey(filename):
                    return int(num_parts(filename)[0])
                self.snapshot_files = sorted(snapshots, key=sortkey)

class ClusterExpansion:
    def __init__(self, mapsrun):
        """
        Container for significant MAPS output.

        Args:
            mapsrun (str): String defining one of the following:
                - Directory location containing lat.in, clusters.out, eci.out,
                  gs_str.out
                - Tar file containing same set of files
                - JSON serialised form of this object                
        """

        self._infiles = (("lat.in", "lat"), ("clusters.out", "clusters"),
                         ("eci.out", "eci"), ("gs_str.out", "gs"))

        if isdir(mapsrun):
            for infile, p in self._infiles:
                with open(join(mapsrun, infile), 'r', encoding="ascii") as f:
                    setattr(self, p, f.read().decode())
        elif tarfile.is_tarfile(mapsrun):
            with tarfile.open(mapsrun, 'r', encoding="ascii") as t:
                for infile, p in self._infiles:
                    with t.extractfile(infile) as f:
                        setattr(self, p, f.read().decode())

        else:
            try:
                with open(mapsrun, 'r') as f:
                    data = json.load(f)
            except ValueError:
                raise Exception("This doesn't seem to be a directory, "
                                "a tar file or a JSON file. I give up.")
            for _, p in self._infiles:
                setattr(self, p, data[p])

    def todict(self):
        """Get a dict of input data files for serialisation"""
        return {"lat": self.lat,
                "clusters": self.clusters,
                "eci": self.eci,
                "gs": self.gs}

    def write_json(self, filename):
        """Serialise cluster expansion parameters to JSON"""
        with open(filename, 'w') as f:
            json.dump(self.todict(), f, sort_keys=True)

    def write_dir(self, dirname):
        """Write cluster expansion files to directory for EMC2 run"""
        if not isdir(dirname):
            os.mkdir(dirname)

        for filename, p in self._infiles:
            with open(join(dirname, filename), 'w') as f:
                f.write(getattr(self, p))

    def write_tar(self, filename):
        """Write cluster expansion files to tar file for archiving"""
        with tarfile.open(filename, 'w') as t:
            for infile, p in self._infiles:
                tarinfo = tarfile.TarInfo(name=infile)
                data = getattr(self, p).encode("ascii")
                tarinfo.size = len(data)
                data_io = io.BytesIO(data)
                t.addfile(tarinfo=tarinfo, fileobj=data_io)
