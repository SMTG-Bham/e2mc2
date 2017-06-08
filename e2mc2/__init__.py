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

from os.path import isdir, join
import io
import tarfile
import json

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

        self._infiles = (("lat.in", "lat"),
                         ("clusters.out", "clusters"),
                         ("eci.out", "eci"),
                         ("gs_str.out", "gs"))


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
