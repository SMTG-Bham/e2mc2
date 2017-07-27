import os
from os.path import join, basename, isfile
import unittest
import tempfile
import tarfile
import filecmp
import shutil
from contextlib import contextmanager

import pkg_resources

json_data = pkg_resources.resource_filename(__name__, 'MgAl2O4.json')

eci_str = [-1.623710, 0.243466, 0.000000,
           0.164667, 0.110744, -0.030886,
           -0.063562, -0.030195, -0.060852]

import e2mc2

class TestClusterExpansion(unittest.TestCase):
    def test_json_read(self):
        ce = e2mc2.ClusterExpansion(json_data)
        self.assertListEqual(ce.eci, (eci_str))

    def test_json_write(self):
        ce = e2mc2.ClusterExpansion(json_data)
        with named_tempfile() as f:
            ce.write_json(f)
            self.assertTrue(filecmp.cmp(json_data, f),
                            "Output JSON does not match input")

    def test_output_cmp_tar_with_files(self):
        """Write cluster expansion data to files and tar, check files match"""
        ce = e2mc2.ClusterExpansion(json_data)
        with named_tempdir() as rundir:
            ce.write_dir(rundir)
            os.mkdir(join(rundir, "tar"))

            ce.write_tar(join(rundir, "tar", "tmp.tar"))
            with tarfile.open(join(rundir, "tar", "tmp.tar"), 'r') as t:
                t.extractall(join(rundir, "tar"))

            for entry in os.listdir(path=rundir):
                if isfile(entry):
                    self.assertTrue(filecmp(join(rundir, entry),
                                    join(rundir, tar, entry)),
                                    "Output file {0} not same in tar "
                                    "output".format(entry))

@contextmanager
def named_tempfile():
    tmpfile, filename = tempfile.mkstemp()
    os.close(tmpfile)
    try:
        yield filename
    finally:
        os.remove(filename)

if __name__ == '__main__':
    unittest.main()

@contextmanager
def named_tempdir():
    tmpdir = tempfile.mkdtemp()
    try:
        yield tmpdir
    finally:
        shutil.rmtree(tmpdir)
