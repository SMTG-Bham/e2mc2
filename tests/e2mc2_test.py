import os
from os.path import join
import unittest
import tempfile
import filecmp
from contextlib import contextmanager

import pkg_resources

json_data = pkg_resources.resource_filename(__name__, 'MgAl2O4.json')

eci_str = (u"-1.623710\n0.243466\n0.000000\n"
           u"0.164667\n0.110744\n-0.030886\n"
           u"-0.063562\n-0.030195\n-0.060852\n")

import e2mc2

class TestClusterExpansion(unittest.TestCase):
    def test_json_read(self):
        ce = e2mc2.ClusterExpansion(json_data)
        self.assertEqual(ce.eci, (eci_str))

    def test_json_write(self):
        ce = e2mc2.ClusterExpansion(json_data)
        with named_tempfile() as f:
            ce.write_json(f)
            self.assertTrue(filecmp.cmp(json_data, f),
                            "Output JSON does not match input")

            
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
