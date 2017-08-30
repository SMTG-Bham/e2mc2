#! /usr/bin/env python3

import ase.io
import glob
import argparse
from os.path import basename

from e2mc2 import atoms_from_sqs

def make_traj(pattern='*snapsh*.out', filename='stitch.traj', files=None):
    """Combine a list of files or pattern match into ASE trajectory format"""
    traj = []
    if files is None or files == []:
        files = glob.glob(pattern)

    for outcar in files:
        if basename(outcar).split('.')[-1]:
            traj += [atoms_from_sqs(outcar)]
        else:
            traj += ase.io.read(outcar, index=':')
    ase.io.write(filename, traj)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('filename', type=str, nargs='?', default='stitch.traj',
                        help="Output file name, including .traj extension")
    parser.add_argument('files', type=str, nargs='*', default=None,
                        help="Ordered list of files to include in trajectory")
    parser.add_argument('--pattern', type=str, default='*/OUTCAR',
                        help="If no file list is provided, use pattern to find"
                             "files; * acts as a wildcard.")

    args = parser.parse_args()

    make_traj(pattern=args.pattern, filename=args.filename, files=args.files)

if __name__ == '__main__':
    main()
