#! /usr/bin/env python3
"""Add atoms to a structure or trajectory"""

import argparse
from os.path import basename

import ase.io

from e2mc2 import atoms_from_sqs

def _read_structure(filename):
    if basename(filename).split('.')[-1] == 'out':
        return atoms_from_sqs(filename)
    else:
        atoms = ase.io.read(filename, index=':')
        if type(atoms) == list and len(atoms) == 1:
            return atoms[0]
        else:
            return atoms

def add_atoms_files(structure, add_atoms, output=None, elements=None,
                    supercell=None):
    """
    Append atoms from one structure file to another, which may be a trajectory

    Args:
        structure (str): Path to main structure or trajectory
        add_atoms (str): Reference file containing atoms to be added
        output (str): Filename for composite output
        elements (set): Element labels to draw from add_atoms. If None, use
            all atoms.
        supercell (3-tuple): Dimensions of supercell expansion to be applied to
            add_atoms. If None, try to calculate this instead.
    """

    traj = _read_structure(structure)
    if type(traj) is not list:
        traj = [traj]
    ref_atoms = _read_structure(add_atoms)

    composite = add_atoms_objects(traj, ref_atoms, elements=elements,
                          supercell=supercell)

    if output is not None:
        if len(composite) == 1:
            ase.io.write(output, composite[0])
        else:
            ase.io.write(output, composite)            
    else:
        print(composite)

def add_atoms_objects(traj, ref_atoms, elements=None, supercell=None):
    """
    Append atoms from one ASE atoms object to a trajectory
    Args:
        traj ([ase.Atoms, ...]): List of ASE Atoms objects.
        ref_atoms (ase.Atoms): Reference structure containing atoms to append
        elements (set or None): If not None, append only these elements
        supercell (3-tuple): Dimensions of ref_atoms supercell. If None, try
            to calculate this instead.

    Returns:
        [ase.Atoms, ...]
        Composite structures as list of Atoms objects

    """

    if elements is None:
        pass
    else:
        cleaned_ref_atoms = ase.Atoms("", cell=ref_atoms.cell)
        [cleaned_ref_atoms.append(atom) for atom in
             filter(lambda x: (x.symbol in elements), ref_atoms)]
        ref_atoms = cleaned_ref_atoms

    if supercell is None:
          sc_dim = (traj[0].cell.diagonal() / ref_atoms.cell.diagonal())
          supercell = sc_dim.astype(int)

    ref_atoms = ref_atoms.repeat(supercell)    
    for atoms in traj:
        atoms.extend(ref_atoms)

    return traj

    
def main():
    parser = argparse.ArgumentParser()

    parser.add_argument('structure', type=str,
                        help="Input structure or trajectory file")
    parser.add_argument('add_atoms', type=str,
                        help="Structure file containing atoms to be added")
    parser.add_argument('-o', '--output', type=str, default=None,
                        help="Output filename")
    parser.add_argument('-s', '--supercell', type=int, nargs=3,
                        default=None)
    parser.add_argument('-e', '--elements', type=str, nargs='?', default=None)


    args = parser.parse_args()

    if args.elements is None:
        elements = None
    else:
        elements = set(args.elements)

    add_atoms_files(args.structure, args.add_atoms, output=args.output,
                    elements=elements, supercell=args.supercell)

if __name__ == '__main__':
    main()

