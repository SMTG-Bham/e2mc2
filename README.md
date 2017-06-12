# E2MC2: "Easy Easy Monte Carlo Code"

This is a little Python API for working with the very useful Monte
Carlo tools included
in
[ATAT](https://www.brown.edu/Departments/Engineering/Labs/avdw/atat/).

## Requirements
- Python3 with its standard library
- A Unix-like filesystem
- ATAT should be on the system PATH

## Installation
In the project directory (containing this README):

    pip3 install --user .

You may need to add your local installation directories to `PATH` and
`PYTHONPATH` if this has not already been done; these will be
something like `~/.local/bin` and `~/.local/python3/site-packages`.

## Goals

The goal is to make it easier to manage sets of many Monte Carlo runs
for convergence testing etc. EMC2 expects to be run in a directory
containing a set of named files, and to freely generate a set of named
files. E2MC2 should be able to serialise and collect the contents of
these files, making it a bit cleaner to manage the inputs and outputs.
JSON and TAR formats are preferred for serialisation.

## Structure

The data management and serialisation should be object-oriented, with
two significant Python classes:

**ClusterExpansion** contains the CE parameterisation and lattice
information needed to initialise a Monte Carlo run. This data is
generally collected from MAPS.

**MonteCarloCalc** collects calculation parameters, and provides the
*run* method which manages the actual calculation. The calculation
results are then attached to the object. To restore the state for
further analysis, it should be possible to initialise this object from
a completed calculation directory.

Further analysis should consist of reasonably pure functions which
operate on the above objects.


## Licensing / Disclaimer

This project is not affiliated with ATAT, Northwestern University,
Brown University or Axel van de Walle.

This project is licensed under the GPLv3
