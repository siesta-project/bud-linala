[![Build Status](https://travis-ci.org/siesta-project/bud-linala.svg?branch=master)](https://travis-ci.org/siesta-project/bud-linala)[![codecov](https://codecov.io/gh/siesta-project/bud-linala/branch/master/graph/badge.svg)](https://codecov.io/gh/siesta-project/bud-linala)


# linala

Library for interfacing linear-algebra libraries to the buds library.

Currently it is recommended to use this as a plugin for the buds library.
Using linala as a plugin will make linala intrinsically usable in the buds library.
To realize the linala library as a plugin you only need to clone/extract this project into
the `plugins` directory using this command:

    git clone git@github.com:siesta-project/bud-linala.git linala

The you need to add this to the `setup.make` file:

    PLUGINS += linala

Then proceed with adding the dependencies on the library linker lines (`LIBS` variable
in the `setup.make` file).

## Requirements

linala depends on a number of linear algebra libraries. Currently these are
required:

- BLAS
- LAPACK
- ScaLAPACK (optional)
- ELPA (optional)


## Documentation

The complete documentation of linala may be found [here](http://siesta-project.github.io/bud-linala/).

## Download

To download the latest linala library please see the [release page](https://github.com/siesta-project/bud-linala/releases).

## Installation

The installation instructions of linala may be found [here](http://siesta-project.github.io/bud-linala/installInstructions.html).



