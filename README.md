# MultiDimInt

C++ library for creating discrete functions on a multi-dimensional grid with various available coordinate spacings

## Installation

### Install dependencies
Before compiling the code, the following dependencies need to be installed:

- [GCC](https://gcc.gnu.org/)
- [GNU Make](https://www.gnu.org/software/make/)
- [GNU Scientific Library](https://www.gnu.org/software/gsl/) (GSL)
- [Cuba](https://feynarts.de/cuba/)
- [Cubature](https://github.com/stevengj/cubature)

### Download and compilation

Clone the repository into the desired location and change into the root directory by running

```bash
git clone https://github.com/rlilow/MultiDimInt.git
cd MultiDimInt
```

If the header or library files of the GSL, Cuba or Cubature cannot be found in the standard search paths, open `Makefile` and set the variables `GSL_INCLUDE_PATH`, `GSL_LIB_PATH`, `CUBA_INCLUDE_PATH`, `CUBA_LIB_PATH`, `CUBATURE_INCLUDE_PATH` and `CUBATURE_LIB_PATH` to the appropriate paths. Do not modify the rest of the file.

Afterwards run

```bash
make
```

## Usage

To use MultiDimInt in your own code, link your program to the static library file `libmultidimint.a` and include the header file `MultiDimInt.hpp`. Both files are in the root directory.

A few small programs demonstrating the usage of MultiDimInt can be found in the directory `demo`. If you modify these, just re-run `make` in the root directory to rebuild them.

## Documentation 

If you have Doxygen (https://www.doxygen.nl/index.html) installed, you can build a detailed documentation of the different classes and functions in CORAS by running

```bash
make doc
```

from within the root directory.
This will create the directory `doc` containing the documentation.
To view it, open the file `doc/documentation.html` in the browser.