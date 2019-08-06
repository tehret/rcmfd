ROBUST COPY-MOVE FORGERY DETECTION BY FALSE ALARM CONTROL
=========================================================

* Author    : EHRET Thibaud <ehret.thibaud@gmail.com>
* Licence   : AGPL v3+, see agpl.txt

OVERVIEW
--------

This source code provides an implementation of "T. Ehret, ROBUST COPY-MOVE FORGERY 
DETECTION BY FALSE ALARM CONTROL, Arxiv 2019" available [here](https://arxiv.org/pdf/1906.00649.pdf).
Please cite it if you use this code as part of your research.
The structure of this code is derived from 
[Fast Image Matching by Affine Simulation](https://github.com/rdguez-mariano/fast_imas_IPOL).

COMPILATION
-----------

The code is compilable on Unix/Linux and hopefully on Mac OS (not tested!). 

**Compilation:** requires the cmake and make programs.

**Dependencies:** 
For image i/o we use [Enric Meinhardt's iio](https://github.com/mnhrdt/iio),
which requires libpng, libtiff and libjpeg.
 
Configure and compile the source code using cmake and make.  It is recommended
that you create a folder for building:

UNIX/LINUX/MAC:
```
$ mkdir build; cd build
$ cmake ..
$ make
```

Binaries will be created in `build`.

USAGE
-----

The following commands have to be run from the `build` folder:

List all available options:</br>
```./main --help```

There is only one mandatory input argument:
* `-i` the input image

Optional arguments are:
* `-o` the path to the output file containting the matches (string)
* `-ps` the patch size of the descriptor (int)
* `-tau` the threshold for the patch matching, also correspond to the NFA threshold when the automatic thershold is activated (float)
* `-gs` convert the the image to grayscale (bool)
* `-auto` Use the threshold defined by NFA (using tau as parameter), otherwise use the provided threshold (bool)

-----

The following command compute matches that could correspond to a copy-move forgery:

```./main -im ../data/forged.tif -auto true```

The image in `data` is from the COVERAGE dataset, available [here](https://github.com/wenbihan/coverage)
