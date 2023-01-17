# Exercise 4: N particle system energy minimiser
Program to find the minimum energy configuration in 3D space for a system of N identical particles for Part II Chemistry Programming Practical Exercise 4.

## Commands
* `minimise <potential> [<strength>] [<max_steps>]` Minimises energy of a system defined by an XYZ file using the specified potential model and the default parameters. This will prompt for a `.xyz` file in the `/xyz_files` subdirectory to use as an input, and then prompt for the name of the file to save the output to.
    * `<potential>`   The model for potential energy between two points to use. See below for a list of available potential functions.
    * `<strength>`  (optional) Value of r_e/σ for the Morse potential. Defaults to 1
    * `<max_steps>` (optional) Maximum number of iterations before program should stop regardless of whether local minimum has been reached. Defaults to 0 (no maximum)
* `minimise_custom <potential> <delta> <lambda> <threshold> [<max_steps>] [<strength>]` Minimises energy of a system defined by an XYZ file using the specified potential model and the specified parameters. This will prompt for a `.xyz` file in the `/xyz_files` subdirectory to use as an input, and then prompt for the name of the file to save the output to.
    * `<potential>`   The model for potential energy between two points to use. See below for a list of available potential functions.
    * `<delta>` Value of delta to use when approximating local gradient of the potential energy
    * `<lambda>` Value of lambda used to determine step size in direction towards local minimum
    * `<threshold>` Value below which a gradient is considered to be 0 and so at a minimum
    * `<strength>`  (optional) Value of r_e/σ for the Morse potential. Defaults to 1
    * `<max_steps>` (optional) Maximum number of iterations before program should stop regardless of whether local minimum has been reached. Defaults to 0 (no maximum)
* `bond_lengths`    Prints out table of inter-particle distances for an XYZ file. This will prompt for a `.xyz` file in the `/xyz_files` subdirectory to use as an input. 
* `energy <potential> [<strength>]` Gets the energy of a system from an XYZ file according to the specified potential.
    * `<potential>`   The model for potential energy between two points to use. See below for a list of available potential functions.
    * `<strength>`  (optional) Value of r_e/σ for the Morse potential. Defaults to 1
* `plot`    Generates a 3D representation of an XYZ file
* `help [<command>]`    Displays a list of available commands. If `<command>` is specified, returns syntax information for specific command
    * `<command>`   (optional) Command to return syntax information for
* `quit`    Exits the program

## Available potential energy functions
### Lennard-Jones
The Lennard-Jones potential between 2 particles a distance r apart is defined as
![alt text](https://github.com/a-lapsley/potential-energy-surfaces/raw/main/img/lennard_jones.PNG "Lennard -Jones formula")
The energies are expressed in terms of the parameter ε (epsilon) and lengths in terms of the parameter σ (sigma).

### Morse
The Morse potential between 2 particles a distance r apart is defined as
![alt text](https://github.com/a-lapsley/potential-energy-surfaces/raw/main/img/morse.PNG "Morse formula")
The energies are expressed in terms of the parameter D_e and lengths in terms of σ (sigma). 
The parameter r_e/σ specifies the strength of the potential and its value can be changed in the program. 

## default_config.json
This file specifies the default parameters to use when minimising the energy of a system with the `minimise` command. Parameters for each potential energy function can be specified separately. 
* `delta` Value of delta to use when approximating local gradient of the potential energy
* `lambda`  Value of lambda used to determine step size in direction towards local minimum
* `threshold`   Value below which a gradient is considered to be 0 and so at a minimum

## Included .xyz files
Some .xyz files are included as inputs and examples of outputs.

Files beginning `input_(N)` contain N arbritrarily placed points as starting points for minimisation. For N > 3, these points all lie on the same plane if the suffix of the file is `_flat` and otherwise do not all lie on the same plane. The planar and non-planar input files may lead to different outputs when minimised as the planar points essentially remain constrained to the plane they start in so the program will generate the minimum energy conformation within that plane, whereas the non-planar input file may lead to a 3-dimensional structure. 

Files beginning with `ex_` are example output files that have already been generated using the provided input files. The suffix gives the potential function used to generated the example file, with `lennard` corresponding to Lennard-Jones, `morse1` a Morse potential with r_e/sigma = 1 and `morse2` a Morse potential with r_e/sigma = 2. 
