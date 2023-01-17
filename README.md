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
The Lennard-Jones potential is defined as
![alt text]https://github.com/a-lapsley/potential-energy-surfaces/raw/main/img/lennard_jones.PNG "Lennard -Jones formula" 


## default_config.json
This file specifies the default parameters to use when minimising the energy of a system with the `minimise` command. Parameters for each potential energy function can be specified separately. 
* `delta` Value of delta to use when approximating local gradient of the potential energy
* `lambda`  Value of lambda used to determine step size in direction towards local minimum
* `threshold`   Value below which a gradient is considered to be 0 and so at a minimum


