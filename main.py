from copy import deepcopy
from Vec3d import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import TABLEAU_COLORS
from random import choice as rand_choice
import json

#File handling

def specify_input_file(name=""):
    #Handles user specifying input XYZ file
    if name != "":
        return "xyz_files\\%s.xyz" % name
    
    valid_file = False
    while not valid_file:
        try:
            print("Please enter name of input file: ")
            name = input().strip().replace(".xyz","")
            dir = "xyz_files\\%s.xyz" % name
            f = open(dir)
            f.close()
            valid_file = True
        except FileNotFoundError:
            print("File not found.")
        except:
            print("Unable to open file.")
    return dir

def specify_output_file(name=""):
    #Handles user specifying data output file

    if name != "":
        return "xyz_files\\%s.xyz" % name

    valid_file = False
    while not valid_file:
        try:
            print("Please enter name of output file: ")
            name = input().strip().replace(".xyz","")
            dir = "xyz_files\\%s.xyz" % name
            f = open(dir, "x")
            f.close()
            valid_file = True
        except FileExistsError:
            print("File already exists")
        except:
            print("Invalid file name")
    return dir 

def read_xyz_file(dir):
    """
    Reads an XYZ file for N atoms and returns an array of vec3d objects 
    containing the coordinates for each atom.
    """
    with open(dir, "r") as f:
        lines = f.readlines()
        n = int(lines[0].strip())
        output = []
        for i, line in enumerate(lines):
            if i > 1:
                coords = []
                for j, value in enumerate(line.split()):
                    if j != 0:
                        value = float(value)
                        coords.append(value)
                output.append(Vec3d(coords))
                        
    return output

def write_xyz_file(dir, info_line, vecs):
    """
    Creates a .XYZ file for the array of 3D vectors 'vecs', using the string 
    'info_line' as the info line of the XYZ file. For the purposes of this 
    program all particles are considered indentical and just labelled 'A'.
    """
    with open(dir, "w") as f:
        f.write("%i\n" % len(vecs))
        f.write(info_line + "\n")
        for vec in vecs:
            xstr = "%.4f" % vec[0]
            ystr = "%.4f" % vec[1]
            zstr = "%.4f" % vec[2]
            line = "A\t"
            line = line + "{0: <8}".format(xstr)
            line = line + "{0: <8}".format(ystr)
            line = line + "{0: <8}\n".format(zstr)
            f.write(line)

def load_command_syntax():
    #Gets information about commands from json file
    global COMMAND_SYNTAX
    with open("commands.json","r") as f:
        COMMAND_SYNTAX = json.load(f)

def load_default_parameters(potential):
    with open("default_config.json","r") as f:
        data = json.load(f)
    
    return data[potential]

#UX handling

def welcome_message():
    print("----------------")
    print("Particle system geometry minimiser")

def commands():
    print("----------------")
    print("Available commands:\n")
    for key in COMMAND_SYNTAX.keys():
        print("{0: <20}".format(key) + COMMAND_SYNTAX[key]["description"])

def command_input():
    #Parse user command
    valid = False
    while not valid:
        inp = input().lower().strip().split(" ")
        command = inp[0]
        args = inp[1:]
        if command == "minimise" or command == "minimise_custom":
            try:
                potential = args[0]
                strength = 1

                if command == "minimise":
                    params = load_default_parameters(potential)
                    delta = params["delta"]
                    scale_factor = params["lambda"]
                    threshold = params["threshold"]
                    arg_index = 1
                
                if command == "minimise_custom":
                    delta = float(args[1])
                    scale_factor = float(args[2])
                    threshold = float(args[3])
                    arg_index = 4
                
                if potential == "morse":
                    try:
                        strength = float(args[arg_index])
                        arg_index += 1
                    except:
                        strength = 1
                        arg_index += 1

                try:
                    max_steps = int(args[arg_index])
                except:
                    max_steps = 0
                
                input_dir = specify_input_file()
                output_dir = specify_output_file()

                start_vecs = read_xyz_file(input_dir)

                finish_vecs = iterate_to_minimum(
                    start_vecs,
                    delta,
                    scale_factor,
                    threshold,
                    potential,
                    max_steps,
                    strength
                )
                
                finish_vecs = centre(finish_vecs)

                info_line = "Generated XYZ file"
                info_line += "\tDelta: %f" % delta
                info_line += "\tLambda: %f" % scale_factor
                info_line += "\tThreshold: %f" % threshold

                write_xyz_file(output_dir, info_line, finish_vecs)

                valid = True
            except:
                print("Invalid command syntax")

        elif command == "plot":
            input_dir = specify_input_file()
            vecs = read_xyz_file(input_dir)
            plot(vecs)
            valid = True
        
        elif command == "energy":
            try:
                potential = args[0]
                potential_func = POTENTIALS[potential]
                if potential == "morse":
                    unit = "D_e"
                    try:
                        strength = float(args[1])
                    except:
                        strength = 1
                elif potential == "lennard_jones":
                    unit = "epsilon"
                    strength = 1

                input_dir = specify_input_file()
                vecs = read_xyz_file(input_dir)

                energy = get_potential(vecs, potential_func, strength)
                print("System energy: %f %s" % (energy, unit))
                valid = True
            except:
                print("Invalid command syntax")
            

        elif command == "bond_lengths":
            input_dir = specify_input_file()
            vecs = read_xyz_file(input_dir)
            bond_lengths(vecs)
            valid = True

        elif command == "help":
            if len(args) == 0:
                commands()
            else:
                try:
                    info  = COMMAND_SYNTAX[args[0]]
                    for key in info.keys():
                        print("{0: <16}".format(key) + info[key])
                except:
                    print("Unknown command")
        elif command == "quit":
            exit()
        else:
            print("Invalid command")

#Potential functions

def lennard_jones(r, strength=1):
    #Function to calculate Lennard Jones potential at distance r.
    #strength parameter is not used in this function but needs to be here
    #so the function arguments are the same as for other potentials
    u = 4 * ( (1 / r)**12 - (1 / r)**6 )
    return u

def morse(r, strength=1):
    #Function to calculate the Morse potential at distance r.
    #strength is the parameter r_e / sigma which specifies the strength of the 
    #potential
    u = (1 - np.exp(strength - r))**2
    return u

POTENTIALS = {
    "lennard_jones":lennard_jones,
    "morse":morse
}

#Computation functions

def add_delta(coord_arr, delta, i, j):
    #Takes an array of coordinates 'coord_arr' and adds value 'delta' to the 
    #jth coordinate of the ith vector. 
    #deepcopy() is used so that the original array is not overwritten
    arr = deepcopy(coord_arr)
    arr[i][j] = arr[i][j] + delta
    return arr

def get_potential(vecs, u, strength):
    """
    Gets the potential energy of a system of coordinates specified by 3D 
    coordinates in array of 3D vectors 'vecs', using the potential function 
    'u'. 'strength' is a parameter of the potential energy function.
    """
    total_energy = 0
    for i in range(len(vecs)):
        for j in range(i+1, len(vecs)):
            rij = vecs[i] - vecs[j]
            total_energy += u(rij.length(), strength)
    return total_energy

def gradient(
     coord_arr,
     delta, 
     vec_index, 
     coord_index, 
     potential_func, 
     strength
     ):
    """
    Gets the gradient of the potential energy at a set of coordinates specified
    by 'coord_arr' in the direction of one of the components (specified by 
    'coord_index') of one of the vectors (specified by 'vec_index').
    'potential' specifies the potential energy function to use, and 'strength'
    is a parameter for this function.

    The gradient is computed numerically, using a small step 'delta'. 
    """
    coord_array_plus = add_delta(coord_arr, delta, vec_index, coord_index)
    coord_array_minus = add_delta(coord_arr, -delta, vec_index, coord_index)
    u_plus = get_potential(coord_array_plus, potential_func, strength)
    u_minus = get_potential(coord_array_minus, potential_func, strength)
    grad = (u_plus - u_minus) / (2 * np.abs(delta))
    return grad

def step(coord_arr, delta, scale_factor, threshold, potential, strength):
    """
    Takes a set of coordinates 'coord_arr' and moves each coordinate towards a 
    local minimum by finding the gradient of the potential with respect to that
    coordinate and taking a small step in the opposite direction.

    Returns the number of coordinates which are not at a stationary point.

    'delta' is a small number used to calculate the gradient numerically.
    'scale_factor' is a small number used to determine how large the step is.
    'threshold' - if the gradient is below this then it is considered to be
    a stationary point. 
    'potential' specifies the potential function to use and 'strength' is a 
    parameter for this function.
    """
    count_over_threshold = 0
    delta_coord_arr = []
    potential_func = POTENTIALS[potential]
    for i, vec in enumerate(coord_arr):
        delta_vec = Vec3d([0.0,0.0,0.0])
        for j in range(len(vec)):
            grad = gradient(coord_arr, delta, i, j, potential_func, strength)
            delta_vec[j] = grad

            if np.abs(grad) > threshold:
                count_over_threshold += 1

        delta_coord_arr.append(delta_vec)
    
    for i in range(len(coord_arr)):
        coord_arr[i] = coord_arr[i] - scale_factor * delta_coord_arr[i]

    return count_over_threshold

def iterate_to_minimum(
                        coord_arr, 
                        delta, 
                        scale_factor, 
                        threshold,
                        potential, 
                        max_steps=0,
                        strength=1
                        ):
    """
    Takes a set of coordinates 'coord_arr' and repeatedly steps this towards 
    a local minimum using the step() function. Continues until all the
    potential in all directions is considered to be at a stationary point (
    i.e., the gradient is less than 'threshold'), or the maximum number of 
    iterations ('max_steps') is reached - if this is 0 it will continue
    indefinitely until a minimum is reached. 
    """
    count = 0
    minimum = False
    while not minimum:
        over_threshold = step(
                                coord_arr, 
                                delta, 
                                scale_factor, 
                                threshold,
                                potential,
                                strength
                            )
        if over_threshold == 0:
            minimum = True
          
        count += 1

        if max_steps != 0:
            if count > max_steps:
                print("Reached maximum number of steps, stopping process")
                return coord_arr

    print("Minimum reached after %i steps" % count)
    return coord_arr

def centre(vecs):
    #Centres the array of vectors 'vecs' so the centre of mass of the system
    #is at the origin

    n = len(vecs)
    sum = Vec3d([0, 0, 0])

    for i in range(n):
        sum += vecs[i]

    centre_of_mass = sum / n

    for i in range(n):
        vecs[i] = vecs[i] - centre_of_mass
    
    return vecs

def plot(coord_arr, range_upper=1, range_lower=-1):
    #Plots the set of 3D coordinates in 'coord_arr'
    #Also draws lines between them and labels the bond lengths

    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    xs = []
    ys = []
    zs = []
    for vec in coord_arr:
        xs.append(vec[0])
        ys.append(vec[1])
        zs.append(vec[2])
    
    for a_x, a_y, a_z in zip(xs, ys, zs):
        for b_x, b_y, b_z, in zip(xs, ys, zs):
            if a_x != b_x or a_y != b_y or a_z != b_z:
                col = TABLEAU_COLORS[rand_choice(list(TABLEAU_COLORS))]
                ax.plot([a_x, b_x], [a_y, b_y], [a_z, b_z], color=col)
                vec = Vec3d([a_x, a_y, a_z]) - Vec3d([b_x, b_y, b_z])
                ax.text(
                    (a_x + b_x) / 2,
                    (a_y + b_y) / 2,
                    (a_z + b_z) / 2,
                    "%.2f" % vec.length(),
                    color=col
                )
    ax.scatter(xs, ys, zs)
    ax.set_xlim3d(range_lower, range_upper)
    ax.set_ylim3d(range_lower, range_upper)
    ax.set_zlim3d(range_lower, range_upper)
    plt.show()

def bond_lengths(coord_arr):
    # Gets and prints the distances between each particle in the system
    for i in range(0, len(coord_arr)):
        for j in range(i+1, len(coord_arr)):
            vec = coord_arr[i] - coord_arr[j]
            print("Particle %i - Particle %i distance: %.3f" % \
                ((i+1, j+1, vec.length()))
                )



####################
#   Main Program   #
####################

#Startup

load_command_syntax()

#Program loop
while True:
    welcome_message()
    commands()
    command_input()
    input("Press enter to continue...")
