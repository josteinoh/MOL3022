# preamble
# region

# this code was written, in its entirety, by Jostein Haraldstad and Frantisek Nentwich, for MOL3022 @ NTNU in the spring semester of 2022
# it should be sourced from https://github.com/josteinoh/MOL3022

# the purpose of this code is to predict the secondary structure of proteins based on their amino acid sequences
# the code is based on a method devised by Chou and Fasman, details can be found at https://archive.ics.uci.edu/ml/datasets/Molecular+Biology+(Protein+Secondary+Structure) and its referenced papers
# more info can be found in the code's attached README file (README.md)

# endregion


# imports
from os import path


# default values for the data sets, based on the Chou-Fasman method
# region
weights = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
helix_formers     = {"E" : 1.37*weights[0], "A" : 1.29*weights[1], "L" : 1.20*weights[2], "H" : 1.11*weights[3], "M" : 1.07*weights[4], "Q" : 1.04*weights[5], "W" : 1.02*weights[6], "V" : 1.02*weights[7], "F" : 1.00*weights[8]}
helix_high_indiff = {"K" : 0.54*weights[9], "I" : 0.50*weights[10]}
helix_breakers    = {"N" : 1.00*weights[11], "Y" : 1.20*weights[12], "P" : 1.24*weights[13], "G" : 1.38*weights[14]}
helix_indiff      = ["K", "I", "D", "T", "S", "R", "C"]
sheet_formers     = {"M" : 1.40*weights[15], "V" : 1.39*weights[16], "I" : 1.34*weights[17], "C" : 1.09*weights[18], "Y" : 1.08*weights[19], "F" : 1.07*weights[20], "Q" : 1.03*weights[21], "L" : 1.02*weights[22], "T" : 1.01*weights[23], "W" : 1.00*weights[24]}
sheet_breakers    = {"K" : 1.00*weights[25], "S" : 1.03*weights[26], "H" : 1.04*weights[27], "N" : 1.14*weights[28], "P" : 1.19*weights[29], "E" : 2.00*weights[30]}
sheet_indiff      = ["A", "R", "G", "D"]
# endregion


# returns the dictionaries used to hold data, updated with the current individual adjustment coefficients
def create_dicts(arg_weight_list):

    h_formers     = {"E" : 1.37*arg_weight_list[0], "A" : 1.29*arg_weight_list[1], "L" : 1.20*arg_weight_list[2], "H" : 1.11*arg_weight_list[3], "M" : 1.07*arg_weight_list[4], "Q" : 1.04*arg_weight_list[5], "W" : 1.02*arg_weight_list[6], "V" : 1.02*arg_weight_list[7], "F" : 1.00*arg_weight_list[8]}
    h_high_indiff = {"K" : 0.54*arg_weight_list[9], "I" : 0.50*arg_weight_list[10]}
    h_breakers    = {"N" : 1.00*arg_weight_list[11], "Y" : 1.20*arg_weight_list[12], "P" : 1.24*arg_weight_list[13], "G" : 1.38*arg_weight_list[14]}
    h_indiff      = ["K", "I", "D", "T", "S", "R", "C"]
    s_formers     = {"M" : 1.40*arg_weight_list[15], "V" : 1.39*arg_weight_list[16], "I" : 1.34*arg_weight_list[17], "C" : 1.09*arg_weight_list[18], "Y" : 1.08*arg_weight_list[19], "F" : 1.07*arg_weight_list[20], "Q" : 1.03*arg_weight_list[21], "L" : 1.02*arg_weight_list[22], "T" : 1.01*arg_weight_list[23], "W" : 1.00*arg_weight_list[24]}
    s_breakers    = {"K" : 1.00*arg_weight_list[25], "S" : 1.03*arg_weight_list[26], "H" : 1.04*arg_weight_list[27], "N" : 1.14*arg_weight_list[28], "P" : 1.19*arg_weight_list[29], "E" : 2.00*arg_weight_list[30]}
    s_indiff      = ["A", "R", "G", "D"]

    return h_formers, h_high_indiff, h_breakers, h_indiff, s_formers, s_breakers, s_indiff


# returns the dictionaries used to hold data, updated with the current individual adjustment coefficients from a file
def create_dicts_file(arg_filename, enable_console=True):

    filename = arg_filename

    if enable_console:
        console_filename = input("Input path + name + extension of the input text file containing optimized weights, or press enter without input to use the default name (\"weights.txt\"): ")
        if console_filename:
            filename = console_filename

    weight_list = []

    if not path.exists(filename):

       write_weights(filename, [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1], False)

    with open(arg_filename, "r") as f:

        for line in f:
            if line != "\n":
                weight_list.append(float(line.strip()))

    h_formers     = {"E" : 1.37*weight_list[0], "A" : 1.29*weight_list[1], "L" : 1.20*weight_list[2], "H" : 1.11*weight_list[3], "M" : 1.07*weight_list[4], "Q" : 1.04*weight_list[5], "W" : 1.02*weight_list[6], "V" : 1.02*weight_list[7], "F" : 1.00*weight_list[8]}
    h_high_indiff = {"K" : 0.54*weight_list[9], "I" : 0.50*weight_list[10]}
    h_breakers    = {"N" : 1.00*weight_list[11], "Y" : 1.20*weight_list[12], "P" : 1.24*weight_list[13], "G" : 1.38*weight_list[14]}
    h_indiff      = ["K", "I", "D", "T", "S", "R", "C"]
    s_formers     = {"M" : 1.40*weight_list[15], "V" : 1.39*weight_list[16], "I" : 1.34*weight_list[17], "C" : 1.09*weight_list[18], "Y" : 1.08*weight_list[19], "F" : 1.07*weight_list[20], "Q" : 1.03*weight_list[21], "L" : 1.02*weight_list[22], "T" : 1.01*weight_list[23], "W" : 1.00*weight_list[24]}
    s_breakers    = {"K" : 1.00*weight_list[25], "S" : 1.03*weight_list[26], "H" : 1.04*weight_list[27], "N" : 1.14*weight_list[28], "P" : 1.19*weight_list[29], "E" : 2.00*weight_list[30]}
    s_indiff      = ["A", "R", "G", "D"]

    return h_formers, h_high_indiff, h_breakers, h_indiff, s_formers, s_breakers, s_indiff


# gets a value from a certain dictionary, returning 0 if the specified key does not exist in the dictionary
def get_value(arg_key, arg_dict):
    
    if arg_key in arg_dict:
        return arg_dict[arg_key]
    return 0


# gets data for training or testing from a file, returning two arrays of strings
# one array holds each test's amino acids, the other holds each test's structures
def read_file(arg_filename, input_string, enable_console=True):

    filename = arg_filename

    if enable_console:
        console_filename = input(input_string)
        if console_filename:
            filename = console_filename

    output_amino_strings = []
    output_structure_strings = []
    current_amino_string = "Z"*6
    current_structure_string = "z"*6
    begin = False

    with open(filename, "r") as f:
        for line in f:

            if line[:2] == "<>":
                begin = True

            elif line[:5] == "<end>" or line[:3] == "end":
                begin = False
                current_amino_string += "Z"*6
                current_structure_string += "z"*6
                output_amino_strings.append(current_amino_string)
                output_structure_strings.append(current_structure_string)
                current_amino_string = "Z"*6
                current_structure_string = "z"*6

            else:
                if begin:
                    current_amino_string += line[0]
                    if len(line) > 2:
                        current_structure_string += line[2]

    return output_amino_strings, output_structure_strings


# writes the list of weights to a file for saving
def write_weights(arg_filename, arg_weight_list, enable_console=True):

    filename = arg_filename

    if enable_console:
        console_filename = input("Input path + name + extension of the target text file for saving weights, or press enter without input to use the default name (\"weights.txt\"): ")
        if console_filename:
            filename = console_filename

    with open(filename, "w") as f:

        for weight in arg_weight_list:
            if 0.9999 < weight < 1.0001:
                f.write(str(1))
            else:
                f.write(str(weight))
            f.write("\n")

    return filename


# writes the list of aminos + structures to a file for saving
def write_structures(arg_filename, arg_aminos, arg_structures, enable_console=True):

    filename = arg_filename

    if enable_console:
        console_filename = input("Input path + name + extension of the target text file for saving predicted structures, or press enter without input to use the default name (\"predicted_structures.txt\"): ")
        if console_filename:
            filename = console_filename

    with open(filename, "w") as f:

        for i in range (len(arg_structures)):
            f.write("\n<>")

            for j in range (len(arg_structures[i])):
                f.write("\n")
                f.write(arg_aminos[i][j+6])
                f.write(" ")
                f.write(arg_structures[i][j])

            f.write("\n<end>")

    return filename


# determines if a helix structure is to be broken, if the current structure is a helix
# the input is a string of length 3: the current amino acid, the previous one, and the next one
def terminate_helix(arg_string):
    
    if (arg_string[1] in helix_breakers) and ( ( (arg_string[0] in helix_breakers) or (arg_string[0] in helix_indiff) ) or (arg_string[2] in helix_breakers) or (arg_string[2] in helix_indiff) ):
        return True
    return False


# determines if a sheet structure is to be broken, if the current structure is a sheet
# the input is a string of length 3: the current amino acid, the previous one, and the next one
def terminate_sheet(arg_string):
    
    if (arg_string[1] in sheet_breakers) and ( ( (arg_string[0] in sheet_breakers) or (arg_string[0] in sheet_indiff) ) or (arg_string[2] in sheet_breakers) or (arg_string[2] in sheet_indiff) ):
        return True
    return False


# determines if a helix structure is to be continued, if the current structure is a helix
# the input is a string of length 3: the current amino acid, the previous one, and the next one
def cont_helix(arg_string):
    
    if arg_string[1] != "P" and ( not terminate_helix(arg_string)):
        return True
    return False


# determines if a sheet structure is to be continued, if the current structure is a sheet
# the input is a string of length 3: the current amino acid, the previous one, and the next one
def cont_sheet(arg_string):
    
    if arg_string[1] != "P" and arg_string[1] != "E" and ( not terminate_sheet(arg_string)):
        return True
    return False


# determines if a helix structure is to be started
# the input is a string of length 13: the five previous, as well as the current, and next five, amino acids
def init_helix(arg_string):
    
    sum_form = 0
    sum_break = 0
    for i in range(len(arg_string)):
        if i == 6:
            continue
        sum_form  += get_value(arg_string[i], helix_formers)
        sum_form  += get_value(arg_string[i], helix_high_indiff)
        sum_break += get_value(arg_string[i], helix_breakers)
    return ((sum_form >= 8) and (sum_break < 4)), sum_form


# determines if a sheet structure is to be started
# the input is a string of length 11: the four previous, as well as the current, and next four, amino acids
def init_sheet(arg_string):

    sum_form = 0
    sum_break = 0
    for i in range(len(arg_string)):
        if i == 5:
            continue
        sum_form  += get_value(arg_string[i], sheet_formers)
        sum_break += get_value(arg_string[i], sheet_breakers)
    return ((sum_form >= 6) and (sum_break < 4)), sum_form


# determines the structure of the current amino acid
# the input is a string of length 13: the five previous, as well as the current, and next five, amino acids
def predict_structure(arg_string, last_structure):

    init_h, score_h = init_helix(arg_string)
    init_e, score_e = init_sheet(arg_string[1:-1])

    if (init_h and init_e):
        if score_h > score_e:
            init_e = False
        else:
            init_h = False
    
    if init_h or (cont_helix(arg_string[5:8]) and last_structure == "h"):
        return "h"
    elif init_e or (cont_sheet(arg_string[5:8]) and last_structure == "e"):
        return "e"
    elif ( last_structure == "h" and terminate_helix(arg_string[5:8]) ) or ( last_structure == "e" and terminate_sheet(arg_string[5:8]) ) or ( last_structure == "_" and ( not init_h ) and ( not init_e ) ):
        return "_"
    else:
        return "_"


# runs a series of structure predictions on a single sequence of amino acids, returning a string containing the structures
# the input is a string of amino acids of any length
def run_sequence(arg_amino):

    return_structure = "_"
    predicted_structure = ""
    for i in range ( 6, len(arg_amino) - 6 ):
        predicted_structure = predict_structure(arg_amino[i-6+1:i+6], return_structure[-1])
        return_structure += predicted_structure
    
    return return_structure[1:]


# compares two structures, returning how many of the structures are the same
# the inputs are two structures
def compare_structures(arg_predicted_structures, arg_correct_structures):

    hits = 0
    for i in range ( min( len(arg_predicted_structures), len(arg_correct_structures) ) ):
        hits += int( arg_predicted_structures[i] == arg_correct_structures[i + 6] )

    return hits


# runs a single test
# returns the score and number of tries
def run_full_test(arg_amino_strings, arg_structure_strings):

    predicted_structures = []
    scores = []
    total_score = 0
    total_tries = 0

    for i in range (len(arg_amino_strings)):
        predicted_structures.append(run_sequence(arg_amino_strings[i]))

    for j in range (len(predicted_structures)):
        scores.append([ compare_structures(predicted_structures[j], arg_structure_strings[j]), len(predicted_structures[j]) ])

    for k in range (len(scores)):
        total_score += scores[k][0]
        total_tries += scores[k][1]
    return total_score, total_tries


# runs a single test
# returns the predicted structures
def run_prediction(arg_amino_strings):

    predicted_structures = []

    for i in range (len(arg_amino_strings)):
        predicted_structures.append(run_sequence(arg_amino_strings[i]))

    return predicted_structures


# runs the optimization algorithm for the weights for the coefficients
# returns the best score, number of tries, and the optimal weights
def run_optimization(filename, enable_console=True):

    file_amino_strings, file_structure_strings = read_file(filename, "Input path + name + extension of the input text file for training, or press enter without input to use the default name (\"train.txt\"): ", enable_console)

    weight_list = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
    global helix_formers
    global helix_high_indiff
    global helix_breakers
    global helix_indiff
    global sheet_formers
    global sheet_breakers
    global sheet_indiff

    current_score = 0
    best_score = 0

    helix_formers, helix_high_indiff, helix_breakers, helix_indiff, sheet_formers, sheet_breakers, sheet_indiff = create_dicts(weight_list)

    for i in range (1, 4):
        step = 10**(-i)

        for j in range (len(weight_list)):

            best_score, tries = run_full_test(file_amino_strings, file_structure_strings)
            current_score = best_score + 1
            while current_score > best_score:
                weight_list[j] += step
                helix_formers, helix_high_indiff, helix_breakers, helix_indiff, sheet_formers, sheet_breakers, sheet_indiff = create_dicts(weight_list)
                current_score, tries = run_full_test(file_amino_strings, file_structure_strings)
                if current_score > best_score:
                    best_score = current_score
                else:
                    weight_list[j] -= step
                    helix_formers, helix_high_indiff, helix_breakers, helix_indiff, sheet_formers, sheet_breakers, sheet_indiff = create_dicts(weight_list)

            current_score = best_score + 1
            while current_score > best_score:
                weight_list[j] -= step
                helix_formers, helix_high_indiff, helix_breakers, helix_indiff, sheet_formers, sheet_breakers, sheet_indiff = create_dicts(weight_list)
                current_score, tries = run_full_test(file_amino_strings, file_structure_strings)
                if current_score > best_score:
                    best_score = current_score
                else:
                    weight_list[j] += step
                    helix_formers, helix_high_indiff, helix_breakers, helix_indiff, sheet_formers, sheet_breakers, sheet_indiff = create_dicts(weight_list)

    return best_score, tries, weight_list[:]


# runs the optimization algorithm for the weights for the coefficients, saving the coefficients to a file
# returns the best score, number of tries, and the optimal weights
def main_opt(arg_filename="train.txt", enable_console=True):

    global helix_formers
    global helix_high_indiff
    global helix_breakers
    global helix_indiff
    global sheet_formers
    global sheet_breakers
    global sheet_indiff

    print()

    filename = arg_filename

    if enable_console:
        console_filename = input("Input path + name + extension of the input text file for training, or press enter without input to use the default name (\"train.txt\"): ")
        if console_filename:
            filename = console_filename

    print()
    print("Optimizing...")
    print()
    opt_best_score, opt_tries, weight_list = run_optimization(filename, False)
    filename_weights = write_weights("weights.txt", weight_list)
    print()
    print("Best optimized score: ",  str(opt_best_score).rjust(8))
    print("Number of amino acids:",  str(opt_tries).rjust(8))
    print("Success rate:         ", (str(round(opt_best_score*100/opt_tries, 2)) + "%").rjust(8))
    print()
    print("Optimized weights successfully written to " + filename_weights)

    return opt_best_score, opt_tries, weight_list


# runs the algorithm against a known sequence + structure, writing the success rate to console
# returns the score and the number of tries
def main_test(arg_filename="test.txt", enable_console=True):

    global helix_formers
    global helix_high_indiff
    global helix_breakers
    global helix_indiff
    global sheet_formers
    global sheet_breakers
    global sheet_indiff

    print()

    helix_formers, helix_high_indiff, helix_breakers, helix_indiff, sheet_formers, sheet_breakers, sheet_indiff = create_dicts_file("weights.txt")

    file_amino_strings, file_structure_strings = read_file(arg_filename, "Input path + name + extension of the input text file for testing, or press enter without input to use the default name (\"test.txt\"): ", enable_console)
    
    total_score, total_tries = run_full_test(file_amino_strings, file_structure_strings)

    print()
    print("Final score in test run:",  str(total_score).rjust(8))
    print("Number of amino acids:  ",  str(total_tries).rjust(8))
    print("Success rate:           ", (str(round(total_score*100/total_tries, 2)) + "%").rjust(8))

    return total_score, total_tries


# runs the algorithm against a sequence with unknown structure, writing the predicted structures to a file
# does not return anything
def main_run(arg_filename="data.txt", enable_console=True):

    global helix_formers
    global helix_high_indiff
    global helix_breakers
    global helix_indiff
    global sheet_formers
    global sheet_breakers
    global sheet_indiff

    print()

    helix_formers, helix_high_indiff, helix_breakers, helix_indiff, sheet_formers, sheet_breakers, sheet_indiff = create_dicts_file("weights.txt")

    file_amino_strings, file_structure_strings = read_file(arg_filename, "Input path + name + extension of the input text file for predicting structures, or press enter without input to use the default name (\"data.txt\"): ", enable_console)
    
    predicted_structures = run_prediction(file_amino_strings)

    filename_structures = write_structures("predicted_structures.txt", file_amino_strings, predicted_structures)

    print()
    print("Predicted structures successfully written to " + filename_structures)


# main program
# gives the user the choice between optimization, testing, and predicting
def main():

    choice = " "

    print()
    print("Type \"opt\" to run the optimization algorithm and save the optimized weights to a file.")
    print("Type \"test\" to to test the optimized structure prediction algorithm on an amino sequence with known structures and test its success rate.")
    print("Type \"run\" to run the optimized structure prediction algorithm and save the structures to a file.")
    print("Hit enter without typing to exit the program.")
    print("For more info, read README.md")

    while choice:

        print()
        choice = input("Choose functionality (opt/test/run): ").lower()

        if choice == "opt":
            main_opt()
        elif choice == "test":
            main_test()
        elif choice == "run":
            main_run()
        elif choice:
            print()
            print("Invalid choice, please try again.")
        else:
            print()
            print("Exiting program.")
            print()
            break

    return 0


main()

