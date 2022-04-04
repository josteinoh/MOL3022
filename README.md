

# README.md

## Ownership and credit

The contents of this README and of `main.py` are, in their entirety, written by Jostein Haraldstad and František Nentwich, for the course MOL3022 at NTNU in the spring semester of 2022.
All files used should be sourced from https://github.com/josteinoh/MOL3022. The repository can be cloned as a whole (recommended), or the files can be downloaded individually. The method used is based on a method by Maclin and Shavlik, which in turn is based on a method by Chou and Fasman. The data used is provided by Qian and Sejnowski. Details on this can be read below.  

## Purpose and functionality

The purpose of this code is to predict the secondary structure of proteins - alpha-helix, beta-sheet or random coil - based on their amino acid sequences. The code is based on a modified implementation of the Chou-Fasman method, devised by Maclin and Shavlik, which can be read about in detail at https://archive.ics.uci.edu/ml/datasets/Molecular+Biology+(Protein+Secondary+Structure) and in its referenced papers. The data used is supplied by Qian and Sejnowski.  

In order to predict the structure for each amino acid in a sequence, the algorithm looks at a window of amino acids around the current amino acid being considered. It considers whether the amino acids in the window have a high chance of starting, ending, or continuing each of the structure types, and assigns a score to each of those options. Based on a comparison of these scores, as well as the structures of the neighboring amino acids, it determines the structure of the amino acid currently being considered. This is then repeated for the full length of every given sequence of amino acids.  

The methods used for considering and comparing the amino acids and their abilities to start/end each structure are used the same way as Maclin and Shavlik (as well as Chou and Fasman) described, with the exception of the size of the window of amino acids being analyzed. Testing revealed that with our data sets, a window of length 11 was most effective (the current amino acid, as well as the five previous and following amino acids).  

Maclin and Shavlik use Qian and Sejnowski's data to assign each amino acid values for how much that acid affects the starting and ending of the three different structures. However, these can be adjusted to attain somewhat better results, and `main.py` contains an optimization algorithm that adjusts these values, in order for the structure prediction algorithm to produce the best results. This optimization algorithm is a simple implementation of a gradient algorithm, as it adjusts each of the values until the results (i.e. the success rate of the structure prediction algorithm) no longer improve.  

## Different modes of usage

When `main.py` is ran, four options will be given:  

- Optimization, chosen by typing "opt" in the console as the program is running. This runs the optimization algorithm, using a specified text file containing amino acids with known structures for training. The program writes the best success rate to console and saves the adjustments for the amino acid values to a text file for later use.  

- Testing, chosen by typing "test" in the console as the program is running. This runs the prediction algorithm on a specified text file containing amino acid sequences *with* known structures. It uses the optimized adjustments for the amino acid values from the corresponding file, or if no such file exists, creates one with default values. It does not save anything to a file as the structures are already known, but it writes the success rate to console.  

- Running, chosen by typing "run" in the console as the program is running. This runs the prediction algorithm on a specified text file containing amino acid sequences *without* known structures. It uses the optimized adjustments for the amino acid values from the corresponding file, or if no such file exists, creates one with default values. It does not write anything to the console as the structures are unknown and success rate cannot be determined, but it saves the structures to a file, using the same formatting as the training and testing files.  

- Exiting the program, chosen by hitting the enter key without typing anything in the console as the program is running.  

## File formatting

The text files used as data for training and testing must be correctly formatted. They can contain any number of amino acid sequences, and these can have any length. The first line in each sequence must be "<>", and the last line must be "\<end>". All lines between these belong to the sequence, and each line represents one amino acid, and possibly its structure. These lines start with a capital letter representing the acid, and if the structure is known, followed by a space and a symbol representing the structure: "h" for helix, "e" for sheet, and "_" for coil. There can be excess lines of text not belonging to any of the sequences before the first "<>" or after the last "\<end>", but not between any of the sequences. Examples of formatting can be seen in `data.txt`, `train.txt`, and `test.txt`, and the names of these files should be used to minimize manual specification of filenames and paths when the program is used.  

If desired, a custom file for adjusting the amino acid values can also be used. The formatting here is simple: 31 lines of text, each consisting exclusively of a number (typically between 0.8 and 1.2) to linearly scale each of the amino acid values by. Examples of this can be seen in `weights.txt`, which also is the default name for this file.  

## Considerations to take before running

Before running, make sure that some version of Python 3 is installed. If the program is to be run from console, make sure that Python is added to your PATH. For more info on this, visit https://www.python.org/downloads/release/python-3100/. To clone the repo from git, git must be installed. For info on how to download and use git, visit https://git-scm.com/. Finally, make sure the repo in its entirety is downloaded, including the text files if needed as data sources.  


## Considerations to take while running

If the program is run via an IDE, open `main.py` **from the main folder of the cloned repo**. If the program is run from a console, **change directory to the main folder of the cloned repo**. On Windows and Linux, this is done by running the command `cd <full path to folder>`. Then, the program can be executed by running it in the IDE, or by typing `py -3 main.py` in the console. From there, the program will be executed, and the functionality as described in this README can be used. Whenever a filename+path is requested, a default name can be used by hitting the enter key without specifying a name. For the data source files, the default names are the names used for the text files in the repo: "train.txt" for training, "test.txt" for testing, "data.txt" for predicting unknown structures, and "weights.txt" for assigning new values to the amino acids without optimizing. There is also a default name for the file that is created to contain predicted structures when the real structures are unknown, "predicted_structures.txt". It should be noted that these default names **will not work** if `main.py` is not properly ran from inside the main repo folder as described earlier. If this is not done, or if errors regarding the location and/or existence of certain files are thrown, the files must be specified with their full path, name, and extension.

