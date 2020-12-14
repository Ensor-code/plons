import os
import os.path
import sys
import datetime                     as dt
import numpy                        as np

# Create file if loc and outputLoc were not specified and save the file in userSettingsFilePath
# There are three quantities that need to be specified:
#   - the prefix of the data files (e.g. prefix_00001, prefix_00002, ...)
#   - loc
#   - outloc
def create(userSettingsFilePath):
    with open(userSettingsFilePath, "w") as file:
        prefix = str(input("Please enter the prefix that has been used to name all the PHANTOM data files \n"
                           "(e.g. the name of the .in and .setup files): "))

        loc = str(input("Enter the path where the PHANTOM output of the models is saved: "))
        while True:
            if not os.path.isdir(loc):
                loc = str(input("Path does not exist, please try again: "))
                continue
            if " " in loc:
                loc = str(input("Path contains spaces, please remove them: "))
                continue
            break

        outputloc = str(input("Enter the path where the pipeline output should be saved: "))
        while True:
            if " " in outputloc:
                outputloc = str(input("Path contains spaces, please remove them: "))
                continue
            break

        file.write("prefix = " + prefix + "\n")
        file.write("data_location = " + loc + "\n")
        file.write("pipeline_output_location = " + outputloc + "\n")
        print("Settings saved at " + userSettingsFilePath)
        print("--------------------------------------------------------------")
        file.close()

# Load userSettingsFile.txt if it exists
def load(userSettingsFilePath):
    dictionary = {}
    splittedFile = []
    print()
    print("--------------------------------------------------------------")
    print("The following user settings will be loaded from %s" % userSettingsFilePath)
    print()
    with open(userSettingsFilePath, "r") as file:
        for line in file.readlines():
            splittedLine = line.split()
            dictionary[splittedLine[0]] = splittedLine[2]
            print(line[:-1])
    print("--------------------------------------------------------------")
    print()

    loc = str(dictionary["data_location"])
    outputloc = str(dictionary["pipeline_output_location"])

    # If the user has changed the directory in userSettingsFile.txt, check if the path still exists
    if not os.path.isdir(loc):
        print(
            "ERROR: The specified directory 'data_location' in %s does not exist anymore, shutting down the pipeline."
            % userSettingsFilePath)
        print('------------------END:', dt.datetime.now(), '---------------------')
        sys.exit()

    if " " in loc:
        print("ERROR: The specified directory 'data_location' in %s contains spaces, shutting down the pipeline."
              % userSettingsFilePath)
        print('------------------END:', dt.datetime.now(), '---------------------')
        sys.exit()

    if not os.path.isdir(outputloc):
        print(
            "ERROR: The specified directory 'pipeline_output_location' in %s does not exist anymore, shutting down the"
            " pipeline." % userSettingsFilePath)
        print('------------------END:', dt.datetime.now(), '---------------------')
        sys.exit()

    if " " in outputloc:
        print(
            "ERROR: The specified directory 'pipeline_output_location' in %s contains spaces, shutting down the pipeline."
            % userSettingsFilePath)
        print('------------------END:', dt.datetime.now(), '---------------------')
        sys.exit()

    return dictionary