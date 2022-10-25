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

        loc = str(input("Enter the path where the PHANTOM outputs of the models are saved: "))
        while True:
            if not os.path.isdir(loc):
                loc = str(input("Path does not exist, please try again: "))
                continue
            if " " in loc:
                loc = str(input("Path contains spaces, please remove them: "))
                continue
            break

        if (loc[-1] == '/' or loc[-1] == '\\'):
            loc = loc[:-1]

        outputloc = str(input("Enter the path where the pipeline output should be saved: "))
        if not os.path.isdir(outputloc):
            os.makedirs(outputloc)
        while True:
            if " " in outputloc:
                outputloc = str(input("Path contains spaces, please remove them: "))
                continue
            break

        if (outputloc[-1] == '/' or outputloc[-1] == '\\'):
            outputloc = outputloc[:-1]

        phantom = str(input("Enter the path where PHANTOM is installed: "))
        while True:
            if not os.path.isdir(phantom):
                phantom = str(input("Path does not exist, please try again: "))
                continue
            if " " in phantom:
                phantom = str(input("Path contains spaces, please remove them: "))
                continue
            break

        if (phantom[-1] == '/' or phantom[-1] == '\\'):
            phantom = phantom[:-1]

        file.write("prefix = " + prefix + "\n")
        file.write("data_location = " + loc + "\n")
        file.write("pipeline_output_location = " + outputloc + "\n")
        file.write("hard_path_to_phantom = " + phantom + "\n")
        print("Settings saved at " + userSettingsFilePath)
        print("--------------------------------------------------------------")
        file.close()

# Load userSettingsFile.txt if it exists
def load(userSettingsFilePath,onlyPathToPhantom=False):
    dictionary = {}
    splittedFile = []
    print()
    print("--------------------------------------------------------------")
    print("The following user settings will be loaded from %s" % userSettingsFilePath)
    print()
    with open(userSettingsFilePath, "r") as file:
        for line in file.readlines():
            splittedLine = line.split()
            if len(splittedLine) == 3:
                dictionary[splittedLine[0]] = splittedLine[2]
            else:
                dictionary[splittedLine[0]] = splittedLine[2:]
            print(line[:-1])
    print("--------------------------------------------------------------")
    print()
    
    if onlyPathToPhantom == False:
      loc = str(dictionary["data_location"])
      outputloc = str(dictionary["pipeline_output_location"])
    phantom = str(dictionary["hard_path_to_phantom"])

    # If the user has changed the directory in userSettingsFile.txt, check if the path still exists
    if onlyPathToPhantom == False:
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

    if not os.path.isdir(phantom):
        print(
            "ERROR: The specified directory 'hard_path_to_phantom' in %s does not exist anymore, shutting down the pipeline."
            % userSettingsFilePath)
        print('------------------END:', dt.datetime.now(), '---------------------')
        sys.exit()

    if " " in phantom:
        print("ERROR: The specified directory 'hard_path_to_phantom' in %s contains spaces, shutting down the pipeline."
              % userSettingsFilePath)
        print('------------------END:', dt.datetime.now(), '---------------------')
        sys.exit()

    return dictionary
