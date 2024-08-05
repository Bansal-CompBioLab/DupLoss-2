# This is a simple Python script for executing multiple runs of DupLoss-2 on the same input dataset using different random seeds.
# The script works for both Linux and macOS and assumes that either DupLoss-2.linux or DupLoss-2.mac, depending on the operating system, is present in the same directory as this script. This can be easily changed below, as needed.   


import subprocess
import sys
import random
import time
import platform


# Get command line arguments
if len(sys.argv) != 5 and len(sys.argv)!=4:
    print("Usage: python MultiRunScript.py <num_executions> <Input file name> <Output file name prefix> <optional random seed>")
    sys.exit(1)

try:
    num_executions = int(sys.argv[1])
except ValueError:
    print("The number of executions must be an integer.")
    sys.exit(1)

# Get the file names from the command line arguments
file1 = sys.argv[2]
file2 = sys.argv[3]



# Get random seed from arguments, if present
if len(sys.argv) == 5:
    try:
        random_seed = int(sys.argv[4])
    except ValueError:
        print("The random seed must be an integer.")
        sys.exit(1)
else:
    # Generate a random number based on system time
    random_seed = int(time.time() * 1000)  # Milliseconds since epoch
    
random.seed(random_seed)
print("Using random seed ", random_seed)   



# Detect operating system

system = platform.system()
if system == "Linux":
    OS = 1
elif system == "Darwin":
    OS = 2
elif system == "Windows":
    OS =3
    print("This script currently only works for Linux and macOS operating systems")
    sys.exit(1)
else:
    print("Unknown OS; Proceeding assuming Linux...\n")  # Unknown OS
    OS = 1




# Execute DupLoss-2 for the specified number of times
for i in range(num_executions):

    random_number = random.randint(1, 10000)  # geerate random number to be used as seed for DupTree-2


    # Constructing command to be executed, depending on operating system
    if OS ==1:
        command = "./DupLoss-2.linux --quiet -i " + file1 + " --seed " + str(random_number) + " -o " + file2 + "_" + str(i+1)
    elif OS == 2:
        command = "./DupLoss-2.mac --quiet -i " + file1 + " --seed " + str(random_number) + " -o " + file2 + "_" + str(i+1)
    
    print("\nExecuting Run #"+ str(i+1) + ": " + command) 

    try:
        subprocess.run([command], shell = True, check=True)
        print("Execution Successful.")
        
        file = file2 +"_" + str(i+1)   # open output file to pull out reconciliation cost

        try:
            with open(file, "r") as file:
                lines = file.readlines()
                print(lines[2].strip())
                file.close()
        except FileNotFoundError:
            print("File " + file + " not found.")



    except subprocess.CalledProcessError as e:
        print("Execution Failed.")
        sys.exit(1)

print("\n" + str(num_executions) + " runs of DupLoss-2 completed.")
