import os
import sys

# check if filename argument is provided
if len(sys.argv) < 2:
    print("Please provide a filename argument without the .cpp extension")
    sys.exit(1)
elif len(sys.argv) == 2:
    filename = sys.argv[1]
    option = ''
elif len(sys.argv) == 3:
    filename = sys.argv[1]
    option = sys.argv[2]

# set compiler flags based on OS
if os.name == "nt":
    # Windows
    compiler_flags = ("-IC:/Users/berko/AppData/Local/Programs/Python/Python311/include "
                      "-LC:/Users/berko/AppData/Local/Programs/Python/Python311/libs "
                      "-lpython311 -fopenmp")
else:
    # Linux or macOS
    compiler_flags = "-I/usr/include/python3.10 -L/usr/lib/x86_64-linux-gnu -lpython3.10 -fopenmp"

# compile the program
os.system(f"g++ {filename}.cpp -o builds/{filename}.exe {compiler_flags}" if os.name == "nt"
          else f"g++ {filename}.cpp -o builds/{filename} {compiler_flags}")
# check if compilation was successful
if os.path.exists(f"builds/{filename}.exe" if os.name == "nt" else os.path.join(os.getcwd(), "builds", filename)):
    # run the program
    os.system(f".\\builds\\{filename}.exe" if os.name == "nt" else f"./builds/{filename}")
    # run the Python script for Lyapunovs and Integrate
    if "Lyapunovs" in filename or "Integrate" in filename or "Energy" in filename:
        os.system("python plottingResults/csvFiles/plot.py {}".format(option) if os.name == "nt" else "python3 plottingResults/csvFiles/plot.py {}".format(option))
else:
    print("Compilation failed")
    sys.exit(1)
