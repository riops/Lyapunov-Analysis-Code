import matplotlib.pyplot as plt
import csv
import os
import glob
import sys

if len(sys.argv) == 2:
    plotted_var = sys.argv[1]
else:
    plotted_var = ''

def transpose(matrix):
    return [[matrix[j][i] for j in range(len(matrix))] for i in range(len(matrix[0]))]

# We will have a list of lists, and each list inside the list will be the spectrum at some time.
csv_files = glob.glob(os.path.join(os.getcwd(), 'plottingResults', 'csvFiles', '*.csv'))
# Sorting the files by last editing time then getting the latest editted csv file.
csv_files.sort(key=os.path.getmtime)
filename = csv_files[-1]

if not 'Energy' in filename:
    with open(filename, 'r') as f:
        reader = csv.reader(f)
        spectrum = transpose([list(map(float, row)) for row in reader])

    cpp_file = filename.split('_')[0].split('\\')[-1] + '.cpp' if os.name == "nt" else filename.split('_')[0].split('/')[-1] + '.cpp'
else:
    Energies = []
    Lyapunovs1 = []
    errors1 = []
    Lyapunovs2 = []
    errors2 = []
    with open(filename, 'r') as f:
        reader = csv.reader(f)
        for i, row in enumerate(reader):
            for element in row:
                if i == 0:
                    Energies.append(float(element))
                elif i == 1:
                    Lyapunovs1.append(float(element))
                elif i == 2:
                    errors1.append(float(element))
                elif i == 3:
                    Lyapunovs2.append(float(element))
                elif i == 4:
                    errors2.append(float(element))

time_list = []
if 'EOM' in filename:
    cpp_file = cpp_file.replace('EOM', 'Integrate')
    y_label = 'Equation of Motion'
    title = 'Equations of Motion with '
    with open(cpp_file, 'r') as generator:
        # We read the file that generated the csv file to get the parameters
        for line in generator:
            if 'double dt' in line:
                dt = float(line.split('=')[-1].split(';')[0])
            elif 'int numIterations' in line:
                numIterations = int(line.split('=')[-1].split(';')[0])
    for i in range(0, numIterations):
        time_list.append((i+1) * dt)
    # We set the title with the necessary information
    title = title + f'dt={dt}, numIterations={numIterations}'
elif 'Lyapunovs' in filename:
    y_label = 'Lyapunov Exponent'
    title = 'Lyapunov Spectrum with '
    with open(cpp_file, 'r') as generator:
        # We read the file that generated the csv file to get the parameters
        for line in generator:
            if 'double dt' in line:
                dt = float(line.split('=')[-1].split(';')[0])
            elif 'int numIterations' in line:
                numIterations = int(line.split('=')[-1].split(';')[0])
            elif 'int numSteps' in line:
                numSteps = int(line.split('=')[-1].split(';')[0])
    for i in range(0, numIterations):
        time_list.append((i+1) * dt * numSteps)
    # We set the title with the necessary information
    title = title + f' dt={dt}, numIterations={numIterations}, numSteps={numSteps}'
elif 'Energy' in filename:
    plt.style.use('seaborn')
    #plt.errorbar(Energies, Lyapunovs_exact, yerr=errors_exact, fmt='.', markersize=4, ecolor='red', markerfacecolor='blue')
    plt.errorbar(Energies, Lyapunovs1, yerr=errors1, fmt='.', markersize=4, ecolor='black', markerfacecolor='grey')
    plt.errorbar(Energies, Lyapunovs2, yerr=errors2, fmt='.', markersize=4, ecolor='red', markerfacecolor='blue')
    title = 'Mean Largest Lyapunov Exponents vs Energy Values'
    plt.xlabel('Energy Values ($E$)')
    plt.ylabel('Mean Largest Lyapunov Exponents ($\lambda$)')

if plotted_var != '':
    if '-' in plotted_var:
        plotted_vars = [int(plotted_var.split('-')[0])-1, int(plotted_var.split('-')[1])-1]
        plt.plot(spectrum[plotted_vars[0]], spectrum[plotted_vars[1]])
        plt.ylabel('var = ' + str(plotted_vars[1]))
        plt.xlabel('var = ' + str(plotted_vars[0]))
    elif ',' in plotted_var:
        plotted_vars = [int(var) for var in plotted_var.split(',')]
        for var in plotted_vars:
            plt.plot(time_list, spectrum[var - 1])
        plt.ylabel(y_label + ' vars = ' + str(plotted_vars))
        plt.xlabel('Time')
    else:
        plt.plot(time_list, spectrum[int(plotted_var) - 1])
        plt.ylabel(y_label + ' var = ' + plotted_var)
        plt.xlabel('Time')
elif 'Energy' not in filename:
    for lyap in spectrum:
        plt.plot(time_list, lyap)
    plt.xlabel('Time')
    plt.ylabel(y_label)

plt.title(title)

figure_name = filename.replace('csvFiles', 'plots').replace('.csv', '.png')
plt.savefig(figure_name)

plt.show()
