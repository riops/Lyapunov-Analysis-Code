#include <iostream>
#include <vector>
#include <fstream>
#include <ctime>
#include <iomanip>
#include <sstream>
#include "cModules/MatrixOperations.h"
#include "cModules/Systems.h"
#include "cModules/LyapunovAnalysis.h"

using namespace std;

int main() {
    int dimension = 3;
    int numIterations = 10000;
    long double dt = 0.001;
    int numSteps = 10;

    vector<long double> point = { 19, 20, 50, 1, 0, 0, 0, 1, 0, 0, 0, 1 };

    auto start = chrono::high_resolution_clock::now();
    cout << endl << "Starting calculation..." << endl << endl;

    vector<vector<long double>> spectrum = LyapunovSpectrum(Lorenz, point, dt, numSteps, numIterations);
    
    auto end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::seconds>(end - start);
    cout << endl << "Time taken by function: " << duration.count() << " seconds" << endl << endl;


    cout << endl << "Done!" << endl << endl << "Writing to file..." << endl << endl;

    // Get the current date and time
    auto now = std::chrono::system_clock::now();
    std::time_t t = std::chrono::system_clock::to_time_t(now);
    std::tm tm = *std::localtime(&t);

    // Format the date and time as a string
    std::ostringstream oss;
    oss << std::put_time(&tm, "%Y-%m-%d_%H-%M-%S");
    std::string datetime = oss.str();

    // Use the formatted date and time as the file name
    std::string filename = "LyapunovsLorenz_" + datetime + ".csv";

    // Write the data to the csv file
    ofstream fs("plottingResults/csvFiles/" + filename);
    for(int j = 0; j < numIterations; j++){
        for(int k = 0; k < dimension; k++){
            if (k == dimension - 1) fs << spectrum[j][k];
            else fs << spectrum[j][k] << ",";
        }
        fs << endl;
    }

    long double maxLyapunov = 0;
    for (int i = 0; i < dimension; i++){
        if (spectrum[numIterations-1][i] > maxLyapunov) maxLyapunov = spectrum[numIterations-1][i];
    }

    cout << "Maximum Lyapunov exponent is: " << maxLyapunov << endl << endl;
    
    long double summed = 0;
    for (int i = 0; i < dimension; i++){
        summed += spectrum[numIterations-1][i];
    }

    cout << "Summation of the Lyapunovs is: " << summed << endl << endl;

    return 0;
}
