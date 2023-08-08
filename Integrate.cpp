#include <iostream>
#include <vector>
#include <fstream>
#include <chrono>
#include <ctime>
#include <iomanip>
#include <sstream>
#include "cModules/MatrixOperations.h"
#include "cModules/Systems.h"
#include "cModules/LyapunovAnalysis.h"

using namespace std;

int main() {
    int dimension = 32;
    int numIterations = 10000;
    long double dt = 0.01;
    vector<long double> initialCondition = {0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

    vector<vector<long double>> spectrum = IntegrateSystem(G, initialCondition, dt, numIterations, true);
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
    std::string filename = "EOM_" + datetime + ".csv";

    // Write the data to the csv file
    ofstream fs("plottingResults/csvFiles/" + filename);
    for(int j = 0; j < numIterations; j++){
        for(int k = 0; k < dimension; k++){
            if (k == dimension - 1) fs << spectrum[j][k];
            else fs << spectrum[j][k] << ",";
        }
        fs << endl;
    }
    return 0;
}
