#pragma clang diagnostic push
#pragma ide diagnostic ignored "openmp-use-default-none"
#include <iostream>
#include <vector>
#include <fstream>
#include <ctime>
#include <iomanip>
#include <chrono>
#include <sstream>

#include "cModules/MatrixOperations.h"
#include "cModules/Systems.h"
#include "cModules/LyapunovAnalysis.h"
#include "cModules/Analysis.h"

#include <omp.h>
#include <Python.h>


using namespace std;

int main() {
    auto now = std::chrono::system_clock::now();
    std::time_t t = std::chrono::system_clock::to_time_t(now);
    std::tm tm = *std::localtime(&t);

    std::ostringstream oss;
    oss << std::put_time(&tm, "%Y-%m-%d_%H-%M-%S");
    std::string datetime = oss.str();

    ofstream toBePlotted("plottingResults/csvFiles/LyapVsEnergy_" + datetime + ".csv");
    ofstream logs("plottingResults/logs/LyapVsEnergy_" + datetime + ".md");

    logs << "# Lyapunov Exponents vs Energy for " << datetime << endl << endl;

    int dimension = 2;
    int eomDimension = 8*dimension*dimension;
    int combinedDimension = eomDimension * eomDimension + eomDimension;
    int numIterations = 6000;
    long double dt = 0.1;
    int numSteps = 1;
    int numThreads = 10;
    int numPoints = 10;

    vector<long double> energies;
    energies = {10, 20, 30, 40, 50, 60, 70, 80, 90, 100};

    vector<long double> averageMaxLyapunovs1;
    averageMaxLyapunovs1.resize(energies.size());
    vector<long double> averageMaxLyapunovs2;
    averageMaxLyapunovs2.resize(energies.size());
    vector<long double> stdErrors1;
    stdErrors1.resize(energies.size());
    vector<long double> stdErrors2;
    stdErrors2.resize(energies.size());

    vector<vector<vector<long double>>> coefficients;

    coefficients.resize(energies.size());
    for (int i = 0; i < energies.size(); i++) {
        coefficients[i].resize(numPoints);
        for (int j = 0; j < numPoints; j++) {
            coefficients[i][j].resize(2);
        }
    }

    int energyCounter = 0;

    for (long double energy : energies) {
        const char *fileName = "systems";
        const char *functionName = "initial_condition_new_model_alternative";

        Py_Initialize();
        PyRun_SimpleString("import sys; sys.path.append('./python/')");

        PyObject * pName, *pModule, *pFunc, *pArgs, *pResult;

        pName = PyUnicode_FromString(fileName);
        pModule = PyImport_Import(pName);

        pFunc = PyObject_GetAttrString(pModule, functionName);
        pArgs = PyTuple_Pack(3, PyFloat_FromDouble(energy), PyFloat_FromDouble(1.0), PyFloat_FromDouble(10.0), PyFloat_FromDouble(dimension));

        cout << "For energy = " << energy << endl;
        vector<vector<long double>> points;
        points.resize(numPoints);

        for (int i = 0; i < numPoints; i++) {
            pResult = PyObject_CallObject(pFunc, pArgs);
            long double alpha = PyFloat_AsDouble(PyList_GetItem(pResult, 0));
            long double beta = PyFloat_AsDouble(PyList_GetItem(pResult, 1));
            coefficients[energyCounter][i][0] = alpha;
            coefficients[energyCounter][i][1] = beta;
            vector<long double> point;
            point = ExtendVector(ExtendVector(ExtendVector(ExtendVector(ScalarVectorMultiplication(alpha, RealPartSU2(dimension)[0]), ScalarVectorMultiplication(alpha, ComplexPartSU2(dimension)[0])), 
                                 ExtendVector(ScalarVectorMultiplication(beta, RealPartSU2(dimension)[1]), ScalarVectorMultiplication(beta, ComplexPartSU2(dimension)[1]))), ZeroMatrix(4)), IdentityMatrix(32));
            points[i] = point;
        }
        cout << "All points have been generated" << endl << endl;

        vector<vector<long double>> result;
        result.resize(numPoints);
        for (int i = 0; i < numPoints; i++) {
            result[i].resize(eomDimension);
        }
        vector<vector<vector<long double>>> spectrum;
        spectrum.resize(numPoints);

        auto start = chrono::high_resolution_clock::now();
        omp_set_num_threads(numThreads);
        cout << "Number of threads: " << numThreads << endl << endl;
        cout << "Starting calculation..." << endl << endl;

#pragma omp parallel
        {
#pragma omp for
            for (int i = 0; i < numPoints; i++) {
                spectrum[i] = LyapunovSpectrum(G, points[i], dt, numSteps, numIterations);
                result[i] = spectrum[i][numIterations - 1];
            }
        };
        auto end = chrono::high_resolution_clock::now();
        auto duration = chrono::duration_cast<chrono::seconds>(end - start);
        cout << endl << "Time taken by function: " << duration.count() << " seconds" << endl << endl;

        vector<long double> maxSpectrum;
        //Finds the maximum value of the vector for each point.
        maxSpectrum.resize(numPoints);
        for (int i = 0; i < numPoints; i++) {
            maxSpectrum[i] = result[i][0];
            for (int j = 1; j < eomDimension; j++) {
                if (result[i][j] > maxSpectrum[i]) {
                    maxSpectrum[i] = result[i][j];
                }
            }
        }

        vector<long double> maxSpectrumSorted  = Sort(maxSpectrum);

        logs << "For Energy " << energy << " the maximum Lyapunov exponents are:<br />" << endl;
        for (int i = 0; i < numPoints; i++) {
            logs << "\tFor point " << i+1 << "<br />\n\t\t$\\alpha$ = " << coefficients[energyCounter][i][0] << "<br />" << endl;
            logs << "\t\t$\\beta$ = " << coefficients[energyCounter][i][1] << "<br />" << endl;
            logs << "\t\tMaximum Lyapunov for this point is " << maxSpectrum[i] << "<br />" << endl;
        }
        
        auto groupingResults = SeparateAsMinStd(maxSpectrumSorted);

        long double average1 = groupingResults[0][0];
        long double stdError1 = groupingResults[0][1];
        
        long double average2 = groupingResults[0][2];
        long double stdError2 = groupingResults[0][3];

        vector<long double> group1 = groupingResults[1];
        vector<long double> group2 = groupingResults[2];

        averageMaxLyapunovs1[energyCounter] = average1;
        stdErrors1[energyCounter] = stdError1;

        averageMaxLyapunovs2[energyCounter] = average2;
        stdErrors2[energyCounter] = stdError2;

        cout << "Group 1 is:" << endl;

        PrintVector(group1);

        cout << endl;

        cout << "Average maximum for group 1 at energy = " << energy << " is: " << average1 << endl;
        cout << "Standard deviation error for group 1 at energy = " << energy << " is: " << stdError1 << endl << endl;

        cout << "Group 2 is:" << endl;

        PrintVector(group2);

        cout << endl;

        cout << "Average maximum for group 2 at energy = " << energy << " is: " << average2 << endl;
        cout << "Standard deviation error for group 2 at energy = " << energy << " is: " << stdError2 << endl << endl; 

        energyCounter++;
    }
    cout << endl << "Done!" << endl << endl << "Writing to file..." << endl << endl;

    for (int i = 0; i < energies.size(); i++){
        if (i == energies.size() - 1){
            toBePlotted << energies[i] << endl;
        }
        else toBePlotted << energies[i] << ",";
    }
    
    for (int i = 0; i < energies.size(); i++){
        if (i == energies.size() - 1){
            toBePlotted << averageMaxLyapunovs1[i] << endl;
        }
        else toBePlotted << averageMaxLyapunovs1[i] << ",";
    }

    for (int i = 0; i < energies.size(); i++){
        if (i == energies.size() - 1){
            toBePlotted << stdErrors1[i] << endl;
        }
        else toBePlotted << stdErrors1[i] << ",";
    }

    for (int i = 0; i < energies.size(); i++){
        if (i == energies.size() - 1){
            toBePlotted << averageMaxLyapunovs2[i] << endl;
        }
        else toBePlotted << averageMaxLyapunovs2[i] << ",";
    }

    for (int i = 0; i < energies.size(); i++){
        if (i == energies.size() - 1){
            toBePlotted << stdErrors2[i] << endl;
        }
        else toBePlotted << stdErrors2[i] << ",";
    }
    
    return 0;

}
#pragma clang diagnostic pop
