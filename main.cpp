#pragma clang diagnostic push
#pragma ide diagnostic ignored "openmp-use-default-none"
#include <iostream>
#include <vector>
#include <fstream>
#include "cModules/MatrixOperations.h"
#include "cModules/Systems.h"
#include "cModules/LyapunovAnalysis.h"
#include <omp.h>
#include <Python.h>

using namespace std;

int main() {
    int dimension = 32;
    int combinedDimension = dimension * dimension + dimension;
    int numIterations = 1500;
    long double dt = 0.001;
    int numSteps = 10;
    int numThreads = 12;
    int numPoints = 10;
    long double energy = 12.0;
    const char* fileName = "systems";
    const char* functionName = "initial_condition_new_model";

    vector<vector<long double>> points;
    points.resize(numPoints);
    for (int i = 0; i < numPoints; i++){
        points[i].resize(combinedDimension);
    }

    Py_Initialize();

    PyRun_SimpleString("import sys; sys.path.append('./python')");
  
    PyObject *pName, *pModule, *pFunc, *pArgs, *pResult;

    pName = PyUnicode_FromString(fileName);
    pModule = PyImport_Import(pName);

    pFunc = PyObject_GetAttrString(pModule, functionName);
    pArgs = PyTuple_Pack(3, PyFloat_FromDouble(energy), PyFloat_FromDouble(2.0), PyFloat_FromDouble(1.0), PyFloat_FromDouble(1.0));

    for (int i = 0; i < numPoints; i++){
        pResult = PyObject_CallObject(pFunc, pArgs);
        for (int j = 0; j < combinedDimension; j++){
            points[i][j] = PyFloat_AsDouble(PyList_GetItem(pResult, j));
        }
        cout << "Point " << i+1 << " generated" << endl;
    }
    cout << endl;
    Py_Finalize();

    vector<vector<long double>> result;
    result.resize(numPoints);
    for (int i = 0; i < numPoints; i++){
        result[i].resize(dimension);
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
        for (int i = 0; i < numPoints; i++){
            spectrum[i] = LyapunovSpectrum(K, points[i], dt, numSteps, numIterations);
            result[i] = spectrum[i][numIterations - 1];
        }
    };
    auto end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::seconds>(end - start);
    cout << endl << "Time taken by function: " << duration.count() << " seconds" << endl << endl;
    /*
    for (int i = 0; i < numPoints; i++){
        for (int j = 0; j < dimension; j++){
            cout << result[i][j] << " ";
        }
        cout << endl;
    }*/


    /*
    cout << endl << "Done!" << endl << endl << "Writing to file..." << endl;

    ofstream fs("Spectrum.txt");
    for(int i = 0; i < numPoints; i++){
        fs << "["; //Start of one iteration
        for(int j = 0; j < numIterations; j++){
            fs << "["; //Start of one point
            for(int k = 0; k < dimension; k++){
                if (k == dimension - 1){
                    fs << spectrum[i][j][k] << "]" << endl;
                }
                else{
                    fs << spectrum[i][j][k] << ", ";
                }
            }
            if (j == numIterations - 1){
                fs << "]" << endl;
            }
            else{
                fs << ", ";
            }
        }
        fs << endl << endl;
    }*/

    /*
    vector<long double> summedSpectrum;
    summedSpectrum.resize(numPoints);
    for (int i = 0; i < numPoints; i++){
        summedSpectrum[i] = 0;
        for (int j = 0; j < dimension; j++){
            summedSpectrum[i] += result[i][j];
        }
    }

    for (int i = 0; i < numPoints; i++){
        cout << "Point " << i+1 << " has the summation: " << summedSpectrum[i] << endl;
    }
    cout << endl;*/

    vector<long double> maxSpectrum;
    //Finds the maximum value of the vector for each point.
    maxSpectrum.resize(numPoints);
    for (int i = 0; i < numPoints; i++){
        maxSpectrum[i] = result[i][0];
        for (int j = 1; j < dimension; j++){
            if (result[i][j] > maxSpectrum[i]){
                maxSpectrum[i] = result[i][j];
            }
        }
    }

    for (int i = 0; i < numPoints; ++i) {
        cout << "Point " << i+1 << " has the maximum: " << maxSpectrum[i] << endl;
    }
    cout << endl;

    long double average = 0;
    for (int i = 0; i < numPoints; i++){
        average += maxSpectrum[i];
    }
    average /= numPoints;
    cout << "Average maximum for energy " << energy << "is: " << average << endl << endl;

    system("pause");
    return 0;

}
#pragma clang diagnostic pop
