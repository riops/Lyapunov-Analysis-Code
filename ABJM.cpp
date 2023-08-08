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
    int dimension = 4;
    int combinedDimension = dimension * dimension + dimension;
    int numIterations = 15000;
    long double dt = 0.001;
    int numSteps = 10;
    int numThreads = 12;
    int numPoints = 40;

    vector<long double> energies;
    energies = {0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100};

    for (long double energy : energies) {
        const char *fileName = "systems";
        const char *functionName = "random_energy_initial";

        vector<vector<long double>> points;
        points.resize(numPoints);
        for (int i = 0; i < numPoints; i++) {
            points[i].resize(combinedDimension);
        }

        Py_Initialize();
        PyRun_SimpleString("import sys; sys.path.append('./python/')");

        PyObject * pName, *pModule, *pFunc, *pArgs, *pResult;

        pName = PyUnicode_FromString(fileName);
        pModule = PyImport_Import(pName);

        pFunc = PyObject_GetAttrString(pModule, functionName);
        pArgs = PyTuple_Pack(1, PyFloat_FromDouble(energy));

        cout << "for energy = " << energy << endl << endl;
        for (int i = 0; i < numPoints; i++) {
            pResult = PyObject_CallObject(pFunc, pArgs);
            for (int j = 0; j < combinedDimension; j++) {
                points[i][j] = PyFloat_AsDouble(PyList_GetItem(pResult, j));
            }
        }
        cout << "All points have been generated for E=" << energy << endl << endl;

        vector<vector<long double>> result;
        result.resize(numPoints);
        for (int i = 0; i < numPoints; i++) {
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
            for (int i = 0; i < numPoints; i++) {
                spectrum[i] = LyapunovSpectrum(K, points[i], dt, numSteps, numIterations);
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
            for (int j = 1; j < dimension; j++) {
                if (result[i][j] > maxSpectrum[i]) {
                    maxSpectrum[i] = result[i][j];
                }
            }
        }
        long double average = 0;
        for (int i = 0; i < numPoints; i++) {
            average += maxSpectrum[i];
            //cout << "Maximum for point " << i + 1 << " is: " << maxSpectrum[i] << endl;
        }
        average /= numPoints;
        cout << "Average maximum for energy " << energy << " is: " << average << endl << endl;
    }
    system("pause");
    return 0;

}
#pragma clang diagnostic pop
