//
// Created by berko on 3/20/2023.
//

#ifndef LYAPUNOVEXPONENTSV_0_3_LYAPUNOVANALYSIS_H
#define LYAPUNOVEXPONENTSV_0_3_LYAPUNOVANALYSIS_H

#include <iostream>
#include <vector>
#include <cmath>
#include <chrono>
#include "MatrixOperations.h"
#include "FunctionalAnalysis.h"

using namespace std;

// This function takes a system as an input (system should be a vector valued function) and a point in the phase space
// and returns the combined system with the Jacobian matrix of the system at that point.
vector<long double> JacobianSystem(vector<long double>(*f)(vector<long double>), const vector<long double>& point) {
    vector<long double> result;

    vector<long double> f1 = f(point);
    int dimension = f1.size();

    vector<long double> initialPhi;
    initialPhi.resize(dimension * dimension);
    for (int i = 0; i < dimension * dimension; i++) {
        initialPhi[i] = point[i + dimension];
    }

    vector<long double> jacobian = Jacobian(f, point);
    vector<long double> phi = MatrixProduct(initialPhi, Transpose(jacobian));
    result = ExtendVector(f1, phi);
    return result;
}
//
vector<long double> RungeKuttaJacobian(vector<long double>(*f)(vector<long double>), const vector<long double>& point, long double dt) {
    vector<long double> k1 = JacobianSystem(f, point);
    vector<long double> k2 = JacobianSystem(f, VectorAddition(point, ScalarVectorMultiplication(dt / 2, k1)));
    vector<long double> k3 = JacobianSystem(f, VectorAddition(point, ScalarVectorMultiplication(dt / 2, k2)));
    vector<long double> k4 = JacobianSystem(f, VectorAddition(point, ScalarVectorMultiplication(dt, k3)));
    vector<long double> result = VectorAddition(
            point, ScalarVectorMultiplication( dt/6,
                                               VectorAddition(
                                                       k1, VectorAddition(
                                                               ScalarVectorMultiplication(2, k2), VectorAddition(
                                                                       ScalarVectorMultiplication(2, k3), k4)))));
    return result;
}


// This function integrates the Jacobian system for a given number of steps and returns the final point.
vector<long double> IntegrateJacobianSystem(vector<long double>(*f)(vector<long double>), const vector<long double>& point, long double dt, int numSteps) {
    vector<long double> currentPoint = point;
    for (int i = 0; i < numSteps; i++) {
        currentPoint = RungeKuttaJacobian(f, currentPoint, dt);
    }
    return currentPoint;
}

// This function takes a point in the phase space and returns the phi matrix.
vector<long double> GetPhi(const vector<long double>& point, int dimension) {
    vector<long double> phi;
    phi.resize(dimension * dimension);
    for (int i = 0; i < dimension * dimension; i++) {
        phi[i] = point[i + dimension];
    }
    return phi;
}

// This function takes in a vector and returns the natural logarithm of the norms of the vectors.
vector<long double> LogNorms(const vector<long double>& vector1) {
    vector<long double> result;
    int dimension = sqrt(vector1.size());
    for (int i = 0; i < dimension; i++) {
        vector<long double> currentVector;
        currentVector.resize(dimension);
        for (int j = 0; j < dimension; j++) {
            currentVector[j] = vector1[i * dimension + j];
        }
        result.push_back(log(Norm(currentVector)));
    }
    return result;
}

// This function takes a vector and normalizes it.
vector<long double> Normalize(const vector<long double>& vector1) {
    vector<long double> result;
    int dimension = sqrt(vector1.size());
    for (int i = 0; i < dimension; i++) {
        vector<long double> currentVector;
        currentVector.resize(dimension);
        for (int j = 0; j < dimension; j++) {
            currentVector[j] = vector1[i * dimension + j];
        }
        currentVector = ScalarVectorMultiplication(1 / Norm(currentVector), currentVector);
        for (int j = 0; j < dimension; j++) {
            result.push_back(currentVector[j]);
        }
    }
    return result;
}

// This function takes a point in the phase space and a vector and replaces the phi matrix that is given as vector with the given vector.
vector<long double> ReplacePhi(const vector<long double>& point, const vector<long double>& newPhi) {
    vector<long double> result;
    int dimension = sqrt(newPhi.size());
    for (int i = 0; i < dimension; i++) {
        result.push_back(point[i]);
    }
    for (int i = 0; i < dimension * dimension; i++) {
        result.push_back(newPhi[i]);
    }
    return result;
}

// This function takes a vector of vectors and sums the vectors.
vector<vector<long double>> SumLogNorms(const vector<vector<long double>>& logNorms, long double dt, int numSteps) {
    vector<vector<long double>> result;
    int dimension = logNorms[0].size();
    result.resize(logNorms.size());
    for (int i = 0; i < logNorms.size(); i++) {
        result[i].resize(dimension);
    }

    for (int i = 0; i < dimension; i++) {
        long double lyapunov = 0.0;
        for (int j = 0; j < logNorms.size(); j++) {
            lyapunov += logNorms[j][i];
            result[j][i] = lyapunov / ((j + 1) * numSteps * dt);
        }
    }
    return result;
}

vector<vector<long double>> LyapunovSpectrum(vector<long double>(*f)(vector<long double>), const vector<long double>& point, long double dt, int numSteps, int numIterations) {
    vector<vector<long double>> result;
    vector<vector<long double>> logNorms;
    int dimension = f(point).size();

    vector<long double> currentPoint = point;
    for (int i = 0; i < numIterations; i++) {
        cout << "\033[A\33[2K\rIteration number: " << i+1 << " / " << numIterations << "\n";
        currentPoint = IntegrateJacobianSystem(f, currentPoint, dt, numSteps);
        vector<long double> phi = GetPhi(currentPoint, dimension);
        phi = Orthogonalize(phi);

        logNorms.push_back(LogNorms(phi));

        phi = Normalize(phi);

        currentPoint = ReplacePhi(currentPoint, phi);
    }

    result = SumLogNorms(logNorms, dt, numSteps);
    return result;
}

#endif //lyapunovexponentsv_0_3_lyapunovanalysis_h
