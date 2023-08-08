//
// Created by berko on 3/20/2023.
//

#ifndef LYAPUNOVEXPONENTSV_0_3_FUNCTIONALANALYSIS_H
#define LYAPUNOVEXPONENTSV_0_3_FUNCTIONALANALYSIS_H

#include <iostream>
#include <vector>
#include <cmath>
#include <chrono>
#include "MatrixOperations.h"

using namespace std;

// This function takes in a vector valued function and a point in the phase space,
// and returns the jacobian matrix of the function at that point as a vector.
/*
vector<long double> Jacobian(vector<long double>(*f)(vector<long double>), const vector<long double>& point) {
    vector<long double> result;
    long double epsilon = 0.001;
    int dimension = f(point).size();
    vector<long double> point1 = point;
    vector<long double> point2 = point;

    for (int i = 0; i < dimension; i++) {
        point1[i] += epsilon;
        point2[i] -= epsilon;
        vector<long double> f1 = f(point1);
        vector<long double> f2 = f(point2);
        for (int j = 0; j < dimension; j++) {
            result.push_back((f1[j] - f2[j]) / (2 * epsilon));
        }
        point1[i] -= epsilon;
        point2[i] += epsilon;
    }

    return result;
}*/
vector<long double> Jacobian(vector<long double>(*f)(vector<long double>), const vector<long double>& point) {
    vector<long double> result;
    long double epsilon = 0.001;
    int dimension = f(point).size();
    result.resize(dimension * dimension);
    vector<long double> point1 = point;
    vector<long double> point2 = point;

    for (int i = 0; i < dimension; i++) {
        point1[i] += epsilon;
        point2[i] -= epsilon;
        vector<long double> f1 = f(point1);
        vector<long double> f2 = f(point2);
        for (int j = 0; j < dimension; j++) {
            result[j*dimension + i] = (f1[j] - f2[j]) / (2 * epsilon);
            }       
        point1[i] -= epsilon;
        point2[i] += epsilon;
        }
    return result;
}

// This function takes in a vector valued function and a point in the phase space, then integrates it using the
// Runge-Kutta method of order 4. The function returns the value of the function at the next time step.
vector<long double> RungeKutta4(vector<long double>(*f)(vector<long double>), const vector<long double>& point, long double dt) {
    vector<long double> k1 = f(point);
    vector<long double> k2 = f(VectorAddition(point, ScalarVectorMultiplication(dt / 2, k1)));
    vector<long double> k3 = f(VectorAddition(point, ScalarVectorMultiplication(dt / 2, k2)));
    vector<long double> k4 = f(VectorAddition(point, ScalarVectorMultiplication(dt, k3)));
    vector<long double> result = VectorAddition(
            point, ScalarVectorMultiplication( dt/6,
                                               VectorAddition(
                                                       k1, VectorAddition(
                                                               ScalarVectorMultiplication(2, k2), VectorAddition(
                                                                       ScalarVectorMultiplication(2, k3), k4)))));
    return result;
}

// This function takes in a vector valued function and an initial point in the phase space and a boolean that decides
// to give all values corresponding to all times. It then integrates the
// function using the Runge-Kutta method of order 4 for a given number of steps. The function returns a vector of vectors
// containing the values of the function at each time step if the boolean is true, otherwise it returns the value of the latest time step.
vector<vector<long double>> IntegrateSystem(vector<long double>(*f)(vector<long double>), const vector<long double>& point, long double dt, int numSteps, bool allValues) {
    int dimension = f(point).size();

    vector<vector<long double>> result;
    vector<long double> currentPoint = point;
    for (int i = 0; i < numSteps; i++) {
        currentPoint = RungeKutta4(f, currentPoint, dt);
        if (allValues) {
            result.push_back(currentPoint);
        }
    }
    if (!allValues) {
        result.push_back(currentPoint);
    }
    return result;
}

#endif //LYAPUNOVEXPONENTSV_0_3_FUNCTIONALANALYSIS_H

