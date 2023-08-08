//
// Created by berk on 3/20/2023.
//

#ifndef LYAPUNOVEXPONENTSV_0_3_SYSTEMS_H
#define LYAPUNOVEXPONENTSV_0_3_SYSTEMS_H

#include <iostream>
#include <vector>
#include <cmath>
#include <chrono>
#include "MatrixOperations.h"

using namespace std;

// This function is the Lorenz system. It takes in a vector of size 3 and returns a vector of size 3.
vector<long double> Lorenz(vector<long double> point) {
    vector<long double> result;
    long double sigma = 16.0;
    long double rho = 45.2;
    long double beta = 4.0;
    result.push_back(sigma * (point[1] - point[0]));
    result.push_back(point[0] * (rho - point[2]) - point[1]);
    result.push_back(point[0] * point[1] - beta * point[2]);
    return result;
}

vector<long double> F(vector<long double> x) {
    vector<long double> result;
    long double m1 = 1.0;
    long double m2 = 60.0;
    result.resize(8);
    result[0] = 2*x[4];
    result[1] = 2*x[5];
    result[2] = 2*x[6];
    result[3] = 2*x[7];
    result[4] = -2*(x[0]*m1 + 2*m2*pow(x[0], 3) + 6*m2*x[0]*(x[1]*x[1] + x[2]*x[2] + x[3]*x[3]));
    result[5] = -2*(2*x[1] + m1*x[1] + m2*x[1]*(x[1]*x[1] + x[2]*x[2] + x[3]*x[3])
            + m2*x[1]*(6*x[0]*x[0] + x[1]*x[1] + x[2]*x[2] + x[3]*x[3]));
    result[6] = -2*(2*x[2] + m1*x[2] + m2*x[2]*(x[1]*x[1] + x[2]*x[2] + x[3]*x[3])
            + m2*x[2]*(6*x[0]*x[0] + x[1]*x[1] + x[2]*x[2] + x[3]*x[3]));
    result[7] = -2*(2*x[3] + m1*x[3] + m2*x[3]*(x[1]*x[1] + x[2]*x[2] + x[3]*x[3])
            + m2*x[3]*(6*x[0]*x[0] + x[1]*x[1] + x[2]*x[2] + x[3]*x[3]));
    return result;
}

long double KroneckerDelta(int i, int j){
    if (i == j){
        return 1.0;
    }
    else{
        return 0.0;
    }
}

vector<vector<long double>> RealPartSU2(int N){
    vector<vector<long double>> result;
    result.resize(3);
    for (int i = 0; i < 3; i++){
        result[i].resize(N*N);
    }
    for (int k = 1; k<=N; k++) {
        for (int l = 1; l<=N; l++){
            result[0][(k-1)*N+(l-1)] = (sqrt((N-k)*k)/2.0*KroneckerDelta(l, k+1)
                    +sqrt((N-l)*l)/2.0*KroneckerDelta(k, l+1));
            result[2][(k-1)*N+(l-1)] = ((N+1)/2.0-k)*KroneckerDelta(k,l);
        }
    }
    return result;
}

vector<vector<long double>> ComplexPartSU2(int N){
    vector<vector<long double>> result;
    result.resize(3);
    for (int i = 0; i < 3; i++){
        result[i].resize(N*N);
    }
    for (int k = 1; k <=N ; k++) {
        for (int l = 1; l<= N; l++){
            result[1][(k-1)*N+(l-1)] = (-sqrt((N-k)*k)/2.0*KroneckerDelta(l, k+1)
                    +sqrt((N-l)*l)/2.0*KroneckerDelta(k, l+1));

        }
    }
    return result;
}

vector<long double> F1(const vector<long double>& X1, const vector<long double>& X2, const vector<long double>& Y1, const vector<long double>& Y2, vector<vector<long double>> C, vector<vector<long double>> R, long double m1, long double m2){
    vector<long double> result;
    int dimension = X1.size();
    result = MatrixProduct(X1, VectorSubtraction(
            MatrixProduct(Y1,Y1), MatrixProduct(Y2,Y2)));
    result = VectorSubtraction(result, MatrixProduct(X2, VectorAddition(
            MatrixProduct(Y1,Y2), MatrixProduct(Y2,Y1))));
    result = VectorSubtraction(result, MatrixProduct(VectorAddition(
            MatrixProduct(Y1,Y2),MatrixProduct(Y2,Y1)), X2));
    result = VectorAddition(result, MatrixProduct(VectorSubtraction(
            MatrixProduct(Y1,Y1), MatrixProduct(Y2,Y2)), X1));
    result = ScalarVectorMultiplication(m2/2, result);

    result = VectorAddition(result, ScalarVectorMultiplication((m1+(dimension-1)/2.0), X1));

    result = VectorAddition(result, ScalarVectorMultiplication(2.0, MatrixProduct(MatrixProduct(R[0], X1), R[0])));
    result = VectorAddition(result, ScalarVectorMultiplication(2.0, MatrixProduct(MatrixProduct(R[2], X1), R[2])));
    result = VectorSubtraction(result, ScalarVectorMultiplication(2.0, MatrixProduct(MatrixProduct(C[1], X1), C[1])));

    result = ScalarVectorMultiplication(-1.0, Transpose(result));
    return result;
}

vector<long double> F2(const vector<long double>& X1, const vector<long double>& X2, const vector<long double>& Y1, const vector<long double>& Y2, vector<vector<long double>> C, vector<vector<long double>> R, long double m1, long double m2){
    vector<long double> result;
    int dimension = X1.size();
    result = MatrixProduct(X2, VectorSubtraction(
            MatrixProduct(Y1,Y1), MatrixProduct(Y2,Y2)));
    result = VectorAddition(result, MatrixProduct(X1, VectorAddition(
            MatrixProduct(Y1,Y2), MatrixProduct(Y2,Y1))));
    result = VectorAddition(result, MatrixProduct(VectorAddition(
            MatrixProduct(Y1,Y2),MatrixProduct(Y2,Y1)), X1));
    result = VectorAddition(result, MatrixProduct(VectorSubtraction(
            MatrixProduct(Y1,Y1), MatrixProduct(Y2,Y2)), X2));
    result = ScalarVectorMultiplication(m2/2, result);

    result = VectorAddition(result, ScalarVectorMultiplication((m1+(dimension-1)/2.0), X2));

    result = VectorAddition(result, ScalarVectorMultiplication(2.0, MatrixProduct(MatrixProduct(R[0], X2), R[0])));
    result = VectorAddition(result, ScalarVectorMultiplication(2.0, MatrixProduct(MatrixProduct(R[2], X2), R[2])));
    result = VectorSubtraction(result, ScalarVectorMultiplication(2.0, MatrixProduct(MatrixProduct(C[1], X2), C[1])));
    result = ScalarVectorMultiplication(-1.0, Transpose(result));
    return result;
}

vector<long double> G(vector<long double> x){
    long double m1 = 1.0;
    long double m2 = 10.0;
    int dimension = 4;
    int N = sqrt(dimension);
    vector<vector<long double>> C = ComplexPartSU2(N);
    vector<vector<long double>> R = RealPartSU2(N);
    vector<long double> result;
    vector<long double> X1;
    for (int i = 0; i < dimension; i++){
        X1.push_back(x[i]);
    }
    vector<long double> X2;
    for (int i = dimension; i < 2*dimension; i++){
        X2.push_back(x[i]);
    }
    vector<long double> Y1;
    for (int i = 2*dimension; i < 3*dimension; i++){
        Y1.push_back(x[i]);
    }
    vector<long double> Y2;
    for (int i = 3*dimension; i < 4*dimension; i++){
        Y2.push_back(x[i]);
    }
    vector<long double> PX1;
    for (int i = 4*dimension; i < 5*dimension; i++){
        PX1.push_back(x[i]);
    }
    vector<long double> PX2;
    for (int i = 5*dimension; i < 6*dimension; i++){
        PX2.push_back(x[i]);
    }
    vector<long double> PY1;
    for (int i = 6*dimension; i < 7*dimension; i++){
        PY1.push_back(x[i]);
    }
    vector<long double> PY2;
    for (int i = 7*dimension; i < 8*dimension; i++){
        PY2.push_back(x[i]);
    }

    vector<long double> px1Dot = F1(X1, X2, Y1, Y2, C, R, m1, m2);
    vector<long double> px2Dot = F2(X1, X2, Y1, Y2, C, R, m1, m2);
    vector<long double> py1Dot = F1(Y1, Y2, X1, X2, C, R, m1, m2);
    vector<long double> py2Dot = F2(Y1, Y2, X1, X2, C, R, m1, m2);
    // We change everything with transposed matrices.
    result = ExtendVector(Transpose(PX1), Transpose(PX2));
    result = ExtendVector(result, Transpose(PY1));
    result = ExtendVector(result, Transpose(PY2));

    result = ExtendVector(result, px1Dot);
    result = ExtendVector(result, px2Dot);
    result = ExtendVector(result, py1Dot);
    result = ExtendVector(result, py2Dot);

    return result;
}

vector<long double> ModifiedDelta(int dim){
    vector<long double> result;
    result.resize(dim*dim*dim*dim);
    for (int i = 0; i < dim; i++){
        for (int j = 0; j < dim; j++){
            for (int k = 0; k < dim; k++){
                for (int l = 0; l < dim; l++){
                    if (j*dim+i == k*dim+l){
                        result[(i*dim+j)*dim*dim + k*dim+l] = 1;
                    }
                    else{
                        result[(i*dim+j)*dim*dim + k*dim+l] = 0;
                    }
                }
            }
        }
    }
    return result;
}

vector<long double> K(vector<long double> x){
    long double k = 1.0;
    long double n = 2.0;
    long double m = 1.0;
    long double pi = 3.14159265358979323846264338327950;
    vector<long double> result;
    result = {x[1]/(2*(n-1)), -(n-1)*x[0]*(k*k*m*m*+4*pi*(12*pi*x[0]*x[0]*x[0]*x[0]-3*pi*x[2]*x[2]*x[2]*x[2]+x[0]*x[0]*(4*k*m-6*pi*x[2]*x[2])))/(k*k), x[3]/(2*(n-1)), 
    -(n-1)*x[2]*(k*k*m*m+4*pi*(-3*pi*x[0]*x[0]*x[0]*x[0]+12*pi*x[2]*x[2]*x[2]*x[2]-2*x[2]*x[2]*(2*k*m+3*pi*x[0]*x[0])))/(k*k)};
    return result;
}

#endif //LYAPUNOVEXPONENTSV_0_3_SYSTEMS_H
