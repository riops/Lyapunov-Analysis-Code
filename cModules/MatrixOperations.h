//
// Created by berk on 3/20/2023.
//

#ifndef LYAPUNOVEXPONENTSV_0_3_MATRIXOPERATIONS_H
#define LYAPUNOVEXPONENTSV_0_3_MATRIXOPERATIONS_H

#include <iostream>
#include <vector>
#include <cmath>
#include <chrono>
#include <fstream>

using namespace std;

// This function prints a matrix that is given as a vector.
void PrintMatrix(vector<long double> matrix){
    int dimension = sqrt(matrix.size());
    for (int i = 0; i < dimension; i++){
        cout << "| ";
        for (int j = 0; j < dimension; j++){
            cout << matrix[i * dimension + j] << " ";
        }
        cout << "|";
        cout << endl;
    }
    cout << endl;
}

// This function writes matrix to the given ofstream.
void WriteMatrix(vector<long double> matrix){
    ofstream file("IntegrationResults.txt", ios_base::app);
    int dimension = sqrt(matrix.size());
    for (int i = 0; i < dimension; i++){
        for (int j = 0; j < dimension; j++){
            file << matrix[i * dimension + j] << " ";
        }
        file << endl;
    }
    file << endl << endl;
}


// This function prints a vector.
void PrintVector(vector<long double> vector){
    int dimension = vector.size();
    for (int i = 0; i < dimension; i++){
        cout << vector[i] << " ";
    }
    cout << endl;
}

vector<long double> ZeroMatrix(int dimension){
    vector<long double> result;
    result.resize(dimension * dimension);
    for (int i = 0; i < dimension; i++){
        for (int j = 0; j < dimension; j++){
            result[i * dimension + j] = 0;
        }
    }
    return result;
}


// This function returns the identity matrix of a given dimension as a vector.
vector<long double> IdentityMatrix(int dimension){
    vector<long double> result;
    result.resize(dimension * dimension);
    for (int i = 0; i < dimension; i++){
        for (int j = 0; j < dimension; j++){
            if (i == j){
                result[i * dimension + j] = 1;
            }
            else{
                result[i * dimension + j] = 0;
            }
        }
    }
    return result;
}

// This function extends a vector with another vector.
vector<long double> ExtendVector(vector<long double> vector1, vector<long double> vector2){
    vector<long double> result;
    int dimension1 = vector1.size();
    int dimension2 = vector2.size();
    result.resize(dimension1 + dimension2);
    for (int i = 0; i < dimension1; i++){
        result[i] = vector1[i];
    }
    for (int i = 0; i < dimension2; i++){
        result[dimension1 + i] = vector2[i];
    }
    return result;
}

// This function adds two vectors/ or matrices that are given as vectors.
vector<long double> VectorAddition(vector<long double> vector1, vector<long double> vector2){
    vector<long double> result;
    int dimension = vector1.size();
    result.resize(dimension);
    for (int i = 0; i < dimension; i++){
        result[i] = vector1[i] + vector2[i];
    }
    return result;
}

// This function subtracts two vectors/ or matrices that are given as vectors.
vector<long double> VectorSubtraction(vector<long double> vector1, vector<long double> vector2){
    vector<long double> result;
    int dimension = vector1.size();
    result.resize(dimension);
    for (int i = 0; i < dimension; i++){
        result[i] = vector1[i] - vector2[i];
    }
    return result;
}

// This function calculates the transpose of a matrix that is given as a vector.
vector<long double> Transpose(vector<long double> matrix){
    vector<long double> result;
    int dimension = sqrt(matrix.size());
    result.resize(dimension * dimension);
    for (int i = 0; i < dimension; i++){
        for (int j = 0; j < dimension; j++){
            result[i * dimension + j] = matrix[j * dimension + i];
        }
    }
    return result;
}

// This function calculates the product of two matrices that are given as vectors.
vector<long double> MatrixProduct(vector<long double> matrix1, vector<long double> matrix2){
    vector<long double> result;
    int dimension = sqrt(matrix1.size());
    result.resize(dimension * dimension);
    for (int i = 0; i < dimension; i++){
        for (int j = 0; j < dimension; j++){
            for (int k = 0; k < dimension; k++){
                result[i * dimension + j] += matrix1[i * dimension + k] * matrix2[k * dimension + j];
            }
        }
    }
    return result;
}

// This function calculates the product of a number and a vector.
vector<long double> ScalarVectorMultiplication(long double scalar, vector<long double> vector1){
    vector<long double> result;
    int dimension = vector1.size();
    result.resize(dimension);

    for (int i = 0; i < dimension; i++){
        result[i] = vector1[i] * scalar;
    }
    return result;
}

// This function calculates the inner product of two vectors.
long double InnerProduct (vector<long double> vector1, vector<long double> vector2){
    long double result = 0.0;
    int dimension = vector1.size();
    for (int i = 0; i < dimension; i++){
        result += vector1[i] * vector2[i];
    }
    return result;
}

// This function calculates the norm of a given vector.
long double Norm (const vector<long double>& vector1){
    long double result;
    result = InnerProduct(vector1, vector1);
    return sqrt(result);
}

// This function calculates the projection of a vector onto another vector.
vector<long double> Projection (vector<long double> vector1, vector<long double> vector2){
    vector<long double> result;
    result.resize(vector1.size());
    if (Norm(vector2) == 0){
        return result;
    }
    else{
        long double scalar = InnerProduct(vector1, vector2) / InnerProduct(vector2, vector2);
        result = ScalarVectorMultiplication(scalar, vector2);
        return result;
    }
}

// This function orthonormalizes given vectors. Vectors are given as a vector, if the given vector is n dimensional
// then there are sqrt(n) vectors that are sqrt(n) dimensional. It uses the Gram-Schmidt process.

vector<long double> Orthogonalize(vector<long double> vectors){
    vector<long double> result;
    int dimension = sqrt(vectors.size());
    result.resize(dimension * dimension);
    for (int i = 0; i <dimension; i++){
        vector<long double> currentVector;
        currentVector.resize(dimension);
        for (int j = 0; j < dimension; j++){
            currentVector[j] = vectors[i * dimension + j];
        }
        for (int j = 0; j < i; j++){
            vector<long double> previousVector;
            previousVector.resize(dimension);
            for (int k = 0; k < dimension; k++){
                previousVector[k] = result[j * dimension + k];
            }
            currentVector = VectorSubtraction(currentVector, Projection(currentVector, previousVector));
        }
        for (int j = 0; j < dimension; j++){
            result[i * dimension + j] = currentVector[j];
        }
    }
    return result;
}


// This function will be used to rewrite a block of a matrix that is given as a vector with another matrix that is also given as a vector.
vector<long double> WriteBlock(vector<long double> matrix, vector<long double> block, int i, int j, int numVariables){
    int dimension = sqrt(matrix.size());
    for(int k = 0; k < dimension/numVariables; k++){
      for(int l = 0; l < dimension/numVariables; l++){
          matrix[(k + i * dimension/numVariables) * dimension + (l + j * dimension/numVariables)] = block[k * dimension/numVariables + l];
      }
    }
    return matrix;
}

#endif //LYAPUNOVEXPONENTSV_0_3_MATRIXOPERATIONS_H
