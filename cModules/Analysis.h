#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

// Function to calculate the average of a vector of numbers
long double CalculateAverage(const std::vector<long double>& numbers) {
    long double sum = 0.0;
    for (const auto& num : numbers) {
        sum += num;
    }
    return static_cast<double>(sum) / numbers.size();
}

// Function to calculate the standard deviation of a vector of numbers
long double CalculateStandardDeviationError(const vector<long double>& values) {
    // Step 1: Calculate the mean
    long double sum = 0.0;
    for (const auto& value : values) {
        sum += value;
    }
    long double mean = sum / values.size();

    // Step 2: Calculate the sum of squared differences
    long double squaredDifferencesSum = 0.0;
    for (const auto& value : values) {
        squaredDifferencesSum += (value - mean) * (value - mean);
    }

    // Step 3: Divide the sum of squared differences by the number of elements
    long double variance = squaredDifferencesSum / values.size();

    // Step 4: Take the square root of the variance to obtain the standard deviation
    long double standardDeviationError = std::sqrt(variance);

    return standardDeviationError;
}


vector<long double> Sort(vector<long double> numbers){
    vector<long double> sortedNumbers;
    sortedNumbers.resize(numbers.size());
    for (int i = 0; i < numbers.size(); i++){
        sortedNumbers[i] = numbers[i];
    }
    for (int i = 0; i < sortedNumbers.size(); i++){
        for (int j = i + 1; j < sortedNumbers.size(); j++){
            if (sortedNumbers[i] > sortedNumbers[j]){
                long double temp = sortedNumbers[i];
                sortedNumbers[i] = sortedNumbers[j];
                sortedNumbers[j] = temp;
            }
        }
    }
    return sortedNumbers;
}

vector<vector<long double>> SeparateAsMinStd(vector<long double> numbers){
    vector<long double> sortedNumbers = Sort(numbers);
    int numberSize = sortedNumbers.size();

    vector<long double> stdErrors;
    stdErrors.resize(numberSize - 1);
    vector<long double> stdErrors1;
    stdErrors1.resize(numberSize - 1);
    vector<long double> stdErrors2;
    stdErrors2.resize(numberSize - 1);
    vector<long double> averages1;
    averages1.resize(numberSize - 1);
    vector<long double> averages2;
    averages2.resize(numberSize - 1);

    for (int separator = 1; separator < numberSize; separator++){
        int firstGroupSize = separator;
        int secondGroupSize = numberSize - separator;

        vector<long double> group1;
        group1.resize(firstGroupSize);
        vector<long double> group2;
        group2.resize(secondGroupSize);

        for (int i = 0; i < firstGroupSize; i++){
            group1[i] = sortedNumbers[i];
        }
        for (int i = 0; i < secondGroupSize; i++){
            group2[i] = sortedNumbers[i + firstGroupSize];
        }
        long double stdError1 = CalculateStandardDeviationError(group1);
        stdErrors1[separator - 1] = stdError1;
        long double stdError2 = CalculateStandardDeviationError(group2);
        stdErrors2[separator - 1] = stdError2;
        long double average1 = CalculateAverage(group1);
        averages1[separator - 1] = average1;
        long double average2 = CalculateAverage(group2); 
        averages2[separator - 1] = average2; 
        stdErrors[separator - 1] = stdError1 + stdError2;
    }

    long double currentMin = stdErrors[0];
    int separationPoint = 1;
    for (int i = 0; i < numberSize - 1; i++){
        if (stdErrors[i] < currentMin){
            currentMin = stdErrors[i];
            separationPoint = i + 1;
        }
    }

    vector<long double> group1;
    group1.resize(separationPoint);
    vector<long double> group2; 
    group2.resize(numberSize - separationPoint);

    for (int i = 0; i < separationPoint; i++){
        group1[i] = sortedNumbers[i];
    }
    for (int i = 0; i < numberSize - separationPoint; i++){
        group2[i] = sortedNumbers[i + separationPoint];
    }

    long double realAverage1 = CalculateAverage(group1);
    long double realAverage2 = CalculateAverage(group2);

    long double realStdError1 = CalculateStandardDeviationError(group1);
    long double realStdError2 = CalculateStandardDeviationError(group2);

    return {{realAverage1, realStdError1, realAverage2, realStdError2}, group1, group2};
} 
