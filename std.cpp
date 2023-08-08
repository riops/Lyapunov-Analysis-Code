#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <tuple>

using namespace std;

// Function to calculate the average of a vector of numbers
double CalculateAverage(const std::vector<long double>& numbers) {
    long double sum = 0.0;
    for (const auto& num : numbers) {
        sum += num;
    }
    return static_cast<double>(sum) / numbers.size();
}

// Function to calculate the standard deviation of a vector of numbers
double CalculateStandardDeviationError(const std::vector<long double>& numbers) {
    long double sum = 0.0;
    long double average = CalculateAverage(numbers);
    for (const auto& num : numbers) {
        sum += std::pow(num - average, 2);
    }
    return sqrt(static_cast<double>(sum) / numbers.size());
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

vector<long double> SeparateAsMinStd(vector<long double> numbers){
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
    for (int i = 0; i < stdErrors.size(); i++){
        cout << stdErrors[i] << " ";
    }
    cout << endl;

    long double currentMin = stdErrors[0];
    int separationPoint = 1;
    for (int i = 0; i < numberSize; i++){
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

    return {realAverage1, realStdError1, realAverage2, realStdError2};
} 


int main() {
    // Input vector of long doubles
    vector<long double> input = {0.647276, 0.686407, 0.195412, 0.192694, 0.69675, 0.625826, 0.614601, 0.608531, 0.581541, 0.633758, 0.651531, 0.678545, 0.194361, 0.611743, 0.4843, 0.19342};

    vector<long double> results = SeparateAsMinStd(input);

    cout << "Group 1 average: " << results[0] << endl;
    cout << "Group 1 standard deviation error: " << results[1] << endl << endl;

    cout << "Group 2 average: " << results[2] << endl;
    cout << "Group 2 standard deviation error: " << results[3] << endl;
    return 0;
}
