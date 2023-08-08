#include "cModules/MatrixOperations.h"
#include "cModules/Systems.h"
#include "cModules/LyapunovAnalysis.h"
#include <vector>
#include <iostream>

using namespace std;

int main(){
    vector<long double> lyaps = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16};
    vector<long double> logs = LogNorms(lyaps);
    PrintVector(logs);
    return 0;
}