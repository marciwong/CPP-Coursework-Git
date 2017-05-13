#include <cmath>
#include <iostream>
#include <string.h>
#include <stdlib.h>
#include <armadillo>
#include <fstream>
#include <sstream>
#include <string>
#include <numeric>

#include "portfolio.h"
#include "csv.h"

using namespace std;

double mean(vector<double> input);
double standardDeviation(vector<double> input , double mean);

Company::Company(){ };

Company::Company(std::vector < std::vector<double> > input, int i)
{
    Days = input.size();
    double returnArray[700];
    for (int j = 0; j <700; j++)
    {
      allReturnVector.push_back (input[i][j]);
    }
    double meanRet = mean(allReturnVector);
    double stdev = standardDeviation(allReturnVector, meanRet);

};

void getCompanyMeanRet()
{
 return meanRet;
};


Portfolio::Portfolio(){ };

Portfolio::Portfolio(std::vector< std::vector<double> > returnVector, int i, int j)
{
    int N = returnVector.sizeof();


};

double mean(vector<double> input)
{      int sum = 0;
       for(int i = 0 ; i < input.size() ; i++ )
               sum += input[i];
       return sum / input.size();
}

double standardDeviation(vector<double> input , double mean)
{
       double sumSq = 0;
       for(int i = 0 ; i < input.size() ; i++)
       {
          sumSq += (input[i] - mean) * (input[i] - mean);
       }

       return sqrt(sumSq / (input.size()-1));
}