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

Company::Company(std::vector < std::vector<double> > input, int i, int timeLength)
{
    int Days = input.size();
    double returnArray[timeLength];
    for (int j = 0; j < timeLength; j++)
    {
      allReturnVector.push_back (input[i][j]);
    }
    double meanRet = mean(allReturnVector);
    double stdev = standardDeviation(allReturnVector, meanRet);

};

double Company::getCompanyMeanRet()
{
 return meanRet;
};


Portfolio::Portfolio(){ };

Portfolio::Portfolio(std::vector< std::vector<double> > returnVector, std::vector<double> vectorOfCompanyMeanRet, double noOfCompany, double time)
{
    std::vector<double> eVector;
    for (int i = 0; i < noOfCompany; i++)
    {
      eVector[i] = -1;
    }
    std::vector<double> meanRetCompany;
    for (int i = 0; i < noOfCompany; i++)
    {
      meanRetCompany[i] = - 1 * vectorOfCompanyMeanRet[i];
    }

    std::vector < std::vector<double> > covarianceMet;
    covarianceMet = covariance(returnVector, noOfCompany, time);
    std::vector< std::vector<double> > inverseMet;
    for (int i = 0; i < noOfCompany; i++)
    {
      for (int j = 0; j < noOfCompany; j++)
      {
        inverseMet[i][j] = covariance[i][j];
      }
    }
   
    inverseMet.push_back(meanRetCompany); 
    inverseMet.push_back(eVector);

};  

//==========================================================================================================
//outside functions

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
std::vector< std::vector<double> > covariance( std::vector< std::vector<double> > returnVector, double size, double timeLength)
{

    std::vector<double> x;
    std::vector<double> y;
    std::vector<std:: vector<double> > cov;
    for(int i = 0 ; i < size ; i++)
    {
      for (int k = 0; i < size ; k++ )
       {  
        for (int j = 0; j< timeLength ; j++)
        {
          x[j] = returnVector[i][j];
          y[j] = returnVector[k][j];
          double xMean = mean(x);
          double yMean = mean(y);
          if (x.size() == size)
            cov[i][j] += (x[i] - xMean) * (y[j] - yMean) / (size - 1);
        }
       }
    }
    return cov;
}
