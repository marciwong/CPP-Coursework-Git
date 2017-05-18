#define ARMA_DONT_USE_WRAPPER
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

using namespace arma;
using namespace std;

double mean(vector<double> input);
double standardDeviation(vector<double> input , double mean);
std::vector< std::vector<double> > getCovariance( std::vector< std::vector<double> > returnVector, double size, double timeLength);
std::vector<std::vector <double> > getPortRetAndStdMat(std::vector< std::vector<double> > Q, std::vector< std::vector<double> > covarianceMet, int noOfCompany, mat tempInv, std::vector<double> vectorOfCompanyMeanRet);
std::vector<std::vector <double> > Multiplication2D (std::vector<std::vector<double> > A, std::vector<std::vector<double> > B);

Company::Company(){ };

Company::Company(std::vector < std::vector<double> > input, int i, int timeLength)
{
    int Days = input.size();
    double returnArray[timeLength];
    for (int j = 0 ; j < timeLength; j++)
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

    covarianceMet = getCovariance(returnVector, noOfCompany, time);

    std:vector<double> negativeRet;
    for (int i = 0; i<83; i++)
    {
      negativeRet[i] = vectorOfCompanyMeanRet[i] * -1;
    }

    for (int i = 0; i < 83; i++ )
    {
      covarianceMet[i].push_back(negativeRet[i]);
      covarianceMet[i].push_back(eVector[i]);
    }

    std:vector<double> negativeRetWithZero = negativeRet;

    negativeRetWithZero.push_back(0);
    negativeRetWithZero.push_back(0);

    std:vector<double> eVectorWithZero = eVector;
    eVectorWithZero.push_back(0);
    eVectorWithZero.push_back(0);
    
    std::vector< std::vector<double> > Q;
    Q = covarianceMet;
    Q.push_back(negativeRetWithZero);
    Q.push_back(eVectorWithZero);

  std::vector<std::vector<double> >portfolioReturnAndStd = getPortRetAndStdMat(Q, covarianceMet, noOfCompany, tempInv, vectorOfCompanyMeanRet);





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
std::vector< std::vector<double> > getCovariance( std::vector< std::vector<double> > returnVector, double size, double timeLength)
{

    std::vector<double> x;
    std::vector<double> y;
    std::vector < std:: vector<double> > cov;
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
            cov[i][k] += (x[i] - xMean) * (y[j] - yMean) / (size - 1);
        }
       }
    }
    return cov;
}


std::vector<std::vector <double> > getPortRetAndStdMat(std::vector< std::vector<double> > Q, std::vector< std::vector<double> > covarianceMet, int noOfCompany, mat tempInv, std::vector<double> vectorOfCompanyMeanRet )
    {
    mat vectorOfOptimisation = randu(3,1);
    // Creating vector for saving portfolio returns and stds depend on different return
    std::vector<double> portfolioRet;
    std::vector<double> portfolioStd;
    double PortMeanRet;
    double PortStd;
    std::vector<double> vectorOfObjectives;
    mat matOfWeightOfCompany = randu(noOfCompany + 2, 1);

    int i = 0;

    while (i <= 0.1)
    {
    // let say portfolio target return is 10%
    vectorOfOptimisation(1,1) = 0;
    vectorOfOptimisation(2,1) = -1 * i;
    vectorOfOptimisation(3,1) = -1;

    //creating a vector of weights and the two langrangian parameters

    //use matrix multiplication
     matOfWeightOfCompany = tempInv * vectorOfOptimisation;

    for (int j = 0; j < noOfCompany + 2; j++)
    {
        vectorOfObjectives[j] = matOfWeightOfCompany(j+1,1);
    }

    // Obtain weight vector from matofWeights
    std::vector<double> vectorOfWeights;
    for (int k = 0; k < noOfCompany; k++ )
    {
      vectorOfObjectives[k] = vectorOfWeights[k];
    }

    // Obtain the lambda and miu
    double lambda = vectorOfObjectives[83];
    double miu = vectorOfObjectives[84];

    // Multiplying the weight with mean company return
    for (int h = 0; h < noOfCompany; h++)
    {
      PortMeanRet += vectorOfCompanyMeanRet[h] * vectorOfWeights[h];
      PortStd += covarianceMet[h][h] * vectorOfWeights[h] * vectorOfWeights[h];
    }

    portfolioRet.push_back(PortMeanRet);
    portfolioStd.push_back(PortStd);
    
    i = i + 0.005;
    }
    std::vector< std::vector<double> > vectorOfReturnAndStd;

    for (int i = 0; i < portfolioRet.size(); i++)
    {
      vectorOfReturnAndStd[0][i] = portfolioRet[i];
      vectorOfReturnAndStd[1][i] = portfolioStd[i];
    }

    return vectorOfReturnAndStd;
    }
std::vector<std::vector <double> > Multiplication2D (std::vector<std::vector<double> > A, std::vector<std::vector<double> > B)
{    
    int columnLengthA = A.size(); //read matrix size horizontally
    int columnLengthB = B.size();
    int rowLengthB = B[0].size(); //read matrix size vertically
    int rowLengthA = A[0].size();
    std::vector< std::vector <double> > multiple;
    std::vector<double> zeros;
        // Initializing elements of matrix mult to 0.
    for (i = 0; i < rowLengthA; i++)
    {
        zeros.push_back(0);
    }
    for(j = 0; j < columnLengthB; ++j)
    {
        multiple.push_back(zeros);
    }

    // Multiplying matrix a and b and storing in array mult.
    for(i = 0; i < rowLengthA; ++i)
        for(j = 0; j < columnLengthB; ++j)
            for(k = 0; k < columnLengthA; ++k)
            {
                multiple[i][j] += A[i][k] * B[k][j];
            }

    return multiple;
}
