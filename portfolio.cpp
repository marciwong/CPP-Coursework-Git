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
std::vector<std::vector <double> > getPortRetAndStdMat(std::vector< std::vector<double> > inverseMet, std::vector< std::vector<double> > covarianceMet, int noOfCompany, mat tempInv, std::vector<double> vectorOfCompanyMeanRet);


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
    
    std::vector< std::vector<double> > inverseMet;
    inverseMet = covarianceMet;
   
    inverseMet.push_back(meanRetCompany); 
    inverseMet.push_back(eVector);

    //inverseMet = join_vert(inverseMet, meanRetCompany, eVector);

    for (int i = 0; i < noOfCompany; i++)
    {
      inverseMet[83][i] = meanRetCompany[i];
      inverseMet[84][i] = eVector[i];
    }

    //create 2 zeros into the matrix as shown in the pdf file
    inverseMet[83][83] = 0;
    inverseMet[83][84] = 0;
    inverseMet[84][83] = 0;
    inverseMet[84][84] = 0;

    mat tempInv = randu(85,85); // using armadillo 

    // assigning matrix tempInv = inverseMet
    
    for (int i = 0; i < noOfCompany + 2; i++)
    {
      for (int j = 0; i < noOfCompany + 2; i++)
      {
        tempInv(i+1,j+1) = inverseMet[i][j];
      }
    }

  std::vector<std::vector<double> >portfolioReturnAndStd = getPortRetAndStdMat(inverseMet, covarianceMet, noOfCompany, tempInv, vectorOfCompanyMeanRet);





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


std::vector<std::vector <double> > getPortRetAndStdMat(std::vector< std::vector<double> > inverseMet, std::vector< std::vector<double> > covarianceMet, int noOfCompany, mat tempInv, std::vector<double> vectorOfCompanyMeanRet )
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
