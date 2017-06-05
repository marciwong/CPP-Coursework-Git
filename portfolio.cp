#include <cmath>
#include <iostream>
#include <string.h>
#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <string>
#include <numeric>

#include "portfolio.h"
#include "csv.h"
#include "matrixOperations.h"
#include "statisticalOperations.h"

using namespace std;
// construction of portfolio
Portfolio::Portfolio(vector< vector<double> > inSampleMat, vector<double> vectorOfCompanyMeanRet, int noOfCompany, int inSampleRollingWindowSize, int numberOfDays, int outOfSampleRollingWindowSize, double noOfTargetReturn, vector<vector<double > > outOfSampleReturn)
{
    
    outOfSampleAverageReturn.resize((1),vector<double> (83));
    vector<vector<double> > tempNegativeRet (1, vector<double> (83));
    
    for (int i = 0; i < 83; i++)
    {
        tempNegativeRet[0][i] = -1 * vectorOfCompanyMeanRet[i];
    }

    for (int k = 0; k < 83; k++)
    {
        for (int i = 0; i < outOfSampleRollingWindowSize; i++)
        {
            outOfSampleAverageReturn[0][k] =  StatisticalOperations::mean(outOfSampleReturn[k]);
        }   
    }
    // creating in sample covariance by calling get covariance function from Statistical Operations class (static)
    inSampleCovariance = StatisticalOperations::getCovariance(inSampleMat, noOfCompany, inSampleRollingWindowSize);
    // creating out of sample covariance by calling get covariance function from Statistical Operations class (static)
    outOfSampleCovariance = StatisticalOperations::getCovariance(outOfSampleReturn, noOfCompany, outOfSampleRollingWindowSize);
    
    // creating Q matrix
    Q.resize((85),vector<double> (85));
    
    for (int j = 0; j < noOfCompany; j++)
    {
        for(int k = 0; k < noOfCompany ; k++)
        {
            Q[j][k] = inSampleCovariance[j][k];
        }
    }

    for (int j = 0; j < noOfCompany + 2; j++)
    {
        Q[j][83] = tempNegativeRet[0][j];
        Q[j][84] = -1;
        Q[83][j] = tempNegativeRet[0][j];
        Q[84][j] = -1;
    }

    Q[83][83] = 0;
    Q[83][84] = 0;
    Q[84][83] = 0;
    Q[84][84] = 0;
    
    //creating temp portfolio weight vector
    vector <double> tempPortfolioWeight(noOfCompany);
    // getWeights function returns the weights after optimisation (Conjugate Gradient Method)
    tempPortfolioWeight = StatisticalOperations::getWeights(Q, noOfCompany, noOfTargetReturn);
    // transforming it into vector of vector since all of the matrix operations are in vector of vector form (i.e. 1x83 or 83x1)
    vector <vector<double> > portfolioWeights;
    portfolioWeights.push_back(tempPortfolioWeight);
    
    // return portfolio varaince and portfolio return
    portfolioCovariance = MatrixOperations::multiple(MatrixOperations::transpose(portfolioWeights), MatrixOperations::multiple(outOfSampleCovariance,portfolioWeights))[0][0];
    actualAverageReturn = MatrixOperations::multiple(MatrixOperations::transpose(outOfSampleAverageReturn),portfolioWeights)[0][0];
};

// all get functions are declared here
vector<vector<double> > Portfolio::getPortfolioWeights()
{
    return portfolioWeight;
};

vector<vector<double> > Portfolio::getPortfolioInSampleCovariance()
{
    return inSampleCovariance;
};

vector<vector<double> > Portfolio::getQ()
{
    return Q;
}

vector<vector<double> > Portfolio::getPortfolioOutOfSampleCovariance()
{
    return outOfSampleCovariance;
}

double Portfolio::getPortfolioCovariance()
{
    return portfolioCovariance;
}

double Portfolio::getPortfolioAverageReturn()
{
    return actualAverageReturn;
}
