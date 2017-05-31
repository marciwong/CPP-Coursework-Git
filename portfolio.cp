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

using namespace std;

Company::Company(){ };

Company::Company(std::vector<std::vector<double> > input, int i, int timeLength)
{
    int Days = timeLength;
    double returnArray[timeLength];
    std::vector<double> allReturnVector;
    for (int j = 0 ; j < timeLength; j++)
        {
            allReturnVector.push_back(input[i][j]);
        }
    meanRet = mean(allReturnVector);
    stdev = standardDeviation(allReturnVector, meanRet);
};

double Company::getCompanyMeanRet()
{
 return meanRet;
};

Portfolio::Portfolio(std::vector< std::vector<double> > inSampleMat, std::vector<double> vectorOfCompanyMeanRet, int noOfCompany, int inSampleRollingWindowSize, int numberOfDays, int outOfSampleRollingWindowSize, double noOfTargetReturn, std::vector<std::vector<double > > outOfSampleReturn)
{   
    outOfSampleAverageReturn.resize((1),std::vector<double> (83));
    std::vector<std::vector<double> > tempNegativeRet (1, std::vector<double> (83));
    
    for (int i = 0; i < 83; i++)
    {
        tempNegativeRet[0][i] = -1 * vectorOfCompanyMeanRet[i];
    }

    for (int k = 0; k < 83; k++)
    {
        for (int i = 0; i < outOfSampleRollingWindowSize; i++)
        {
            outOfSampleAverageReturn[0][k] =  mean(outOfSampleReturn[k]);
        }   
    }

    inSampleCovariance = getCovariance(inSampleMat, noOfCompany, inSampleRollingWindowSize);
    outOfSampleCovariance = getCovariance(outOfSampleReturn, noOfCompany, outOfSampleRollingWindowSize);

    Q.resize((85),std::vector<double> (85));
    
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

    std::vector <double> tempPortfolioWeight(noOfCompany);
    tempPortfolioWeight = getWeights(Q, noOfCompany, noOfTargetReturn);
    std::vector <std::vector<double> > portfolioWeights;
    portfolioWeights.push_back(tempPortfolioWeight);
    portfolioCovariance = Multiplication(transpose(portfolioWeights),Multiplication(outOfSampleCovariance,portfolioWeights))[0][0];
    actualAverageReturn = Multiplication(transpose(outOfSampleAverageReturn),portfolioWeights)[0][0];
};  

std::vector<std::vector<double> > Portfolio::getPortfolioWeights()
{
    return portfolioWeight;
};

std::vector<std::vector<double> > Portfolio::getPortfolioInSampleCovariance()
{
    return inSampleCovariance;
};

std::vector<std::vector<double> > Portfolio::getQ()
{
    return Q;
}

std::vector<std::vector<double> > Portfolio::getPortfolioOutOfSampleCovariance()
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

//==========================================================================================================
//outside functions

double mean(std::vector<double> input)
{
   double sum = 0.0;
   for(int i = 0 ; i < input.size() ; i++ )
   {
      sum += input[i];
   }  
   double average = (sum / input.size());
   return average;
}

double standardDeviation(std::vector<double> input , double mean)
{
       double sumSQ = 0.0;
       for(int i = 0 ; i < input.size() ; i++)
       {
          sumSQ += (input[i] - mean) * (input[i] - mean);
       }
       double std = (sqrt(sumSQ / (input.size() - 1 )));
       return std;
}

std::vector<std::vector <double> > Multiplication(std::vector<std::vector<double> > A, std::vector<std::vector<double> > B)
{
    std::vector< std::vector <double> > multiple(B.size(), std::vector<double>(A[0].size()));
    // Multiplying matrix a and b and storing in array mult.
    for (int i = 0; i < A[0].size(); i++)
    {
        for (int j = 0; j < B.size(); j++)
        {
            for (int k = 0; k < A.size(); k++)
            {
                multiple[j][i] += A[k][i] * B[j][k];
            }
        }
    }
    return multiple;
}

std::vector<std::vector<double> > Minus(std::vector<std::vector<double> > A, std::vector<std::vector<double> > B)
{
    std::vector<std::vector<double> > sum(B.size(), std::vector<double>(A[0].size()));
    for (int i = 0; i < B.size(); i++)
    {
        for (int j = 0; j < A[0].size(); j++)
        {
            sum[i][j] = A[i][j] - B[i][j];
        }
    }
    return sum;
}

std::vector<std::vector<double> > Plus(std::vector<std::vector<double> > A, std::vector<std::vector<double> > B)
{
    std::vector<std::vector<double> > sum(B.size(), std::vector<double>(A[0].size()));
    for (int i = 0; i < B.size(); i++)
    {
        for (int j = 0; j < A[0].size(); j++)
        {
            sum[i][j] = A[i][j] + B[i][j];    
        }
    }
    return sum;
}


std::vector<std::vector<double> > transpose(std::vector<std::vector<double> > A)
{
    std::vector<std::vector<double> > tempA (A[0].size(),vector <double>(A.size()));

    for (int i = 0; i < A.size(); i++)
    {
        for (int j = 0; j < A[0].size();j++)
        {
            tempA[j][i] = A[i][j];
        }
    }
    return tempA;
};

std::vector<std::vector<double> > scalarMultiplication(double alpha, std::vector<std::vector<double> > A)
{
    std::vector<std::vector<double> > tempA (A.size(), std::vector<double>(A[0].size()));
    for (int i = 0; i < A.size(); i++) {
        for (int j = 0; j < A[0].size(); j++) {
            tempA[i][j] = alpha * A[i][j];
        }
    }
    return tempA;
};



std::vector< std::vector<double> > getCovariance(std::vector< std::vector<double> > returnVector, int numberOfCompany, int timeLength)
{

    std::vector< std::vector<double> > cov(numberOfCompany, std::vector<double>(numberOfCompany));

    std::vector<double> firstCompany(timeLength);
    std::vector<double> secondCompany(timeLength);

    for (int i = 0; i < numberOfCompany; i++)
    {
        for (int k = 0; k < numberOfCompany; k++)
        {  
            for (int j = 0; j < timeLength; j++)
            {
                firstCompany[j] = returnVector[i][j];
                secondCompany[j] = returnVector[k][j];
            }

            double firstCompanyMean = mean(firstCompany);
            double secondCompanyMean = mean(secondCompany);           

            for (int j = 0; j < timeLength; j++)
            {
                cov[i][k] += (firstCompany[j] - firstCompanyMean) * (secondCompany[j] - secondCompanyMean) / (timeLength - 1);    
            }
        }
    }
    return cov;
}


std::vector<double> getWeights(std::vector< std::vector<double> > Q, double numberOfCompany, double noOfTargetReturn) {
    double tolerence = 0.000006;

    // Set up x
    std::vector< std::vector<double> > x(1, std::vector<double>(numberOfCompany + 2));
    for (int i = 0; i < numberOfCompany; i++) 
    {
        x[0][i] = 1.0 / numberOfCompany;
    }
    x[0][numberOfCompany] = 1.0; // lambda
    x[0][numberOfCompany + 1] = 1.0; // mu

    // Set up b
    std::vector< std::vector<double> > b(1, std::vector<double>(numberOfCompany + 2));
    for (int i = 0; i < numberOfCompany; i++) 
    {
        b[0][i] = 0.0;
    }
    b[0][numberOfCompany] = -1.0 * noOfTargetReturn; // -r_p
    b[0][numberOfCompany + 1] = -1.0;

    std::vector< std::vector<double> > s = Minus(b, Multiplication(Q, x));
    std::vector< std::vector<double> > p(s);

    double sTs = Multiplication(transpose(s), s)[0][0];
    while (sTs > tolerence) 
    {
        double alpha = sTs / (Multiplication(Multiplication(transpose(p), Q), p)[0][0]);

        x = Plus(x, (scalarMultiplication(alpha, p)));

        std::vector< std::vector<double> > s_plus1 = Minus(s, (scalarMultiplication(alpha, (Multiplication(Q, p)))));

        sTs = Multiplication(transpose(s_plus1), s_plus1)[0][0];

        double beta = (sTs) / (Multiplication(transpose(s), s)[0][0]);

        p = Plus(s_plus1, (scalarMultiplication(beta, p)));

        s = s_plus1;
    }

    std::vector<double> weights (numberOfCompany);
    for (int i = 0; i < weights.size(); i++)
    {
        weights[i] = x[0][i];
    }

    return weights;
}

void printMatrix(std::vector< std::vector<double> > input) {
    for (int i = 0; i < input.size(); i++) {
        for (int j = 0; j < input[0].size(); j++) {
            cout << input[i][j] << " ";
        }
        cout << "" << endl;
    }
}
