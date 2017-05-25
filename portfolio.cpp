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

double mean(vector<double> input);
double meanArray(double input[]);
double standardDeviation(vector<double> input , double mean);
std::vector< std::vector<double> > getCovariance( std::vector< std::vector<double> > returnVector, int size, int timeLength);
std::vector<double> getWeights(std::vector< std::vector<double> > Q, double numberOfCompany, double noOfTargetReturn);
std::vector<std::vector <double> > Multiplication(std::vector<std::vector<double> > A, std::vector<std::vector<double> > B);
std::vector<std::vector<double> > Minus(std::vector<std::vector<double> > A, std::vector<std::vector<double> > B);
std::vector<std::vector<double> > Plus (std::vector<std::vector<double> > A, std::vector<std::vector<double> > B);
double ATransposeA(std::vector<std::vector<double> > A);
std::vector<double> Minus1D(std::vector<double> A, std::vector<double> B);
std::vector<double> Plus1D(std::vector<double> A, std::vector<double> B);

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

Portfolio::Portfolio(std::vector< std::vector<double> > inSampleMat, std::vector<double> vectorOfCompanyMeanRet, int noOfCompany, int inSampleRollingWindowSize, int numberOfDays, int outOfSampleRollingWindowSize, double noOfTargetReturn)
{   
    std::vector<double> tempEVector;
    std::vector<double> tempNegativeRet;

    for (int j = 0; j < noOfCompany ; j++)
    {
        tempEVector.push_back(-1);
        tempNegativeRet.push_back(vectorOfCompanyMeanRet[j] * -1);
    }
    tempEVector.push_back(0);
    tempEVector.push_back(0);
    tempNegativeRet.push_back(0);
    tempNegativeRet.push_back(0);

    std::vector< std::vector<double> > covarainceMat;
    covarainceMat = getCovariance(inSampleMat, noOfCompany, inSampleRollingWindowSize);

    std::vector< std::vector<double> > Q;
    std::vector<double> oneDzeros;
    for (int j = 0; j < noOfCompany; j++)
    {   
        Q.push_back(tempEVector);
    }

 
    for (int j = 0; j < noOfCompany; j++)
    {
        for(int k = 0; k < noOfCompany; k++)
        {
            Q[j][k] = covarainceMat[j][k];   
        }
        Q[j].push_back(vectorOfCompanyMeanRet[j]);
        Q[j].push_back(-1);
    }

    Q.push_back(tempNegativeRet);
    Q.push_back(tempEVector);

    // cout << Q[0][0] << endl;
    // cout << Q[84][84] << endl;

    portfolioWeight = getWeights(Q, noOfCompany, noOfTargetReturn);

    // cout << portfolioWeight[0].size() << endl;
    // cout << portfolioWeight.size() << endl;

};  

std::vector<double> Portfolio::getPortfolioWeights()
{
    return portfolioWeight;
};


//==========================================================================================================
//outside functions

double mean(vector<double> input)
{
       double sum = 0;
       for(int i = 0 ; i < input.size() ; i++ )
       {
          sum += input[i];
       }  
       double average = (sum / input.size());
   return average;
}

double standardDeviation(vector<double> input , double mean)
{
       double sumSq = 0;
       for(int i = 0 ; i < input.size() ; i++)
       {
          sumSq += (input[i] - mean) * (input[i] - mean);
       }
       double std = (sqrt(sumSq / (input.size() - 1 )));
       return std;
}

std::vector<std::vector <double> > Multiplication(std::vector<std::vector<double> > A, std::vector<std::vector<double> > B)
{    
    std::vector< std::vector <double> > multiple;
    std::vector<double> zeros;
    double sum;
        // Initializing elements of matrix mult to 0.

        for (int i = 0; i < A[0].size(); i++)
        {
            zeros.push_back(0.0);
        }
        for(int j = 0; j < B.size(); j++)
        {
            multiple.push_back(zeros);
        }

    // Multiplying matrix a and b and storing in array mult.
        for(int i = 0; i < A[0].size(); i++)
        {
            for(int j = 0; j < B.size(); j++)
            {
                for(int k = 0; k < A.size(); k++)
                {
                    multiple[j][i] += A[k][i] * B[j][k];
                }
            }
        }
    return multiple;
}

std::vector<std::vector<double> > Minus(std::vector<std::vector<double> > A, std::vector<std::vector<double> > B)
{
    std::vector<std::vector<double> > sum;
    std::vector<double> zeros;

        for(int j = 0; j < B.size(); ++j)
        {
        for (int i = 0; i < A[0].size(); i++)
        {
            zeros.push_back(0);
        }
            sum.push_back(zeros);
        }

        for(int i = 0; i < B.size(); ++i)
        {
            for(int j = 0; j < A[0].size(); ++j)
            {
                sum[i][j] += A[i][j] - B[i][j];    
            }
        }
    return sum;
}

std::vector<std::vector<double> > Plus(std::vector<std::vector<double> > A, std::vector<std::vector<double> > B)
{
    std::vector<std::vector<double> > sum;
    std::vector<double> zeros;

        for(int j = 0; j < B.size(); ++j)
        {
        for (int i = 0; i < A[0].size(); i++)
        {
            zeros.push_back(0);
        }
            sum.push_back(zeros);
        }

        for(int i = 0; i < B.size(); ++i)
        {
            for(int j = 0; j < A[0].size(); ++j)
            {
                sum[i][j] += A[i][j] + B[i][j];    
            }
        }
    return sum;
}


double ATransposeA(std::vector<std::vector<double> > A)
{
    double sum = 0;
    for (int i = 0; i < A.size(); i++)
    {
        sum = sum + A[0][i]*A[0][i];
    }
    return sum;
}

std::vector< std::vector<double> > getCovariance(std::vector< std::vector<double> > returnVector, int numberOfCompany, int timeLength)
{
    // Initialise cov matrix
    std::vector< std::vector<double> > cov;
    for (int i = 0; i < numberOfCompany; i++) 
    {
      std::vector<double> column;
      for (int j = 0; j < numberOfCompany; j++) 
      {
        column.push_back(0);
      }
      cov.push_back(column);
    }

    std::vector<double> firstCompany;
    std::vector<double> secondCompany;

    for (int j = 0; j < timeLength; j++)
    {
      firstCompany.push_back(0);
      secondCompany.push_back(0);
    }

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


std::vector<double> getWeights(std::vector< std::vector<double> > Q, double numberOfCompany, double noOfTargetReturn)
{   
    std::vector<double> weights;
    std::vector<double> eightyThreezeros;
    for (int i = 0; i < 83; i++)
    {
        weights.push_back(0);
        eightyThreezeros.push_back(0);
    }

        std::vector<std::vector<double> > s;
        std::vector<std::vector<double> > b;
        std::vector<double> bZeros;
        std::vector<std::vector<double> > x;
        std::vector<std::vector<double> > Qx;
        std::vector<double> zeros;
        std::vector<double> xVector;

        for (int i = 0; i < 83; i++)
        {
            zeros.push_back(0.0);
            xVector.push_back(0.01204819277);
            bZeros.push_back(0.0);

        }
        zeros.push_back(0.0); //portfolio return
        zeros.push_back(-1.0);

        xVector.push_back(1.0); // adding lagrangian multipliers into the vector
        xVector.push_back(1.0);

        double negativeTargetReturn = -1 * noOfTargetReturn;
        bZeros.push_back(negativeTargetReturn);
        bZeros.push_back(-1.0);
        b.push_back(bZeros);
        s.push_back(zeros);
        x.push_back(xVector);
        Qx.push_back(zeros);

        double alpha;
        double sTs0 = 1.0;
        double sTs1;
        double beta;
        std::vector<std::vector<double> > p;
        std::vector<std::vector<double> > pT;
        std::vector<std::vector<double> > Qp;
        std::vector<std::vector<double> > alphaQp;

        p.push_back(zeros);
        for (int i = 0; i < pT.size(); i++)
        {
            pT[i].push_back(0.0);
        }
        double pQpTtemp;
        // // //initialise 
        int k = 0;
        int i = 0;
        Qx = Multiplication(Q,x);
        for (int i = 0; i < 85; i++)
        {
            cout << Qx[0][i] << endl;
        }
        cout << endl;
        s = Minus(b,Qx);
        p = s;

        while (sTs0 <= 0.000006)
        {   
            std::vector<double> temp;
            for (int i = 0; i < p[0].size(); i++)
            {
                temp.push_back(p[0][i]);
            }
            pT.push_back(temp);
            Qp = Multiplication(Q,p);
            pQpTtemp = Multiplication(p, Qp)[0][0];
            alpha = ATransposeA(s) / pQpTtemp;
            for (int i = 0; i < pT.size(); i++)
            {
                pT[i][0] = alpha * pT[i][0];
            }
            x = Plus(x, pT);
            for (int i = 0; i < Qp.size(); i++)
            {
                alphaQp[0][i] = alpha * Qp[0][i];
            }
            sTs0 = ATransposeA(s);
            s = Minus(s, alphaQp);
            sTs1 = ATransposeA(s);
            beta = sTs1 / sTs0;
        }

        int noOfWeights = 0;
        for (int i = 0; i < (x[0].size()-2); i++)
        {
            weights[i] = x[0][i];
        }

return weights;

}




