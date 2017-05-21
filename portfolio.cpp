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
double meanArray(double input[]);
double standardDeviation(vector<double> input , double mean);
std::vector< std::vector<double> > getCovariance( std::vector< std::vector<double> > returnVector, int size, int timeLength);
std::vector<std::vector <double> > getWeights(std::vector< std::vector<double> > Q, double numberOfCompany);
std::vector<std::vector <double> > Multiplication(std::vector<std::vector<double> > A, std::vector<std::vector<double> > B);
std::vector<std::vector<double> > Minus(std::vector<std::vector<double> > A, std::vector<std::vector<double> > B);
std::vector<std::vector<double> > Plus (std::vector<std::vector<double> > A, std::vector<std::vector<double> > B);
double ATransposeA(std::vector<double> A);
std::vector<double> Minus1D(std::vector<double> A, std::vector<double> B);
std::vector<double> Plus1D(std::vector<double> A, std::vector<double> B);

Company::Company(){ };

Company::Company(std::vector < std::vector<double> > input, int i, int timeLength)
{
    int Days = input.size();
    double returnArray[timeLength];
    std::vector<double> allReturnVector;
    for (int j = 0 ; j < timeLength; j++)
    {
      allReturnVector.push_back(input[i][j]);
    }
    double meanRet = mean(allReturnVector);
    double stdev = standardDeviation(allReturnVector, meanRet);
};

double Company::getCompanyMeanRet()
{
 return meanRet;
};


Portfolio::Portfolio(){ };

Portfolio::Portfolio(std::vector< std::vector<double> > returnVector, std::vector<double> vectorOfCompanyMeanRet, int noOfCompany, int numberOfDays)
{
    std::vector<double> eVector;
    std::vector<double> negativeRet;
    for (int i = 0; i < noOfCompany; i++)
    {
      eVector.push_back(-1);
      negativeRet.push_back(vectorOfCompanyMeanRet[i] * -1);
    }

    std::vector < std::vector<double> > covarianceMet = getCovariance(returnVector, noOfCompany, numberOfDays);

    for (int i = 0; i < 83; i++ )
    {
      covarianceMet[i].push_back(negativeRet[i]);
      covarianceMet[i].push_back(eVector[i]);
    }

    std::vector<double> negativeRetWithZero = negativeRet;

    negativeRetWithZero.push_back(0);
    negativeRetWithZero.push_back(0);

    std::vector<double> eVectorWithZero = eVector;
    eVectorWithZero.push_back(0);
    eVectorWithZero.push_back(0);
    
    std::vector< std::vector<double> > Q;
    Q = covarianceMet;
    Q.push_back(negativeRetWithZero);
    Q.push_back(eVectorWithZero);

    std::vector<std::vector<double> >portfolioWeight = getWeights(Q, noOfCompany);
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
    int columnLengthA = A.size(); //read matrix size horizontally
    int columnLengthB = B.size();
    int rowLengthB = B[0].size(); //read matrix size vertically
    int rowLengthA = A[0].size();
    std::vector< std::vector <double> > multiple;
    std::vector<double> zeros;
    double sum;
        // Initializing elements of matrix mult to 0.

        for (int i = 0; i < rowLengthA; i++)
        {
            zeros.push_back(0);
        }
        for(int j = 0; j < columnLengthB; j++)
        {
            multiple.push_back(zeros);
        }

    //Multiplying matrix a and b and storing in array mult.
        for(int i = 0; i < rowLengthA; i++)
        {
            for(int j = 0; j < columnLengthB; j++)
            {
                for(int k = 0; k < columnLengthA; k++)
                {
                    sum += A[k][i] * B[j][k];
                }

                multiple[j][i] = sum;
            }
        }
    return multiple;
}

std::vector<std::vector<double> > Minus(std::vector<std::vector<double> > A, std::vector<std::vector<double> > B)
{
    int columnLengthA = A.size(); //read matrix size horizontally
    int columnLengthB = B.size();
    int rowLengthB = B[0].size(); //read matrix size vertically
    int rowLengthA = A[0].size();
    std::vector<std::vector<double> > sum;
    std::vector<double> zeros;

        for(int j = 0; j < columnLengthA; ++j)
        {
        for (int i = 0; i < rowLengthA; i++)
        {
            zeros.push_back(0);
        }
            sum.push_back(zeros);
        }

        for(int i = 0; i < rowLengthA; ++i)
        {
            for(int j = 0; j < columnLengthB; ++j)
            {
                sum[i][j] += A[i][j] - B[i][j];    
            }
        }
    return sum;
}

std::vector<std::vector<double> > Plus(std::vector<std::vector<double> > A, std::vector<std::vector<double> > B)
{
    int columnLengthA = A.size(); //read matrix size horizontally
    int columnLengthB = B.size();
    int rowLengthB = B[0].size(); //read matrix size vertically
    int rowLengthA = A[0].size();
    std::vector<std::vector<double> > sum;
    std::vector<double> zeros;

        for (int i = 0; i < rowLengthA; i++)
        {
          zeros.push_back(0);
        }
        for(int j = 0; j < columnLengthA; ++j)
        {
            sum.push_back(zeros);
        }
        for(int i = 0; i < rowLengthA; ++i)
        {
            for(int j = 0; j < columnLengthB; ++j)
            {

                sum[i][j] += A[i][j] + B[i][j];   
            }
        }
    return sum;
}


double ATransposeA(std::vector<double> A)
{
    double sum = 0;
    for (int i = 0; i < A.size(); i++)
    {
        sum = sum + A[i]*A[i];
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
          cov[i][k] += (firstCompany[j] - firstCompanyMean) * (secondCompany[j] - secondCompanyMean) / (numberOfCompany - 1);    
        }
      }
    }
    return cov;
}


std::vector<std::vector <double> > getWeights(std::vector< std::vector<double> > Q, double numberOfCompany)
    {
        std::vector<std::vector<double> > s;
        std::vector<std::vector<double> > b;
        std::vector<std::vector<double> > x;
        std::vector<std::vector<double> > Qx;
        std::vector<double> zeros;
        std::vector<double> xVector;

        for (int i = 0; i < 83; i++)
        {
            zeros.push_back(0);
            xVector.push_back(1/83);

        }
        zeros.push_back(0); //portfolio return
        zeros.push_back(-1);

        xVector.push_back(0); // adding lagrangian multipliers into the vector
        xVector.push_back(0);

        b.push_back(zeros);
        s.push_back(zeros);
        x.push_back(xVector);
        Qx.push_back(zeros);

        double alpha;
        double sTs;
        double beta;
        std::vector<std::vector<double> > p;
        std::vector<std::vector<double> > pT;

        // // //initialise 
        int k = 0;
        int i = 0;
        Qx = Multiplication(Q,x);
        
        // s = Minus(b,Qx);
        // p = s;
        // while (sTs < 0.000006)
        // {   
        //     if (k == 0)
        //     {
        //         pT.push_back(p);
        //         alpha = ATransposeA(s) / Multiplication(p, Multiplication(Q,pT));
        //         x = Plus(x, Multiplication(alpha,pT));
        //         s = Minus(s, Multiplication(alpha, Multiplication(Q,pT)));
        //         sVector.push_back(s);
        //         pVector.push_back(p);c
        //     }

        //     else
        //     {
        //         pT.push_back(p);
        //         alpha = ATransposeA(s) / Multiplication(p, Multiplication(Q,pT));
        //         x = Plus(x, Multiplication(alpha,pT));
        //         s = Minus(s, Multiplication(alpha, Multiplication(Q,pT)));
        //         sVector.push_back(s);
        //         beta = ATransposeA(sVector);
        //         p = Plus(s,beta * p);
        //         pVector.push_back(p);

        //     }
        //     sTs = ATransposeA(sVector);
        //     k += 1;
        // }

        // while (i <= 0.1)
        //     {

        //         b[0][83] = i * -1;

        //         i += 0.05;
        //     }

    return b;
    }

