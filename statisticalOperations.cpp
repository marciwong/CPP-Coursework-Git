#include <cmath>
#include <vector>

#include "statisticalOperations.h"
#include "matrixOperations.h"

using namespace std;
	
double StatisticalOperations::mean(vector<double> input)
{
    double sum = 0.0;
    for (int i = 0; i < input.size(); i++)
    {
        sum += input[i];
    }  
    return (sum / input.size());
}

double StatisticalOperations::standardDeviation(vector<double> input , double mean)
{
    double sumSQ = 0.0;
    for (int i = 0; i < input.size(); i++)
    {
        sumSQ += (input[i] - mean) * (input[i] - mean);
    }    
    return (sqrt(sumSQ / (input.size() - 1 )));
}

vector< vector<double> > StatisticalOperations::getCovariance(vector< vector<double> > returnVector, int numberOfCompany, int timeLength)
{
    vector< vector<double> > cov(numberOfCompany, vector<double>(numberOfCompany));

    vector<double> firstCompany(timeLength);
    vector<double> secondCompany(timeLength);

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

// this includes conjugate gradient method

vector<double> StatisticalOperations::getWeights(vector< vector<double> > Q, double numberOfCompany, double noOfTargetReturn)
{
    double tolerence = 0.000001;

    // Set up x
    vector< vector<double> > x(1, vector<double>(numberOfCompany + 2));
    for (int i = 0; i < numberOfCompany; i++) 
    {
        x[0][i] = 1.0 / numberOfCompany;
    }
    x[0][numberOfCompany] = 1.0; // lambda
    x[0][numberOfCompany + 1] = 1.0; // mu

    // Set up b
    vector< vector<double> > b(1, vector<double>(numberOfCompany + 2));
    for (int i = 0; i < numberOfCompany; i++) 
    {
        b[0][i] = 0.0;
    }
    b[0][numberOfCompany] = -1.0 * noOfTargetReturn; // -r_p
    b[0][numberOfCompany + 1] = -1.0;

    vector< vector<double> > s = MatrixOperations::minus(b, MatrixOperations::multiple(Q, x));
    vector< vector<double> > p(s);

    double sTs = MatrixOperations::multiple(MatrixOperations::transpose(s), s)[0][0];
    while (sTs > tolerence) 
    {
        double alpha = sTs / (MatrixOperations::multiple(MatrixOperations::multiple(MatrixOperations::transpose(p), Q), p)[0][0]);

        x = MatrixOperations::plus(x, (MatrixOperations::scalarMultiple(alpha, p)));

        vector< vector<double> > s_plus1 = MatrixOperations::minus(s, (MatrixOperations::scalarMultiple(alpha, (MatrixOperations::multiple(Q, p)))));

        sTs = MatrixOperations::multiple(MatrixOperations::transpose(s_plus1), s_plus1)[0][0];

        double beta = (sTs) / (MatrixOperations::multiple(MatrixOperations::transpose(s), s)[0][0]);

        p = MatrixOperations::plus(s_plus1, (MatrixOperations::scalarMultiple(beta, p)));

        s = s_plus1;
    }

    vector<double> weights (numberOfCompany);
    for (int i = 0; i < weights.size(); i++)
    {
        weights[i] = x[0][i];
    }

    return weights;
}


