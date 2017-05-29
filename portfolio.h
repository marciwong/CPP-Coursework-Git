#ifndef Portfolio_h
#define Portfolio_h
#include <armadillo>
#include <cmath>
#include <iostream>
#include <string.h>
#include <stdlib.h>
#include <armadillo>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

using namespace std;

double mean(vector<double> input);
double meanArray(double input[]);
double standardDeviation(vector<double> input , double mean);
std::vector< std::vector<double> > getCovariance( std::vector< std::vector<double> > returnVector, int size, int timeLength);
std::vector<double> getWeights(std::vector< std::vector<double> > Q, double numberOfCompany, double noOfTargetReturn);
std::vector<std::vector <double> > Multiplication(std::vector<std::vector<double> > A, std::vector<std::vector<double> > B);
std::vector<std::vector<double> > Minus(std::vector<std::vector<double> > A, std::vector<std::vector<double> > B);
std::vector<std::vector<double> > scalarMultiplication(double alpha, std::vector<std::vector<double> > A);
std::vector<std::vector<double> > Plus (std::vector<std::vector<double> > A, std::vector<std::vector<double> > B);
std::vector<std::vector<double> > transpose(std::vector<std::vector<double> > A);
std::vector<double> Minus1D(std::vector<double> A, std::vector<double> B);
std::vector<double> Plus1D(std::vector<double> A, std::vector<double> B);

class Portfolio
{
	private:
		std::vector< std::vector<double> > inSampleCovariance;
		std::vector< std::vector<double> > outOfSampleCovariance;
		std::vector< std::vector<double> > outOfSampleAverageReturn;
		double actualAverageReturn;
		std::vector<std::vector<double > > portfolioWeight;
		std::vector< std::vector<double> > Q;
		double portfolioCovariance;

	public:
		Portfolio(std::vector< std::vector<double> > inSampleMat, std::vector< double > matrixOfCompanyMeanRet, int noOfCompany, int inSampleRollingWindowSize, int numberOfDays, int outOfSampleRollingWindowSize, double noOfTargetReturn, std::vector<std::vector<double> > outOfSampleReturn);
		std::vector<std::vector<double> > getPortfolioWeights();
		std::vector<std::vector<double> > getPortfolioInSampleCovariance();
		std::vector<std::vector<double> > getPortfolioOutOfSampleCovariance();
		std::vector<std::vector<double> > getQ();
		double getPortfolioCovariance();
		double getPortfolioAverageReturn();
};

class Company
{
	private:
		int Days;
		double meanRet;
		double stdev;
		std::vector<double> allReturnVector;
		// double weight;   //**crucial add it in later


	public:
		Company();
		Company(std::vector< vector<double> > input, int i, int days );
		double getCompanyMeanRet();
};

#endif
