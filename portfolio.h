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

class timeSeriesPortfolio
{
	private:
		std::vector<std::vector<double> > weights;
		double portfolioRet;

	public:
		double getTimeSeriesPortfolioRet();
};

class Portfolio: public timeSeriesPortfolio
{
	private:
		// int N;
		// std::vector<double> ret;
		std::vector< std::vector<double> > covariance;
		// double lambda;
		// double miu;
		std::vector<double> portfolioWeight;
		std::vector< std::vector<double> > Q;
		// double portRet;

	public:
		Portfolio(std::vector< std::vector<double> > inSampleMat, std::vector< double > matrixOfCompanyMeanRet, int noOfCompany, int inSampleRollingWindowSize, int numberOfDays, int outOfSampleRollingWindowSize, double noOfTargetReturn);
		std::vector<double> getPortfolioWeights();
		std::vector<std::vector<double> > getPortfolioCovariance();
		std::vector<std::vector<double> > getQ();

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
