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

class Portfolio
{
	private:
		int N;
		std::vector<double> ret;
		std::vector< std::vector<double> > covariance;
		double lambda;
		double miu;
		std::vector<double> weight;
		double portRet;

	public:
		Portfolio();
		Portfolio(std::vector< std::vector<double> > returnVector, std::vector<double> vectorOfCompanyMeanRet, int noOfCompany, int time);
		std::vector<double> getPortfolioWeights;
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
