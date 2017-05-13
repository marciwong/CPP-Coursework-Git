#ifndef Portfolio_h
#define Portfolio_h

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

	public:
		void getCovariance();

};

class Company
{
	private:
		int Days;
		double meanRet;
		double stdev;
		std::vector<double> allReturnVector;
		// std::vector<double> inSampleVector;
		// std::vector<double> outSampleVector;

	public:
		Company();
		Company(std::vector< vector<double> > input, int i );
		// void getCompanyReturn(std::vector<vector<double> > input)
		// void getInSampleArray();
		// void getOutSampleArray();
};

#endif