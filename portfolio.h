#ifndef Portfolio_h
#define Portfolio_h

#include <cmath>
#include <iostream>
#include <string.h>
#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

using namespace std;

class Portfolio
{
	private:
		vector< vector<double> > inSampleCovariance;
		vector< vector<double> > outOfSampleCovariance;
		vector< vector<double> > outOfSampleAverageReturn;
		double actualAverageReturn;
		vector<vector<double > > portfolioWeight;
		vector< vector<double> > Q;
		double portfolioCovariance;

	public:
		Portfolio(vector< vector<double> > inSampleReturn, vector< double > matrixOfCompanyMeanRet, int noOfCompany, int inSampleRollingWindowSize, int numberOfDays, int outOfSampleRollingWindowSize, double noOfTargetReturn, vector<vector<double> > outOfSampleReturn);
		vector<vector<double> > getPortfolioWeights();
		vector<vector<double> > getPortfolioInSampleCovariance();
		vector<vector<double> > getPortfolioOutOfSampleCovariance();
		vector<vector<double> > getQ();
		double getPortfolioCovariance();
		double getPortfolioAverageReturn();
};
#endif
