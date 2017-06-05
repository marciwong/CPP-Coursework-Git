#ifndef StatisticalOperations_h
#define StatisticalOperations_h

using namespace std;

class StatisticalOperations
{
	public:
		static double mean(vector<double> input);
		static double meanArray(double input[]);
		static double standardDeviation(vector<double> input , double mean);
		static vector< vector<double> > getCovariance(vector< vector<double> > returnVector, int size, int timeLength);
		static vector<double> getWeights(vector< vector<double> > Q, double numberOfCompany, double noOfTargetReturn);
};

#endif
