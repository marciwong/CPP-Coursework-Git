#ifndef MatrixOperations_h
#define MatrixOperations_h

using namespace std;

class MatrixOperations
{
	public:
		static vector<double> plus(vector<double> matrix1, vector<double> matrix2);
		static vector< vector<double> > plus(vector< vector<double> > matrix1, vector< vector<double> > matrix2);
		static vector<double> minus(vector<double> matrix1, vector<double> matrix2);
		static vector< vector<double> > minus(vector< vector<double> > matrix1, vector< vector<double> > matrix2);
		static vector< vector<double> > multiple(vector< vector<double> > matrix1, vector< vector<double> > matrix2);
		static vector< vector<double> > scalarMultiple(double scalar, vector< vector<double> > matrix);
		static vector< vector<double> > transpose(vector< vector<double> > matrix);
};

#endif
