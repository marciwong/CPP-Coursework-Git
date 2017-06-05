#include <cmath>
#include <vector>

#include "matrixOperations.h"

using namespace std;
	
vector< vector<double> > MatrixOperations::plus(vector< vector<double> > matrix1, vector< vector<double> > matrix2)
{
	int width = matrix1.size();
	int height = matrix1[0].size();
	vector<vector<double> > result(width, vector<double>(height));
    for (int i = 0; i < width; i++)
    {
        for (int j = 0; j < height; j++)
        {
            result[i][j] = matrix1[i][j] + matrix2[i][j];
        }
    }
    return result;
}
	
vector< vector<double> > MatrixOperations::minus(vector< vector<double> > matrix1, vector< vector<double> > matrix2)
{
	int width = matrix1.size();
	int height = matrix1[0].size();
	vector<vector<double> > result(width, vector<double>(height));
    for (int i = 0; i < width; i++)
    {
        for (int j = 0; j < height; j++)
        {
            result[i][j] = matrix1[i][j] - matrix2[i][j];
        }
    }
    return result;
}

// this is a function to multiply two matrices
vector< vector<double> > MatrixOperations::multiple(vector< vector<double> > matrix1, vector< vector<double> > matrix2)
{
	vector< vector <double> > result(matrix2.size(), vector<double>(matrix1[0].size()));

    // Multiplying matrix a and b and storing in array mult.
    for (int i = 0; i < matrix1[0].size(); i++)
    {
        for (int j = 0; j < matrix2.size(); j++)
        {
            for (int k = 0; k < matrix1.size(); k++)
            {
                result[j][i] += matrix1[k][i] * matrix2[j][k];
            }
        }
    }
    return result;
}

// this is a function to multiply a scalar with a matrix
vector< vector<double> > MatrixOperations::scalarMultiple(double scalar, vector< vector<double> > matrix)
{
	int width = matrix.size();
	int height = matrix[0].size();
	vector<vector<double> > result(width, vector<double>(height));
    for (int i = 0; i < width; i++) 
    {
        for (int j = 0; j < height; j++)
        {
            result[i][j] = scalar * matrix[i][j];
        }
    }
    return result;
}
	
vector< vector<double> > MatrixOperations::transpose(vector< vector<double> > matrix)
{
	int width = matrix.size();
	int height = matrix[0].size();
	vector<vector<double> > result(height, vector<double>(width));
    for (int i = 0; i < width; i++)
    {
        for (int j = 0; j < height; j++)
        {
            result[j][i] = matrix[i][j];
        }
    }
    return result;
}


