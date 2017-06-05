#include <cmath>
#include <iostream>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

// git add --all
// git commit -m "changed"
// git push
//g++ -o portfolio mainwo.cpp portfolio.cp csv.cp matrixOperations.cpp statisticalOperations.cpp


#include "portfolio.h"
#include "statisticalOperations.h"
#include "csv.h"

using namespace std;

double string_to_double( const string& s );
void readData(double **data,string fileName);
vector<vector<vector<double> > > inSampleRollingWindow (int inSampleRollingWindowSize, int outOfSampleRollingWindowSize, int numberOfAssets, int numberOfDays, vector<vector<double> > returnVector);
vector<vector<vector<double> > > outOfSampleRollingWindow (int inSampleRollingWindowSize, int outOfSampleRollingWindowSize, int numberOfAssets, int numberOfDays, vector<vector<double> > returnVector);


int main()
{
    // declaring all the variables that I need to use
  int numberOfAssets = 83;
  int numberOfDays = 700;
  int inSampleRollingWindowSize = 100;
  int outOfSampleRollingWindowSize = 12;
  int numberOfRollingWindows = 50;
  int numberOfPortfolioReturns = 21;

  double **returnMatrix = new double*[numberOfAssets]; //matrix to store the return data by allocating memroy for return data
  for (int i =0; i< numberOfAssets; i++)
    returnMatrix[i] = new double[numberOfDays];
  string fileName = "asset_returns.csv";
  readData(returnMatrix,fileName);
  
  vector< vector<double> > returnVector (numberOfAssets, vector<double>(numberOfDays));

    // transforming array to vector
  for (int i = 0; i < numberOfAssets; i++)
  {
     for (int j = 0; j < numberOfDays; j++)
    {
     returnVector[i][j] = returnMatrix[i][j];
   }
   }
    // construction of 3D vectors for in sample rolling windows and out of sample rolling windows
  vector<vector<vector<double> > > inSampleReturn = inSampleRollingWindow (inSampleRollingWindowSize, outOfSampleRollingWindowSize, numberOfAssets, numberOfDays, returnVector);
  vector<vector<vector<double> > > outOfSampleReturn = outOfSampleRollingWindow (inSampleRollingWindowSize, outOfSampleRollingWindowSize, numberOfAssets, numberOfDays, returnVector);
    // creating 2D vector for each company mean return over 50 rolling windows
    // 1D vector for pushing into 2D vector
  vector<double> VectorOfcompanyMeanRet (numberOfAssets);
  vector<vector<double> > matrixOfCompanyMeanReturn;

  for  (int j = 0; j < numberOfRollingWindows; j++)
  {
    for (int i = 0 ; i < numberOfAssets; i++)
    {
        VectorOfcompanyMeanRet[i] = StatisticalOperations::mean(inSampleReturn[j][i]);
    }
    matrixOfCompanyMeanReturn.push_back(VectorOfcompanyMeanRet);
  }
    // creating 2D vectors for saving all portfolios out of sample portfolio returns and out of sample portfolio variance
  vector<vector<double> > oosAverageReturn (numberOfPortfolioReturns, vector<double>(numberOfRollingWindows));
  vector<vector<double> > oosCovariance (numberOfPortfolioReturns, vector<double>(numberOfRollingWindows));
    // setting target return to be 0 to initialise different portfolio's target returns (in total 21 portfolios)
  double targetReturn = 0.0;
    //loop through number of portfolio returns, which is 21
  for (int j = 0; j < numberOfPortfolioReturns; j++)
  { //loop through number of portfolio rolling windows, which is 50
    for (int i = 0; i < numberOfRollingWindows; i++)
    {
        //constructing different portfolios
      Portfolio portfolio(inSampleReturn[i], matrixOfCompanyMeanReturn[i], numberOfAssets, inSampleRollingWindowSize, numberOfDays, outOfSampleRollingWindowSize, targetReturn, outOfSampleReturn[i]);
        // getting different return and variance from different portfolios
      oosAverageReturn[j][i] = portfolio.getPortfolioAverageReturn();
      oosCovariance[j][i] = portfolio.getPortfolioCovariance();
    }
    targetReturn += 0.005; //increment by 0.5% each time to create 21 portfolios
  }

    // output of csv files (these two files contain information for portfolio return and portfolio variance
    
  ofstream myfile2;
    myfile2.open ("oosAverageReturn.csv");
    for (int j = 0 ; j < numberOfPortfolioReturns; j++)
        {for(int i = 0 ; i < numberOfRollingWindows; i++)
            {myfile2 << oosAverageReturn[j][i] << ",";}
            myfile2 << "\n";
        }
    myfile2.close();

  ofstream myfile1;
    myfile1.open ("oosCovariance.csv");
    for (int j = 0 ; j < numberOfPortfolioReturns; j++)
        {for(int i = 0 ; i < numberOfRollingWindows; i++)
            {myfile1 << oosCovariance[j][i] << ",";}
            myfile1 << "\n";
        }
    myfile1.close();


}
// these are codes, which are provided by lecturer
double string_to_double( const string& s )
{
	istringstream i(s);
	double x;
	if (!(i >> x))
		return 0;
	return x;
} 

void readData(double **data,string fileName)
{
	char tmp[20];
	ifstream file (strcpy(tmp, fileName.c_str()));
	Csv csv(file);
	string line;
	if (file.is_open()) 
	{
		int i=0;
		while (csv.getline(line) != 0) {
         	for (int j = 0; j < csv.getnfield(); j++)
            {
               double temp=string_to_double(csv.getfield(j));
               //cout << "Asset " << j << ", Return "<<i<<"="<< temp<<"\n";
               data[j][i]=temp;
            }
            i++;
		}		
		file.close();
	}
	else 
	{
		cout <<fileName <<" missing\n";exit(0);
	}
}

// in sample rolling window function (3D matrix)
vector<vector<vector<double> > > inSampleRollingWindow (int inSampleRollingWindowSize, int outOfSampleRollingWindowSize, int numberOfAssets, int numberOfDays, vector<vector<double> > returnVector)
{
  vector<vector<vector<double> > > tempBacktest;
  //(50, vector<vector<double> >(numberOfAssets, vector<double>(inSampleRollingWindowSize)));
  vector<vector<double> > tempReturnVector (numberOfAssets, vector<double> (inSampleRollingWindowSize));
    for (int j = 0; j < numberOfDays - inSampleRollingWindowSize; j += 12)
    {
      for (int k = 0; k < numberOfAssets; k++)
      {
        for (int i = 0; i < 100; i++)
        {   
            tempReturnVector[k][i] = returnVector[k][(i+j)];
        }
      }
    tempBacktest.push_back(tempReturnVector);
    }
  return tempBacktest;
}

// out of sample rolling window function (3D matrix)
vector<vector<vector<double> > > outOfSampleRollingWindow (int inSampleRollingWindowSize, int outOfSampleRollingWindowSize, int numberOfAssets, int numberOfDays, vector<vector<double> > returnVector)
{
  vector<vector<vector<double> > > tempBacktest;
  vector<vector<double> > tempReturnVector (numberOfAssets, vector<double> (outOfSampleRollingWindowSize));
  for (int j = 100; j < numberOfDays; j += 12)
  {
    for (int k = 0; k < numberOfAssets; k++)
    {
      for (int i = 0; i < outOfSampleRollingWindowSize; i++)
      {   
        tempReturnVector[k][i] = returnVector[k][(i+j)];
      }
    }
    tempBacktest.push_back(tempReturnVector);
  }

  return tempBacktest;
}









