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
//g++ -o portfolio mainwo.cpp portfolio.cp csv.cp


#include "portfolio.h"
#include "csv.h"

using namespace std;

double string_to_double( const std::string& s );
void readData(double **data,string fileName);
std::vector<std::vector<std::vector<double> > > inSampleRollingWindow (int inSampleRollingWindowSize, int outOfSampleRollingWindowSize, int numberOfAssets, int numberOfDays, std::vector<vector<double> > returnVector);
std::vector<std::vector<std::vector<double> > > outOfSampleRollingWindow (int inSampleRollingWindowSize, int outOfSampleRollingWindowSize, int numberOfAssets, int numberOfDays, std::vector<vector<double> > returnVector);


int main()
{
  int numberOfAssets = 83;
  int numberOfDays = 700;
  int inSampleRollingWindowSize = 100;
  int outOfSampleRollingWindowSize = 12;
  int numberOfRollingWindows = 49;
  int numberOfPortfolioReturns = 20; 

  double **returnMatrix = new double*[numberOfAssets]; //matrix to store the return data by allocating memroy for return data
  for (int i =0; i< numberOfAssets; i++)
    returnMatrix[i] = new double[numberOfDays];
  string fileName = "asset_returns.csv";
  readData(returnMatrix,fileName);
  
  std::vector< std::vector<double> > returnVector (numberOfAssets, std::vector<double>(numberOfDays));

  for (int i = 0; i < numberOfAssets; i++)
  {
     for (int j = 0; j < numberOfDays; j++)
    {
     returnVector[i][j] = returnMatrix[i][j]; //transforming 2d array into 2d vector
   }
   }

  std::vector<std::vector<std::vector<double> > > inSampleReturn = inSampleRollingWindow (inSampleRollingWindowSize, outOfSampleRollingWindowSize, numberOfAssets, numberOfDays, returnVector);
  std::vector<std::vector<std::vector<double> > > outOfSampleReturn = outOfSampleRollingWindow (inSampleRollingWindowSize, outOfSampleRollingWindowSize, numberOfAssets, numberOfDays, returnVector);

  std::vector<double> VectorOfcompanyMeanRet (numberOfAssets);
  std::vector<std::vector<double> > matrixOfCompanyMeanReturn;

  for  (int j = 0; j < numberOfRollingWindows; j++)
  {
    for (int i = 0 ; i < numberOfAssets; i++)
    {
  		  Company company(inSampleReturn[j], i, inSampleRollingWindowSize);
  		  VectorOfcompanyMeanRet[i] = (company.getCompanyMeanRet());
    }
    matrixOfCompanyMeanReturn.push_back(VectorOfcompanyMeanRet);
  }

  std::vector<std::vector<double> > oosAverageReturn (numberOfPortfolioReturns, std::vector<double>(numberOfRollingWindows));
  std::vector<std::vector<double> > oosCovariance (numberOfPortfolioReturns, std::vector<double>(numberOfRollingWindows));
  std::vector<std::vector<Portfolio> > portfolioMat;
  std::vector<Portfolio> tempPortfolioVector;
  double targetReturn = 0.005;

  for (int j = 0; j < numberOfPortfolioReturns; j++)
  {
    for (int i = 0; i < numberOfRollingWindows; i++)
    {
      Portfolio portfolio(inSampleReturn[i], matrixOfCompanyMeanReturn[i], numberOfAssets, inSampleRollingWindowSize, numberOfDays, outOfSampleRollingWindowSize, targetReturn, outOfSampleReturn[i]);
      oosAverageReturn[j][i] = portfolio.getPortfolioAverageReturn();
      oosCovariance[j][i] = portfolio.getPortfolioCovariance();
    }
    targetReturn += 0.005;
  }

  ofstream myfile2;
    myfile2.open ("oosAverageReturn.csv");
    for (int j = 0 ; j < numberOfPortfolioReturns; j++)
        {for(int i = 0 ; i < numberOfRollingWindows; i++)
            {myfile2 << oosAverageReturn[j][i] << ",";}
            myfile2 << "\n";
        }
    myfile2.close();

  ofstream myfile;
    myfile.open ("oosReturn.csv");
    for (int j = 0 ; j < 83; j++)
        {for(int i = 0 ; i < 12; i++)
            {myfile << outOfSampleReturn[48][j][i] << ",";}
            myfile << "\n";
        }
    myfile.close();

  ofstream myfile1;
    myfile1.open ("oosCovariance.csv");
    for (int j = 0 ; j < numberOfPortfolioReturns; j++)
        {for(int i = 0 ; i < numberOfRollingWindows; i++)
            {myfile1 << oosCovariance[j][i] << ",";}
            myfile1 << "\n";
        }
    myfile1.close();


}



double string_to_double( const std::string& s )
{
	std::istringstream i(s);
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

std::vector<std::vector<std::vector<double> > > inSampleRollingWindow (int inSampleRollingWindowSize, int outOfSampleRollingWindowSize, int numberOfAssets, int numberOfDays, std::vector<vector<double> > returnVector)
{
  std::vector<std::vector<std::vector<double> > > tempBacktest;
  //(50, vector<vector<double> >(numberOfAssets, vector<double>(inSampleRollingWindowSize)));
  std::vector<std::vector<double> > tempReturnVector (numberOfAssets, std::vector<double> (inSampleRollingWindowSize));
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

std::vector<std::vector<std::vector<double> > > outOfSampleRollingWindow (int inSampleRollingWindowSize, int outOfSampleRollingWindowSize, int numberOfAssets, int numberOfDays, std::vector<vector<double> > returnVector)
{
  std::vector<std::vector<std::vector<double> > > tempBacktest;
  std::vector<std::vector<double> > tempReturnVector (numberOfAssets, vector<double> (outOfSampleRollingWindowSize));
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









