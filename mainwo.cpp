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
//g++ -o portfolio mainwo.cpp portfolio.cpp csv.cp


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
  int numberOfRollingWindows = (numberOfDays- inSampleRollingWindowSize)/outOfSampleRollingWindowSize;
  int numberOfPortfolioReturns = 20; 

  double **returnMatrix = new double*[numberOfAssets]; //matrix to store the return data by allocating memroy for return data
  for (int i =0; i< numberOfAssets; i++)
    returnMatrix[i] = new double[numberOfDays];
  string fileName = "asset_returns.csv";
  readData(returnMatrix,fileName);
  
  std::vector<vector<double> > returnVector (numberOfAssets,vector<double>(numberOfDays));

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
  		  Company company(inSampleReturn[j],i, inSampleRollingWindowSize);
  		  VectorOfcompanyMeanRet[i] = (company.getCompanyMeanRet());
    }
    matrixOfCompanyMeanReturn.push_back(VectorOfcompanyMeanRet);
  }

  std::vector<std::vector<double> > oosAverageReturn (numberOfPortfolioReturns, std::vector<double>(numberOfRollingWindows));
  std::vector<std::vector<double> > oosCovariance (numberOfPortfolioReturns, std::vector<double>(numberOfRollingWindows));

  for (int j = 0; j < numberOfPortfolioReturns; j++)
  {
    for (int i = 0; i < numberOfRollingWindows; i++)
    {
      for (double targetReturn = 0.0; targetReturn < 0.100000; targetReturn +=0.005)
      {
        Portfolio portfolio(inSampleReturn[i], matrixOfCompanyMeanReturn[i], numberOfAssets, inSampleRollingWindowSize, numberOfDays, outOfSampleRollingWindowSize, targetReturn, outOfSampleReturn[i]);
        // portfolioMatrix[j][i] = portfolio;
      }
    }
  }

  // ofstream myfile;
  //   myfile.open ("portfolios.csv");
  //   for (int i = 0; i < 83; i++)
  //       {for(int j = 0;j < 83;j++)
  //           {myfile << portfolioMatrix[0][0].getPortfolioOutOfSampleCovariance()[i][j] << ",";}
  //           myfile << "\n";
  //       }
  //   myfile.close();
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
  std::vector<std::vector<double> > tempReturnVector;
  std::vector<double> hundredZeros;

  for (int j = 0; j < inSampleRollingWindowSize; j++)
  {
    hundredZeros.push_back(0);
  }
  for (int i = 0; i < numberOfAssets; i++)
  {
    tempReturnVector.push_back(hundredZeros);
  }
    for (int j = 0; j < numberOfDays - inSampleRollingWindowSize + 1; j += 12)
    {
      for (int h = 0; h < 50; h++)
      {
      for (int k =0; k < numberOfAssets; k++)
      {
        for (int i = 0; i < 100; i++)
        {   
            tempReturnVector[k][i] = returnVector[k][(i+j)];
        }
      }
    }
  tempBacktest.push_back(tempReturnVector);
  }
  return tempBacktest;
}

std::vector<std::vector<std::vector<double> > > outOfSampleRollingWindow (int inSampleRollingWindowSize, int outOfSampleRollingWindowSize, int numberOfAssets, int numberOfDays, std::vector<vector<double> > returnVector)
{
  std::vector<std::vector<std::vector<double> > > tempBacktest;
  std::vector<std::vector<double> > tempReturnVector;
  std::vector<double> zeros;

  for (int j = 0; j < outOfSampleRollingWindowSize; j++)
  {
    zeros.push_back(0);
  }
  for (int i = 0; i < numberOfAssets; i++)
  {
    tempReturnVector.push_back(zeros);
  }

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









