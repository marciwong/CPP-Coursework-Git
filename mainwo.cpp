#include <cmath>
#include <iostream>
#include <string.h>
#include <stdlib.h>
#include <armadillo>
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

int main()
{
  int numberOfAssets = 83;
  int numberOfDays = 700;
  double **returnMatrix = new double*[numberOfAssets]; //matrix to store the return data by allocating memroy for return data
  for (int i =0; i< numberOfAssets; i++)
    returnMatrix[i] = new double[numberOfDays];
  string fileName = "asset_returns.csv";
  readData(returnMatrix,fileName);
  
  std::vector<vector<double> > returnVector;//(83,vector<double>(83));
  std::vector<double> zeroVector;
  // std::vector<Company> companyVector;

  int inSampleRollingWindowSize = 100;
  int outOfSampleRollingWindowSize = 12;

  for (int i = 0; i < numberOfDays; i++)
  {
    zeroVector.push_back(0);
  }

  for (int i = 0; i < numberOfAssets; i++)
  {
    returnVector.push_back(zeroVector);;
  }

  for (int i = 0; i < numberOfAssets; i++)
   {
     for (int j = 0; j < numberOfDays; j++)
   {
     returnVector[i][j] = returnMatrix[i][j]; //transforming 2d array into 2d vector
   }
   }

  std::vector<std::vector<std::vector<double> > > inSampleMat = inSampleRollingWindow(inSampleRollingWindowSize, outOfSampleRollingWindowSize, numberOfAssets, numberOfDays, returnVector);
  std::vector<double> VectorOfcompanyMeanRet;
  std::vector<std::vector<double> > MatrixOfcompanyMeanRet;
  std::vector<Company> vectorOfCompany;
  std::vector<std::vector<Company> > rollingCompanyMat;

  for  (int j = 0; j < ((numberOfDays- inSampleRollingWindowSize)/outOfSampleRollingWindowSize); j++)
  {
    for (int i = 0 ; i < numberOfAssets; i++)
    {
  		  Company company(inSampleMat[j],i, inSampleRollingWindowSize);
  		  VectorOfcompanyMeanRet.push_back(company.getCompanyMeanRet());
    }
    MatrixOfcompanyMeanRet.push_back(VectorOfcompanyMeanRet);
  }

  std::vector<std::vector<std::vector<double> > > portfoliosWeightsMatrix;
  std::vector<std::vector<double> > portfolioTwoDWeights;
  double noOfTargetReturn = 0.0;
  for (int i = 0; i < ((numberOfDays- inSampleRollingWindowSize)/outOfSampleRollingWindowSize); i++)
  {
    for (int noOfTargetReturn = 0; noOfTargetReturn < 0.1000001; noOfTargetReturn +=0.005)
    {
      Portfolio portfolio(std::vector< std::vector<double> > inSampleMat, std::vector<double> vectorOfCompanyMeanRet, int noOfCompany, int inSampleRollingWindowSize, int numberOfDays, int outOfSampleRollingWindowSize, double noOfTargetReturn);  
      portfolioTwoDWeights.push_back(portfolio.getPortfolioWeights());
    }
    portfoliosWeightsMatrix.push_back(portfolioTwoDWeights);
  }

  // cout << portfoliosWeightsMatrix[0].size() << endl;
  // cout << portfoliosWeightsMatrix.size() << endl;
  // cout << portfoliosWeightsMatrix[0][0].size() << endl;

//   for (int i = 0; i < 50; i++)
//   {
//     for (int j = 0; j < 20; j++ )
//     {
//       for (int h = 0; h <83; h++)
//       {
//         cout << portfoliosWeightsMatrix[i][j][h];
//       }
//     }
//   }
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







