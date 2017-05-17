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
//g++ -o portfolio main.cpp portfolio.cpp csv.cp


#include "portfolio.h"
#include "csv.h"

using namespace std;

double string_to_double( const std::string& s );
void readData(double **data,string fileName);

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

  for (int i = 0; i < 83; i++)
  	{
  		for (int j = 0; j < 83; j++)
  	{
  		returnVector[i][j] = returnMatrix[i][j]; //transforming 2x2 array into 2x2 vector
  	}

  	}
  std::vector<Company> companyVector;

  // Company.getCompanyReturn(companyVector)
  std::vector<Company> vectorOfCompanyRet;
  for (int i = 0 ; i < 83; i++)
  {
  		Company company(returnVector,i, numberOfDays);
  		vectorOfCompanyRet.push_back(company);
  }

  std::vector<double> vectorOfCompanyMeanRet;
  for (int i = 0; i < 83; i++)
  {

      vectorOfCompanyMeanRet.push_back(vectorOfCompanyRet[i].getCompanyMeanRet());
  }

  //Portfolio fullSamplePort(vectorOfCompanyMeanRet, vectorOfCompanyMeanRet, numberOfAssets,numberOfDays);
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



