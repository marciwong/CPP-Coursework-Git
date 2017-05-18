#include <cmath>
#include <iostream>
#include <vector>

using namespace std;

std::vector<std::vector <double> > Multiplication2D (std::vector<std::vector<double> > A, std::vector<std::vector<double> > B);

int main ()
{   
    // double n = 100;
    std::vector<std::vector<double> > A;
    std::vector<double> B;
    B.push_back(10);
    B.push_back(11);
    A.push_back(B);
    A.push_back(B);
    A.push_back(B);
    A.push_back(B);
    // std::vector<double> C;
    // C.push_back(10);
    // C.push_back(11);
    // A.push_back(B);
    // A.push_back(C);
    // A[0].push_back(110000);
    // std::vector<std::vector<double> > D;
    // D = A;
    // std::vector<std::vector<double> > E;
    // E = Multiplication2D(A,D);
    // std::vector<std::vector <double> > B;

    // A[0].push_back(2);
    // A[1].push_back(2);    
    // B[0][0] = 5;
    // B[0][1] = 4;
    // B[1][0] = 3;
    // B[1][1] = 2;


    // std::vector<std::vector <double> > D = Multiplication2D(A,B);  
    double columnLengthB = A.size();

    cout << As.size() << endl;  


};

std::vector<std::vector <double> > Multiplication2D (std::vector<std::vector<double> > A, std::vector<std::vector<double> > B)
 {   
    double columnLengthA = A.size(); //read matrix size horizontally
    double columnLengthB = B.size();
    double rowLengthB = B[0].size(); //read matrix size vertically
    double rowLengthA = A[0].size();
    std::vector< std::vector <double> > multiple;
    
    // if (columnLengthA == rowLengthB)
    // {
        // Initializing elements of matrix mult to 0.
    for(double i = 0; i < rowLengthA; ++i)
        for(double j = 0; j < columnLengthB; ++j)
        {
            multiple[i][j]=0;
        }

    // Multiplying matrix a and b and storing in array mult.
    for(double i = 0; i < rowLengthA; ++i)
        for(double j = 0; j < columnLengthB; ++j)
            for(double k = 0; k < columnLengthA; ++k)
            {
                multiple[i][j] += A[i][k] * B[k][j];
            }
    // }
  return multiple;
  }
