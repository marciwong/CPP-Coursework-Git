std::vector<std::vector <double> > getPortRetAndStdMat(std::vector< std::vector<double> > Q)
    {
        std::vector<double> b;
        for (int i = 0; i < 83; i++)
        {
            b.push_back(0);
        }
        b.push_back(0); //portfolio return
        b.push_back(-1);

        std::vector<double> s;

        std::vector<double> x;
        for (int i = 0; i < 85; i++)
        {
            x.push_back(0);
        }

        std::vector<double> p;
        std::vector<vector<double> > pT;
        std::vector<vector<double> > pVector;
        std::vector<double> Qx;
        std::vector<double> beta;
        std::vector<double> Qx;
        long sTs;
        std::vector<vector<double> > sVector;

        int i = 0;
        b.push_back(0);
        b.push_back(0.1);
        b.push_back(-1);

        //initialise 
        int k = 0;
        int i = 0;
        Qx = Multiplication(Q,x);
        s = Minus(b,Qx);
        p = s;
       while (sTs < 0.000006)
        {   
            if (k == 0)
            {
                pT.push_back(p);
                a = ATransposeA(s) / Multiplication(p, Multiplication(Q,pT));
                x = Plus(x, Multiplication(a,pT));
                s = Minus(s, Multiplication(a, Multiplication(Q,pT)));
                sVector.push_back(s);
                pVector.push_back(p);
            }

            else
            {
                pT.push_back(p);
                a = ATransposeA(s) / Multiplication(p, Multiplication(Q,pT));
                x = Plus(x, Multiplication(a,pT));
                s = Minus(s, Multiplication(a, Multiplication(Q,pT)));
                sVector.push_back(s);
                beta = ATransposeA(sVector(k));
                p = Plus(s,beta * p);
                pVector.push_back(p);

            }
            sTs = ATransposeA(sVector(k));
            k += 1;
        }

        while (i <= 0.1)
            {

                b[1] = i * -1;

                i += 0.05;
            }

    }

std::vector<std::vector <double> > Multiplication (std::vector<std::vector<double> > A, std::vector<std::vector<double> > B)
{    
    int columnLengthA = A.size(); //read matrix size horizontally
    int columnLengthB = B.size();
    int rowLengthB = B[0].size(); //read matrix size vertically
    int rowLengthA = A[0].size();
    std::vector< std::vector <double> > multiple;
    std::vector<double> zeros;
        // Initializing elements of matrix mult to 0.
    if (columnLengthA == rowLengthB)
    {
        for (i = 0; i < rowLengthA; i++)
        {
            zeros.push_back(0);
        }
        for(j = 0; j < columnLengthB; ++j)
        {
            multiple.push_back(zeros);
        }

    // Multiplying matrix a and b and storing in array mult.
        for(i = 0; i < rowLengthA; ++i)
        {
            for(j = 0; j < columnLengthB; ++j)
            {
                for(k = 0; k < columnLengthA; ++k)
                {
                multiple[i][j] += A[i][k] * B[k][j];
                }
            }
        }
        return multiple;
    }
    else
    {
        return 0;
    }

}

std::vector<std::vector<double> > Minus (std::vector<std::vector<double> > A, std::vector<std::vector<double> > B)
{
    int columnLengthA = A.size(); //read matrix size horizontally
    int columnLengthB = B.size();
    int rowLengthB = B[0].size(); //read matrix size vertically
    int rowLengthA = A[0].size();
    std::vector<std::vector<double> > sum;

    if (columnLengthA = columnLengthB) && (rowLengthB = rowLengthA)
    {
        for (i = 0; i < rowLengthA; i++)
        {
            zeros.push_back(0);
        }
        for(j = 0; j < columnLengthA; ++j)
        {
            sum.push_back(zeros);
        }
        for(i = 0; i < rowLengthA; ++i)
        {
            for(j = 0; j < columnLengthB; ++j)
            {

                sum[i][j] += A[i][j] - B[i][j];
                
            }
        }
    return sum;
    }

    else
    {
        return 0;
    }
}

std::vector<std::vector<double> > Plus (std::vector<std::vector<double> > A, std::vector<std::vector<double> > B)
{
    int columnLengthA = A.size(); //read matrix size horizontally
    int columnLengthB = B.size();
    int rowLengthB = B[0].size(); //read matrix size vertically
    int rowLengthA = A[0].size();
    std::vector<std::vector<double> > sum;

    if (columnLengthA = columnLengthB) && (rowLengthB = rowLengthA)
    {
        for (i = 0; i < rowLengthA; i++)
        {
            zeros.push_back(0);
        }
        for(j = 0; j < columnLengthA; ++j)
        {
            sum.push_back(zeros);
        }
        for(i = 0; i < rowLengthA; ++i)
        {
            for(j = 0; j < columnLengthB; ++j)
            {

                sum[i][j] += A[i][j] + B[i][j];   
            }
        }
    return sum;
    }

    else
    {
        return 0;
    }
}

double ATransposeA(std::vector<std::vector<double> > A)
{
    double sum;
    for (int i = 0; i < A.size(); i++)
    {
        sum += A[i] * A[i];
    }
}