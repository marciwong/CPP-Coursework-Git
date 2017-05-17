std::vector<std::vector <double> > getPortRetAndStdMat(std::vector< std::vector<double> > inverseMet, std::vector< std::vector<double> > covarianceMet, int noOfCompany, mat tempInv, std::vector<double> vectorOfCompanyMeanRet )
    {
        std::vector<std::vector <double> > b;
        std::vector<double> s0;
        std::vector<double> x0;
        std::vector<double> p0;

        int i = 0;
            while (i <= 0.1)
                {
                    b[0] = 0;
                    b[1] = i * -1; //loop through different portfolio target
                    b[2] = -1;
                    for (int i = 0; i < noOfCompany; i++)
                    {
                        s0[i] = b - inverseMet * x0[i]
                    }



                }

    }
