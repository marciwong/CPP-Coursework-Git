std::vector<std::vector <double> > getPortRetAndStdMat(std::vector< std::vector<double> > inverseMet, std::vector< std::vector<double> > covarianceMet, int noOfCompany, mat tempInv, std::vector<double> vectorOfCompanyMeanRet )
    {
        mat b = randu(3,1);
        mat s = randu(noOfCompany,1);
        mat x = randu(noOfCompany,1);
        mat p = randu(noOfCompany,1);
        mat a = randu(noOfCompany,1);
        mat sT = 

        int i = 0;
            while (i <= 0.1)
                {
                    b(1,1) = 0;
                    b(2,1) = i * -1; //loop through different portfolio target
                    b(3,1) = -1;
                    s(1,1) = b - inverseMet * x0;
                    p(1,1) = s(1,1);
                    for (int k = 0; st*s <= 0.000001; k++) //may use while loop instead
                        {   
                            mat sT = s.t();
                            a(k,1) = (sT * s)/(pT * Q * p)
                        }
                }

    }
