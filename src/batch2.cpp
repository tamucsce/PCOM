#include<iostream>
#include<cmath>

using namespace std;

// Probs that b bins overflow where r is size of each bin and n is total number of elements

mpf_class C(int n, int r, int b) 
{
        std::cout.precision(50);
        mpf_class err(1.00,50);
        mpf_class dummy(1.00,50);
        
        cout<<err<<endl;
        
        for (int i=n; i>n-r; i--)
                err = err * ((i*dummy)/(b*(i-(n-r))));
        for (int i=0; i<n-i; i++)
                err = err * ((b-1)*dummy/b);
        return err;

}

int main()
{
     std::cout.precision(50);

    //======================== debug ================================
    int n,tmp;
        //double error= pow(2.0, -40);
        mpf_class error(pow(2.0, -40),1024);
        mpf_class err(1.00,1024), total;

        n=8;
        for (int i=0;i<16;i++)
        {
                n = n*2;
                for (int j=2; j<n; j++) //Counter for number of bins
                {
                    err = pow((1.000/j),n);
                    total = j*err;
                    for (int k=n-1; k>=2; k--) //Counter for bin size
                    {
            
                        tmp = n-k;  
                        err = (err * (n-tmp-1) * (k-1))/tmp;
                        total = total + (j*err);  
                        cout<<err<<endl;            
                        if (err < error)
                        {
                                // cout << "RUNNING :" << " Input size = " << n << ", #BIN=" << j << " ,BIN-SIZE=" << k << endl;
                        }
                        else
                        {
                            cout << "SUCCESS :" << " Input size = " << n << ", #BIN=" << j << " ,BIN-SIZE=" << k+1 << endl;
                            break;
                        }   
                    }

                }
        }
}
    //======================== debug ================================
