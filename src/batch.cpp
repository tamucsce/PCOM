#include<iostream>
#include<cmath>

using namespace std;

// Probs that b bins overflow where r is size of each bin and n is total number of elements

double C(int n, int r, int b) 
{
        double err=1.00;
        double dummy=1.00;
        for (int i=n; i>n-r; i--)
                err = err * ((i*dummy)/(b*(i-(n-r))));
        for (int i=0; i<n-i; i++)
                err = err * ((b-1)*dummy/b);
        return err;

}

int main()
{       
        int n;
        double error= pow(2.0, -40);
        cout<<error<<endl;
        double err=0.000;
        n=8;
        for (int i=0;i<14;i++)
        {
                cout<<"i: "<<i<<endl;
                n = n*2;
                for (int j=2; j<n; j++) //Counter for number of bins
                {
                        cout<<"j: "<<j<<endl;
                        for (int k=2; k<n; k++) //Counter for bin size
                        {
                                for (int l = k; l<=n; l++) //For summing the probabilities of a bin having more than l elements in a bin
                                        err = err + C(n,l,j);
                                err = j*err;
                                cout<<"k: "<<k<<endl;
                                cout<<err<<endl;
                                if (err < error)
                                {
                                        cout << "SUCCESS :" << " Input size = " << n << ", #BIN=" << j << " ,BIN-SIZE=" << k << endl;
                                        break;
                                }
                        }

                }
        }
        return 0;
}
