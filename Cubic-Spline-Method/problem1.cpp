// HW6 Problem1
// Huachen Qin

#include<iostream>
#include<cmath>
#include <array>

using namespace std;

double bond_price(double coupon, double * spot_rate, double T);
void cubic_splines(double *x, double *y, int size, double *a_coeff, double *b_coeff, double *c_coeff, double *d_coeff);
double linear_exterpolation(double target, double start, double start_value, double end, double end_value);
void print(double *arr,int size);
double spline1();
double *spline2();
double *spline3();
double *spline4();
double *spline5();
double *spline6();
double forward_rate(double T, double delta);

const double face_value = 100;

int main()
{
    double T, delta;
    
    cout<<"The results of first extropolation is: "<<endl;
    cout<<spline1()<<endl;
    cout<<endl;
    
    cout<<"The results of second extropolation is: "<<endl;
    print(spline2(),2);
    cout<<endl;
    
    cout<<"The results of third extropolation is: "<<endl;
    print(spline3(),5);
    cout<<endl;
    
    cout<<"The results of fourth extropolation is: "<<endl;
    print(spline4(),8);
    cout<<endl;
    
    cout<<"The results of fifth extropolation is: "<<endl;
    print(spline5(),13);
    cout<<endl;
    
    cout<<"The results of sixth extropolation is: "<<endl;
    print(spline6(),52);
    cout<<endl;
    
    cout<<"Now Let's use forward rate functions to calculate forward rate."<<endl;
    cout<<"Please enter future starting time T: ";
    cin>>T;
    cout<<"Please enter future time period delta: ";
    cin>>delta;
    cout<<"The forward rate between "<<T<<" and "<<delta+T<<" is:"<<forward_rate(T,delta)<<endl;
    cout<<endl;
    
    cout<<"The entire program ends, please see the interest rate curve in the attached excel file."<<endl;
}

double spline1()
{
    int size = 6; double target = 1.5; double rate;
    double * a = new double[size - 1];
    double * b = new double[size - 1];
    double * c = new double[size - 1];
    double * d = new double[size - 1];
    double * x = new double[size]{0,0.0795,0.2521,0.5014,0.9425,1.9836};
    double * y = new double[size]{0,0,0,0.000797926,0.002443129,0.006259908};
    
    cubic_splines(x, y, size, a, b, c, d);
    
    int i = 0;
    while (!(target > x[i] && target <= x[i + 1])) {
        i++;
    }
    
    rate = a[i] + b[i] * (target - x[i]) + c[i] * pow(target - x[i], 2) + d[i] * pow(target - x[i], 3);
    
    delete[]a;
    delete[]b;
    delete[]c;
    delete[]d;
    delete[]x;
    delete[]y;
    return rate;
    
}

double *spline2()
{
    int size = 7; static double rate[2];
    double * target = new double[2]{1.5,2.5};
    double * a = new double[size - 1];
    double * b = new double[size - 1];
    double * c = new double[size - 1];
    double * d = new double[size - 1];
    double * x = new double[size]{0,0.0795,0.2521,0.5014,0.9425,1.9836,3.0247};
    double * y = new double[size]{0,0,0,0.000797926,0.002443129,0.006259908,0.007792};
    
    cubic_splines(x, y, size, a, b, c, d);
    
    for(int j=0;j<2;j++)
    {
        int i = 0;
        while (!(target[j] > x[i] && target[j] <= x[i + 1])) i++;
        rate[j] = a[i] + b[i] * (target[j] - x[i]) + c[i] * pow(target[j] - x[i], 2) + d[i] * pow(target[j] - x[i], 3);
    }
    
    delete[]a;
    delete[]b;
    delete[]c;
    delete[]d;
    delete[]x;
    delete[]y;
    delete[]target;
    return rate;
    
}

double *spline3()
{
    int size = 8; static double rate2[5];
    double * target = new double[5]{1.5,2.5,3.5,4,4.5};
    double * a = new double[size - 1];
    double * b = new double[size - 1];
    double * c = new double[size - 1];
    double * d = new double[size - 1];
    double * x = new double[size]{0,0.0795,0.2521,0.5014,0.9425,1.9836,3.0247,4.9863};
    double * y = new double[size]{0,0,0,0.000797926,0.002443129,0.006259908,0.007792,0.013741};
    
    cubic_splines(x, y, size, a, b, c, d);
    
    for(int j=0;j<5;j++)
    {
        int i = 0;
        while (!(target[j] > x[i] && target[j] <= x[i + 1])) i++;
        rate2[j] = a[i] + b[i] * (target[j] - x[i]) + c[i] * pow(target[j] - x[i], 2) + d[i] * pow(target[j] - x[i], 3);
    }
    
    delete[]a;
    delete[]b;
    delete[]c;
    delete[]d;
    delete[]x;
    delete[]y;
    delete[]target;
    return rate2;
}

double *spline4()
{
    int size = 9; static double rate3[8];
    double * target = new double[8]{1.5,2.5,3.5,4,4.5,5.5,6,6.5};
    double * a = new double[size - 1];
    double * b = new double[size - 1];
    double * c = new double[size - 1];
    double * d = new double[size - 1];
    double * x = new double[size]{0,0.0795,0.2521,0.5014,0.9425,1.9836,3.0247,4.9863,6.9863};
    double * y = new double[size]{0,0,0,0.000797926,0.002443129,0.006259908,0.007792,0.013741,0.017704};
    
    cubic_splines(x, y, size, a, b, c, d);
    
    for(int j=0;j<8;j++)
    {
        int i = 0;
        while (!(target[j] > x[i] && target[j] <= x[i + 1])) i++;
        rate3[j] = a[i] + b[i] * (target[j] - x[i]) + c[i] * pow(target[j] - x[i], 2) + d[i] * pow(target[j] - x[i], 3);
    }
    
    delete[]a;
    delete[]b;
    delete[]c;
    delete[]d;
    delete[]x;
    delete[]y;
    delete[]target;
    return rate3;
}

double *spline5()
{
    int size = 10; static double rate4[13];
    double * target = new double[13]{1.5,2.5,3.5,4,4.5,5.5,6,6.5,7.5,8,8.5,9,9.5};
    double * a = new double[size - 1];
    double * b = new double[size - 1];
    double * c = new double[size - 1];
    double * d = new double[size - 1];
    double * x = new double[size]{0,0.0795,0.2521,0.5014,0.9425,1.9836,3.0247,4.9863,6.9863,9.863};
    double * y = new double[size]{0,0,0,0.000797926,0.002443129,0.006259908,0.007792,0.013741,0.017704,0.020575};
    
    cubic_splines(x, y, size, a, b, c, d);
    
    for(int j=0;j<13;j++)
    {
        int i = 0;
        while (!(target[j] > x[i] && target[j] <= x[i + 1])) i++;
        rate4[j] = a[i] + b[i] * (target[j] - x[i]) + c[i] * pow(target[j] - x[i], 2) + d[i] * pow(target[j] - x[i], 3);
    }
    
    delete[]a;
    delete[]b;
    delete[]c;
    delete[]d;
    delete[]x;
    delete[]y;
    delete[]target;
    return rate4;
}

double *spline6()
{
    int size = 11; static double rate5[52];
    double * target = new double[52]{1.5,2.5,3.5,4,4.5,5.5,6,6.5,7.5,8,8.5,9,9.5,10.5,11,11.5,12,12.5,13,13.5,14,14.5,15,15.5,16,16.5,17,17.5,18,18.5,19,19.5,20,20.5,21,21.5,22,22.5,23,23.5,24,24.5,25,25.5,26,26.5,27,27.5,28,28.5,29,29.5};
    double * a = new double[size - 1];
    double * b = new double[size - 1];
    double * c = new double[size - 1];
    double * d = new double[size - 1];
    double * x = new double[size]{0,0.0795,0.2521,0.5014,0.9425,1.9836,3.0247,4.9863,6.9863,9.863,29.877};
    double * y = new double[size]{0,0,0,0.000797926,0.002443129,0.006259908,0.007792,0.013741,0.017704,0.020575,0.027214};
    
    cubic_splines(x, y, size, a, b, c, d);
    
    for(int j=0;j<52;j++)
    {
        int i = 0;
        while (!(target[j] > x[i] && target[j] <= x[i + 1])) i++;
        rate5[j] = a[i] + b[i] * (target[j] - x[i]) + c[i] * pow(target[j] - x[i], 2) + d[i] * pow(target[j] - x[i], 3);
    }
    
    delete[]a;
    delete[]b;
    delete[]c;
    delete[]d;
    delete[]x;
    delete[]y;
    delete[]target;
    return rate5;
}

double forward_rate(double T, double delta)
{
    if(T+delta>30) {cout<<"Error Input! Please enter future time and delta again."<<endl; return 0;}
    else
    {int size = 11; double rate1, rate2, rate;
        
        double * a = new double[size - 1];
        double * b = new double[size - 1];
        double * c = new double[size - 1];
        double * d = new double[size - 1];
        double * x = new double[size]{0,0.0795,0.2521,0.5014,0.9425,1.9836,3.0247,4.9863,6.9863,9.863,29.877};
        double * y = new double[size]{0,0,0,0.000797926,0.002443129,0.006259908,0.007792,0.013741,0.017704,0.020575,0.027214};
        
        cubic_splines(x, y, size, a, b, c, d);
        
        int i = 0;
        while (!(T > x[i] && T <= x[i + 1])) i++;
        rate1 = a[i] + b[i] * (T - x[i]) + c[i] * pow(T - x[i], 2) + d[i] * pow(T - x[i], 3);
        
        int j = 0;
        while(!((T+delta) > x[j] && (T+delta) <=x[j+1])) j++;
        rate2 = a[j] + b[j] * (T - x[j]) + c[j] * pow(T - x[j], 2) + d[j] * pow(T - x[j], 3);
        
        rate = ((T +delta) * rate2 - T *rate1) / delta;
        
        delete[]a;
        delete[]b;
        delete[]c;
        delete[]d;
        delete[]x;
        delete[]y;
        return rate;}
    
}

void cubic_splines(double *x, double *y, int size, double *a_coeff, double *b_coeff, double *c_coeff, double *d_coeff)
{
    double *h = new double[size-1];
    double *d = new double[size-2];
    double *a = new double[size-2];
    double *b = new double[size-2];
    double *c = new double[size-2];
    double *c_hat = new double[size-2];
    double *d_hat = new double[size-2];
    double *M = new double[size];
    
    for(int i=0;i<size-1;i++) h[i] = x[i+1] - x[i];
    
    for(int k=0;k<size-2;k++) d[k] = (y[k+2] - y[k+1]) / h[k+1] - (y[k+1] - y[k]) / h[k];
    
    for(int m=1;m<size-2;m++) a[m] = h[m];a[0] = 0;
    
    for(int n=0;n<size-2;n++) b[n] = 2 * (h[n] + h[n+1]);
    
    for(int l=0;l<size-3;l++) c[l] = h[l+1]; c[size-3] = 0;
    
    c_hat[0] = c[0] / b[0]; for(int p=1;p<size-2;p++) c_hat[p] = c[p] / (b[p] - c_hat[p-1] * a[p]);
    
    d_hat[0] = d[0] / b[0]; for(int q=1;q<size-2;q++) d_hat[q] = (d[q] - d_hat[q-1] * a[q]) / (b[q] - c_hat[q-1] * a[q]);
    
    M[0] = M[size-1] = 0; M[size-2] = d_hat[size-3];
    for(int r=size-3;r>0;--r) M[r] = d_hat[r-1] - c_hat[r-1] * M[r+1];
    
    for(int i=0;i<size-1;i++)
    {   a_coeff[i] = y[i];
        b_coeff[i] = (y[i+1] - y[i]) / h[i] - h[i] * M[i] / 2 - h[i] * (M[i+1] - M[i]) / 6;
        c_coeff[i] = M[i] / 2;
        d_coeff[i] = (M[i+1] - M[i]) / (6 * h[i]);
    }
    
    delete []h; delete []d; delete []a; delete []b; delete []c; delete []c_hat; delete []d_hat;
    delete []M;
}

double linear_exterpolation(double target, double start, double start_value, double end, double end_value)
{
    return ((target - start) * end_value - (target - end) * start_value) / (end - start);
}

double bond_price(double coupon, double * spot_rate, double T)
{
    int m = 0;
    while(abs(0.5*m-T)>0.2) m++;
    
    if(T<0.5)
        return face_value * exp(-spot_rate[m-1]*T);
    else
    {
        double price = face_value * exp(-spot_rate[m-1] * T);
        for(int i=0;i<m;i++) price += coupon/2 * face_value * exp(-spot_rate[m-i-1] * (T - 0.5 * i));
        return price;
    }
}

void print(double *arr,int size)
{
    int len = size;
    for(int i=0;i<len-1;i++) cout<<arr[i]<<", ";
    cout<<arr[len-1]<<endl;
}

