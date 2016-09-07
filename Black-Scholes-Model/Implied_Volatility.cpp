// HW5 Problem2
// Huachen Qin
// HW3 Problem4
//Huachen Qin

#include<iostream>
#include<cmath>

using namespace std;

// declare functions
double auxillary_function(double x);
double auxillary_function2(double x, double K);
double BlackScholesCall(double S, double K, double sigma, double tau, double r);
double BlackScholesPut(double S, double K, double sigma, double tau, double r);
double interval_bisection(double target,double left, double right, double error, double K);
double newton_raphson(double target, double initial, double error, double K);
double secant_method(double target, double initial1, double initial2, double error, double K);

// define global constant variables
const double pi = 3.1415926535897;
const double S = 28;
const int size = 7;
const double K[size] = {17.5,20,22.5,25,27.5,30,32.5};
const double tau = 1;
const double r = 0;
const double put_price[size] = {0.28,0.48,0.81,1.41,2.42,3.8,5.56};



int main()
{
    double sigma[size];
    double sigma2[size];
    double sigma3[size];
    
    double left = 0;
    double right = 0.5;
    double error = 0.001;
    double initial = 0.5;
    double initial1 = 0.1;
    double initial2 = 0.2;
    
    cout<<"The result by Bisection Method is:"<<endl;
    for(int i=0;i<size;i++)
        sigma[i] = interval_bisection(put_price[i], left, right, error, K[i]);
    for(int j=0;j<size;j++)
        cout<<"When the strike price K = "<<K[j]<<" , implied volatility = "<<sigma[j]<<endl;
    cout<<endl;
    
    cout<<"The result by Newton Raphson Method is:"<<endl;
    for(int m=0;m<size;m++)
        sigma2[m] = newton_raphson(put_price[m], initial, error, K[m]);
    for(int n=0;n<size;n++)
        cout<<"When the strike price K = "<<K[n]<<" , implied volatility = "<<sigma2[n]<<endl;
    cout<<endl;
    
    cout<<"The result by Secant Method is:"<<endl;
    for(int a=0;a<size;a++)
        sigma3[a] = secant_method(put_price[a], initial1, initial2, error, K[a]);
    for(int b=0;b<size;b++)
        cout<<"When the strike price K = "<<K[b]<<" , implied volatility = "<<sigma3[b]<<endl;
    cout<<endl;
        
}

double newton_raphson(double target, double initial, double error, double K)
{
    double estimate = initial;
    double temp;
    
    do
    {
        temp = estimate;
        estimate = estimate - (BlackScholesPut(S, K,estimate, tau, r) - target) /(S * auxillary_function2(estimate, K) * sqrt(tau));
    } while(abs(temp-estimate)>=error);
    
    return estimate;
}

double secant_method(double target, double initial1, double initial2, double error, double K)
{
    double estimate;
    double temp1 = initial1;
    double temp2 = initial2;
    
    do
    {
        estimate = temp2 - (BlackScholesPut(S, K, temp2, tau, r) - target) * ( temp2 - temp1) / ((BlackScholesPut(S, K, temp2, tau, r) - target) - (BlackScholesPut(S, K, temp1, tau, r) - target));
        temp1 = temp2;
        temp2 = estimate;
    } while(abs(estimate-temp1)>=error);
    
    return estimate;
}

double interval_bisection(double target,double left, double right, double error, double K)   // define the fuction that uses bisection method to approximate implied volatility
{
    double mid = 0.5 * (left + right);
    
    if((BlackScholesPut(S, K, left, tau, r)-target) * (BlackScholesPut(S, K, right, tau, r) - target) < 0)
    {
        
        do
        {
            if((BlackScholesPut(S, K, mid, tau, r)-target) * (BlackScholesPut(S, K, right, tau, r) - target) < 0)
                left = mid;
            
            if((BlackScholesPut(S, K, mid, tau, r)-target) * (BlackScholesPut(S, K, left, tau, r) - target) < 0)
                right = mid;
            
            mid = 0.5 * (left + right);
            
        } while(abs(left - right) > error);
        
        return mid;
    }
    
    else                                                    // if initial bounds [a,b] do not contain the function root, the function would return error
    {
        cout<<"The initial bounds are wrong!"<<endl;
        return 0;
    }
}


double BlackScholesCall(double S, double K, double sigma, double tau, double r) // define the function _ Call
{
    double d_positive = (log(S/K) + (r + sigma * sigma /2)*tau) / (sigma * pow(tau,0.5));
    double d_negative = (log(S/K) + (r - sigma * sigma /2)*tau) / (sigma * pow(tau,0.5));
    
    double call_price = 0.0;
    
    if(d_positive>=0 && d_negative>=0)                                        // consider the situation that d_positive and d_negative could be negative
        call_price = S * auxillary_function(d_positive) - K * exp(-r * tau) * auxillary_function(d_negative);
    if(d_positive>=0 && d_negative<=0)
        call_price = S * auxillary_function(d_positive) - K * exp(-r * tau) * (1 - auxillary_function(-d_negative));
    if(d_positive<=0 && d_negative>=0)
        call_price = S * (1 - auxillary_function(-d_positive)) - K * exp(-r * tau) * auxillary_function(d_negative);
    if(d_positive<=0 && d_negative<=0)
        call_price = S * (1 - auxillary_function(-d_positive)) - K * exp(-r * tau) * (1 - auxillary_function(d_negative));
    return call_price;
}

double BlackScholesPut(double S, double K, double sigma, double tau, double r)  // define the function _ Put
{
    double d_positive = (log(S/K) + (r + sigma * sigma /2)*tau) / (sigma * pow(tau,0.5));
    double d_negative = (log(S/K) + (r - sigma * sigma /2)*tau) / (sigma * pow(tau,0.5));
    
    double put_price = 0.0;
    
    if(d_positive<=0 && d_negative<=0)                                        // consider the situation that d_positive and d_negative could be negative
        put_price = K * exp(-r * tau) * auxillary_function(-d_negative) - S * auxillary_function(-d_positive);
    if(d_positive>=0 && d_negative<=0)
        put_price = K * exp(-r * tau) * auxillary_function(-d_negative) - S * (1 - auxillary_function(d_positive));
    if(d_positive<=0 && d_negative>=0)
        put_price = K * exp(-r * tau) * (1 - auxillary_function(d_negative)) - S * auxillary_function(-d_positive);
    if(d_positive>=0 && d_negative>=0)
        put_price = K * exp(-r * tau) * (1 - auxillary_function(d_negative)) - S * (1-auxillary_function(d_positive));
    return put_price;
}

double auxillary_function(double x)   // define the auxiliary function to calculate Î¦(x)
{
    double b0 = 0.2316419, b1 = 0.3193815300, b2 = -0.356563782, b3 = 1.7814779370, b4 = -1.821255978, b5 = 1.3302744290;
    double phi = ( 1 / pow(2*pi,0.5)) * exp(-x*x/2);
    double t = 1/(1+b0*x);
    return 1 - phi * (b1*t+b2*pow(t,2)+b3*pow(t,3)+b4*pow(t,4)+b5*pow(t,5));
}

double auxillary_function2(double x, double K)     // define the auxiliary function to calculate phi(x)
{
    double d_positive = (log(S/K) + (r + x * x /2)*tau) / (x * pow(tau,0.5));
    return ( 1 / sqrt(2*pi)) * exp(-d_positive * d_positive / 2);
}
