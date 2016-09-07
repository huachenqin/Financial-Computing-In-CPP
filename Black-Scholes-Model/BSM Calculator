// HW5 Problem1
// Huachen Qin

#include<iostream>
#include <cstdlib>
#include <ctime>
#include<cmath>

using namespace std;
double call_payoff(double St, double K);
double call_price_BM(double S0, double K, double r, double sigma, double T);
double call_price_AJ(double S0, double K, double r, double sigma, double T);
double call_price_Marsaglia(double S0, double K, double r, double sigma, double T);

const double pi = 3.14159265359;

int main()
{
    double T = 1;
    double S0 = 100;
    double K = 100;
    double r = 0.05;
    double sigma = 0.2;
    cout<<"The call option price calculated by Box-Muller Mtheod is:"<<call_price_BM(S0, K, r, sigma, T)<<endl;
    cout<<endl;
    
    cout<<"The call option price calculated by Accept-Reject Mtheod is:"<<call_price_AJ(S0, K, r, sigma, T)<<endl;
    cout<<endl;
    
    cout<<"The call option price calculated by Marsaglia Mtheod is:"<<call_price_Marsaglia(S0, K, r, sigma, T)<<endl;
    cout<<endl;
}

double call_payoff(double St, double K)                    // question (a), define payoff function of European call price
{
    if(St>=K)
        return St - K;
    else
        return 0;
}

double call_price_BM(double S0, double K, double r, double sigma, double T)     // question (b), Box-Muller Method
{
    const int n = 10000;                           // number of the experiments is 10000
    double call_price[n];
    double sum = 0;
    double temp1, temp2, Z, ST;
    
    for(int i=0;i<n;i++)
    {
        temp1 = rand()/(static_cast<double>(RAND_MAX));
        temp2 = rand()/(static_cast<double>(RAND_MAX));
        Z = sqrt(-2 * log(temp1)) * cos(2 * pi * temp2);
        ST = S0 * exp((r - sigma * sigma / 2) * T + sigma * sqrt(T) * Z);
        call_price[i] = exp(-r * T) * call_payoff(ST, K);
    }
    
    for(int j=0;j<n;j++)
        sum = sum + call_price[j];
    
    return sum/n;
}

double call_price_AJ(double S0, double K, double r, double sigma, double T)           // question (c), Accept-Reject Method
{
    const int n = 10000;                           // number of the experiments is 10000
    double call_price[n];
    double sum = 0;
    
    double temp1, temp2, Y, Z, ST;
    
    for(int i=0;i<n;i++)
    {
        do
        {
            temp1 = rand()/(static_cast<double>(RAND_MAX));
            temp2 = rand()/(static_cast<double>(RAND_MAX));
            Y = - log(1 - temp2);                            // generate exponential random variable by Inverse Transform Method (lambda = 1)
        } while(temp1>exp(-(Y-1)*(Y-1)/2));
        
        double temp3 = rand()/(static_cast<double>(RAND_MAX));
        temp3 <= 0.5? Z = abs(Y) : Z = - abs(Y);
        ST = S0 * exp((r - sigma * sigma / 2) * T + sigma * sqrt(T) * Z);
        call_price[i] = exp(-r * T) * call_payoff(ST, K);
    }
    
    for(int j=0;j<n;j++)
        sum = sum + call_price[j];
    
    return sum/n;
}

double call_price_Marsaglia(double S0, double K, double r, double sigma, double T)        // question (c), Marsaglia Method
{
    const int n = 10000;                           // number of the experiments is 10000
    double call_price[n];
    double sum = 0;
    double temp1, temp2, S, Z, ST;
    
    for(int i=0;i<n;i++)
    {
        do
        {
            temp1 = rand()/(static_cast<double>(RAND_MAX)) * 2.0 - 1.0;
            temp2 = rand()/(static_cast<double>(RAND_MAX)) * 2.0 - 1.0;
            S = temp1 * temp1 + temp2 * temp2;
        } while(S>=1);
        
        Z = temp2 * sqrt(-2 * log(S)/S);
        ST = S0 * exp((r - sigma * sigma / 2) * T + sigma * sqrt(T) * Z);
        call_price[i] = exp(-r * T) * call_payoff(ST, K);
    }
    
    for(int j=0;j<n;j++)
        sum = sum + call_price[j];
    
    return sum/n;
}
