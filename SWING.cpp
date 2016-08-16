//
//  SWING.cpp
//  Final_Problem1
//
//  Created by Huachen Qin on 12/17/15.
//  Copyright © 2015 Huachen Qin. All rights reserved.
//

#include "SWING.h"

// constructor
 SwingOption::SwingOption(double strike_, double expiry_, int num_of_exercises_, double r_, double sigma_, vector<double>& exercises_)
 {
     strike = strike_;
     expiry = expiry_;
     r = r_;
     sigma = sigma_;
     num_of_exercises = num_of_exercises_;
     exercises.resize(num_of_exercises);
     exercises = exercises_;
 }
 
 // destructor
 SwingOption::~SwingOption()
 {
     exercises.clear();
 }

// payoff calculator of call options
double Call_PayOff(double spot, double strike)
{
    return spot>strike?(spot-strike):0.0;
}

// payoff calculator of put options
double Put_PayOff(double spot, double strike)
{
    return strike>spot?(strike-spot):0.0;
}

//the A_n_m, B_n_m, C_n_m coefficients in Crank Nicolson Method
vector<double> SwingOption::coeff_ABC(double m, double n, double delta_t,double delta_s)
{
    double a = pow(sigma*m*delta_s,2) / 2.0;
    double b = r * m * delta_s;
    double c = -r;
    double v_1 = delta_t / pow(delta_s,2);
    double v_2 = delta_t / delta_s;
    
    vector<double> coeff = {(-v_1 * a / 2.0 + v_2 * b / 4.0), (v_1 * a- delta_t * c / 2.0), (-v_1 * a / 2.0 - v_2 * b / 4.0)};
    return coeff;
}
 
 // get time to maturity
 double SwingOption::GetExpiry() const
 {
     return expiry;
 }

// get number of remaining exercise opportunities
 int SwingOption::GetNumExpiry() const
 {
     return num_of_exercises;
 }

// pricing option given the exercise date
 double SwingOption::pricing(double time, double Spot, int k_, double exercise_time1, double exercise_time2)
 {
     assert((expiry-time)<=expiry);
     
     double M = 40;        // divide Stock price S in M steps
     double N = pow(3,4);  // divide time to maturity T in 3^4 steps
     vector<vector<double>> grid;     // finite difference grid, stores all the values of v
     grid.resize(N+1);
     for(int i=0;i<=N;i++)
         grid[i].resize(M+1);
 
     double S_max = 4 * strike;      // S_max
     double delta_s = S_max / M;     // delta_s
     double delta_t = expiry / N;    // delta_t
     
     if((exercise_time1==1.0) || exercise_time2==1.0)  // if exercise time includes terminal time T, or not
     {
         for(int i=0;i<=M;i++)
             grid[N][i] = Call_PayOff(i*delta_s, strike);
     }
     else
     {
         for(int i=0;i<=M;i++)
             grid[N][i] = 0;
     }
     
     for(int j=0;j<=N;j++)
     {
         grid[j][0] = 0;
         grid[j][M] = Call_PayOff(S_max, strike);
     }
     
     // construct matrices of linear equation from Crank Nicolson method and solve by LU decomposion and Guassian Method
     for(int k=N-1;k>=0;k--)
     {
         Matrix M_L(M-1,M-1);                  // M_L matrix
         M_L.setelements(0, 0, 1 + coeff_ABC(1, k, delta_t, delta_s)[1]);
         M_L.setelements(0, 1, coeff_ABC(1, k, delta_t, delta_s)[2]);
         for(int i=1;i<(M-2);i++)
         {
             M_L.setelements(i, i-1, coeff_ABC(i+1, k, delta_t, delta_s)[0]);
             M_L.setelements(i, i, 1 + coeff_ABC(i+1, k, delta_t, delta_s)[1]);
             M_L.setelements(i, i+1, coeff_ABC(i+1, k, delta_t, delta_s)[2]);
         }
         M_L.setelements(M-2, M-3, coeff_ABC(M-1, k, delta_t, delta_s)[0] - coeff_ABC(M-1, k, delta_t, delta_s)[2]);
         M_L.setelements(M-2, M-2, 1 + coeff_ABC(M-1, k, delta_t, delta_s)[1] + 2 * coeff_ABC(M-1, k, delta_t, delta_s)[2]);
         
         Matrix M_R(M-1,M-1);                 // M_R matrix
         M_R.setelements(0, 0, 1 - coeff_ABC(1, k+1, delta_t, delta_s)[1]);
         M_R.setelements(0, 1, -coeff_ABC(1, k+1, delta_t, delta_s)[2]);
         for(int i=1;i<(M-2);i++)
         {
             M_R.setelements(i, i-1, -coeff_ABC(i+1, k+1, delta_t, delta_s)[0]);
             M_R.setelements(i, i, 1 - coeff_ABC(i+1, k+1, delta_t, delta_s)[1]);
             M_R.setelements(i, i+1, -coeff_ABC(i+1, k+1, delta_t, delta_s)[2]);
         }
         M_R.setelements(M-2, M-3, -coeff_ABC(M-1, k+1, delta_t, delta_s)[0] + coeff_ABC(M-1, k+1, delta_t, delta_s)[2]);
         M_R.setelements(M-2, M-2, 1 - coeff_ABC(M-1, k+1, delta_t, delta_s)[1] - 2 *  coeff_ABC(M-1, k+1, delta_t, delta_s)[2]);
         
         Matrix v_N(M-1,1);                 // v_N matrix
         for(int i=1;i<M;i++)
             v_N.setelements(i-1, 0, grid[k+1][i]);
        
         Matrix r_(M-1,1);                 // r_ matrix
         r_.setelements(0, 0, (-coeff_ABC(1, k+1, delta_t, delta_s)[0] - coeff_ABC(1, k, delta_t, delta_s)[1]) * 0.0);
         
         Matrix b(M-1,1);
         b = ((M_R * v_N) + r_);
         
         Matrix result(M-1,1);        // solve v such that M_L * v = b by LU decomposion and Gaussian Method
         result = GSmethod(M_L, b);
         
         for(int i=1;i<M;i++)
             grid[k][i] =  result.getelements(i-1,0);
         
         if(k == (exercise_time1 / delta_t) || k == (exercise_time2 / delta_t))   // calcualte if early exercise is optimal
         {
             for(int j=1;j<M;j++)
             {
                 if(grid[k][j]<Call_PayOff(j*delta_s, strike))
                     grid[k][j] += Call_PayOff(j*delta_s, strike);
             }
         }
         
     }
     
     return grid[0][(Spot/delta_s)];     // return option price when S = 50
}

// general pricing of all possible scenarios
double SwingOption::general_pricing(double T, double Spot, int k)
{
    assert((expiry-T)<=expiry);
    double result_1 = pricing(T, Spot, k, exercises[0], exercises[1])>pricing(T, Spot, k, exercises[0], exercises[2])?pricing(T, Spot, k, exercises[0], exercises[1]):pricing(T, Spot, k, exercises[0], exercises[2]);
    double result_2 = pricing(T, Spot, k, exercises[1], exercises[2])>result_1?pricing(T, Spot, k, exercises[1], exercises[2]):result_1;
    return result_2;
}

void SwingOption::print_grid(double time, double Spot, int k_, double exercise_time1, double exercise_time2)
{
    assert((expiry-time)<=expiry);
    
    double M = 40;        // divide Stock price S in M steps
    double N = pow(3,4);  // divide time to maturity T in 3^4 steps
    vector<vector<double>> grid;     // finite difference grid, stores all the values of v
    grid.resize(N+1);
    for(int i=0;i<=N;i++)
        grid[i].resize(M+1);
    
    double S_max = 4 * strike;      // S_max
    double delta_s = S_max / M;     // delta_s
    double delta_t = expiry / N;    // delta_t
    
    if((exercise_time1==1.0) || exercise_time2==1.0)  // if exercise time includes terminal time T, or not
    {
        for(int i=0;i<=M;i++)
            grid[N][i] += Call_PayOff(i*delta_s, strike);
    }
    else
    {
        for(int i=0;i<=M;i++)
            grid[N][i] = 0;
    }
    
    for(int j=0;j<=N;j++)
    {
        grid[j][0] = 0;
        grid[j][M] = Call_PayOff(S_max, strike);
    }
    
    // construct matrices of linear equation from Crank Nicolson method and solve by LU decomposion and Guassian Method
    for(int k=N-1;k>=0;k--)
    {
        Matrix M_L(M-1,M-1);                  // M_L matrix
        M_L.setelements(0, 0, 1 + coeff_ABC(1, k, delta_t, delta_s)[1]);
        M_L.setelements(0, 1, coeff_ABC(1, k, delta_t, delta_s)[2]);
        for(int i=1;i<(M-2);i++)
        {
            M_L.setelements(i, i-1, coeff_ABC(i+1, k, delta_t, delta_s)[0]);
            M_L.setelements(i, i, 1 + coeff_ABC(i+1, k, delta_t, delta_s)[1]);
            M_L.setelements(i, i+1, coeff_ABC(i+1, k, delta_t, delta_s)[2]);
        }
        M_L.setelements(M-2, M-3, coeff_ABC(M-1, k, delta_t, delta_s)[0] - coeff_ABC(M-1, k, delta_t, delta_s)[2]);
        M_L.setelements(M-2, M-2, 1 + coeff_ABC(M-1, k, delta_t, delta_s)[1] + 2 * coeff_ABC(M-1, k, delta_t, delta_s)[2]);
        
        Matrix M_R(M-1,M-1);                 // M_R matrix
        M_R.setelements(0, 0, 1 - coeff_ABC(1, k+1, delta_t, delta_s)[1]);
        M_R.setelements(0, 1, -coeff_ABC(1, k+1, delta_t, delta_s)[2]);
        for(int i=1;i<(M-2);i++)
        {
            M_R.setelements(i, i-1, -coeff_ABC(i+1, k+1, delta_t, delta_s)[0]);
            M_R.setelements(i, i, 1 - coeff_ABC(i+1, k+1, delta_t, delta_s)[1]);
            M_R.setelements(i, i+1, -coeff_ABC(i+1, k+1, delta_t, delta_s)[2]);
        }
        M_R.setelements(M-2, M-3, -coeff_ABC(M-1, k+1, delta_t, delta_s)[0] + coeff_ABC(M-1, k+1, delta_t, delta_s)[2]);
        M_R.setelements(M-2, M-2, 1 - coeff_ABC(M-1, k+1, delta_t, delta_s)[1] - 2 *  coeff_ABC(M-1, k+1, delta_t, delta_s)[2]);
        
        Matrix v_N(M-1,1);                 // v_N matrix
        for(int i=1;i<M;i++)
            v_N.setelements(i-1, 0, grid[k+1][i]);
        
        Matrix r_(M-1,1);                 // r_ matrix
        r_.setelements(0, 0, (-coeff_ABC(1, k+1, delta_t, delta_s)[0] - coeff_ABC(1, k, delta_t, delta_s)[1]) * 0.0);
        
        Matrix b(M-1,1);
        b = ((M_R * v_N) + r_);
        
        Matrix result(M-1,1);        // solve v such that M_L * v = b by LU decomposion and Gaussian Method
        result = GSmethod(M_L, b);
        
        for(int i=1;i<M;i++)
            grid[k][i] =  result.getelements(i-1,0);
        
        if(k == (exercise_time1 / delta_t) || k == (exercise_time2 / delta_t))   // calcualte if early exercise is optimal
        {
            for(int j=1;j<M;j++)
            {
                if(grid[k][j]<Call_PayOff(j*delta_s, strike))
                    grid[k][j] += Call_PayOff(j*delta_s, strike);
            }
        }
        
    }
    
    for(int i=0;i<=M;i++)
        cout<<"when S = "<<(i*delta_s)<<", v = "<<grid[0][i]<<endl;

}






