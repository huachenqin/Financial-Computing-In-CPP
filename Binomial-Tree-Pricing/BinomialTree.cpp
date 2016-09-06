// HW10 problem1
// Huachen Qin


#include "BinomialTree.h"

/*---------------------------------------CRR Binomial Tree Model-------------------------------------------*/

CRRBTree::CRRBTree()                   // default constructor function
{
    option_price = S = K = r = sigma = T = u = pu = d = pd = 0;
    N =0;
    cout<< "The binomial tree is NULL!"<<endl;
}

CRRBTree::~CRRBTree()                 // destructor function
{
    stockprice.clear();
    optionvalue.clear();
}

CRRBTree::CRRBTree(double S_, double K_, double T_, double r_, double sigma_, int N_)   // constructor function
{
    S = S_; K = K_; T = T_; r = r_; sigma = sigma_; N = N_;
    assert(N>=1);
    delta_t = T / N;
    discount_factor = exp(-r * delta_t);
    
    u = exp(sigma * sqrt(delta_t));
    d = 1 / u;
    pu = (exp(r * delta_t) - d) / (u - d);
    pd = 1 - pu;
    
    stockprice.resize(N+1);
    for(int i=0;i<N+1;i++)
        stockprice[i].resize(i+1);
    optionvalue.resize(N+1);
    for(int i=0;i<N+1;i++)
        optionvalue[i].resize(i+1);
    
    for(int i=0;i<N+1;i++)
        for(int j=0;j<=i;j++)
            stockprice[i][j] = pow(u,i-j) * pow(d,j) * S;
}

double CRRBTree::optionpricing()        // option pricing function, return the option value
{
    if(N<1)
    {
        cout<<"Error! The step N must be greater than or equal to 1"<<endl;
        exit(1);
    }
    else
    {
        for(int i=0;i<N+1;i++)
            optionvalue[N][i] = stockprice[N][i]>K?(stockprice[N][i]-K):0;
        for(int i=N-1;i>=0;i--)
            for(int j=0;j<=i;j++)
                optionvalue[i][j] = discount_factor * (pu * optionvalue[i+1][j] + pd * optionvalue[i+1][j+1]);
        
        option_price = optionvalue[0][0];
        return option_price;
    }
}

double CRRBTree::pruning_1()            // pruning method 1, return option value after pruning
{
    double S_u = S * exp((r - sigma * sigma / 2) * T + 6 * sigma * sqrt(T));
    double S_d = S * exp((r - sigma * sigma / 2) * T - 6 * sigma * sqrt(T));
    
    if(N<1)
    {
        cout<<"Error! The step N must be greater than or equal to 1"<<endl;
        exit(1);
    }
    else
    {
        for(int i=0;i<N+1;i++)
        {
            if(stockprice[N][i]<=S_u && stockprice[N][i]>=S_d)
                optionvalue[N][i] = stockprice[N][i]>K?(stockprice[N][i]-K):0;
            else
                optionvalue[N][i] = 0;
        }
        for(int i=N-1;i>=0;i--)
            for(int j=0;j<=i;j++)
            {
                if(stockprice[i][j]<=S_u && stockprice[i][j]>=S_d)
                    optionvalue[i][j] = discount_factor * (pu * optionvalue[i+1][j] + pd * optionvalue[i+1][j+1]);
                else
                    optionvalue[i][j] = 0;
            }
        
        double option_price = optionvalue[0][0];
        return option_price;
    }
}

double CRRBTree::pruning_2()                  // pruning method 2, returns option value after pruning
{
    double S_u = S * exp((r - sigma * sigma / 2) * T + 6 * sigma * sqrt(T));
    double S_d = S * exp((r - sigma * sigma / 2) * T - 6 * sigma * sqrt(T));
    
    if(N<1)
    {
        cout<<"Error! The step N must be greater than or equal to 1"<<endl;
        exit(1);
    }
    else
    {
        for(int i=0;i<N+1;i++)
            optionvalue[N][i] = stockprice[N][i]>K?(stockprice[N][i]-K):0;
        for(int i=N-1;i>=0;i--)
            for(int j=0;j<=i;j++)
            {
                if(stockprice[i][j]<=S_u && stockprice[i][j]>=S_d)
                    optionvalue[i][j] = discount_factor * (pu * optionvalue[i+1][j] + pd * optionvalue[i+1][j+1]);
                else
                    optionvalue[i][j] = S > K*exp(-r*(T-i*delta_t))?(S-K*exp(-r*(T-i*delta_t))):0;
            }
        
        option_price = optionvalue[0][0];
        return option_price;
    }
}

double CRRBTree::pruning_3()                            // pruning method 3, return option value after pruning
{
    double S_u = S * exp((r - sigma * sigma / 2) * T + 6 * sigma * sqrt(T));
    
    if(N<1)
    {
        cout<<"Error! The step N must be greater than or equal to 1"<<endl;
        exit(1);
    }
    else
    {
        int tmp = 0;
        while(stockprice[N][tmp]>S_u) tmp++;
        for(int j=tmp;j<=N-tmp;j++)
            optionvalue[N][j] = stockprice[N][j]>K?(stockprice[N][j]-K):0;
        for(int j=0;j<tmp;j++)
            optionvalue[N][j] = 1/2 * (optionvalue[N][tmp] + optionvalue[N][tmp]);
        for(int j=N-tmp+1;j<N+1;j++)
            optionvalue[N][j] = 1/2 * (optionvalue[N][N-tmp] + optionvalue[N][N-tmp-1]);

        for(int i=N-1;i>=0;i--)
        {
            int temp = 0;
            while(stockprice[i][temp]>S_u) temp++;
            for(int j=temp;j<=i-temp;j++)
                optionvalue[i][j] = discount_factor * (pu * optionvalue[i+1][j] + pd * optionvalue[i+1][j+1]);
            for(int j=0;j<temp;j++)
                optionvalue[i][j] = 1/2 * (optionvalue[i][temp] + optionvalue[i][temp+1]);
            for(int j=i-temp+1;j<=i;j++)
                optionvalue[i][j] = 1/2 * (optionvalue[i][i-temp] + optionvalue[i][i-temp-1]);
        }
        
        option_price = optionvalue[0][0];
        return option_price;
    }
}

double CRRBTree::error_analysis(int N_)
{
    CRRBTree  p1(S,K,T,r,sigma,N_);
    return p1.optionpricing() - BlackScholesCall(S, K, sigma, T, r);
}

double CRRBTree::ratio_analysis(int N_)            // error analysis
{
    
    CRRBTree p1(S,K,T,r,sigma,N_);
    double P_2deltaT = p1.optionpricing();
    
    CRRBTree p2(S,K,T,r,sigma,2*N_);
    double P_deltaT = p2.optionpricing();
    
    double P_true = BlackScholesCall(S, K, sigma, T, r);
    
    double ratio = log2((P_2deltaT - P_true) / (P_deltaT - P_true));
    return ratio;
}

double CRRBTree::Richardson(int N_)
{
    CRRBTree p1(S,K,T,r,sigma,N_);
    double P_2deltaT = p1.optionpricing();
    
    CRRBTree p2(S,K,T,r,sigma,2*N_);
    double P_deltaT = p2.optionpricing();
    
    double coeff = ratio_analysis(N_);

    return (pow(2,coeff) * P_deltaT - P_2deltaT) / (pow(2,coeff) - 1);
}

/*---------------------------------------------Tian Binomial Tree Modle----------------------------------------*/


TianBTree::TianBTree()
{
    option_price = S = K = r = sigma = T = u = pu = d = pd = 0;
    N =0;
    cout<< "The binomial tree is NULL!"<<endl;
}


TianBTree::~TianBTree()
{
    stockprice.clear();
    optionvalue.clear();
}

TianBTree::TianBTree(double S_, double K_, double T_, double r_, double sigma_, int N_)
{
    S = S_; K = K_; T = T_; r = r_; sigma = sigma_; N = N_;
    
    delta_t = T / N;
    discount_factor = exp(-r * delta_t);
    
    double r_N = exp(r * delta_t);
    double v_N = exp(sigma * sigma * delta_t);
    u = r_N * v_N / 2 * (v_N + 1 + sqrt(v_N * v_N + 2 * v_N - 3));
    d = r_N * v_N / 2 * (v_N + 1 - sqrt(v_N * v_N + 2 * v_N - 3));
    pu = (exp(r * delta_t) - d) / (u - d);
    pd = 1 - pu;
    
    stockprice.resize(N+1);
    for(int i=0;i<N+1;i++)
        stockprice[i].resize(i+1);
    optionvalue.resize(N+1);
    for(int i=0;i<N+1;i++)
        optionvalue[i].resize(i+1);
    
    for(int i=0;i<N+1;i++)
        for(int j=0;j<=i;j++)
            stockprice[i][j] = pow(u,i-j) * pow(d,j) * S;
}

double TianBTree::optionpricing()
{
    if(N<1)
    {
        cout<<"Error! The step N must be greater than or equal to 1"<<endl;
        exit(1);
    }
    else
    {
        for(int i=0;i<N+1;i++)
            optionvalue[N][i] = stockprice[N][i]>K?(stockprice[N][i]-K):0;
        for(int i=N-1;i>=0;i--)
            for(int j=0;j<=i;j++)
                optionvalue[i][j] = discount_factor * (pu * optionvalue[i+1][j] + pd * optionvalue[i+1][j+1]);
        
        option_price = optionvalue[0][0];
        return option_price;
    }
}

double TianBTree::pruning_1()
{
    double S_u = S * exp((r - sigma * sigma / 2) * T + 6 * sigma * sqrt(T));
    double S_d = S * exp((r - sigma * sigma / 2) * T - 6 * sigma * sqrt(T));
    
    if(N<1)
    {
        cout<<"Error! The step N must be greater than or equal to 1"<<endl;
        exit(1);
    }
    else
    {
        for(int i=0;i<N+1;i++)
        {
            if(stockprice[N][i]<=S_u && stockprice[N][i]>=S_d)
                optionvalue[N][i] = stockprice[N][i]>K?(stockprice[N][i]-K):0;
            else
                optionvalue[N][i] = 0;
        }
        for(int i=N-1;i>=0;i--)
            for(int j=0;j<=i;j++)
            {
                if(stockprice[i][j]<=S_u && stockprice[i][j]>=S_d)
                    optionvalue[i][j] = discount_factor * (pu * optionvalue[i+1][j] + pd * optionvalue[i+1][j+1]);
                else
                    optionvalue[i][j] = 0;
            }
        
        double option_price = optionvalue[0][0];
        return option_price;
    }
}

double TianBTree::pruning_2()
{
    double S_u = S * exp((r - sigma * sigma / 2) * T + 6 * sigma * sqrt(T));
    double S_d = S * exp((r - sigma * sigma / 2) * T - 6 * sigma * sqrt(T));
    
    if(N<1)
    {
        cout<<"Error! The step N must be greater than or equal to 1"<<endl;
        exit(1);
    }
    else
    {
        for(int i=0;i<N+1;i++)
            optionvalue[N][i] = stockprice[N][i]>K?(stockprice[N][i]-K):0;
        for(int i=N-1;i>=0;i--)
            for(int j=0;j<=i;j++)
            {
                if(stockprice[i][j]<=S_u && stockprice[i][j]>=S_d)
                    optionvalue[i][j] = discount_factor * (pu * optionvalue[i+1][j] + pd * optionvalue[i+1][j+1]);
                else
                    optionvalue[i][j] = S > K*exp(-r*(T-i*delta_t))?(S-K*exp(-r*(T-i*delta_t))):0;
            }
        
        option_price = optionvalue[0][0];
        return option_price;
    }
}

double TianBTree::pruning_3()
{
    double S_u = S * exp((r - sigma * sigma / 2) * T + 6 * sigma * sqrt(T));
    
    if(N<1)
    {
        cout<<"Error! The step N must be greater than or equal to 1"<<endl;
        exit(1);
    }
    else
    {
        int tmp = 0;
        while(stockprice[N][tmp]>S_u) tmp++;
        for(int j=tmp;j<=N-tmp;j++)
            optionvalue[N][j] = stockprice[N][j]>K?(stockprice[N][j]-K):0;
        for(int j=0;j<tmp;j++)
            optionvalue[N][j] = 1/2 * (optionvalue[N][tmp] + optionvalue[N][tmp]);
        for(int j=N-tmp+1;j<N+1;j++)
            optionvalue[N][j] = 1/2 * (optionvalue[N][N-tmp] + optionvalue[N][N-tmp-1]);
        
        for(int i=N-1;i>=0;i--)
        {
            int temp = 0;
            while(stockprice[i][temp]>S_u) temp++;
            for(int j=temp;j<=i-temp;j++)
                optionvalue[i][j] = discount_factor * (pu * optionvalue[i+1][j] + pd * optionvalue[i+1][j+1]);
            for(int j=0;j<temp;j++)
                optionvalue[i][j] = 1/2 * (optionvalue[i][temp] + optionvalue[i][temp+1]);
            for(int j=i-temp+1;j<=i;j++)
                optionvalue[i][j] = 1/2 * (optionvalue[i][i-temp] + optionvalue[i][i-temp-1]);
        }
        
        option_price = optionvalue[0][0];
        return option_price;
    }
}

double TianBTree::error_analysis(int N_)
{
    TianBTree  p1(S,K,T,r,sigma,N_);
    return p1.optionpricing() - BlackScholesCall(S, K, sigma, T, r);
}

double TianBTree::ratio_analysis(int N_)
{
    
    TianBTree p1(S,K,T,r,sigma,N_);
    double P_2deltaT = p1.optionpricing();
    
    TianBTree p2(S,K,T,r,sigma,2*N_);
    double P_deltaT = p2.optionpricing();
    
    double P_true = BlackScholesCall(S, K, sigma, T, r);
    
    double ratio = log2(abs((P_2deltaT - P_true) / (P_deltaT - P_true)));
    
    return ratio;
}

double TianBTree::Richardson(int N_)
{
    TianBTree p1(S,K,T,r,sigma,N_);
    double P_2deltaT = p1.optionpricing();
    
    TianBTree p2(S,K,T,r,sigma,2*N_);
    double P_deltaT = p2.optionpricing();
    
    double coeff = ratio_analysis(N_);
    
    return (pow(2,coeff) * P_deltaT - P_2deltaT) / (pow(2,coeff) - 1);
}



/*---------------------------------------------Trinomial Tree----------------------------------------------------*/

TriTree::TriTree()
{
    option_price = S = K = r = sigma = T = u = pu = d = pd = m = pm = 0;
    N =0;
    cout<< "The trinomial tree is NULL!"<<endl;
}

TriTree::~TriTree()
{
    stockprice.clear();
    optionvalue.clear();
}

TriTree::TriTree(double S_, double K_, double T_, double r_, double sigma_, int N_)
{
    S = S_; K = K_; T = T_; r = r_; sigma = sigma_; N = N_;
    delta_t = T / N;
    discount_factor = exp(-r * delta_t);
    
    double delta_x = sigma * sqrt(2 * delta_t);
    double coeff1 = sigma * sigma * delta_t / pow(delta_x,2);
    double coeff2 = (r - sigma*sigma/2) * delta_t / delta_x;
    
    u = exp(delta_x);
    d = exp(-delta_x);
    m=1;
    pu = (coeff1 + pow(coeff2,2) + coeff2) / 2;
    pm = 1 - coeff1 - pow(coeff2,2);
    pd = (coeff1 + pow(coeff2,2) - coeff2) / 2;
    
    stockprice.resize(N+1);
    for(int i=0;i<N+1;i++)
        stockprice[i].resize(2*i+1);
    optionvalue.resize(N+1);
    for(int i=0;i<N+1;i++)
        optionvalue[i].resize(2*i+1);
    
    for(int i=0;i<N+1;i++)
    {
        for(int j=0;j<i;j++)
            stockprice[i][j] = S * pow(u,i-j);
        stockprice[i][i] = S;
        for(int j=i+1;j<2*i+1;j++)
            stockprice[i][j] = S * pow(d,j-i);
    }
}

double TriTree::optionpricing()
{
    assert(N>=1);
    
    for(int i=0;i<2*N+1;i++)
        optionvalue[N][i] = stockprice[N][i]>K?(stockprice[N][i]-K):0;
    for(int i=N-1;i>=0;i--)
        for(int j=0;j<2*i+1;j++)
            optionvalue[i][j] = discount_factor * (pu * optionvalue[i+1][j] + pm * optionvalue[i+1][j+1] + pd * optionvalue[i+1][j+2]);
    
    option_price = optionvalue[0][0];
    return option_price;
}

double TriTree::pruning_1()
{
    double S_u = S * exp((r - sigma * sigma / 2) * T + 6 * sigma * sqrt(T));
    double S_d = S * exp((r - sigma * sigma / 2) * T - 6 * sigma * sqrt(T));
    
    if(N<1)
    {
        cout<<"Error! The step N must be greater than or equal to 1"<<endl;
        exit(1);
    }
    else
    {
        for(int i=0;i<2*N+1;i++)
        {
            if(stockprice[N][i]<=S_u && stockprice[N][i]>=S_d)
                optionvalue[N][i] = stockprice[N][i]>K?(stockprice[N][i]-K):0;
            else
                optionvalue[N][i] = 0;
        }
        for(int i=N-1;i>=0;i--)
            for(int j=0;j<2*i+1;j++)
            {
                if(stockprice[i][j]<=S_u && stockprice[i][j]>=S_d)
                    optionvalue[i][j] = discount_factor * (pu * optionvalue[i+1][j] + pm * optionvalue[i+1][j+1] + pd * optionvalue[i+1][j+2]);

                else
                    optionvalue[i][j] = 0;
            }
        
        double option_price = optionvalue[0][0];
        return option_price;
    }

}

double TriTree::pruning_2()
{
    double S_u = S * exp((r - sigma * sigma / 2) * T + 6 * sigma * sqrt(T));
    double S_d = S * exp((r - sigma * sigma / 2) * T - 6 * sigma * sqrt(T));
    
    if(N<1)
    {
        cout<<"Error! The step N must be greater than or equal to 1"<<endl;
        exit(1);
    }
    else
    {
        for(int i=0;i<2*N+1;i++)
            optionvalue[N][i] = stockprice[N][i]>K?(stockprice[N][i]-K):0;
        for(int i=N-1;i>=0;i--)
            for(int j=0;j<2*i+1;j++)
            {
                if(stockprice[i][j]<=S_u && stockprice[i][j]>=S_d)
                    optionvalue[i][j] = discount_factor * (pu * optionvalue[i+1][j] + pm * optionvalue[i+1][j+1] + pd * optionvalue[i+1][j+2]);
                else
                    optionvalue[i][j] = S > K*exp(-r*(T-i*delta_t))?(S-K*exp(-r*(T-i*delta_t))):0;
            }
        
        option_price = optionvalue[0][0];
        return option_price;
    }

}

double TriTree::pruning_3()
{
    double S_u = S * exp((r - sigma * sigma / 2) * T + 6 * sigma * sqrt(T));
    
    if(N<1)
    {
        cout<<"Error! The step N must be greater than or equal to 1"<<endl;
        exit(1);
    }
    else
    {
        int tmp = 0;
        while(stockprice[N][tmp]>S_u) tmp++;
        for(int j=tmp;j<2*N+1-tmp;j++)
            optionvalue[N][j] = stockprice[N][j]>K?(stockprice[N][j]-K):0;
        for(int j=0;j<tmp;j++)
            optionvalue[N][j] = 1/2 * (optionvalue[N][tmp] + optionvalue[N][tmp]);
        for(int j=2*N+1-tmp;j<2*N+1;j++)
            optionvalue[N][j] = 1/2 * (optionvalue[N][N-tmp] + optionvalue[N][N-tmp-1]);
        
        for(int i=N-1;i>=0;i--)
        {
            int temp = 0;
            while(stockprice[i][temp]>S_u) temp++;
            for(int j=temp;j<2*i+1-temp;j++)
                optionvalue[i][j] = discount_factor * (pu * optionvalue[i+1][j] + pm * optionvalue[i+1][j+1] + pd * optionvalue[i+1][j+2]);
            for(int j=0;j<temp;j++)
                optionvalue[i][j] = 1/2 * (optionvalue[i][temp] + optionvalue[i][temp+1]);
            for(int j=2*i+1-temp;j<2*i+1;j++)
                optionvalue[i][j] = 1/2 * (optionvalue[i][2*i-temp] + optionvalue[i][2*i-temp-1]);
        }
        
        option_price = optionvalue[0][0];
        return option_price;
    }

}

double TriTree::error_analysis(int N_)
{
    TriTree  p1(S,K,T,r,sigma,N_);
    return p1.optionpricing() - BlackScholesCall(S, K, sigma, T, r);
}

double TriTree::ratio_analysis(int N_)
{
    
    TriTree p1(S,K,T,r,sigma,N_);
    double P_2deltaT = p1.optionpricing();
    
    TriTree p2(S,K,T,r,sigma,2*N_);
    double P_deltaT = p2.optionpricing();
    
    double P_true = BlackScholesCall(S, K, sigma, T, r);
    
    double ratio = log2((P_2deltaT - P_true) / (P_deltaT - P_true));
    return ratio;
}

double TriTree::Richardson(int N_)
{
    TriTree p1(S,K,T,r,sigma,N_);
    double P_2deltaT = p1.optionpricing();
    
    TriTree p2(S,K,T,r,sigma,2*N_);
    double P_deltaT = p2.optionpricing();
    
    double coeff = ratio_analysis(N_);
    
    return (pow(2,coeff) * P_deltaT - P_2deltaT) / (pow(2,coeff) - 1);
}


/*--------------------------------------Black Scholes Pricing Function-------------------------------------------*/

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


double auxillary_function(double x)   // define the auxiliary function to calculate Φ(x)
{
    double b0 = 0.2316419, b1 = 0.3193815300, b2 = -0.356563782, b3 = 1.7814779370, b4 = -1.821255978, b5 = 1.3302744290;
    double phi = ( 1 / pow(2*M_PI,0.5)) * exp(-x*x/2);
    double t = 1/(1+b0*x);
    return 1 - phi * (b1*t+b2*pow(t,2)+b3*pow(t,3)+b4*pow(t,4)+b5*pow(t,5));
}

