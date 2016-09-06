// HW11 problem1
//Huachen Qin

#include "HullWhiteTree.h"

// constructor function
HWTree::HWTree()
{
    alpha = 1;
    sigma = 0.2;
    T = 30;
    N = 30 * 12;
    delta_t = T / N;
    delta_x = sigma * sqrt(3*delta_t);
    jmax = ceil(0.184/(alpha*delta_t));
    
    tree.resize(N+1);
    for(int i=0;i<=jmax;i++)
        tree[i].resize(2*i+1);
    for(int i=jmax+1;i<N+1;i++)
        tree[i].resize(2*jmax+1);
    
    for(int i=0;i<=jmax;i++)
        for(int j=0;j<2*i+1;j++)
            tree[i][j] = (i-j) * delta_x;
    for(int i=jmax+1;i<N+1;i++)
        for(int j=0;j<2*jmax+1;j++)
            tree[i][j] = (jmax-j) * delta_x;
    
    rt.resize(N+1);
    for(int i=0;i<=jmax;i++)
        rt[i].resize(2*i+1);
    for(int i=jmax+1;i<N+1;i++)
        rt[i].resize(2*jmax+1);
    
    vector<double> bond = {0.0,0.000798,0.002443,0.004517,0.006260,0.007073,0.007792,0.009172,0.010695,0.012247,0.013741,0.014819,0.015820,0.016787,0.017704,0.018238,0.018746,0.019246,0.019738,0.020224,0.0575,0.020833,0.021032,0.021228,0.021422,0.021613,0.021802,0.021988,0.022172,0.022354,0.022533,0.022711,0.02886,0.023059,0.023231,0.023400,0.023568,0.023734,0.023899,0.024062,0.024223,0.024383,0.024542,0.024699,0.024855,0.025010,0.025164,0.025318,0.025470,0.025621,0.025771,0.025921,0.026070,0.026219,0.026367,0.026515,0.026662,0.026810,0.026957,0.027103,0.027214};
    
    r_tudi.resize(N+1);
    for(int i=0;i<=jmax;i++)
        r_tudi[i].resize(2*i+1);
    for(int i=jmax+1;i<N+1;i++)
        r_tudi[i].resize(2*jmax+1);
    
    bond_p.resize(N+1);
    gamma.resize(N+1);
    mt.resize(N+1);
    
    for(int i=0;i<bond.size();i++)
        for(int j=0;j<6;j++)
            bond_p[i*6+j] = (6-j) * bond[i] / 6 + j * bond[i+1] / 6;
    bond_p[N] = bond[bond.size()-1];
}


// deconstructor function
HWTree::~HWTree()
{
    tree.clear();
    r_tudi.clear();
    gamma.clear();
    mt.clear();
    bond_p.clear();
}

// get the process of gamma
void HWTree::get_gamma()
{
    gamma[0] = bond_p[1];
    r_tudi[0][0] = exp(-bond_p[0]*delta_t);
    
    for(int step=1;step<=jmax;step++)
    {
        for(int j=0;j<step;j++)
        {
            r_tudi[step][j] = 0;
            for(int k=0;k<=j;k++)
                r_tudi[step][j] += r_tudi[step-1][k] * normal_prob(step-1,k)[j-k] * exp(-(gamma[step-1]+tree[step-1][k])*delta_t);
        }
        
        if(step==1)
            r_tudi[1][1] = r_tudi[0][0] * normal_prob(0,0)[1] * exp(-(gamma[0] + tree[0][0]) * delta_t);
        else
        {
            r_tudi[step][step] = 0;
            for(int k=0;k<3;k++)
                r_tudi[step][step] += r_tudi[step-1][step-2+k] * normal_prob(step-1,step-2+k)[2-k] * exp(-(gamma[step-1] + tree[step-1][step-2+k])*delta_t);
        }
        
        for(int j=2*step;j>step;j--)
        {
            r_tudi[step][j] = 0;
            for(int k=0;k<=2*step-j;k++)
                r_tudi[step][j] += r_tudi[step-1][j-2+k] * normal_prob(step-1,j-2+k)[2-k] * exp(-(gamma[step-1]+tree[step-1][j-2+k])*delta_t);
        }
        
        double temp = 0;
        for(int i=0;i<2*step+1;i++)
            temp += r_tudi[step][i] * exp(-tree[step][i]*delta_t);
        
        gamma[step] = log(temp)/delta_t + bond_p[step+1] * (step+1);
        
        
    }
    
    for(int step=jmax+1;step<N;step++)
    {
        
        r_tudi[step][0] = r_tudi[step-1][0] * top_prob(0)[0] * exp(-(gamma[step-1] + tree[step-1][0]) * delta_t);
        r_tudi[step][0] += r_tudi[step-1][1] * normal_prob(jmax,1)[0] * exp(-(gamma[step-1] + tree[step-1][1]) * delta_t);
        
        r_tudi[step][1] = r_tudi[step-1][0] * top_prob(0)[1] * exp(-(gamma[step-1] + tree[step-1][0]) * delta_t);
        r_tudi[step][1] += r_tudi[step-1][1] * normal_prob(jmax,1)[1] * exp(-(gamma[step-1] + tree[step-1][1]) * delta_t);
        r_tudi[step][1] += r_tudi[step-1][2] * normal_prob(jmax,2)[0] * exp(-(gamma[step-1] + tree[step-1][2]) * delta_t);
        
        
        r_tudi[step][2] = r_tudi[step-1][0] * top_prob(0)[2] * exp(-(gamma[step-1] + tree[step-1][0]) * delta_t);
        for(int i=1;i<4;i++)
            r_tudi[step][2] += r_tudi[step-1][i] * normal_prob(jmax,i)[3-i] * exp(-(gamma[step-1] + tree[step-1][i]) * delta_t);
        
        r_tudi[step][3] = 0;
        for(int i=0;i<3;i++)
            r_tudi[step][3] += r_tudi[step-1][2+i] * normal_prob(jmax,2+i)[2-i] * exp(-(gamma[step-1] + tree[step-1][2+i]) * delta_t);
        
        r_tudi[step][4] = r_tudi[step-1][6] * bot_prob(6)[0] * exp(-(gamma[step-1] + tree[step-1][6]) * delta_t);
        for(int i=0;i<3;i++)
            r_tudi[step][4] += r_tudi[step-1][3+i] * normal_prob(jmax,3+i)[2-i] * exp(-(gamma[step-1] + tree[step-1][3+i]) * delta_t);
        
        r_tudi[step][5] = r_tudi[step-1][4] * normal_prob(jmax,4)[2] * exp(-(gamma[step-1] + tree[step-1][4]) * delta_t);
        r_tudi[step][5] += r_tudi[step-1][5] * normal_prob(jmax,5)[1] * exp(-(gamma[step-1] + tree[step-1][5]) * delta_t);
        r_tudi[step][5] += r_tudi[step-1][6] * bot_prob(6)[1] * exp(-(gamma[step-1] + tree[step-1][6]) * delta_t);
        
        r_tudi[step][6] = r_tudi[step-1][5] * normal_prob(jmax,5)[2] * exp(-(gamma[step-1] + tree[step-1][5]) * delta_t);
        r_tudi[step][6] += r_tudi[step-1][6] * bot_prob(2*jmax)[2] * exp(-(gamma[step-1] + tree[step-1][6]) * delta_t);
        
        double temp = 0;
        for(int i=0;i<2*jmax+1;i++)
            temp += r_tudi[step][i] * exp(- tree[step][i] * delta_t);
        
        gamma[step] = log(temp)/delta_t + (step+1) * bond_p[step+1];
    }
}

// get process of m_delta_t
void HWTree::get_mt()
{
    mt[0] = 0;
    for(int i=1;i<N;i++)
        mt[i] = (gamma[i+1] - gamma[i]) / (alpha * delta_t) + gamma[i];
}

// get the calibrated interate tree rt
void HWTree::get_rt()
{
    for(int i=0;i<=jmax;i++)
        for(int j=0;j<2*i+1;j++)
            rt[i][j] = tree[i][j] + gamma[i];
    for(int i=jmax+1;i<N;i++)
        for(int j=0;j<2*jmax+1;j++)
            rt[i][j] = tree[i][j] + gamma[i];
}

// print gamma
void HWTree::print_gamma()
{
    cout<<"Here are the results of gamma:"<<endl;
    for(int i=0;i<N;i++)
        cout<<"i = "<<i<<", gamma = "<<gamma[i]<<endl;
    cout<<"------------------------------------"<<endl;
    cout<<endl;
}


// print the calibrated tree
void HWTree::print_tree()
{
    cout<<"After we construct the Hull White Tree and Calculate Gamma, we can calibrate interest rate tree:"<<endl;
    cout<<endl;
    for(int i=0;i<=jmax;i++)
    {
        cout<<"------------------------------------"<<endl;
        cout<<"step = "<<i<<endl;
        for(int j=0;j<2*i+1;j++)
            cout<<"Node "<<j<<" : "<<rt[i][j]<<endl;
    }
    
    for(int i=jmax+1;i<N;i++)
    {
        cout<<"------------------------------------"<<endl;
        cout<<"step = "<<i<<endl;
        for(int j=0;j<2*jmax+1;j++)
            cout<<"Node "<<j<<" : "<<rt[i][j]<<endl;
    }
}


vector<double> HWTree::normal_prob(int step, int j_)
{
    int j = step - j_;
    vector<double> prob;
    prob.resize(3);
    prob[0] = 1.0/6 + (pow(alpha*j*delta_t,2) - alpha * j * delta_t) / 2;
    prob[1] = 2.0/3 - pow(alpha*j*delta_t,2);
    prob[2] = 1.0/6 + (pow(alpha*j*delta_t,2) + alpha * j * delta_t) / 2;
    return prob;
}

vector<double> HWTree::top_prob(int j_)
{
    int j = jmax - j_;
    vector<double> prob;
    prob.resize(3);
    prob[0] = 7.0/6 + (pow(alpha*j*delta_t,2) - 3 * alpha * j * delta_t) / 2;
    prob[1] = -1.0/3 - pow(alpha*j*delta_t,2) + 2 * alpha * j * delta_t;
    prob[2] = 1.0/6 + (pow(alpha*j*delta_t,2) - alpha * j * delta_t) / 2;
    return prob;
}

vector<double> HWTree::bot_prob(int j_)
{
    int j= jmax - j_;
    vector<double> prob;
    prob.resize(3);
    prob[0] = 1.0/6 + (pow(alpha*j*delta_t,2) + alpha * j * delta_t) / 2;
    prob[1] = -1.0/3 - pow(alpha*j*delta_t,2) - 2 * alpha * j * delta_t;
    prob[2] = 7.0/6 + (pow(alpha*j*delta_t,2) + 3 * alpha * j * delta_t) / 2;
    return prob;
}


