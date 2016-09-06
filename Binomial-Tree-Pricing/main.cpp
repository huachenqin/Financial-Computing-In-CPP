// HW10 problem1
// Huachen Qin

#include "BinomialTree.h"

int main()
{
    CRRBTree v1(100,100,1,0.01,0.2,100);
    TianBTree v2(100,100,1,0.01,0.2,100);
    TriTree v3(100,100,1,0.01,0.2,100);
    
    cout<<"The option prices from CRR binomial tree model are:"<<endl;
    cout<<"Case 1, without pruning:"<<endl;
    cout<<v1.optionpricing()<<endl;
    cout<<"Case 2, with pruning method 1:"<<endl;
    cout<<v1.pruning_1()<<endl;
    cout<<"Case 3, with pruning method 2:"<<endl;
    cout<<v1.pruning_2()<<endl;
    cout<<"Case 4, with pruning method 3"<<endl;
    cout<<v1.pruning_3()<<endl;
    cout<<"Case 5, use Richardson Extrapolation Method:"<<endl;
    cout<<v1.Richardson(100)<<endl;
    cout<<"Now let us run the error analysis, N = 100, 200, 300, 400, 500"<<endl;
    for(int i=1;i<6;i++)
        cout<<"When N= "<<i*100<<" , the error = "<<v1.error_analysis(i*100)<<", the speed of convergence is: 2^k = "<<pow(2,v1.ratio_analysis(i*100))<<endl;
    cout<<endl;
    
    cout<<"The option prices from Tian binomial tree model are:"<<endl;
    cout<<"Case 1, without pruning:"<<endl;
    cout<<v2.optionpricing()<<endl;
    cout<<"Case 2, with pruning method 1:"<<endl;
    cout<<v2.pruning_1()<<endl;
    cout<<"Case 3, with pruning method 2:"<<endl;
    cout<<v2.pruning_2()<<endl;
    cout<<"Case 4, with pruning method 3"<<endl;
    cout<<v2.pruning_3()<<endl;
    cout<<"Case 5, use Richardson Extrapolation Method:"<<endl;
    cout<<v2.Richardson(50)<<endl;
    cout<<"Now let us run the error analysis, N=100, 200, 300, 400, 500"<<endl;
    for(int i=1;i<6;i++)
        cout<<"When N= "<<i*100<<" , the error = "<<v2.error_analysis(i*100)<<", the speed of convergence is: 2^k = "<<pow(2,v2.ratio_analysis(i*100))<<endl;
    cout<<endl;
    
    cout<<"The option prices from trinomial tree model are:"<<endl;
    cout<<"Case 1, without pruning:"<<endl;
    cout<<v3.optionpricing()<<endl;
    cout<<"Case 2, with pruning method 1:"<<endl;
    cout<<v3.pruning_1()<<endl;
    cout<<"Case 3, with pruning method 2:"<<endl;
    cout<<v3.pruning_2()<<endl;
    cout<<"Case 4, with pruning method 3"<<endl;
    cout<<v3.pruning_3()<<endl;
    cout<<"Case 5, use Richardson Extrapolation Method:"<<endl;
    cout<<v3.Richardson(200)<<endl;
    cout<<"Now let us run the error analysis, N=100, 200, 300, 400, 500"<<endl;
    for(int i=1;i<6;i++)
        cout<<"When N= "<<i*100<<" , the error = "<<v3.error_analysis(i*100)<<", the speed of convergence is: 2^k = "<<pow(2,v3.ratio_analysis(i*100))<<endl;
    cout<<endl;
    
    cout<<"From the output above, we can observe that:"<<endl;
    cout<<"1. As delta T decreases (or N increases), the error decreases for the CRR, Tian and Trinomial Tree Model"<<endl;
    cout<<"2. For the speed of convergence, we use 2^k to measure the speed. We can observe that CRR model and Trinomial Model has stable speed, which are around 1.98~2.00; For Tian Model, the convergence speed is not stable (around 0.5~2.5)"<<endl;
    cout<<"3. All the three methods can efficient approximate the value of European Call Option."<<endl;
}