// HW8 problem1
// Huachen Qin

#include "matrix.h"

int main()
{
    cout<<"Now Let's test our Matrix Class:"<<endl;
    cout<<endl;
    
    cout<<"create matrix m:"<<endl;
    Matrix m(3,3);
    vector<double> v1 = {1,2,3,0,1,4,5,6,0};
    m.setmatrix(v1);
    cout<<m<<endl;
    
    cout<<"create matrix n:"<<endl;
    Matrix n(3,3);
    vector<double> v2 = {25,15,-5,15,18,0,-5,0,11};
    n.setmatrix(v2);
    cout<<n<<endl;
    
    cout<<"Let's test some matrix calculations:"<<endl;
    cout<<"m + n = "<<endl; Matrix test1 = m+n;
    cout<<test1<<endl; m.setmatrix(v1);
    cout<<"m * n = "<<endl; Matrix test2 = m*n;
    cout<<test2<<endl;
    cout<<"the transpose of m is:"<<endl; Matrix test3 = m.transpose();
    cout<<test3<<endl;
    cout<<"the cholesky decomposition of n is:"<<endl;
    cout<<n.choleskyDecompose()<<endl;
    cout<<"the inverse of m is:"<<endl;
    m.inverse();cout<<m<<endl;
    
    cout<<"Now Let's test the LU Decomposition and Gauss-Seidel method: "<<endl;
    vector<double> v3 = {1,1,0,3,2,1,-1,1,3,-1,-1,2,-1,2,3,-1};
    Matrix E(4,4); E.setmatrix(v3);
    cout<<"create matrix E:"<<endl;
    cout<<E<<endl;
    E.luDecompose();

    cout<<"Now we have linear equations A * x = b, "<<endl;
    vector<double> v4 = {3,1,0,0,5,2,2,3,5}, v5 = {9,17,18};
    Matrix A(3,3), b(3,1);
    A.setmatrix(v4); b.setmatrix(v5);
    cout<<"A = "<<endl; cout<<A<<endl;
    cout<<"b = "<<endl; cout<<b<<endl;
    cout<<"the root of the linear equations is x = "<<endl;
    cout<<GSmethod(A,b)<<endl;
    
    cout<<"END of PROGEAM!"<<endl;
    
}

