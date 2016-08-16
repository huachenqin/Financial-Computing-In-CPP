// Final Project
// Huachen Qin

#include "matrix.h"
#include <iostream>

using namespace std;

Matrix::Matrix() {}

Matrix::Matrix(const Matrix &m)
{
    numRows = m.numRows;
    numColumns = m.numColumns;
    elements.resize(numRows);
    for(int i=0;i<numRows;i++)
        elements[i].resize(numColumns);
    for(int i=0;i<numRows;i++)
        for(int j=0;j<numColumns;j++)
            elements[i][j] = m.elements[i][j];
}

Matrix::Matrix(int rows, int columns)
{
    numRows = rows;
    numColumns = columns;
    elements.resize(numRows);
    for(int i=0;i<numRows;i++)
        elements[i].resize(numColumns);
    
    for(int i=0;i<numRows;i++)
        for(int j=0;j<numColumns;j++)
            elements[i][j] = 0;
}


void Matrix::subtractRow(int i, int j, double d)
{
    for(int k=0;k<numColumns;k++)
        elements[i][k] -=  elements[j][k] * d;
}

void Matrix::swapRows(int i, int j)
{
    for(int k=0;k<numColumns;k++)
    {
        elements[i][k] = elements[i][k] + elements[j][k];
        elements[j][k] = elements[i][k] - elements[j][k];
        elements[i][k] = elements[i][k] - elements[j][k];
    }
    
}

void Matrix::multiplyRow(int i, double d)
{
    for(int k=0;k<numColumns;k++)
        elements[i][k] *= d;
}

Matrix &Matrix::operator=(const Matrix&m)
{
    if(numColumns==m.numColumns && numRows==m.numRows)
    {
        for(int i=0;i<numRows;i++)
            for(int j=0;j<numColumns;j++)
                elements[i][j] = m.elements[i][j];
        numRows = m.numRows;
        numColumns = m.numColumns;
    }
    else
        cout<<"Error!The dimension does not match!"<<endl;
    
    return *this;
}

Matrix Matrix::operator+(const Matrix &m)
{
    assert(numRows==m.numRows && numColumns==m.numColumns);
    Matrix n(*this);
    {
        for(int i=0;i<numRows;i++)
            for(int j=0;j<numColumns;j++)
                n.elements[i][j] += m.elements[i][j];
    }
    return n;
}

Matrix Matrix::operator-(const Matrix &m)
{
    assert(numRows==m.numRows && numColumns==m.numColumns);
    Matrix n(*this);
    {
        for(int i=0;i<numRows;i++)
            for(int j=0;j<numColumns;j++)
                n.elements[i][j] -= m.elements[i][j];
    }
    return n;
}

Matrix Matrix::operator*(const Matrix &m)
{
    assert(numColumns==m.numRows);
    Matrix n(numRows,m.numColumns);
    
    for(int i=0;i<n.numRows;i++)
        for(int j=0;j<n.numColumns;j++)
        {
            for(int k=0;k<n.numRows;k++)
            {
                n.elements[i][j] += elements[i][k] * m.elements[k][j];
            }
        }
    return n;
}

Matrix Matrix::operator*(double d)
{
    Matrix m = *this;
    for(int i=0;i<numRows;i++)
        for(int j=0;j<numColumns;j++)
            m.elements[i][j] *= d;
    return m;
}

vector<double> Matrix::operator[](int i)
{
    return elements[i];
}

Matrix Matrix::transpose()
{
    Matrix m(numColumns,numRows);
    for(int i=0;i<m.numRows;i++)
        for(int j=0;j<m.numColumns;j++)
        {
            if(i!=j)
            {
                m.elements[i][j] = elements[j][i];
            }
            else
                m.elements[i][j] = elements[i][j];
        }
    return m;
}

Matrix Matrix::choleskyDecompose()
{
    Matrix m(numRows,numColumns);
    if(numRows == numColumns)
    {
        double * L = new double [numRows*numColumns];
        for(int m=0;m<numRows;m++)
            for(int n=0;n<numColumns;n++)
                L[m*numRows+n] = elements[m][n];
        for(int i=0;i<numRows;i++)
            for(int j=0;j<(i+1);j++)
            {
                double temp = 0;
                for(int k=0;k<j;k++)
                    temp += m.elements[i][k] * m.elements[j][k];
                m.elements[i][j] = (i==j)?sqrt(L[i*numRows+i]-temp):(1.0/m.elements[j][j]*(L[i*numRows+j]-temp));
            }
        delete [] L;
    }
    else
        cout<<"Error! The dimension of the matrix does not match!"<<endl;
    return m;
}

void Matrix::luDecompose()
{
    assert(numRows==numColumns);
    int n = numColumns;
    Matrix a(*this);
    Matrix L(numRows,numColumns), U(numRows,numColumns);
    
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<n;j++)
        {
            if(j<i)
                U.elements[i][j] = 0.0;
            else
            {
                U.elements[i][j] = a.elements[i][j];
                for(int k=0;k<i;k++)
                    U.elements[i][j] -= U.elements[k][j] * L.elements[i][k];
            }
        }
        
        for(int j=0;j<n;j++)
        {
            if(j<i)
                L.elements[j][i] = 0.0;
            else if(j==i)
                L.elements[j][i] = 1.0;
            else
            {
                L.elements[j][i] = a.elements[j][i] / U.elements[i][i];
                for(int k=0;k<i;k++)
                    L.elements[j][i] -= L.elements[j][k] * U.elements[k][i] / U.elements[i][i];
            }
        }

        
    }
    
    cout<<"The LU Decomposition of the matrix is:"<<endl;
    cout<<"L = "<<endl; cout<<L<<endl;
    cout<<"U = "<<endl; cout<<U<<endl;
}

void Matrix::inverse()									//inverts the matrix
{
    if(numRows==numColumns)
    {
        Matrix temp(numRows, 2*numColumns);
        for (int i =0;i<numRows; i++)
        {
            temp.elements[i][i+numColumns] = 1;
            for (int j = 0; j < numColumns; j++)
                temp.elements[i][j] = elements[i][j];
        }
        for (int j=0; j < numColumns; j++)
        {
            temp.multiplyRow(j, 1/temp.elements[j][j]);
            for (int i = j+1; i < numRows; i++)
                temp.subtractRow(i, j, temp.elements[i][j]);
        }
        
        for(int j=numColumns-1;j>=0;j--)
        {
            for (int i = j-1;i>= 0;i--)
                temp.subtractRow(i, j, temp.elements[i][j]);
        }
        
        for (int i =0;i<numRows;i++)
        {
            for (int j = numColumns; j < 2*numColumns; j++)
                elements[i][j-numColumns] = temp.elements[i][j];
        }}
    else
        cout<<"Error! The dimension does not match!"<<endl;
}

Matrix GSmethod(Matrix &A,Matrix &b)
{
    double error = 0.0001; int m = A.numRows, n = A.numColumns;
    assert(m==n);
    
    Matrix x(n,1), L(m,n), U(m, n);
    for(int i=0;i<m;i++)
    {
        for(int j=0;j<n;j++)
        {
            if(j<=i)
                L.elements[i][j] = A.elements[i][j];
            else
                U.elements[i][j] = A.elements[i][j];
        }
    }
    
    for(int i=0;i<n;i++) x.elements[i][0] = 1;
    L.inverse(); Matrix T(m,n), C(m,1),temp(x), errormatrix(x);
    C = L * b;
    T = L * U;  T = T * (-1);
    
    do
    {
        x = temp;
        temp = T * x;
        temp = temp + C;
        errormatrix = temp -x;
    } while(maxelement(errormatrix)>=error);
    
    return x;
}

ostream& operator<<(ostream& os,const Matrix& m)
{
    const char separator = ' ';               // initialize table
    const int tablewidth = 10;
    for(int i=0;i<m.numRows;i++)
    {
        for(int j=0;j<m.numColumns;j++)
            os<<left<<setw(tablewidth)<<setfill(separator)<<m.elements[i][j];
        os<<endl;
    }
    return os;
}

void Matrix::setmatrix(vector<double> arr)
{
    for(int i=0;i<numRows;i++)
        for(int j=0;j<numColumns;j++)
            elements[i][j] = arr[i*numColumns+j];
    
}

void Matrix::setelements(int i, int j,double value)
{
    elements[i][j] = value;
}

double maxelement(Matrix &m)
{
    double temp = abs(m.elements[0][0]);
    for(int i=0;i<m.numRows;i++)
        for(int j=0;j<m.numColumns;j++)
            temp = abs(m.elements[i][j])>temp?abs(m.elements[i][j]):temp;
    return temp;
}

double Matrix::getelements(int i, int j)
{
    return elements[i][j];
}