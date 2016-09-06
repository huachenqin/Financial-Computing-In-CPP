#ifndef MATRIX_H
#define MATRIX_H


#include <iostream>
#include <ctime>
#include<cmath>
#include<iomanip>
#include<vector>
#include<assert.h>

using namespace std;

class Matrix{
private:
    int numRows;
    int numColumns;
    //helper functions - useful for finding the inverse
    void subtractRow(int i, int j, double d);					//Subtracts row j from row i
    void swapRows(int i, int j);					//Swaps rows i and j
    void multiplyRow(int i, double d);				//Multiplies row i by constant d
    
    vector<vector<double>>  elements;								//the array
public:
    Matrix();
    Matrix(const Matrix &m);						//copy constructor
    Matrix(int rows, int columns);					//constructor
    ~Matrix() {}								//destructor
    
    Matrix& operator=(const Matrix& m);				//assignment operator
    Matrix operator + (const Matrix& m);			//addition operator
    Matrix operator-(const Matrix &m);
    Matrix operator *(const Matrix& m);				//matrix multiplication
    Matrix operator *(double d);					//multiplication by const
    vector<double> operator[](int i);						//overloaded operator [] - note the return type
    
    void setmatrix(vector<double> v);
    void setelements(int i,int j,double value);
    
    void inverse();									//inverts the matrix
    Matrix transpose();								//returns a tranposed matrix
    Matrix choleskyDecompose();						//finds and returns a Cholecky Decomposition of the matrix
    void luDecompose();          // finds and returns a LU Decomposition of the matrix
    
    friend Matrix GSmethod(Matrix &A, Matrix &b);   //solves linear equations using Gauss-Seidel method
    friend double maxelement(Matrix &m);
    friend ostream& operator<<(ostream& os,const Matrix& m);//Stream "<<" overloading
};

const Matrix EMPTY_MATRIX(0, 0);					//empty matrix- in case of error

#endif