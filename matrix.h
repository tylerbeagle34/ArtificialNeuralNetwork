/***************************************************************************
 *   Copyright (C) 2005 by Dr.Homayoun Valafar                             *
 *   homayoun@cse.sc.edu                                                   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#ifndef MATRIX_H
#define MATRIX_H
#include <iostream>
#include "matrix.h"
#include <math.h>
#include <stdlib.h>


#define SIGN(a,b) ((b) > 0.0 ? fabs(a) : -fabs(a))
#define SQR(a) (((a)) == 0.0 ? 0.0 : a*a)
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define PYTH(a,b) sqrt(a*a+b*b)
#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);\
			    a[k][l] = h + s * (g - h * tau);

using namespace std;
class Matrix;

/**
  This class provides the basic matrix data type and operations.

  @author Homayoun Valafar
 */
class Matrix {
    friend ostream &operator<<(ostream &, const Matrix &);
    friend istream &operator>>(istream &, Matrix &);
public:

    static string svn_id() {
        return "$Id: matrix.h 418 2012-10-25 20:48:10Z siminm $";
    }

    Matrix(int = 3, int = 3, double = 0.0);
    Matrix(const Matrix &);
    ~Matrix();

    void Set_Mij(int i, int j, double x) {
        M[i][j] = x;
    }

    double Get_Mij(int i, int j) const {
        return M[i][j];
    }

    int GetRow() const {
        return row;
    }

    int GetCol() const {
        return col;
    }

    const Matrix & operator=(const Matrix &);
    Matrix operator+(const Matrix &);
    Matrix operator*(const Matrix &);
    Matrix Transpose();
    bool operator==(const Matrix &);
    bool operator!=(const Matrix &);
    int Diag(Matrix &, Matrix &);
    double Det();
    bool isSquare();
    bool isSym();
    double Trace();

    void makeIdentity();
    bool isIdentity();
    void AppendRows(Matrix &N);
    void Copy(Matrix &N);
    Matrix * canon();

    void Print();

    double **M;
protected:
    int row, col;
};

#endif
