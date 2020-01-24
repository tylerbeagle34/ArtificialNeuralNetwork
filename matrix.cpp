/***************************************************************************
 *   Copyright (C) 2005 by Dr. Homayoun Valafar                            *
 *   homayoun@cse.sc.edu                                                   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILI.MTY or FITNESS FOR A PARTICULAR PURPOSE.  See the       *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#include "matrix.h"
using namespace std;

#define SMALL_NUMBER 1e-6

Matrix::Matrix(int row_size, int col_size, double Init_value)
: row(row_size), col(col_size) {
	M = new double *[row];

	for (int i = 0; i < row; i++) {
		M[i] = new double [col];
		for (int j = 0; j < col; j++)
			M[i][j] = Init_value;
	}
}

Matrix::Matrix(const Matrix & m1)
: row(m1.row), col(m1.col) {
	M = new double *[row];

	for (int i = 0; i < row; i++) {
		M[i] = new double [col];
		for (int j = 0; j < col; j++)
			M[i][j] = m1.M[i][j];
	}
}

Matrix::~Matrix() {
	for (int i = 0; i < row; i++) {
		delete [] M[i];
	}

	delete [] M;
}


ostream & operator<<(ostream & s, const Matrix & rhs) {

	//s.setf(ios::scientific);
	//s.precision(10);
	for (int i = 0; i < rhs.row; i++) {
		for (int j = 0; j < rhs.col; j++)
			s << rhs.M[i][j] << " ";
		s << endl;
	}
	return s;
}

istream & operator>>(istream & s, Matrix & rhs) {

	for (int i = 0; i < rhs.row; i++) {
		for (int j = 0; j < rhs.col; j++)
			s >> rhs.M[i][j];
	}
	return s;
}

Matrix Matrix::operator+(const Matrix & rhs) {
	//Matrix *temp_ptr = new Matrix(*this);
	Matrix temp(*this);

	if ((row != rhs.row) || (col != rhs.col)) {
		cout << "Miss-matched dimensions in + operator.\n";
		exit(1);
	}

	for (int i = 0; i < row; i++)
		for (int j = 0; j < col; j++)
			temp.M[i][j] = M[i][j] + rhs.M[i][j];

	return temp;
}

const Matrix & Matrix::operator=(const Matrix & rhs) {

	if ((row != rhs.row) || (col != rhs.col)) {
		cout << "Miss-matched dimensions in = operator. " << row << " vs. " << rhs.row << " and " << col << " vs. " << rhs.col << endl;
		exit(1);
	}

	for (int i = 0; i < row; i++) {
		for (int j = 0; j < col; j++)
			M[i][j] = rhs.M[i][j];
	}
	return *this;
}


Matrix Matrix::operator*(const Matrix & rhs) {

	Matrix temp(row, rhs.col);

	if (col != rhs.row) {
		cerr << "Miss-matched dimensions in * operator: (" << row << "x" << col << ") * ("<<rhs.row << "x" << rhs.col << ")" << endl;
		exit(1);
	}

	for (int i = 0; i < row; i++) {
		for (int j = 0; j < rhs.col; j++) {
			double s = 0.0;
			for (int k = 0; k < col; k++)
				s += M[i][k] * rhs.M[k][j]; //s.m[i][j] = 0 at initiation
			temp.M[i][j] = s;
		}
	}

	return temp;
}

bool Matrix::operator==(const Matrix & rhs) {

	if ((row != rhs.row) || (col != rhs.col)) {
		cout << "Miss-matched dimensions in == operator.\n";
		exit(1);
	}

	for (int i = 0; i < row; i++)
		for (int j = 0; j < col; j++)
			if (M[i][j] != rhs.M[i][j])
				return (false);
	return (true);
}

bool Matrix::operator !=(const Matrix & rhs) {
	return (!(*this == rhs));
}

Matrix Matrix::Transpose() {
	//Matrix *temp;
	//temp = new Matrix(col,row,0);

	Matrix temp(col, row, 0);

	for (int i = 0; i < temp.row; i++)
		for (int j = 0; j < temp.col; j++)
			temp.M[i][j] = M[j][i];

	return temp;
}

int Matrix::Diag(Matrix &d, Matrix &v) {
	int j, iq, ip, i, nrot, n;
	double tresh, theta, tau, t, sm, s, h, g, c, *b, *z;
	double a[row + 1][col + 1];

	if (row != col) {
		cerr << "Error: Matrix is not square\n";
		return -1;
	}

	for (i = 0; i < row; i++)
		for (j = 0; j < i; j++)
			if (fabs(M[i][j] - M[j][i]) > SMALL_NUMBER) {
				cerr << "Error: Matrix not symmetric" << endl;
				cerr << M[i][j] << " != " << M[j][i] << "diff: " << fabs(M[i][j] - M[j][i]) << endl;
				cerr << (*this) << endl;
				exit(0);
				return -1;
			}

	n = row;
	for (i = 0; i < n; i++)
		for (j = 0; j < n; j++)
			a[i + 1][j + 1] = M[i][j];

	b = (double *) malloc((n + 1) * sizeof (double));
	if (b == NULL) {
		cerr << "Error allocating b." << endl;
		exit(1);
	}

	z = (double *) malloc((n + 1) * sizeof (double));
	if (z == NULL) {
		cout << "Error allocating z." << endl;
		exit(1);
	}

	for (ip = 1; ip <= n; ip++) {
		for (iq = 1; iq <= n; iq++)
			v.M[ip - 1][iq - 1] = 0.0;
		v.M[ip - 1][ip - 1] = 1.0;
	}

	for (ip = 1; ip <= n; ip++) {
		b[ip] = d.M[ip - 1][ip - 1] = a[ip][ip];
		z[ip] = 0.0;
	}

	nrot = 0;
	for (i = 1; i <= 50; i++) {
		sm = 0.0;
		for (ip = 1; ip <= n - 1; ip++) {
			for (iq = ip + 1; iq <= n; iq++)
				sm += fabs(a[ip][iq]);
		}
		if (sm == 0.0) {
			free(z);
			free(b);
			return nrot;
		}
		if (i < 4)
			tresh = 0.2 * sm / (n * n);

		else
			tresh = 0.0;
		for (ip = 1; ip <= n - 1; ip++) {
			for (iq = ip + 1; iq <= n; iq++) {
				g = 100.0 * fabs(a[ip][iq]);
				if (i > 4
					&& (double) (fabs(d.M[ip - 1][ip - 1]) + g) == (double) fabs(d.M[ip - 1][ip - 1])
					&& (double) (fabs(d.M[iq - 1][iq - 1]) + g) == (double) fabs(d.M[iq - 1][iq - 1]))
					a[ip][iq] = 0.0;

				else if (fabs(a[ip][iq]) > tresh) {
					h = d.M[iq - 1][iq - 1] - d.M[ip - 1][ip - 1];
					if ((double) (fabs(h) + g) == (double) fabs(h))
						t = (a[ip][iq]) / h;
					else {
						theta = 0.5 * h / (a[ip][iq]);
						t =
							1.0 / (fabs(theta) +
							sqrt(1.0 + theta * theta));
						if (theta < 0.0)
							t = -t;
					}
					c = 1.0 / sqrt(1 + t * t);
					s = t * c;
					tau = s / (1.0 + c);
					h = t * a[ip][iq];
					z[ip] -= h;
					z[iq] += h;
					d.M[ip - 1][ip - 1] -= h;
					d.M[iq - 1][iq - 1] += h;
					a[ip][iq] = 0.0;
					for (j = 1; j <= ip - 1; j++) {
						ROTATE(a, j, ip, j, iq);
					}
					for (j = ip + 1; j <= iq - 1; j++) {
						ROTATE(a, ip, j, j, iq);
					}
					for (j = iq + 1; j <= n; j++) {
						ROTATE(a, ip, j, iq, j);
					}
					for (j = 1; j <= n; j++) {
						ROTATE(v.M, j - 1, ip - 1, j - 1, iq - 1);
					}
					++nrot;
				}
			}
		}
		for (ip = 1; ip <= n; ip++) {
			b[ip] += z[ip];
			d.M[ip - 1][ip - 1] = b[ip];
			z[ip] = 0.0;
		}
	}
	cout << "Too many iterations in routine Diag without convergence" << endl;
	exit(1);
}

double Matrix::Det() {
	if (!isSquare()) {
		cerr << "Taking a Determinant of a non square matrix." << endl;
		exit(1);
	}
	if (row == 1) return M[0][0];
	if (row == 2) return M[0][0] * M[1][1] - M[1][0] * M[0][1];
	if (row == 3) {
		return M[0][0]*(M[1][1] * M[2][2] - M[2][1] * M[1][2]) -
			M[1][0]*(M[0][1] * M[2][2] - M[2][1] * M[0][2]) +
			M[2][0]*(M[0][1] * M[1][2] - M[1][1] * M[0][2]);
	}
	cerr << "Taking a Determinant of a matrix greater than 3x3." << endl;
	exit(1);
}

bool Matrix::isSquare() {
	if (row != col)
		return false;
	else
		return true;
}

bool Matrix::isSym() {
	if (!(this->isSquare()))
		return false;

	for (int i = 0; i < row; i++)
		for (int j = 0; j < col; j++)
			if (M[i][j] != M[j][i])
				return false;

	return true;
}

double Matrix::Trace() {
	double trace = 0;

	if (!(this->isSquare())) {
		cout << "Error: Matrix is not square\n";
		return -1;
	}

	for (int i = 0; i < row; i++)
		trace += M[i][i];

	return trace;
}

void Matrix::makeIdentity() {
	for (int i = 0; i < row; i++)
		for (int j = 0; j < row; j++)
			if (i != j) M[i][j] = 0;
			else M[i][j] = 1;

}

bool Matrix::isIdentity() {
	for (int i = 0; i < row; i++)
		for (int j = 0; j < row; j++)
			if (i != j && fabs(M[i][j]) > 1e-30) return false;
			else if (i == j && fabs(1 - M[i][j]) > 1e-30) return false;

	return true;
}

void Matrix::AppendRows(Matrix &N) {
	if (col != N.col) {
		cout << "Error: Array dimensions are incompatible\n";
		return;
	}
	double **oldM = M;
	M = new double *[row + N.row];

	for (int i = 0; i < row; i++) M[i] = oldM[i];
	delete [] oldM;

	for (int i = row; i < row + N.row; i++) M[i] = N.M[i - row];
	delete [] N.M;
	N.M = NULL;
	row += N.row;
	N.row = 0;
}

void Matrix::Copy(Matrix &N) {
	for (int i = 0; i < MIN(row, N.row); i++) {
		for (int j = 0; j < MIN(col, N.col); j++) {
			M[i][j] = N.M[i][j];
		}
	}
}

void Matrix::Print() {
	cout << "(" << row << ", " << col << ")\n";
	for (int i = 0; i < row; i++) {
		cout << "[ ";
		for (int j = 0; j < col; j++) {
			cout << M[i][j] << " ";
		}
		cout << "]\n";
	}
}

Matrix * Matrix::canon() { // Returns the rotation to be applied to canonicalize this matrix
	Matrix *R = new Matrix(3, 3);
	R->makeIdentity();

	double Sxy = M[0][1];
	double Sxz = M[0][2];
	double Syz = M[1][2];
				
	// dominating by absolute value (assume all are, then find the smallest one);
	bool dXY = true;
	bool dXZ = true;
	bool dYZ = true;

	// only one of these can be true
	     if (fabs(Sxy) < fabs(Sxz) && fabs(Sxy) < fabs(Syz)) dXY = false;
	else if (fabs(Sxz) < fabs(Sxy) && fabs(Sxz) < fabs(Syz)) dXZ = false;
	else if (fabs(Syz) < fabs(Sxy) && fabs(Syz) < fabs(Sxz)) dYZ = false;

	bool nXY = (dXY && Sxy < 0);
	bool nXZ = (dXZ && Sxz < 0);
	bool nYZ = (dYZ && Syz < 0);

	if (nXY) { // we have to negate XY
		// here we rotate such that the smallest by abs() gets rotated too
		// and the other largest is untouched

		if (!dXZ) // XZ can take one for the team
			R->M[0][0] = -1;
		else	  // in this case XZ is dominant, so we negate XY and YZ
			R->M[1][1] = -1;
	}
	if (nXZ) {
		if (!dXY)
			R->M[0][0] = -1;
		else
			R->M[2][2] = -1;
	}
	if (nYZ) {
		if (!dXY)
			R->M[1][1] = -1;
		else
			R->M[2][2] = -1;
	}

	return R;
}

