/*
 * File: ANN.h
 * Author: Warren Beagle
 * 
 * Created on December 12 2018, 2:43PM
 *
 * Simulates an artificial neural network. Utilizes funcionalities like File I/O and pointers.
 *           Incorporates different functions to make the neural network multi-staged
 *           Class inheritance, composition, operator overloading, functions, arrays and pointers
 */

#ifndef ANN_H
#define ANN_H
#include <iostream>
#include "matrix.h"
#include <math.h>
#include <stdlib.h>

using namespace std;
class ANN;

class ANN {

public:
	//Constructors
	ANN(); //Default
	ANN(int _P, int _L, int _M, int _N); //Parameterized
	//Deconstructor
	~ANN();
	
	//Other functions
	Matrix forwardPropogate(Matrix *input, Matrix *weight); //Multiplies two matrices to move forward in the neural network
	void train(int epochs); //Trains the neural network with a number of epochs
	void loadData(string file); //Loads the data from an input file
	double sig(double input); //Calculates the sigmoid of a value
	double sd(Matrix *input, Matrix *output); //Calculates the error and sums the squares
	void updateWeights(); //Updates the weights to random values
	void evaluate(); //Evaluates the neural network; prints
	
private:
	//Member variables
	//Matrices
	Matrix *I;
	Matrix *H;
	Matrix *O;
	Matrix *D;
	Matrix *W1;
	Matrix *W2;
	Matrix *bH;
	Matrix *bO;
	Matrix *tW1;
	Matrix *tW2;
	Matrix *tbH;
	Matrix *tbO;
	Matrix *tO;
	//Sizes
	int P;
	int L;
	int M;
	int N;
	//Score
	double score;
	
};

#endif
