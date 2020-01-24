/*
 *  File: ANN.cpp
 *  Author: Warren Beagle
 *  Final Project
 *
 *  Created on December 12 2018, 2:43 PM
 */
#include "ANN.h"
#include <fstream>
#include <string>
#include <sstream>


using namespace std;

//Constructors
ANN::ANN() : P(0), L(0), M(0), N(0) //Default
{
	//Sets values to 0
}
ANN::ANN(int _L, int _M, int _N, int _P) : P(_P), L(_L), M(_M), N(_N) //Parameterized
{
	//Sets values to input values
	//Creates new matrices using the values
	I = new Matrix(_P, _L, 0);
	H = new Matrix(_P, _M, 0);
	O = new Matrix(_P, _N, 0);
	D = new Matrix(_P, _N, 0);
	W1 = new Matrix(_L, _M, 0);
	W2 = new Matrix(_M, _N, 0);
	bH = new Matrix(1, _M, 0);
	bO = new Matrix(1, _N, 0);
	tW1 = new Matrix(_L, _M, 0);
	tW2 = new Matrix(_M, _N, 0);
	tbH = new Matrix(1, _M, 0);
	tbO = new Matrix(1, _N, 0);
	tO = new Matrix(_P, _N, 0);
}
//Deconstructor
ANN::~ANN()
{
	//Deletes the matrices
	delete I;
	delete H;
	delete O;
	delete D;
	delete W1;
	delete W2;
	delete bH;
	delete bO;
	delete tW1;
	delete tW2;
	delete tbH;
	delete tbO;
	delete tO;
}

//Other functions
Matrix ANN::forwardPropogate(Matrix *input, Matrix *input2) //Takes in 2 matrices
{
	//Multiplies the matrices together to get the next matrix in the neural network
	int row = input->GetRow();
	int column = input2->GetCol();
	Matrix temp(row, column, 0);
	temp = *input * *input2;
	return temp; //Returns a matrix
}
void ANN::train(int epochs) //Trains the neural network over a input number of epochs
{
	for(int i = 0; i < epochs; i++)
	{
		updateWeights();
		*H = forwardPropogate(I, W1);
		//Account for bias and sigmoid for the Hidden layer
		for(int j = 0; j < H->GetRow(); j++)
		{
			for(int k = 0; k < H->GetCol(); k++)
			{
				H->Set_Mij(j, k, (H->Get_Mij(j, k) + bH->Get_Mij(0, k)));
				H->Set_Mij(j, k, sig(H->Get_Mij(j, k)));
			}
		}
		*O = forwardPropogate(H, W2);
		//Account for bias and sigmoid for the Output laye
		for(int j = 0; j < O->GetRow(); j++)
		{
			for(int k = 0; k < O->GetCol(); k++)
			{
				O->Set_Mij(j, k, (O->Get_Mij(j, k) + bO->Get_Mij(0, k)));
				O->Set_Mij(j, k, sig(O->Get_Mij(j, k)));
			}
		}
		//Tests the score 
		if(sd(O, D) < score || i == 0)
		{
			//Sets best matrices
			*tW1 = *W1;
			*tW2 = *W2;
			*tbH = *bH;
			*tbO = *bO;
			*tO = *O;
			score = sd(O, D);
			//Prints data
			cout << "Epoch: " << i << endl;
			cout << "W1: " << endl; 
			W1->Print();
			cout << "W2: " << endl; 
			W2->Print();
			cout << "bH: " << endl; 
			bH->Print();
			cout << "bO: " << endl; 
			bO->Print();
			cout << "Actual Output: " << endl;
			O->Print();
			cout << "Score: " << score << endl;
			cout << "---------------------------------------" << endl;
		}
	}
}
void ANN::loadData(string file) //Loads the data from an input file
{
	ifstream inFile;
	inFile.open(file.c_str());
	if(inFile.fail())
	{
		cout << "Unable to load file" << endl;
		exit(0);
	}
	else
	{
		double input = 9;
		for(int i = 0; i < I->GetRow(); i++)
		{
			for(int j = 0; j <= I->GetCol(); j++)
			{
				//Inputs
				if(j < I->GetCol())
				{
					inFile >> input;
					I->Set_Mij(i, j, input);
				}
				//Desired Output
				else
				{
					inFile >> input;
					D->Set_Mij(i, (j - I->GetCol()), input);
				}
			}
		}
	}
}
double ANN::sig(double input) //Gets the sigmoid of an input
{
	return 1/(1 + (exp (-1 * input)));
}
double ANN::sd(Matrix *desiredOutput, Matrix *output) //Gets the error using two matrices
{
	double sum(0);
	for(int i = 0; i < desiredOutput->GetRow(); i++)
	{
		for(int j = 0; j < output->GetCol(); j++)
		{
			sum += ((output->Get_Mij(i, j) - desiredOutput->Get_Mij(i, j)) * (output->Get_Mij(i, j) - desiredOutput->Get_Mij(i, j)));
		}
	}
	return sum;
}
void ANN::updateWeights() //Updates the matrices with random values 
{
	for(int i = 0; i < W1->GetRow(); i++)
	{
		for(int j = 0; j < W1->GetCol(); j++)
		{
			//Sets the values to random numbers between -2 and 2
			W1->Set_Mij(i, j, (4 * drand48() - 2));
		}
	}
	for(int i = 0; i < W2->GetRow(); i++)
	{
		W2->Set_Mij(i, 0, (4 * drand48() - 2));
	}
	for(int i = 0; i < bH->GetCol(); i++)
	{
		bH->Set_Mij(0, i, (4 * drand48() - 2));
	}
	bO->Set_Mij(0, 0, (4 * drand48() - 2));
}
void ANN::evaluate() //Evaluates the neural network
{
	//Prints 
	cout << "---------------------------------------" << endl;
	cout << "Input: " << endl; 
	I->Print();
	cout << "Desired Output: " << endl; 
	D->Print();
	cout << "Actual Ouput: " << endl;
	tO->Print();
	cout << "Best Weights: " << endl; 
	cout << "W1: " << endl;
	tW1->Print();
	cout  << "W2: " << endl; 
	tW2->Print();
	cout << "Best H Biases: " << endl;
	tbH->Print();
	cout << "Best O Biases: " << endl;
	tbO->Print();
	cout << "Best Score: " << score << endl;
}

