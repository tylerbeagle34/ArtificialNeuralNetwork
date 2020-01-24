#include <cstdlib>
#include <iostream>
#include "ANN.h"
#include "matrix.h"

using namespace std;

int main(int argc, char** argv) {
    
    	int iSize = 0;
    	int hSize = 0;
    	int oSize = 0;
    	int pSize = 0;
    	int epochs = 0;
	
	if (argc != 7) 
	{
        cout << "Usage: " << argv[0] << " <L> <M> <N> <Epochs> <P> <TrainingFile>" << endl;
        exit(0);
    	}
    	iSize = atoi(argv[1]);
    	hSize = atoi(argv[2]);
    	oSize = atoi(argv[3]);
    	epochs = atoi(argv[4]);
    	pSize = atoi(argv[5]);
	
	ANN myAnn(iSize, hSize, oSize, pSize);
	myAnn.loadData(argv[6]);
	myAnn.train(epochs);
	myAnn.evaluate();
  
    return 0;
    
}

