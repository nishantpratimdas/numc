//Note: f(c)==0 => c is the root

#include <iostream>
#include <iomanip>
#include <math.h>
#include <fstream>
#include <algorithm>
#include <stdarg.h>
#include <ios>
#include "myfunctions.cpp"

float logMinusSin(float x)
{
	return log(x)-sin(x);
}

float minusXMinusCos(float x)
{
	return -x-cos(x);
}


int main()
{
	float root;
	string path;
	ofstream file;

	cout << "Error: 0.000001\n\n";

	cout << "Q1A. Root of f(x) = ln(x)-sin(x)\n\n" << endl;
	path = "data1_a.txt";
	
	file.open(path, ios::app);
	file << "Quetion 1A: f(x) = ln(x)-sin(x)\n\n";
	file.close();
	if(bisection(path, &root,1.5, 2.5, logMinusSin))
		cout << "Using Bisection: " << fixed << setprecision(6) <<  root  << endl;
	if(regulaFalsi(path, &root,1.5, 2.5, logMinusSin))
		cout << "Using Regula Falsi: " << fixed << setprecision(6) <<  root  << endl;
	if(newtonRaphson(path, &root,2, logMinusSin))
		cout << "Using Newton Raphson: " << fixed << setprecision(6) <<  root  << endl;

	cout << "\n\nQ1B. Root of f(x) = ln(x)-sin(x)\n\n" << endl;
	path = "data2_a.txt";

	file.open(path, ios::app);
	file << "Quetion 1B: f(x) = -x-cos(x)\n\n";
	file.close();	

	float bracketStart = -10;
	float bracketEnd = 10;
	enhanceBracket(&bracketStart, &bracketEnd, minusXMinusCos);	


	if(bisection(path, &root,bracketStart, bracketEnd, minusXMinusCos))
		cout << "Using Bisection: " << fixed << setprecision(6) <<  root  << endl;
	if(regulaFalsi(path, &root, bracketStart, bracketEnd, minusXMinusCos))
		cout << "Using Regula Falsi: " << fixed << setprecision(6) <<  root  << endl;
	if(newtonRaphson(path, &root, bracketEnd, minusXMinusCos))
		cout << "Using Newton Raphson: " << fixed << setprecision(6) <<  root  << endl;

	cout << "\n\nQ3: roots of the polynomial:\n\n";
	Polynomials p1(5 ,1.0, -3.0, -7.0, 27.0, -18.0); 
	p1.printPolynomial();
	
	cout << "real roots are: ";
	for (int i = 0; i<4; i++){
		root = laguerreMethod(p1, 1.5);
		p1 = p1.syntheticDivision(root);
	}
	cout << "\b\b. " << endl;
}

// Error: 0.000001

// Q1A. Root of f(x) = ln(x)-sin(x)


// Using Bisection: 2.219108
// Using Regula Falsi: 2.219107
// Using Newton Raphson: 2.219107


// Q1B. Root of f(x) = ln(x)-sin(x)


// Using Bisection: -0.739085
// Using Regula Falsi: -0.739085
// Using Newton Raphson: -0.739085


// Q3: roots of the polynomial:

// P(x) = +1.000000(x^4) -3.000000(x^3) -7.000000(x^2) +27.000000(x^1) -18.000000

// real roots are: 1.000000, -3.000000, 3.000000, 2.000000. 
