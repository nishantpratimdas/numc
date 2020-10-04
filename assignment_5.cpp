//Note: f(c)==0 => c is the root

#include <iostream>
#include <iomanip>
#include <math.h>
#include <fstream>
#include <algorithm>
#include <stdarg.h>
#include <ios>




#define epsilon 1/pow(10,6)

using namespace std;

float derivative(float f(float), float x)
{
	float h = 1/pow(10,6);
	return (f(x+h)-f(x-h))/(2*h);
}

float logMinusSin(float x)
{
	return log(x)-sin(x);
}

float minusXMinusCos(float x)
{
	return -x-cos(x);
}

void enhanceBracket(float* a, float* b, float f(float), float factor = 0.5, int iteration = 0){
	if(iteration == 20)
		return;
	iteration++;
//	cout << "here we go again " << abs(*b - *a) << " " << *a << "  " << *b << endl;
	if(abs(*b-*a)<1)
		return;
	while(f(*a + factor*abs(*b-*a))*f(*b)<0)
		*a += factor*abs(*b-*a);
	while(f(*b - factor*abs(*b-*a))*f(*a)<0)
		*b -= factor*abs(*b-*a);
	return f(*a + factor*abs(*b-*a))*f(*b)<0
	?
		enhanceBracket(a, b, f, factor,iteration)
	:
		enhanceBracket(a, b, f, factor/2,iteration);
}


bool bisection(string path, float *root_store, float a, float b, float f(float), int iteration = 0)
{
	if(iteration == 200)
		return false;
	
	ofstream file;
	file.open(path, ios::app);

	if(iteration==0)
		file << "Bisection \niteration " << "       " << " error\n\n";
	
	file <<setw(5) << iteration << "            " << fixed << setprecision(6) << abs(a-b) << endl;	
	iteration++;

	if (abs(a-b)<epsilon)
		return true;
		
	float c = (a+b)/2;
	*root_store = c;
	
	return f(c)*f(a)>0
	?
		bisection(path, root_store, c, b, f, iteration)
	:
		bisection(path, root_store, a, c, f, iteration);
}

bool regulaFalsi(string path, float *root_store, float a, float b, float f(float), int iteration = 0)
{
	if(iteration == 200)
		return false;
	ofstream file;
	file.open(path, ios::app);

	if(iteration==0)
		file << "\n\nRegula Falsi\niteration " << "       " << " error\n\n";

	file << setw(5) << iteration << "           " << fixed << setprecision(6) << abs(b-a) << endl;
	iteration++;

	float c =  (a*f(b)-b*f(a))/(f(b)-f(a));
	if(abs(c - *root_store)<epsilon)
		return true;

	*root_store = c;

	return f(c)*f(a)>0
		?
		regulaFalsi(path, root_store, c, b, f, iteration)
		:
		regulaFalsi(path, root_store, a, c, f, iteration);

}

bool newtonRaphson(string path, float *root_store, float a, float f(float), int iteration = 1)
{
	if(iteration == 200)
		return false;
	ofstream file;
	file.open(path, ios::app);

	if(iteration==1)
		file << "\n\nNewton-Raphson\n iteration " << "       " << " error\n\n";

	float a_next = a - (f(a)/derivative(f,a));

	file << setw(5) << iteration << "           " << fixed << setprecision(6) << abs(a_next-a) << endl;
	iteration++;

	if(abs(a_next-a)<epsilon)
		return true;
	*root_store = a_next;

	return newtonRaphson(path, root_store, a_next, f, iteration);
}

class Polynomials {

	private:
		int degree;
		double *coeffiecients;

	public:
		double * coeffiecientsCopy;
		Polynomials(int degree_incremented, ...){
			degree = degree_incremented -1;
			coeffiecients = new double[degree_incremented];

			va_list coeffiecientsList;
			reset();
			va_start(coeffiecientsList, degree_incremented);

			for (int i =degree; i>=0; i--){
				double argument = va_arg(coeffiecientsList, double);
				//				cout << argument << endl;
				*coeffiecientsCopy = argument;
				//				cout << *coeffiecientsCopy << endl;
				coeffiecientsCopy++;
			}
		}
		Polynomials(double * coeffiecientsList, int degree_incremented){
			degree = degree_incremented - 1;

			coeffiecients = new double[degree_incremented];
			reset();
			for (int i = degree; i>=0; i--){
				*coeffiecientsCopy = *coeffiecientsList;
				coeffiecientsList++;
				coeffiecientsCopy++;
			}

		}

		void printPolynomial(){
			reset();
			cout << "P(x) = ";
			for (int i =degree; i>0; i--){
				cout << showpos << *coeffiecientsCopy << noshowpos << "(x^" << i << ") ";
				coeffiecientsCopy++;
			}
			cout << showpos << *coeffiecientsCopy << noshowpos; 
			cout << endl << endl;
		}

		void reset(){
			coeffiecientsCopy = coeffiecients;
		}

		float getPolynomial (float x){
			float value = 0;
			reset();
			for (int i = degree; i>=0; i--) {
				value += float(*coeffiecientsCopy) * pow(x, i); 			
				//				cout << float(*coeffiecientsCopy) << " x " << pow(x, i) << endl; 			
				coeffiecientsCopy++;
			}


			return value;
		}

		float getFirstDerivative (float x){

			float value = 0;
			reset();

			for(int i = degree; i>0; i--){
				value += float(*coeffiecientsCopy)*i * pow(x, i-1);
				//				cout << *coeffiecientsCopy << " x " << pow(x, i-1) << endl;
				coeffiecientsCopy++;
			}

			return value;
		}

		float getSecondDerivative (float x){

			float value = 0;
			reset();

			for(int i = degree; i>1; i--){
				value += float(*coeffiecientsCopy)*i*(i-1) * pow(x, i-2);
				//				cout << *coeffiecientsCopy << " x " << pow(x, i-2) << endl;
				coeffiecientsCopy++;
			}

			return value;
		}
		int getDegree(){
			return float(degree);
		}

		double * getCoefficients(){
			return coeffiecientsCopy;
		}


		Polynomials syntheticDivision(double root){
			double* finalCoefficients = new double[degree];
			reset();
			double *finalCoefficientsCopy = finalCoefficients;
			*finalCoefficientsCopy = *coeffiecientsCopy;
			for(int i = 1; i<degree; i++){
				finalCoefficientsCopy++;
				coeffiecientsCopy++;
				*finalCoefficientsCopy = *coeffiecientsCopy + root* (*(finalCoefficientsCopy-1));
			}
			Polynomials p3(finalCoefficients, degree);
			return p3;
		}
};

float laguerreMethod (Polynomials p, float x, int iteration = 0)
{
//      Uncomment these lines to see the iterations in a tabular format
//	if (iteration == 0)
//		cout << "  iteration        estimate        P(estimate) \n\n";
//	cout << setw(6) << iteration << "            " << setw(7) << x << "          " <<  setw(7) << p.getPolynomial(x) << "\n";
	if (abs(p.getPolynomial(x)) < epsilon){
		cout << x << ", ";
		return x;
	}
	iteration++;

	float px = p.getPolynomial(x);
	float dp1 = p.getFirstDerivative(x);
	float dp2 = p.getSecondDerivative(x);
	float n = p.getDegree();

	float G = dp1/px;
	float H = pow(G,2) - dp2/px;

	float denominator = max(G + sqrt((n-1)*(n*H - pow(G,2))),
			G -  sqrt((n-1)*(n*H - pow(G,2))));
	//	cout << denominator << endl;
	float x_next = x - n/denominator;

	return laguerreMethod(p, x_next, iteration);

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
