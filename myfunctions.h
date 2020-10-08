#ifndef MATRIX_OPERATIONS
#define MATRIX_OPERATIONS

#include <string>
#include <math.h>
#include <iomanip>
#include <iostream>

float derivative(float f(float), float x);
void enhanceBracket(float* a, float* b, float f(float), float factor = 0.5, int iteration = 0);
bool bisection(std::string path, float *root_store, float a, float b, float f(float), int iteration = 0);
bool regulaFalsi(std::string path, float *root_store, float a, float b, float f(float), int iteration = 0);
bool newtonRaphson(std::string path, float *root_store, float a, float f(float), int iteration = 1);


void matrixMultiply(float *matrix_1, float *matrix_2, float *product, int n = 3);
void matrixMultiply(float *matrix_1, float *matrix_2, float *product, int rows_1, int common, int columns_2);
void partialPivot(float *arr, int position, bool* solutions, int n, int m);
void matrixPrint(std::string message, float *arr, int rows = 3, int columns = -1, bool augmented = false, int augwhere = 3);
void matrixRead(std::string path, float *arr, int columns = 3, int rows = 3, bool semiRead = false, int columnRead = 3);
void matrixRead(std::string path, int rows, int columns, int num, ...);
void matrixSum(float *matrix_1, float *matrix_2, float *sum, int n = 3);
void matrixSum(float *matrix_1, float *matrix_2, float *sum, int rows, int columns);
bool gaussJordanElimination(float *arr, int n, int m);
void extractVector(float *matrix, float *vector, int n, int m, int position);
bool matrixCompare(float *matrix_1, float *matrix_2, int n);
void matrixSetIdentity(float *matrix, int n, bool semiFill = false, int totalColumns = 6);
void matrixSetZero(float *matrix, int n);
bool decomposeLU(int n, int m, float *mat, float *L, float *U);
void substitution(float *L, float*U, float *B, float *X, int n = 3, bool XColumn = true);

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
			std::cout << "P(x) = ";
			for (int i =degree; i>0; i--){
				std::cout << std::showpos << *coeffiecientsCopy << std::noshowpos << "(x^" << i << ") ";
				coeffiecientsCopy++;
			}
			std::cout << std::showpos << *coeffiecientsCopy << std::noshowpos; 
			std::cout << std::endl << std::endl;
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

float laguerreMethod (Polynomials p, float x, int iteration = 0);

#endif