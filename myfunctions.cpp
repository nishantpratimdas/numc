#ifndef MATRIX_OPERATIONS
#define MATRIX_OPERATIONS
#include <iostream>
#include <string>
#include <iomanip>
#include <fstream>
#include <stdarg.h>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <fstream>
#include <algorithm>
#include <stdarg.h>
#include <ios>

#define epsilon 1/pow(10,6)


using namespace std;

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

// float laguerreMethod (Polynomials p, float x, int iteration = 0);

// ******************************* ROOTS ************************
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
	//	cout << p.getPolynomial(x) << endl;
	//	cout << p.getFirstDerivative(x) << endl;
	//	cout << p.getSecondDerivative(x) << endl;

	//	return 1;
}

float derivative(float f(float), float x)
{
	float h = 1/pow(10,6);
	return (f(x+h)-f(x-h))/(2*h);
}



void enhanceBracket(float* a, float* b, float f(float), float factor = 0.5, int iteration=0){
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

bool newtonRaphson(string path, float *root_store, float a, float f(float), int iteration = 0)
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




#endif

//********************MATRIX OPERATIONS **********************************

// This print function prints a 3x3 matrix if no dimensions are specified
// If a dimension is specified, it assumes we want to print a square matrix
// By default, augmented is set to false. 
// If we want to print a augmented matrix, bool augmented has to be set to true
// Where the lines (|) should appear can also be specified 
// defaults to between 3rd and 4th column
void matrixPrint(string message, float *arr, int rows = 3, int columns = -1, bool augmented = false, int augwhere = 3)
{ 
    cout << message << "\n\n";

    if(columns == -1)
        columns = rows;
    
    for (int i=0; i<rows; i++)
    {
        for(int j = 0; j<columns; j++)
        {
            if(augmented && j == augwhere) 
                cout << "  |" << setw(7) << setprecision(2) << *arr;
            else
                cout << setw(7) << setprecision(2) << *arr;

            arr++;
        }
        cout << "\n";
    }
    cout << "\n\n";
}
// Although the code is leaner this way, i.e, if the augmented check is done during the loop is running
// It is more efficient (I believe), if the check is done once and then we have two different loops
// Note to self - confirm and change. 

// This function interchanges two rows of a matrix(partial pivot)
// iff the diagonal element at the ith position is 0
// If a all zero row is found, it doesn't interchange
// and set bool* solutions to false
void partialPivot(float *arr, int position, bool* solutions, int n, int m)
{
    int i = position;
    if(*(arr+(m*i)+i) == 0)
    {
        // cout << "reach  " << i << endl;
        // matrixPrint("A", arr, 4, 8, true, 4);

        *solutions = false;

        bool solutionsFlag = true;
        for(int j = i+1; solutionsFlag && j<n ;j++)
        {
            // cout << "           " << *(arr+(m*j)+i);
            if(*(arr+(m*j)+i) != 0)
            {
                *solutions = true;
                solutionsFlag = false;
            }
        }
            

        if(*solutions == false)
            cout << "Error. Failed to resolve row " << i+1 << endl;
        //Improve this message


        //this loop interchanges the pivot row with the first possible interchange-able row
        bool interchanged = false;
        int j = i+1;
        while(j<n && !interchanged && *solutions)
        {
            if(*(arr+(m*j)+i) != 0)
            {
                // cout << "exchange " << i << " " << j << endl;
                for(int k = 0; k<m; k++)
                {
                    float holder = *(arr+(m*i)+k);
                    *(arr+(m*i)+k) = *(arr+(m*j)+k);
                    *(arr+(m*j)+k) = holder;
                }
                interchanged = true;
            }
            j++;
        }
        // matrixPrint("A", arr, 4, 8, true, 4);

    }
}

//this function first performs Gauss Jordan Elimination
//and then returns true if the it was possible,
//that is, no row row was found in the first matrix
//returns false otherwise.
//Note to self - default: 3xn
bool gaussJordanElimination(float *arr, int n, int m)
{
    bool uniqueSolutions = true;

    for(int i = 0; uniqueSolutions && i<n; i++)
    {
        partialPivot(arr, i, &uniqueSolutions, n, m);

        float scalef = *(arr+(m*i)+i);
        for (int j = 0; j<m; j++)
            *(arr+(m*i)+j) /= scalef;

        for(int j = 0; j<n; j++)
            if(i != j)
            {
                float factor = *(arr+(m*j)+i);
                for(int k=0; k<m; k++)
                    *(arr+(m*j)+k) -= (factor*(*(arr+(m*i)+k)));
            }

    }

    return uniqueSolutions;
}

//This function stores the product of two square matrices, matrix_1 and matrix_2,
//in product. The default dimensions are 3x3.
//The matrunner runs through the matrix as 'cursors' and also
//are used to run through rows and columns to calculate individual elements of the product
void matrixMultiply(float*matrix_1, float*matrix_2, float*product, int n=3)
{
    float *matrunner_1, *matrunner_2;
    for(int i = 0; i<n; i++)
        for(int j = 0; j<n; j++)
        {
            matrunner_1 = matrix_1+(i*n);
            matrunner_2 = matrix_2+j;
            *product = 0;
            for(int k = 0; k<n; k++)
            {
                *product += (*matrunner_1)*(*matrunner_2);
                matrunner_1++;
                matrunner_2+=n;
            }
            product++;
        }
}
//This is a more general version. Multiplies (rows_1xcommon) x (commonxcolumns_2) and stores it in 
//product which will be (rows_1xcolumns_2)
void matrixMultiply(float*matrix_1, float*matrix_2, float*product, int rows_1, int common, int columns_2)
{
    float *matrunner_1, *matrunner_2;
    for(int i = 0; i<rows_1; i++)
        for(int j = 0; j<columns_2; j++)
        {
            matrunner_1 = matrix_1+(i*common);
            matrunner_2 = matrix_2+j;
            for(int k = 0; k<common; k++)
            {
                *product += (*matrunner_1)*(*matrunner_2);
                matrunner_1++;
                matrunner_2+=columns_2;
            }
            product++;
        }
}

//Unfortunately, as of now, there are no checks to know if an 'invalid array' has been passed.
//Note to self - try to implement. 

//This function read 1 matrix from path
//a 3x3 matrix is read, if no dimensions are specified
//a 3xn matrix is read, if only one dimension is specified
//a m(2nd argument) x n(1st arguement) is read, if both dimensions are specified
//although 2nd argument x 1st argument is counter-intuitive, to the best of my knowledge,
//it has to be done this way since C++ does not support named arguments
bool matrixCompare(float *matrix_1, float *matrix_2, int n)
{
    bool equal = true, flag=true;
    for (int i = 0; i<n && flag; i++)
    {
        for(int j=0; j<n && flag; j++)
        {
            if( *(matrix_1+(n*i)+j) != *(matrix_2+(n*i)+j))
            {
                equal = false;
                flag = false;
            }
        }
    }

    return equal;
}

//Extracts a vector from the given position and stores it in vector
//Then it deletes the column from the matrix
void extractVector(float *matrix, float *vector, int n, int m, int position)
{
    vector+=n-1;
    for(int i=n-1; i>=0; i--)
    {
        for(int j=0; j<m; j++)
            if(j == position-1)
                *vector-- = *(matrix+(m*i)+j);
 
        for(int j=0; j<((n-1-i)*m)+(m-position); j++)
            *(matrix+(m*i)+position-1+j) = *(matrix+(m*i)+position+j);
    }

}

//Reads 1 matrix from a file
//if semiRead is true the function is capable of reading a nxm matrix and storing it
//in a axb matrix when a>=n, b>=m.
void matrixRead(string path, float *arr, int columns=3, int rows=3, bool semiRead=false, int columnRead=3)
{
    if (semiRead == false)
        columnRead = columns;

    ifstream file;
    file.open(path);

    float *matrunner;
    for (int i = 0; i < rows; i++) 
    {
        matrunner = arr + (columns*i);
        for (int j = 0; j < columnRead; j++) 
        {
            file >> *matrunner;
            matrunner++;
        }
    }
     
		
    file.close();
} 


//This is a more general function which reads num number of matrices from path
//and stores it in the variable argument list provided
//here it wasn't necessary to accept the columns first and then the row
//nonetheless, it has been done that way so that it is coherent with 
//the simpler matrixRead defined above
//This needs improvement.
void matrixRead(string path, int columns, int rows, int num, ...)
{

    ifstream file;
    file.open(path);

    va_list listOfMatrix;
    va_start(listOfMatrix, num);

    for(int n = 0; n<num; n++)
    { 
        float *arr = va_arg(listOfMatrix, float *);
        for (int i = 0; i < rows; i++) 
            for (int j = 0; j < columns; j++) 
            {
                file >> *arr;
                arr++;
            }
    }
		
    file.close();
} 

//Adds two square matrices
void matrixSum(float *matrix_1, float *matrix_2, float *sum, int n=3)
{
    for(int i=0; i<n; i++)
        for(int j=0; j<n; j++)
        {
            *sum = 0;
            *sum = (*matrix_1)+(*matrix_2);
             
            sum++;
            matrix_1++;
            matrix_2++;
        }
}

//Adds two nxm matrices
void matrixSum(float *matrix_1, float *matrix_2, float *sum, int rows, int columns)
{
    for(int i=0; i<rows; i++)
        for(int j=0; j<columns; j++)
        {
            *sum = 0;
            *sum = (*matrix_1)+(*matrix_2);

            sum++;
            matrix_1++;
            matrix_2++;
        }
}

//Sets the given matrix as I_4
//If semiFill is true, it is cabable of setting a part of a matrix as I
//This is useful when we augment a matrix with I to find it's inverse
void matrixSetIdentity(float *matrix, int n, bool semiFill=false, int totalColumns=6)
{
    if(semiFill == false)
        totalColumns = n;

    for(int i=0; i<n; i++)
        for(int j=0; j<n; j++)
        {
            if(i==j)
                *(matrix+(totalColumns*i)+j)=1;
            else 
                *(matrix+(totalColumns*i)+j)=0;
        }
}

//sets the given matrix to 0 
void matrixSetZero(float *matrix, int n)
{
    for(int i=0; i<n; i++)
        for(int j=0; j<n; j++)
            *(matrix+(n*i)+j)=0;
        
}

//Function to decompose a matrix mat, and store its's L part in L, U part in U.
//returns true if the decomposition was possible, false otherwise.
bool decomposeLU(int n, int m, float *mat, float *L, float *U)
{
    bool decomposePossible = true;
    for(int i=0; i<n && decomposePossible; i++)
    {
        partialPivot(mat, i, &decomposePossible, n, m);
        for(int j=0; j<n; j++)
        {
            float sum = 0;
            if(j>i)
            {
                for(int k=0; k<i; k++)
                    sum+= *(L+k+(n*j))*(*(U+i+(n*k)));
                *(L+(n*j)+i)=(*(mat+(m*j)+i)-sum)/(*(U+(n*i)+i));
            }
            else
            {
                for(int k=0; k<j; k++)
                    sum+= *(L+k+(n*j))*(*(U+i+(n*k)));
                *(U+(n*j)+i)=*(mat+(m*j)+i)-sum;
            }           
        }
    }

    return decomposePossible;
}

//Function to do the forward-backward substitution.
//The result is stored in X
void substitution(float *L, float*U, float *B, float *X, int n=3, bool XColumn=true)
{
    int fillX;
    if(XColumn == true)
        fillX = n;
    else
        fillX = 1;
    

    //use L and B to find Y, LY = B
    //forward subsitution

    float Y[n];

	for(int i =0; i<n; i++)
		*(Y+i) = 0;

    for(int i=0; i<n; i++)
    {
        float sum = 0;

        for(int j=0; j<i; j++)
            sum += *(L+(n*i)+j)*Y[j];

        Y[i] = *(B+i) - sum;
    }
    
    //use U and Y to find X. UX = Y
    //backward subsitution

    for(int i=n-1; i>=0; i--)
    {
        float sum = 0;

        for(int j=i+1; j<n; j++)
            sum += *(U+(n*i)+j)*(*(X+(fillX*j)));
        
        *(X+(fillX*i)) = (Y[i] - sum)/(*(U+(n*i)+i));
    }

    //Using the following instead of the printing function in matrix_operations to make to output look cleaner
    cout << "For B=       Y=      X= \n\n";
    for(int i=0; i<n; i++)
        cout << setw(7) << setprecision(2) << *(B+i) << setw(8) << setprecision(2) << *(Y+i) << setw(9) << setprecision(2) << *(X+(fillX*i)) << endl;
    cout << "\n\n";
}
