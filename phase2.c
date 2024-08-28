#include<stdio.h>
#include<math.h>
#include<stdlib.h>

/////////////////////////////////////////////////////////////Function prototype///////////////////////////////////////////
double sumSquare(int n, double x[]);
void Sgradient(int n, double x[], double grad[], double h);                               // n = total variables, x[] = array with initial values of the variables,  h = very small value to find partial derivative (0.00001)

double rosenBrock(int n, double x[]);
void Rgradient(int n, double x[], double grad[], double h);          //grad[] = array to store gardient values

double dixonPrice(int n, double x[]);
void DPgradient(int n, double x[], double grad[], double h);

double trid(int n, double x[]);
void Tgradient(int n, double x[], double grad[], double h);

double zakharov(int n, double x[]);
void Zgradient(int n, double x[], double grad[], double h);

void SfindHessian(int n, double x[], double Hessian[][n], double h);            // Hessian[][n], a 2D array will store the values of hessian matrix
void RfindHessian(int n, double x[], double Hessian[][n], double h);
void DPfindHessian(int n, double x[], double Hessian[][n], double h);
void TfindHessian(int n, double x[], double Hessian[][n], double h);
void ZfindHessian(int n, double x[], double Hessian[][n], double h);

void check(int n, double Hessain[][n], double matrix_inv[][n]);                   // matrix_inv = will store the value of hessian invers. Will check wheather Hessian matrix is positive definite/ semidefinite

void inverseUpperTriangular(int n, double UpperTriangular[][n], double matrix_inv[][n]);      //use to find the inverse of a UPPER TRIANGULAR FUNCTION
double newton(int n, double x[], double grad[],double Hessian[][n], double matrix_inv[][n], double h,double product[]);
double norm(int n, double vector[]);
double matrixMult(int n, double A[][n], double B[],double product[]);                             //will do mutiplication between a matrix A (hessian inverse) and a column vector B (gradient), store the results in 1D array product[]
//////////////////////////////////////////////////////////////////////main() function////////////////////////////////////
FILE* fp;                     //file pointer to read the input data from "input.txt"

int main(){
	fp = fopen("input.txt","r");
	int n;                           //total number of variables
	double h =0.00001;             //small value to find derivatives
	fscanf(fp,"%d",&n);            //reading 
	//printf("%d\n",n);
	double x[n];                //1D array to store the initial point
	for (int i =1;i<=n;i++){            //loop reading the file and storing initial point
		fscanf(fp,"%lf",&x[i-1]);           
		printf("%lf\n",x[i-1]);
	}

	double grad[n];                  //1D array to store the gradient values
	double matrix_inv[n][n];             //2D to store the inverse of HESSIAN MATRIX in particular
	double product[n];              //1D array to store MULTIPLICATION results of HESSIAN INVERSE and GRADIENT
	double Hessian[n][n];                //2D array to store the elements of hessian matrix
	
	//Rgradient(n,x,grad,0.00001);
	Zgradient(n,x,grad,h);                     //calling gradient and hessian functions
	ZfindHessian(n,x,Hessian,h);
	check(n,Hessian,matrix_inv);               //calling check() to check hessian is positive definite/semidefinite or not
	//printf("%lf\n",norm(n,grad));
	newton(n,x,grad,Hessian,matrix_inv,h,product);
	matrixMult(n,matrix_inv,grad,product);
	
return 0;
}

///////////////////////////////////////////////////////////////function definition////////////////////////////////////////////
/*double sumSquare(int n, double x[]) {
    double result = 0.0;
    for (int i = 1; i <= n; i++) {
        result += i*pow(x[i-1],2);
    }
    return result;
}*/
/*
double rosenBrock(int n, double x[]) {
    double result = 0.0;
    for (int i = 0; i < n - 1; i++) {
        double firstTerm = pow((1 - x[i]),2);
        double secondTerm = 100 * pow((x[i + 1] - x[i] * x[i]),2);
        result += firstTerm + secondTerm;
    }
    return result;
}
*/
/*
double dixonPrice(int n, double x[]) {
    double result = pow((x[0] - 1), 2);
    for (int i = 1; i < n; i++) {
        double firstTerm = i + 1;
        double secondTerm = pow((2 * x[i] * x[i] - x[i - 1]),2);
        result += firstTerm * secondTerm;
    }
    return result;
}
*/
/*
double trid(int n, double x[]) {
    double result = 0.0;
    
    // Calculate the first part of the function
    double firstSum = 0.0;
    for (int i = 0; i < n; i++) {
        firstSum += pow((x[i] - 1),2);
    }
    
    // Calculate the second part of the function
    double secondSum = 0.0;
    for (int i = 1; i < n; i++) {
        secondSum += x[i] * x[i - 1];
    }
    
    result = firstSum - secondSum;
    
    return result;
}
*/

double zakharov(int n, double x[]) {
    double firstSum = 0.0;
    double secondSum = 0.0;
    
    for (int i = 0; i < n; i++) {
        firstSum += x[i] * x[i];
        secondSum += 0.5 * (i + 1) * x[i];
    }
    
    double result = firstSum + pow(secondSum,2) + pow(secondSum,4);
    
    return result;
}

/*double f(int n, double x[]) {
    double result = 0.0;
    for (int i = 1; i <= n; i++) {
        result += x[i-1] * x[i-1];
    }
    return result;
}*/
 
/////////////////////////////////////////////////////////////////Gradient functions///////////////////////////////////////////////////////////////
/*
void Sgradient(int n, double x[], double grad[], double h) {
    for (int i = 1; i <= n; i++) {
        double original = x[i-1];
        x[i-1] += h;
        double f_plus_h = sumSquare(n, x);
        x[i-1] = original - h;
        double f_minus_h = sumSquare(n, x);
        x[i-1] = original;
        grad[i-1] = (f_plus_h - f_minus_h) / (2 * h);
        printf("Gradient of f with respect x[%d] = %lf\n",i,grad[i-1]);
    }
    printf("\n");
    
} */
/*
void Rgradient(int n, double x[], double grad[], double h) {
    for (int i = 1; i <= n; i++) {
        double original = x[i-1];
        x[i-1] += h;
        double f_plus_h = rosenBrock(n, x);
        x[i-1] = original - h;
        double f_minus_h = rosenBrock(n, x);
        x[i-1] = original;
        grad[i-1] = (f_plus_h - f_minus_h) / (2 * h);
        printf("Gradient of f with respect x[%d] = %lf\n",i,grad[i-1]);
    }
    printf("\n");
    
}
*/
/*
void DPgradient(int n, double x[], double grad[], double h) {
    for (int i = 1; i <= n; i++) {
        double original = x[i-1];
        x[i-1] += h;
        double f_plus_h = dixonPrice(n, x);
        x[i-1] = original - h;
        double f_minus_h = dixonPrice(n, x);
        x[i-1] = original;
        grad[i-1] = (f_plus_h - f_minus_h) / (2 * h);
        printf("Gradient of f with respect x[%d] = %lf\n",i,grad[i-1]);
    }
    printf("\n");
    
}
*/
/*
void Tgradient(int n, double x[], double grad[], double h) {
    for (int i = 1; i <= n; i++) {
        double original = x[i-1];
        x[i-1] += h;
        double f_plus_h = trid(n, x);
        x[i-1] = original - h;
        double f_minus_h = trid(n, x);
        x[i-1] = original;
        grad[i-1] = (f_plus_h - f_minus_h) / (2 * h);
        printf("Gradient of f with respect x[%d] = %lf\n",i,grad[i-1]);
    }
    printf("\n");
    
}
*/
void Zgradient(int n, double x[], double grad[], double h) {
    for (int i = 1; i <= n; i++) {
        double original = x[i-1];     //storing original value of initial point
        x[i-1] += h;                           //updating x[i-1]
        double f_plus_h = zakharov(n, x);                //x updtaed array
        x[i-1] = original - h;
        double f_minus_h = zakharov(n, x);
        x[i-1] = original;
        grad[i-1] = (f_plus_h - f_minus_h) / (2 * h);
        printf("Gradient of f with respect x[%d] = %lf\n",i,grad[i-1]);
    }
    printf("\n");
    
}

///////////////////////////////////////////////////////////Hessian matrix////////////////////////////////////////////////////////
/*
void SfindHessian(int n, double x[], double Hessian[][n], double h) {
     printf("//////////////////////////////////////////////////////Required Hessian Matrix/////////////////////////////////////////////\n");

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            double original_i = x[i];
            double original_j = x[j];
	    //printf("I am Ist oriI = %lf oriJ = %lf\n",original_i,original_j);
            if (i==j){                                              //will find diagonal elements
            	double f0 = 2*sumSquare(n,x);
            	x[i] += h;
            	//printf("x[%d] = %lf\n",i,x[i]);
            	double f1 = sumSquare(n,x);
            	x[i] = original_i;
            	x[i] -= h;
            	//printf("x[%d] = %lf\n",i,x[i]);
            	double f2 = sumSquare(n,x);
            	Hessian[i][j] = ((f1 - f0 + f2) / (h * h));
            	//printf("f0 = %lf f1 = %lf f2 = %lf Hessian[%d][%d] = %lf\n",f0,f1,f2,i,j,Hessian[i][j]);
            }
            
            else{                                                   //will find off-diagonal elements
            // Compute f(x + h, x_i)
            	x[i] += h;
            	x[j] += h;
            	double f1 = sumSquare(n, x);

            // Compute f(x - h, x_i)
            	x[i] = original_i;
            	x[i] += h;
            	x[j] = original_j;
            	x[j] -= h;
            	double f2 = sumSquare(n, x);

            // Compute f(x, x_j + h)
            	x[i] = original_i;
            	x[i] -= h;
            	x[j] = original_j;
            	x[j] += h;
            	double f3 = sumSquare(n, x);

            // Compute f(x, x_j - h)
            	x[i] = original_i;
            	x[i] -=h;
            	x[j] = original_j;
            	x[j] -=h;
            	double f4 = sumSquare(n, x);

            // Calculate the second partial derivative
            	Hessian[i][j] = ((f1 - f2 - f3 + f4) / (4 * h * h));
	    }
            // Restore the original values
            x[i] = original_i;
            x[j] = original_j;
            printf("%10.3lf ", Hessian[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}
*/
/*
void RfindHessian(int n, double x[], double Hessian[][n], double h) {
     printf("//////////////////////////////////////////////////////Required Hessian Matrix/////////////////////////////////////////////\n");

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            double original_i = x[i];
            double original_j = x[j];
	    //printf("I am Ist oriI = %lf oriJ = %lf\n",original_i,original_j);
            if (i==j){                                                    //will find diagonal elements
            	double f0 = 2*rosenBrock(n,x);
            	x[i] += h;
            	//printf("x[%d] = %lf\n",i,x[i]);
            	double f1 = rosenBrock(n,x);
            	x[i] = original_i;
            	x[i] -= h;
            	//printf("x[%d] = %lf\n",i,x[i]);
            	double f2 = rosenBrock(n,x);
            	Hessian[i][j] = ((f1 - f0 + f2) / (h * h));
            	//printf("f0 = %lf f1 = %lf f2 = %lf Hessian[%d][%d] = %lf\n",f0,f1,f2,i,j,Hessian[i][j]);
            }
            
            else{                                                                  // will find off-diagonal elements
            // Compute f(x + h, x_i)
            	x[i] += h;
            	x[j] += h;
            	double f1 = rosenBrock(n, x);

            // Compute f(x - h, x_i)
            	x[i] = original_i;
            	x[i] += h;
            	x[j] = original_j;
            	x[j] -= h;
            	double f2 = rosenBrock(n, x);

            // Compute f(x, x_j + h)
            	x[i] = original_i;
            	x[i] -= h;
            	x[j] = original_j;
            	x[j] += h;
            	double f3 = rosenBrock(n, x);

            // Compute f(x, x_j - h)
            	x[i] = original_i;
            	x[i] -=h;
            	x[j] = original_j;
            	x[j] -=h;
            	double f4 = rosenBrock(n, x);

            // Calculate the second partial derivative
            	Hessian[i][j] = ((f1 - f2 - f3 + f4) / (4 * h * h));
	    }
            // Restore the original values
            x[i] = original_i;
            x[j] = original_j;
            printf("%10.3lf ", Hessian[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}
*/
/*
void DPfindHessian(int n, double x[], double Hessian[][n], double h) {
     printf("//////////////////////////////////////////////////////Required Hessian Matrix/////////////////////////////////////////////\n");

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            double original_i = x[i];
            double original_j = x[j];
	    //printf("I am Ist oriI = %lf oriJ = %lf\n",original_i,original_j);
            if (i==j){                                                                                  //will find diagonal elements
            	double f0 = 2*dixonPrice(n,x);
            	x[i] += h;
            	//printf("x[%d] = %lf\n",i,x[i]);
            	double f1 = dixonPrice(n,x);
            	x[i] = original_i;
            	x[i] -= h;
            	//printf("x[%d] = %lf\n",i,x[i]);
            	double f2 = dixonPrice(n,x);
            	Hessian[i][j] = ((f1 - f0 + f2) / (h * h));
            	//printf("f0 = %lf f1 = %lf f2 = %lf Hessian[%d][%d] = %lf\n",f0,f1,f2,i,j,Hessian[i][j]);
            }
            
            else{                                                                                       //will find off-diagonal elements
            // Compute f(x + h, x_i)
            	x[i] += h;
            	x[j] += h;
            	double f1 = dixonPrice(n, x);

            // Compute f(x - h, x_i)
            	x[i] = original_i;
            	x[i] += h;
            	x[j] = original_j;
            	x[j] -= h;
            	double f2 = dixonPrice(n, x);

            // Compute f(x, x_j + h)
            	x[i] = original_i;
            	x[i] -= h;
            	x[j] = original_j;
            	x[j] += h;
            	double f3 = dixonPrice(n, x);

            // Compute f(x, x_j - h)
            	x[i] = original_i;
            	x[i] -=h;
            	x[j] = original_j;
            	x[j] -=h;
            	double f4 = dixonPrice(n, x);

            // Calculate the second partial derivative
            	Hessian[i][j] = ((f1 - f2 - f3 + f4) / (4 * h * h));
	    }
            // Restore the original values
            x[i] = original_i;
            x[j] = original_j;
            printf("%10.3lf ", Hessian[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}
*/
/*
void TfindHessian(int n, double x[], double Hessian[][n], double h) {
     printf("//////////////////////////////////////////////////////Required Hessian Matrix/////////////////////////////////////////////\n");

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            double original_i = x[i];
            double original_j = x[j];
	    //printf("I am Ist oriI = %lf oriJ = %lf\n",original_i,original_j);
            if (i==j){                                                             //will find diagonal elements
            	double f0 = 2*trid(n,x);
            	x[i] += h;
            	//printf("x[%d] = %lf\n",i,x[i]);
            	double f1 = trid(n,x);
            	x[i] = original_i;
            	x[i] -= h;
            	//printf("x[%d] = %lf\n",i,x[i]);
            	double f2 = trid(n,x);
            	Hessian[i][j] = ((f1 - f0 + f2) / (h * h));
            	//printf("f0 = %lf f1 = %lf f2 = %lf Hessian[%d][%d] = %lf\n",f0,f1,f2,i,j,Hessian[i][j]);
            }
            
            else{                                                                  //will find off-diagonal elements
            // Compute f(x + h, x_i)
            	x[i] += h;
            	x[j] += h;
            	double f1 = trid(n, x);

            // Compute f(x - h, x_i)
            	x[i] = original_i;
            	x[i] += h;
            	x[j] = original_j;
            	x[j] -= h;
            	double f2 = trid(n, x);

            // Compute f(x, x_j + h)
            	x[i] = original_i;
            	x[i] -= h;
            	x[j] = original_j;
            	x[j] += h;
            	double f3 = trid(n, x);

            // Compute f(x, x_j - h)
            	x[i] = original_i;
            	x[i] -=h;
            	x[j] = original_j;
            	x[j] -=h;
            	double f4 = trid(n, x);

            // Calculate the second partial derivative
            	Hessian[i][j] = ((f1 - f2 - f3 + f4) / (4 * h * h));
	    }
            // Restore the original values
            x[i] = original_i;
            x[j] = original_j;
            printf("%10.3lf ", Hessian[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}
*/

void ZfindHessian(int n, double x[], double Hessian[][n], double h) {
    printf("//////////////////////////////////////////////////////Required Hessian Matrix/////////////////////////////////////////////\n");

    for (int i = 0; i < n; i++) {                        //capturing the row of hessian matrix
        for (int j = 0; j < n; j++) {                   //capturing the column hessian matrix
            double original_i = x[i];        //saving the the values of initial point
            double original_j = x[j];
	    //printf("I am Ist oriI = %lf oriJ = %lf\n",original_i,original_j);
            
            ////////////////////////////////////////////////////////////////////////////////Numerically finding the each element of Hessian matrix///////////////////////////////////////////////
            if (i==j){                                                                      //will find diagonal elements
            	double f0 = 2*zakharov(n,x);
            	x[i] += h;
            	//printf("x[%d] = %lf\n",i,x[i]);
            	double f1 = zakharov(n,x);
            	x[i] = original_i;
            	x[i] -= h;
            	//printf("x[%d] = %lf\n",i,x[i]);
            	double f2 = zakharov(n,x);
            	Hessian[i][j] = ((f1 - f0 + f2) / (h * h));
            	//printf("f0 = %lf f1 = %lf f2 = %lf Hessian[%d][%d] = %lf\n",f0,f1,f2,i,j,Hessian[i][j]);
            }
            
            else{                                                                                //will find off-diagonal elements
            // Compute f(x_i + h, x_j + h)
            	x[i] += h;
            	x[j] += h;
            	double f1 = zakharov(n, x);

            // Compute f(x_i + h, x_j - h)
            	x[i] = original_i;
            	x[i] += h;
            	x[j] = original_j;
            	x[j] -= h;
            	double f2 = zakharov(n, x);

            // Compute f(x_i - h, x_j + h)
            	x[i] = original_i;
            	x[i] -= h;
            	x[j] = original_j;
            	x[j] += h;
            	double f3 = zakharov(n, x);

            // Compute f(x_i - h, x_j - h)
            	x[i] = original_i;
            	x[i] -=h;
            	x[j] = original_j;
            	x[j] -=h;
            	double f4 = zakharov(n, x);

            // Calculate the second partial derivative
            	Hessian[i][j] = ((f1 - f2 - f3 + f4) / (4 * h * h));
	    }
            // Restore the original values
            x[i] = original_i;
            x[j] = original_j;
            printf("%10.3lf ", Hessian[i][j]);                                           //will print the HESSIAN MATRIX
        }
        printf("\n");
    }
    printf("\n");
}
/////////////////////////////////////////////////will do elementary row operations (WITHOUT PARTIAL PIVOTING) to convert the hessian into UPPER TRIANGULAR MATRIX///////////////////////////////////////////////

void check(int n, double Hessian[][n],double matrix_inv[][n]){                       
	
	for (int j = 1;j<=n;j++){                            
		
		for (int i=1;i<=n;i++){
			if (i>j){
				double constant = (Hessian[i-1][j-1]) / (Hessian[j-1][j-1]);
				for (int k=1;k<=n;k++){
					Hessian[i-1][k-1] = (Hessian[i-1][k-1]) - (constant*Hessian[j-1][k-1]);            //updating elements of (i-1)th row 
					//printf("%10lf ",Hessian[i-1][k-1]);
					//fprintf(fp1,"%10lf ",Hessian[i-1][k-1]);
				}
				//printf("\n");
				//fprintf(fp1,"/n");
			}
			else{                  //row won't change therefore we have to print it
				for (int a=i;a<=i;a++){              //for n
					for (int b=1;b<=n;b++){           //for columns
						double ele = Hessian[a-1][b-1];
						//printf("%10lf ",ele);
						//fprintf(fp1,"%10lf ",ele);
					}
					//printf("\n");
					//fprintf(fp1,"\n");
				}	
			}
			
		}
		//printf("\n");
		//fprintf(fp1,"\n");
		
	}
	///////////////////////////////////////////print the Upper triangular form of Hessian Matrix///////////////////////////////////////////////////
	printf("/////////////////////////////////////////////////Upper triangular form of Hessian Matrix/////////////////////////////////////////////\n");
	for(int i=0;i<n;i++){
		for(int j=0;j<n;j++){
			printf("%10.3lf ",Hessian[i][j]);
		}
		printf("\n");
	
	}
	printf("\n");
	/////////////////////////////////////////Checking the diagonal elements of Hessian matrix//////////////////////////////////////////////
	double hessianDiagonal[n];                 //1D array will store the DIAGONAL ELEMENTS OF Hessian Matrix
	for(int i=0;i<n;i++){
		hessianDiagonal[i] = Hessian[i][i];	
			
	}
	double min = hessianDiagonal[0];                              //variable will store minimum diagonal element of Hessian Matrix
	for(int j =1;j<n;j++){                                   
		if (min> hessianDiagonal[j]){
			min = hessianDiagonal[j];
		}
		
	}	
	if (min>=0){                                                              //if minimum diagonal element >=0 => all diagonal elements of Hessian matrix >=0
		printf("Hessian Matrix is Positive Definite/semidefinte\n");
		inverseUpperTriangular(n,Hessian,matrix_inv);
	}
	else{                                                                          
		printf("Hessian Matrix is Negative Definite/semidefinte\n");       //if minimum diagonal element < 0 => AT LEAST ONE diagonal elements of Hessian matrix < 0
	}

}
////////////////////////////////////////////////////////////////////////////////////finding the inverse of Hessian Matrix//////////////////////////////////////////////////////
void inverseUpperTriangular(int n,double UpperTriangular[][n], double matrix_inv[][n]) {
    for (int i = 0; i < n; i++) {
        // Initialize the inverse matrix column-wise
        for (int j = 0; j < n; j++) {
            matrix_inv[i][j] = 0.0;
        }

        // Diagonal elements of the inverse are 1 over the diagonal elements of U
        matrix_inv[i][i] = 1.0 / UpperTriangular[i][i];

        // Compute the rest of the elements in the current column
        for (int j = i - 1; j >= 0; j--) {
            double sum = 0.0;
            for (int k = j + 1; k < n; k++) {
                sum += UpperTriangular[j][k] * matrix_inv[k][i];
            }
            matrix_inv[j][i] = -sum / UpperTriangular[j][j];
        }
    }
    printf("//////////////////////////Inverse of Hessian//////////////////////////////////////////////\n");
    
    for(int i = 0;i<n;i++){                   //printing the elements of hessian inverse
    		for(int j= 0;j<n;j++){
    			printf("%10.3lf ", matrix_inv[i][j]);
    	}
    	printf("\n");
    
    }
}

double newton(int n, double x[], double grad[],double Hessian[][n], double matrix_inv[][n],double h,double product[]){
	double s[n];                 //1D array,  will SEARCH DIRECTION
	double m =1000;  // max number of iterations
	int k = 0;
	double E1 , E2 = 0.001;
	////step-2/////
	Zgradient(n,x,grad,h);
	/////step-3/////
	s[k] = -matrixMult(n,matrix_inv,grad,product);
	printf("%lf %lf\n",s[0],s[1]);
	while((norm(n,grad) > E1) && k<m){
			 //s[k] = -matrixMult(n,matrix_inv,grad,product);
	
	
	}

}

double norm(int n, double vector[]){                         //vector of size n, will find its length/norm
	double sum =0;
	for(int i=0;i<n;i++){
		sum = sum + pow(vector[i],2);
	}
	double length = pow(sum,0.5);
	return length;

}


double matrixMult(int n,double A[][n], double B[],double product[]){                   //will multiple a 2d array-A (size = n i.e. hessian inverse) with column vector - B (size =n, i.e. gradient) 
	for (int i=0;i<n;i++){
		double sum=0;
		for (int j=0;j<n;j++){
			sum = sum + A[i][j]*B[j];
		
		}
		product[i] = sum;
	
	}
	for(int i=0;i<n;i++){
		printf("%.3lf\n",product[i]);
	}
}
