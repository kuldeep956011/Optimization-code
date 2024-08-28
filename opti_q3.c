




/////submitted by Kuldeep singh(234103005) & Abhijeet(234103001)////////////////





#include<stdio.h>
#include<stdlib.h>
#include<math.h>
double hB(int n, double x[]);               // n= no of variables x = vector to store initial values of the variable
double sumSquare(int n, double x[]);
double rosenBrock(int n, double x[]);
double dixonPrice(int n, double x[]);
double trid(int n, double x[]);
double zakharov(int n, double x[]);


void gradient(int n, double x[], double grad[], double h);          // n= no of variables x = vector to store initial values of the variable, grad = vectore to store gradient values
void Sgradient(int n, double x[], double grad[], double h);
void Rgradient(int n, double x[], double grad[], double h);
void DPgradient(int n, double x[], double grad[], double h);
void Tgradient(int n, double x[], double grad[], double h);
void Zgradient(int n, double x[], double grad[], double h);

void findHessian(int n, double x[], double Hessian[][n], double h);   //no need of gradient, only x[],h and function definition needed
void SfindHessian(int n, double x[], double Hessian[][n], double h);
void RfindHessian(int n, double x[], double Hessian[][n], double h); // n= no of variables x = vector to store initial values of the variable, Hessian = matrix to store elements of Hessian
void DPfindHessian(int n, double x[], double Hessian[][n], double h);
void TfindHessian(int n, double x[], double Hessian[][n], double h);
void ZfindHessian(int n, double x[], double Hessian[][n], double h);



double UpperTriangular(int n,double A[][n],double B[][n]);                //A = known matrix(hessian) B=upper traigular form of known matrix(A)

double matrixInverse(int n , double A[][n],double B[][n],double C[][n]);           //A = known matrix(hessian) C =copy matrix(of known matrix i.e. A) B=Inverse of known matrix(A = hessian)
void swapRows(int n,double A[][n], int row1, int row2);
void scaleRow(int n,double A[][n], int row, double factor);
void subtractRows(int n,double A[][n], int destRow, int srcRow, double factor);

double length(int n, double vector[]);                       //finding the length/magniture of a vector
double matrixMult(int n,double A[][n], double B[],double re[]);          //multiplying A matrix(hessian inverse) with B vector(gradient) and storing result in r vector ( search direction) 
double dotProduct(int n, double Vec1[],double Vec2[]);                  //finding DOT PRODUCT of 2 vectors Vec1 and Vec2
double newton(int n, double a[],double g[],double H[][n],double u[][n],double mv[][n],double c[][n],double re[],double x_new[],double a1,double b1,double h,double r);             //a=vector(initial guess) , g =vector(gradient vector),H=Hessian matrix,u=upper matrix to store UPPER TRIANGULAR for of Hessian, mv=matrix to store Matrix Inverse(hessian Inverse), c = copy matrix(copy the Hessian Matrix), re = vector will storre RESULT OF MATRIX MULTIPLICATION(hessian inverse & gradient/SERACH DIRECTION),x_new=vector(to store new point),s=vector = SERACH DIIRECTION, (a1,b1)=range for variable(alpha),  (h=0.0001=small value to find grad,hessian) 
double randNum();               //will return random float between [0,1]
double randInterval(double a1, double b1);                 // //will return random float between [a1,b1]
double BoundingPhase(int n, double a1, double b1,double x[],double x_new[],double s[],double r);       // n = number of variables, (a1,b1)=range of variable, x=vector(initial point), x_new=vector(to store new point),s=vector = SERACH DIIRECTION
double intervalHalf(int n,double a1,double b1,double x[],double x_new[],double s[],double r);           // n = number of variables, (a1,b1)=range of variable obtained from BOUNDING PHASE METHOD, x=vector(initial point), x_new=vector(to store new point),s=vector = SERACH DIIRECTION
double objectiveFunction(int n, double alp, double x[], double x_new[],double s[],double r);        //n = number of variables, alp =value of alpha ,x=vector(initial point), x_new=vector(to store new point),s=vector = SERACH DIIRECTION



FILE* fp;               //file pointer to read the input
FILE* fp1;                     //file pointer to store the results


int fevalu;

double Fx(int n, double x[]) {
    double result = 0.0;
    result =  x[0] + x[1] + x[2];
    printf("Fx = %lf\n",result);
    return result;
}

double Gx1(int n, double x[]) {
    double result = 0.0;
    
    result = 1.0 - 0.0025*(x[3] + x[5]);
        
   // printf("Gx1 = %lf\n",result);
    if(result>=0.0){
    return 0;
	}
    else{
    	return result;
    }
}

double Gx2(int n, double x[]) {
    double result = 0.0;
    result = 1.0 - 0.0025*(-x[3] + x[4] + x[6]);
   // printf("Gx2 = %lf\n",result);
    if(result>=0.0){
    return 0;
	}
    else{
    	return result;
    }
    
}

double Gx3(int n, double x[]) {
    double result = 0.0;
    
    result = 1.0 - 0.01*(-x[5] + x[7]);
    
   // printf("Gx1 = %lf\n",result);
    if(result>=0.0){
    return 0;
	}
    else{
    	return result;
    }
}

double Gx4(int n, double x[]) {
    double result = 0.0;
    
    result = -100.0*x[0] + x[0]*x[5] - 833.33252*x[3] + 83333.33;
    
    //printf("Gx1 = %lf\n",result);
    if(result>=0.0){
    return 0;
	}
    else{
    	return result;
    }
}

double Gx5(int n, double x[]) {
    double result = 0.0;
    
    result = - x[1]*x[3] + x[1]*x[6] + 1250.0*x[3] - 1250.0*x[4];
    
   // printf("Gx1 = %lf\n",result);
    if(result>=0.0){
    return 0;
	}
    else{
    	return result;
    }
}

double Gx6(int n, double x[]) {
    double result = 0.0;
    
    result = - x[2]*x[4] + x[3]*x[7] + 2500.0*x[4] - 1250000;
    
    //printf("Gx1 = %lf\n",result);
    if(result>=0.0){
    return 0;
	}
    else{
    	return result;
    }
}


double omega(int n, double x[], double r) {
    double result = 0.0;
    double result1,result2;
    
	result = r*(pow(Gx1(n,x),2) + pow(Gx2(n,x),2)+ pow(Gx3(n,x),2)+ pow(Gx4(n,x),2)+ pow(Gx5(n,x),2)+ pow(Gx6(n,x),2));
	printf("omega = %lf\n",result);
	return result;
}

double Px(int n, double x[], double r) {

	fevalu +=1;
    double result = 0.0;
    result= Fx(n,x) + omega(n,x,r);
    
   printf("Px = %lf\n",result);
	return result;
}

void Pgradient(int n, double x[], double grad[], double h, double r) {
    for (int i = 1; i <= n; i++) {
        double original = x[i-1];
        x[i-1] += h;
        double f_plus_h = Px(n,x,r);
        x[i-1] = original - h;
        double f_minus_h = Px(n,x,r);
        x[i-1] = original;
        grad[i-1] = (f_plus_h - f_minus_h) / (2 * h);
        fprintf(fp1,"Gradient of f with respect x[%d] = %lf\n",i,grad[i-1]);
    }
    fevalu = fevalu + 2*n;
    fprintf(fp1,"\n");
    
} 


void PfindHessian(int n, double x[], double Hessian[][n], double h, double r) {
     printf("//////////////////////////////////////////////////////Required Hessian Matrix/////////////////////////////////////////////\n");

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            double original_i = x[i];
            double original_j = x[j];
	    //printf("I am Ist oriI = %lf oriJ = %lf\n",original_i,original_j);
            if (i==j){                                              //will find diagonal elements
            	double f0 = 2*Px(n,x,r);
            	x[i] += h;
            	//printf("x[%d] = %lf\n",i,x[i]);
            	double f1 = Px(n,x,r);
            	x[i] = original_i;
            	x[i] -= h;
            	//printf("x[%d] = %lf\n",i,x[i]);
            	double f2 = Px(n,x,r);
            	Hessian[i][j] = ((f1 - f0 + f2) / (h * h));
            	//printf("f0 = %lf f1 = %lf f2 = %lf Hessian[%d][%d] = %lf\n",f0,f1,f2,i,j,Hessian[i][j]);
            }
            
            else{                                                   //will find off-diagonal elements
            // Compute f(x + h, x_i)
            	x[i] += h;
            	x[j] += h;
            	double f1 = Px(n, x, r);

            // Compute f(x - h, x_i)
            	x[i] = original_i;
            	x[i] += h;
            	x[j] = original_j;
            	x[j] -= h;
            	double f2 = Px(n, x, r);

            // Compute f(x, x_j + h)
            	x[i] = original_i;
            	x[i] -= h;
            	x[j] = original_j;
            	x[j] += h;
            	double f3 = Px(n, x, r);

            // Compute f(x, x_j - h)
            	x[i] = original_i;
            	x[i] -=h;
            	x[j] = original_j;
            	x[j] -=h;
            	double f4 = Px(n, x, r);

            // Calculate the second partial derivative
            	Hessian[i][j] = ((f1 - f2 - f3 + f4) / (4 * h * h));
	    }
            // Restore the original values
            x[i] = original_i;
            x[j] = original_j;
            fprintf(fp1,"%10.3lf ", Hessian[i][j]);
        }
        fprintf(fp1,"\n");
    }
    fprintf(fp1,"\n");
    fevalu = fevalu + 2*n + 1; 
}



double UpperTriangular(int n,double A[][n],double B[][n]){                                //A = known matrix(hessian), B = upper triangular form of KNOWN MATRIX (A=HESSIAN)
double g;                                                  //useful in checking positive definiteness
fprintf(fp1,"-----------Hessian Matrix------------\n");
//printf("-----------Hessian Matrix------------\n");
for(int i=0;i<n;i++){
	for(int j=0;j<n;j++){
		
		fprintf(fp1,"%10.3lf ",A[i][j]);
		//printf("%10.3lf ",A[i][j]);
	}
	fprintf(fp1,"\n");
	//printf("\n");
}
fprintf(fp1,"\n");
//printf("\n");
//printf("---------------------Copying of Hessian in another matrix ----------------\n");
for(int i=0;i<n;i++){
	for(int j=0;j<n;j++){
		B[i][j] = A[i][j];
		//printf("%.3lf ",B[i][j]);
	}
	//printf("\n");

}
//printf("\n");
////////////////////////////////////////upperTriangular Conversion using elementary row operations//////////////////////
for (int j = 1;j<=n;j++){
		
		for (int i=1;i<=n;i++){
			if (i>j){
				double constant = (B[i-1][j-1]) / (B[j-1][j-1]);
				for (int k=1;k<=n;k++){
					B[i-1][k-1] = (B[i-1][k-1]) - (constant*B[j-1][k-1]);            //updating elements of (i-1)th row 
					
				}
				
			}
			else{                  //row won't change therefore we have to print it
				for (int a=i;a<=i;a++){              //for n
					for (int b=1;b<=n;b++){           //for columns
						double ele = B[a-1][b-1];
					}
					
				}	
			}
			
		}
		
		
}
fprintf(fp1,"----------------Upper Traigular form of Matrix(Hessian)------------------\n");
//printf("----------------Upper Traigular form of Matrix(Hessian)------------------\n");
for(int i=0;i<n;i++){
	for(int j=0;j<n;j++){
		
		fprintf(fp1,"%10.3lf ",B[i][j]);
		//printf("%10.3lf ",B[i][j]);

	}
	fprintf(fp1,"\n");
	//printf("\n");

}
fprintf(fp1,"\n");
//printf("\n");
/////////////////////////Checking Positve definiteness/////////////////////////////////////
double hessianDiagonal[n];                 //1D array will store the DIAGONAL ELEMENTS OF Upper triagular  Matrix
	for(int i=0;i<n;i++){
		hessianDiagonal[i] = B[i][i];	
			
	}
	double min = hessianDiagonal[0];                              //variable will store minimum diagonal element of Hessian Matrix
	for(int j =1;j<n;j++){                                   
		if (min> hessianDiagonal[j]){
			min = hessianDiagonal[j];
		}
		
	}	
	if (min>=0){                                                              //if minimum diagonal element >=0 => all diagonal elements of Hessian matrix >=0
		fprintf(fp1,"Hessian Matrix is Positive Definite/semidefinte\n");
		//printf("Hessian Matrix is Positive Definite/semidefinte\n");
		g=1.0;
		
	}
	else{                                                                          
		fprintf(fp1,"Hessian Matrix is Negative Definite/semidefinte\n");       //if minimum diagonal element < 0 => AT LEAST ONE diagonal elements of Hessian matrix < 0
		//printf("Hessian Matrix is Negative Definite/semidefinte\n");
		g=0.0;
	}
/*
printf("--------------------Hessian remain same after upper triangular operations----------------\n");
for(int i=0;i<n;i++){
	for(int j=0;j<n;j++){
		
		printf("%.3lf ",A[i][j]);
	}
	printf("\n");

}
printf("\n");
*/
return g;
}




double matrixInverse(int n,double A[][n],double B[][n],double C[][n]){                //A = known matrix(hessian) , B=will store elements of INVERSE, C = copy matrix(of A)
//printf("---------------------Copying of Hessian in another matrix C ----------------\n");
for(int i=0;i<n;i++){
	for(int j=0;j<n;j++){
		C[i][j] = A[i][j];
		//printf("%.3lf ",C[i][j]);
	}
	//printf("\n");

}
//printf("\n");

///////////////making B Idntity matrix////////////////
	for(int i=0;i<n;i++){
		for(int j=0;j<n;j++){
			if(i==j){
				B[i][j] = 1.0;
			}
			else{
				B[i][j] = 0.0;
			}
		
		}
	
	}
	
	// Perform Gauss-Jordan elimination
      for (int i = 0; i < n; i++) {
        // Find the pivot element and swap rows
        int pivotRow = i;
        while (C[pivotRow][i] == 0) {
            pivotRow++;
            if (pivotRow == n) {
                fprintf(fp1,"Matrix is singular. Cannot find the inverse.\n");
                //printf("Matrix is singular. Cannot find the inverse.\n");
                return 1.0;
            }
        }
        swapRows(n,C, i, pivotRow);
        swapRows(n,B, i, pivotRow);

        // Scale the pivot row
        float pivot = C[i][i];
        scaleRow(n,C, i, 1.0 / pivot);
        scaleRow(n,B, i, 1.0 / pivot);

        // Eliminate other rows
        for (int j = 0; j < n; j++) {
            if (j != i) {
                float factor = C[j][i];
                subtractRows(n,C, j, i, factor);
                subtractRows(n,B, j, i, factor);
            }
        }
    }
    	fprintf(fp1,"------------Inverse--------------------\n");
    	//printf("------------Inverse--------------------\n");
	for(int i=0;i<n;i++){
		for(int j=0;j<n;j++){	
			fprintf(fp1,"%10.3lf ",B[i][j]);
			//printf("%10.3lf ",B[i][j]);
		}
		fprintf(fp1,"\n");
		//printf("\n");

	}
	fprintf(fp1,"\n");
	//printf("\n");
	/*
	printf("------------copy matrix after gauss jordan-------\n");
	printf("\n");
	for(int i=0;i<n;i++){
		for(int j=0;j<n;j++){	
			printf("%.3lf ",C[i][j]);
		}
		printf("\n");

	}
	printf("\n");
	*/
	/*
	printf("-----------------Original Matrix----------\n");
	for(int i=0;i<n;i++){
		for(int j=0;j<n;j++){	
			printf("%.3lf ",A[i][j]);
		}
		printf("\n");

	}
	printf("\n");
	*/
}


void swapRows(int n,double A[][n], int row1, int row2) {
    for (int i = 0; i < n; i++) {
        float temp = A[row1][i];
        A[row1][i] = A[row2][i];
        A[row2][i] = temp;
    }
}

void scaleRow(int n,double A[][n], int row, double factor) {
    for (int i = 0; i < n; i++) {
        A[row][i] *= factor;
    }
}

void subtractRows(int n,double A[][n], int destRow, int srcRow, double factor) {
    for (int i = 0; i < n; i++) {
        A[destRow][i] -= factor * A[srcRow][i];
    }
}

double newton(int n, double a[],double g[],double H[][n],double u[][n],double mv[][n],double c[][n],double re[],double x_new[],double a1,double b1,double h,double r){            //a=vector(initial guess) , g =vector(gradient vector),H=Hessian matrix,u=upper matrix to store UPPER TRIANGULAR for of Hessian, mv=matrix to store Matrix Inverse(hessian Inverse), c = copy matrix(copy the Hessian Matrix), re = vector will storre RESULT OF MATRIX MULTIPLICATION(hessian inverse & gradient or SEARCH DIRECTION), x_new = vector will STORE NEW POINT , (a1,b1) = range for variable (alpha),(h=0.0001=small value to find grad,hessian) 
    
    int k=0;            //counter for iteration number
    printf("r=%lf\n",r);
	int m=1000;       //setting MAXIMUM NUMBER OF ITERATIONS
	double E = 0.001;  //TERMINATION PARAMETER
	int terminate;   //use to BREAK THE LOOP CORRECTLY
	int count=0;          //helpful in checking linear dependency between 2 search directions
	double angle;                     // find the angle between 2 consecutive search directions
	double functionValue;    //store the value of function at a particular point
	
	double sPrevious[n];              //1D array will store PREVIOUS SEARCH DIRECTION, useful in checking LINEAR INDEPENDENCY between 2 directions
	double descent;      //check descency of direction when hessian isn't positive definite
	double magnitude;     //will store the value of magniture of vector(particulaly GRADIENT)
	double difference[n];             // 1D array to store the DIsTANCE BETWEEN CURRENT POINT AND NEW POINT, useful in terminating the program	
	
	//gradient(n,a,g,h);
	//Sgradient(n,a,g,h);
	Pgradient(n,a,g,h,r);
	//Rgradient(n,a,g,h);
	//DPgradient(n,a,g,h);
	//Tgradient(n,a,g,h);
	//Zgradient(n,a,g,h);
	magnitude = length(n,g);
	
while(magnitude>E && k<m){
	fevalu =0;
	r = 10.0*r;	
        printf("r = ////////////////////////%lf\n",r);
	fprintf(fp1,"----------------------k = %d and Iteration Number = %d--------------------------------\n",k,k+1);
	//printf("----------------------k = %d and Iteration Number = %d--------------------------------\n",k,k+1);
	//functionValue = hB(n,a);
	//functionValue = sumSquare(n,a);
	functionValue = Px(n,a,r);
	//functionValue = rosenBrock(n,a);
	//functionValue = dixonPrice(n,a);
	//functionValue = trid(n,a);
	//functionValue = zakharov(n,a);
	fprintf(fp1,"Function Value = %lf\n",functionValue);
	//printf("Function Value = %lf\n",functionValue);
	
	//gradient(n,a,g,h);
	//Sgradient(n,a,g,h);
	Pgradient(n,a,g,h,r);
	//Rgradient(n,a,g,h);
	//DPgradient(n,a,g,h);
	//Tgradient(n,a,g,h);
	//Zgradient(n,a,g,h);
	magnitude = length(n,g);
	//////////Printing previous Serach Direction//////////////////
	for(int i=0;i<n;i++){
		fprintf(fp1,"sPrevious[%d] = %lf\n",i,sPrevious[i]);
		//printf("sPrevious[%d] = %lf\n",i,sPrevious[i]);

	}
	
	//findHessian(n,a,H,h);
	//SfindHessian(n,a,H,h);
	PfindHessian(n,a,H,h,r);
	//RfindHessian(n,a,H,h);
	//DPfindHessian(n,a,H,h);
	//TfindHessian(n,a,H,h);
	//ZfindHessian(n,a,H,h);
	
	double out = UpperTriangular(n,H,u);      //output 1 or 0,  1, if hessian is positve definite and 0, if hessian isn't positive definite/semidefinite
	matrixInverse(n,H,mv,c);
	//printf("out = %lf\n",out);
	if(out==0){         //out = 0, hessian isn't positve definite
		matrixMult(n,mv,g,re);            //multiplying Matrix Inverse=mv=(hessian inverse) & a Vector=g=(gradient), storing result in re(vector)	
		descent = dotProduct(n,re,g);       //dot product of re(hessian inverse*gradient)*g(gradient)
		if(descent<=0){
			fprintf(fp1,"Even though the Hessian Matrix isn't Positive Definite, the DIRECTION is DESCENT with value = %lf\n",descent);
			//printf("Even though the Hessian Matrix isn't Positive Definite, the DIRECTION is DESCENT with value = %lf\n",descent);

			for(int i=0;i<n;i++){
				fprintf(fp1,"New Direction[%d] = %lf\n",i,re[i]);
				//printf("New Direction[%d] = %lf\n",i,re[i]);

			}
			double alphaOptimal = BoundingPhase(n,a1,b1,a,x_new,re,r);                  //CALLING BOUNDING PHASE and withing bounding phase we are calling INTERVAL HALF METHOD
			//printf("alphaOptimal = %lf\n",alphaOptimal);
			fprintf(fp1,"alphaOptimal = %lf\n",alphaOptimal);              //storing alphaOptimal in OUTPUT FILE
			//printf("alphaOptimal = %lf\n",alphaOptimal);
			//////////////////////finding the angle-LINEAR INDEPENDENCY CHECK/////////////////////////////////////////////////////
			if(count>=1){           //skiping 1st iteration, because after 1st iteration(when count =0) we'll be having values of "s" & "sPrevious"
				double directionProduct = dotProduct(n,re,sPrevious);
				double lengthProduct = length(n,re)*length(n,sPrevious);
				double value = acos(directionProduct / lengthProduct);   
				angle = (value*180) / 3.14159;						//changed from RADIANS TO DEGREES 
				fprintf(fp1,"directionProduct = %lf lengthProduct = %lf\n",directionProduct,lengthProduct);
				//printf("directionProduct = %lf lengthProduct = %lf\n",directionProduct,lengthProduct);	  
				fprintf(fp1,"Angle = %.3lf\n",angle);          
				//printf("Angle = %.3lf\n",angle);          

				if(angle<=5.0){
					terminate = 4;
					fprintf(fp1,"Search Directions are Linearly Dependent\n");
					//printf("Search Directions are Linearly Dependent\n");

					break;	
				}
			}
			count += 1;
			
			for (int j=0;j<n;j++){                                
				sPrevious[j] =  re[j];      //updating the SEARCH DIRECTION
			}
			///////////////////////////////////////finding the new point using previous point and search direction//////////////////////////////
			for(int i=0;i<n;i++){
				//printf("x[%d] = %lf alphaOptimal = %lf s[%d] = %lf\n",i,a[i],alphaOptimal,i,re[i]);
				x_new[i] = a[i] + alphaOptimal*re[i];
				//printf("x_new[%d] = %lf\n",i,x_new[i]);	
				fprintf(fp1,"x_new[%d] = %lf\n",i,x_new[i]);                //storing x_new point in OUTPUT FILE
				//printf("x_new[%d] = %lf\n",i,x_new[i]); 
				difference[i] = x_new[i] - a[i];               //checking distnace between 2 points
			
			}
			
			double distance = (length(n,difference) / length(n,a));
			fprintf(fp1,"Distance = %lf\n",distance);
			//printf("Distance = %lf\n",distance);
			if (distance > E){
				k = k+1;
				//updating the values of a[i] to x_new[i] for next iteration
				for(int i=0;i<n;i++){
					a[i] = x_new[i];
					
				}
				//functionValue = hB(n,a);                   //finding function value at UPDATED POINT
				//functionValue = sumSquare(n,a);
				functionValue = Px(n,a,r);
				//functionValue = rosenBrock(n,a);
				//functionValue = dixonPrice(n,a);
				//functionValue = trid(n,a);
				//functionValue = zakharov(n,a);
				fprintf(fp1,"Function Value = %lf\n",functionValue);
				//printf("Function Value = %lf\n",functionValue);

				for(int i=0;i<n;i++){
					fprintf(fp1,"UpdatedPoint[%d] = %lf\n",i,a[i]);
					//printf("UpdatedPoint[%d] = %lf\n",i,a[i]);

				}
				//gradient(n,a,g,h);
				//Sgradient(n,a,g,h);
				Pgradient(n,a,g,h,r);
				//Rgradient(n,a,g,h);
				//DPgradient(n,a,g,h);
				//Tgradient(n,a,g,h);
				//Zgradient(n,a,g,h);
				magnitude = length(n,g);
				terminate = 6 ;                   //when gradient magnitue<E, we will come out of loop
			}
			else{               //distancs <E, terminate
				terminate = 2;
				break;
			}
		}
		else{
			terminate = 1;
			fprintf(fp1,"Direction isn't Desecent, with value = %lf\n",descent);
			//printf("Direction isn't Desecent, with value = %lf\n",descent);
			fprintf(fp1,"Change the INITIAL Point in the INPUT FILE\n");
			//printf("Change the INITIAL Point in the INPUT FILE\n");

			break;
		
		}
		
	
	}
	else{                     //out = 1, hessian is positve definite
		matrixMult(n,mv,g,re);            //multiplying Matrix Inverse=mv=(hessian inverse) & a Vector=g=(gradient), stroing result in re(vector)	
		for(int i=0;i<n;i++){
			fprintf(fp1,"New Direction[%d] = %lf\n",i,re[i]);
			//printf("New Direction[%d] = %lf\n",i,re[i]);
		}
		double alphaOptimal = BoundingPhase(n,a1,b1,a,x_new,re,r);                  //CALLING BOUNDING PHASE and withing bounding phase we are calling INTERVAL HALF METHOD
		//printf("alphaOptimal = %lf\n",alphaOptimal);
		fprintf(fp1,"alphaOptimal = %lf\n",alphaOptimal);              //storing alphaOptimal in OUTPUT FILE
		//printf("alphaOptimal = %lf\n",alphaOptimal);
		
		//////////////////////finding the angle-LINEAR INDEPENDENCY CHECK/////////////////////////////////////////////////////
		if(count>=1){           //skiping 1st iteration, because after 1st iteration(when count =0) we'll be having values of "s" & "sPrevious"
			double directionProduct = dotProduct(n,re,sPrevious);
			double lengthProduct = length(n,re)*length(n,sPrevious);
			double value = acos(directionProduct / lengthProduct);   
			angle = (value*180) / 3.14159;						//changed from RADIANS TO DEGREES 
			fprintf(fp1,"directionProduct = %lf lengthProduct = %lf\n",directionProduct,lengthProduct);	  
			//printf("directionProduct = %lf lengthProduct = %lf\n",directionProduct,lengthProduct);	  

			fprintf(fp1,"Angle = %.3lf\n",angle);          
			//printf("Angle = %.3lf\n",angle);          

			if(angle<=5.0){
				terminate = 5;
				fprintf(fp1,"Search Directions are Linearly Dependent\n");
				//printf("Search Directions are Linearly Dependent\n");

				break;	
			}
		}
			count += 1;
			
		for (int j=0;j<n;j++){                                
			sPrevious[j] =  re[j];      //updating the SEARCH DIRECTION
		}
		///////////////////////////////////////finding the new point using previous point and search direction//////////////////////////////
		for(int i=0;i<n;i++){
			//printf("x[%d] = %lf alphaOptimal = %lf s[%d] = %lf\n",i,a[i],alphaOptimal,i,re[i]);
			x_new[i] = a[i] + alphaOptimal*re[i];
			
			//printf("x_new[%d] = %lf\n",i,x_new[i]);	
			fprintf(fp1,"x_new[%d] = %lf\n",i,x_new[i]);                //storing x_new point in OUTPUT FILE
			//printf("x_new[%d] = %lf\n",i,x_new[i]);
			difference[i] = x_new[i] - a[i];               //checking distnace between 2 points
			
		}
		
		//////////////////////////////////////////////////////////new point calulation over//////////////////////////////////////////////////	
		double distance = (length(n,difference) / length(n,a));
		fprintf(fp1,"Distance = %lf\n",distance);
		//printf("Distance = %lf\n",distance);
		
		if (distance > E){
			k = k+1;
			//updating the values of a[i] to x_new[i] for next iteration
			for(int i=0;i<n;i++){
				a[i] = x_new[i];
				
			}
			//functionValue = hB(n,a);                   //finding function value at UPDATED POINT
			//functionValue = sumSquare(n,a);
			functionValue = Px(n,a,r);
			//functionValue = rosenBrock(n,a);
			//functionValue = dixonPrice(n,a);
			//functionValue = trid(n,a);
			//functionValue = zakharov(n,a);
			fprintf(fp1,"Function Value = %lf\n",functionValue);
			//printf("Function Value = %lf\n",functionValue);
			for(int i=0;i<n;i++){
				fprintf(fp1,"UpdatedPoint[%d] = %lf\n",i,a[i]);
				//printf("UpdatedPoint[%d] = %lf\n",i,a[i]);

			}
			//gradient(n,a,g,h);
			//Sgradient(n,a,g,h);
			Pgradient(n,a,g,h,r);
			//Rgradient(n,a,g,h);
			//DPgradient(n,a,g,h);
			//Tgradient(n,a,g,h);
			//Zgradient(n,a,g,h);
			magnitude = length(n,g);
			terminate = 7;              //when gradient magnitue<E, we will come out of loop
		}
		else{               //distancs <E, terminate
			terminate = 3;
			break;
		}
	}
	fprintf(fp1,"Number of Function Evalutions = %d\n",fevalu);

}
/////////////////////////////if we don't enter the loop///////////////////
		fprintf(fp1,"The Optimal Point and Function Value are:\n");
		//printf("The Optimal Point and Function Value are:\n");

		for(int i=0;i<n;i++){
			fprintf(fp1,"x[%d] = %lf\n",i,a[i]);
			//printf("x[%d] = %lf\n",i,x_new[i]);
		}
		//functionValue = hB(n,a);
		//functionValue = sumSquare(n,a);
		functionValue = Px(n,a,r);
		fprintf(fp1,"Number of Function Evalutions = %d\n",fevalu);
		//functionValue = rosenBrock(n,a);
		//functionValue = dixonPrice(n,a);
		//functionValue = trid(n,a);
		//functionValue = zakharov(n,a);
		fprintf(fp1,"Function Value = %lf\n",functionValue);
		//printf("Function Value = %lf\n",functionValue);
/////////////////////////////////////////////////////////////////////////

	if(terminate == 2){
		fprintf(fp1,"The Optimal Point and Function Value are:\n");
		//printf("The Optimal Point and Function Value are:\n");

		for(int i=0;i<n;i++){
			fprintf(fp1,"x[%d] = %lf\n",i,x_new[i]);
			//printf("x[%d] = %lf\n",i,x_new[i]);
		}
		//functionValue = hB(n,a);
		//functionValue = sumSquare(n,a);
		functionValue = Px(n,a,r);
		fprintf(fp1,"Number of Function Evalutions = %d\n",fevalu);
		//functionValue = rosenBrock(n,a);
		//functionValue = dixonPrice(n,a);
		//functionValue = trid(n,a);
		//functionValue = zakharov(n,a);
		fprintf(fp1,"Function Value = %lf\n",functionValue);
		//printf("Function Value = %lf\n",functionValue);
	}

	else if (terminate==3){
		fprintf(fp1,"The Optimal Point and Function Value are:\n");
		//printf("The Optimal Point and Function Value are:\n");
		for(int i=0;i<n;i++){
			fprintf(fp1,"x[%d] = %lf\n",i,x_new[i]);
			//printf("x[%d] = %lf\n",i,x_new[i]);
		}
		//functionValue = hB(n,a);
		//functionValue = sumSquare(n,a);
		functionValue = Px(n,a,r);
		fprintf(fp1,"Number of Function Evalutions = %d\n",fevalu);
		//functionValue = rosenBrock(n,a);
		//functionValue = dixonPrice(n,a);
		//functionValue = trid(n,a);
		//functionValue = zakharov(n,a);
		fprintf(fp1,"Function Value = %lf\n",functionValue);
		//printf("Function Value = %lf\n",functionValue);
	
	}

	///when SEARCH DIRECTIONS ARE DEPENDENT/////////
	else if(terminate == 4){
		fprintf(fp1,"Change the Initial Points in the Input File and RE_RUN the program:\n");
		//printf("Change the Initial Points in the Input File and RE_RUN the program:");

		for(int i=0;i<n;i++){
			fprintf(fp1,"x_new[%d] = %lf\n",i,x_new[i]);
			//printf("x_new[%d] = %lf\n",i,x_new[i]);
		}
		//functionValue = hB(n,a);
		//functionValue = sumSquare(n,a);
		functionValue = Px(n,a,r);
		fprintf(fp1,"Number of Function Evalutions = %d\n",fevalu);
		//functionValue = rosenBrock(n,a);
		//functionValue = dixonPrice(n,a);
		//functionValue = trid(n,a);
		//functionValue = zakharov(n,a);
		fprintf(fp1,"Function Value = %lf\n",functionValue);
	}
	
	else if(terminate == 5){
		fprintf(fp1,"Change the Initial Points in the Input File and RE_RUN the program:\n");
		//printf("Change the Initial Points in the Input File and RE_RUN the program:");

		for(int i=0;i<n;i++){
			fprintf(fp1,"x_new[%d] = %lf\n",i,x_new[i]);
			//printf("x_new[%d] = %lf\n",i,x_new[i]);
		}
		//functionValue = hB(n,a);
		//functionValue = sumSquare(n,a);
		functionValue = Px(n,a,r);
		fprintf(fp1,"Number of Function Evalutions = %d\n",fevalu);
		//functionValue = rosenBrock(n,a);
		//functionValue = dixonPrice(n,a);
		//functionValue = trid(n,a);
		//functionValue = zakharov(n,a);
		fprintf(fp1,"Function Value = %lf\n",functionValue);
	}

	else if(terminate == 6 || terminate == 7 || k>=m){
		fprintf(fp1,"The Optimal Point and Function Value are:\n");
		//printf("The Optimal Point and Function Value are:\n");
		for(int i=0;i<n;i++){
			fprintf(fp1,"x[%d] = %lf\n",i,x_new[i]);
			//printf("x[%d] = %lf\n",i,x_new[i]);

		}
		//functionValue = hB(n,a);
		//functionValue = sumSquare(n,a);
		functionValue = Px(n,a,r);
		fprintf(fp1,"Number of Function Evalutions = %d\n",fevalu);
		//functionValue = rosenBrock(n,a);
		//functionValue = dixonPrice(n,a);
		//functionValue = trid(n,a);
		//functionValue = zakharov(n,a);
		fprintf(fp1,"Function Value = %lf\n",functionValue);
		//printf("Function Value = %lf\n",functionValue);
	}	
	
}





//////////////////////////////////////////////////////////////////////EXTRA USEFUL FUNCTIONS/////////////////////////////////////////////////////////////////////////////////////

double length(int n, double vector[]){                         //vector of size n, will find its length/norm
	double sum =0;
	for(int i=0;i<n;i++){
		sum = sum + pow(vector[i],2);
	}
	double magnitude = pow(sum,0.5);
	//printf("norm = %lf\n",length);
	return magnitude;

}


double matrixMult(int n,double A[][n], double B[],double re[]){          //will multiple a 2d array-A (size = n i.e. hessian inverse) with column vector - B (size =n, i.e. gradient) 
	for (int i=0;i<n;i++){
		double sum=0;
		for (int j=0;j<n;j++){
			sum = sum + A[i][j]*B[j];
		
		}
		re[i] = -sum;               // - negative sign added to satify the SEARCH direction requirement ( direction = -hessianInverse(A)*gradient(B))
	
	}
	/*for(int i=0;i<n;i++){
		printf("%.3lf\n",r[i]);
	}*/
}

double dotProduct(int n, double Vec1[],double Vec2[]){                               //will perform dot product between 2 vectors Veec1 & Vec2 of size n
	double sum = 0;
	for (int i=0;i<n;i++){
		sum += Vec1[i]*Vec2[i];
	
	}
	//printf("dot product = %lf\n",sum);
	return sum;
}




double BoundingPhase(int n, double a1, double b1,double x[],double x_new[],double s[],double r){       //(a1,b1)=range of variable, x=vector(initial point), x_new=vector(to store new point),s=vector = SERACH DIIRECTION
	double optimal;
	//FILE* fp;
	//fp = fopen("result.txt","w");
	//if(fp!=NULL){
	for(int i=0; i<1;i++){                            // will generate 10 initial guesses
		double xvalue[1000]={};                       //will store value of ALPHA at differnt iterations, size is kept large otherwise overflow error may occur
		//fprintf(fp,"\n//////////////////////////////Guess Number = %d////////////////////////////////////////////////\n\n",i+1);
		double delta = 0.001;
		double x0=randInterval(a1,b1);                         /// will generate initial guess for ALPHA between [a,b]
		xvalue[0] = x0;
		int k = 0;
		///////////////END OF STEP--1///////////////////////////////
		//printf("%lf\n",x0);
	
		double fl,f0,fm;
		fl = objectiveFunction(n,x0-fabs(delta),x,x_new,s,r);                       // f at x0-delta
		f0 = objectiveFunction(n,x0,x,x_new,s,r);                                      // f at x0
		fm = objectiveFunction(n,x0+fabs(delta),x,x_new,s,r);	                     // f at x0+ delta
		
		while(!((fl >= f0 && f0 >= fm) || (fl <= f0 && f0 <= fm))){                            //checking condition 
			
			
				x0 = randInterval(a1,b1);
				fl = objectiveFunction(n,x0-fabs(delta),x,x_new,s,r);                       
				f0 = objectiveFunction(n,x0,x,x_new,s,r);                                       
				fm = objectiveFunction(n,x0+fabs(delta),x,x_new,s,r);                       	
				//printf("x0 = %lf fl = %lf f0 = %lf fm = %lf delta = %lf\n",x0,fl,f0,fm,delta);
				//fprintf(fp, "x0 = %lf fl = %lf f0 = %lf fm = %lf delta = %lf\n",x0,fl,f0,fm,delta);
		}
			
			
			
		if (fl >= f0 && f0 >= fm){
			delta = delta;
			//printf("x0 = %lf fl = %lf f0 = %lf fm = %lf delta = %lf\n",x0,fl,f0,fm,delta);
			//fprintf(fp, "x0 = %lf fl = %lf f0 = %lf fm = %lf delta = %lf\n",x0,fl,f0,fm,delta);
			
				
		}
		else if (fl <= f0 && f0 <= fm){
			delta = (-1)*delta;
			//printf("x0 = %lf fl = %lf f0 = %lf fm = %lf delta = %lf\n",x0,fl,f0,fm,delta);
			//fprintf(fp, "x0 = %lf fl = %lf f0 = %lf fm = %lf delta = %lf\n",x0,fl,f0,fm,delta);
		}
		//////////	END OF STEP-2 //////////////////////////////////////////////////////////////
		double fprev = f0;                      //variable to store previous value of f
		int xstep=1;                           // counter for ALPHA at a particular iteration
		x0 = x0 + pow(2,k)*delta;             // updating value of x0
		
		
		////end of STEP-3/////////////////////////
		f0 = objectiveFunction(n,x0,x,x_new,s,r);
		while(f0 < fprev){
			//printf("\nf%d = %lf f%d = %lf x%d = %lf k = %d xvalue = %lf\n",xstep-1,fprev,xstep,f0,xstep,x0,k,xvalue[k]);
            		//fprintf(fp, "\nf%d = %lf f%d = %lf x%d = %lf k = %d xvalue = %lf\n",xstep-1,fprev,xstep,f0,xstep,x0,k,xvalue[k]);
				
			k=k+1;
			fprev = f0;
			xvalue[k]=x0;
			
			x0 = x0 + pow(2,k)*delta;
			xstep ++;
			f0 = objectiveFunction(n,x0,x,x_new,s,r);
			//printf("f%d = %lf f%d = %lf x%d = %lf k = %d xvalue = %lf\n",xstep-1,fprev,xstep,f0,xstep,x0,k,xvalue[k]);
			//fprintf(fp, "f%d = %lf f%d = %lf x%d = %lf k = %d xvalue = %lf\n",xstep-1,fprev,xstep,f0,xstep,x0,k,xvalue[k]);
			
			fevalu = fevalu +2;
		}
		//printf("Minima point lies (x%d,x%d) = (%lf, %lf)\n",k-1,k+1,xvalue[k-1],x0);
		//fprintf(fp, "Minima point lies (x%d,x%d) = (%lf, %lf)\n",k-1,k+1,xvalue[k-1],x0);
		double newx,newy;                           //interval for interval half method
		if (xvalue[k-1]>x0){                         //[newx,newy] , ensuring newx>newy
			newx = x0;
			newy = xvalue[k-1];
			//printf("\nnewx = %lf,newy = %lf\n",newx,newy);
			//fprintf(fp1,"newx = %lf,newy = %lf\n",newx,newy);
		}
		else if (xvalue[k-1]<x0){
			newx = xvalue[k-1];
			newy = x0;
			//printf("\nnewx = %lf,newy = %lf\n",newx,newy);
			//fprintf(fp1, "newx = %lf,newy = %lf\n",newx,newy);
		}
		optimal = intervalHalf(n,newx,newy,x,x_new,s,r);                             //passing the obtained bounding interval [a=newx,b=newy] in interval half to get ACCURATE BOUNDED INTERVAL for E=0.001 (termination condition)
		
	}
	
	return optimal;
//}
//fclose(fp);
}

double intervalHalf(int n,double a1,double b1,double x[],double x_new[],double s[],double r){          // n = number of variables, (a1,b1)=range of variable obtained from BOUNDING PHASE METHOD, x=vector(initial point), x_new=vector(to store new point),s=vector = SERACH DIIRECTION
	
	double E = 0.001;                 //E = 0.001 for termination of the loop
	double xm,L,f,x1,x2,fx1,fx2,fxm;                 // xm = midpoint of (x1,x2),  L = size of (x1,x2)
	int iteration = 0;  /* count iteration number*/
	
	//FILE* fp1;                //storing the results in an output file
	//fp1 = fopen("result.txt","a");
	
	//if(fp1!=NULL){
	L = b1 - a1;                                          //fx1, fx2, fxm are the function values at x1,x2 and xm respectively
	
	while(fabs(L) > E){       //E = 0.001 for termination of the loop
		iteration ++;
		//printf("\n********Iteration number =  %d ********\n",iteration);
		//fprintf(fp, "\n********Iteration number =  %d ********\n",iteration); 
		xm = (a1 + b1) / 2;
		//printf("xm = %lf\n",xm);
		//fprintf(fp, "xm = %lf\n",xm);
		fxm = objectiveFunction(n,xm,x,x_new,s,r);
		
		//printf("fxm = %lf\n",fxm);
		//fprintf(fp, "fxm = %lf\n",fxm);
		x1 = a1 + (L/4);
		
		//printf("x1 = %lf\n",x1);
		//fprintf(fp, "x1 = %lf\n",x1);
		
		
		x2 = b1 - (L/4);
		
		//printf("x2 = %lf\n",x2);
		//fprintf(fp, "x2 = %lf\n",x2);
		
		fx1 = objectiveFunction(n,x1,x,x_new,s,r);
		
		//printf("fx1 = %lf\n",fx1);
		//fprintf(fp, "fx1 = %lf\n",fx1);
		
		fx2 = objectiveFunction(n,x2,x,x_new,s,r);
		
		//printf("fx2 = %lf\n",fx2);
		//fprintf(fp, "fx2 = %lf\n",fx2);
		
		/* Eliminate region based on the following Conditions
		  Correspondily updating the values */
		if (fx1 < fxm){
			b1 = xm;     
			xm = x1;
			//printf("I am in 1st if\n");
		
		}
		
		else if (fx2 < fxm){         
			a1 = xm;
			xm = x2;
			//printf("I am in 2nd if\n");
		}
		else{
			a1=x1;
			b1=x2;
			xm = (a1+b1)/2;
			//printf("I am in else\n");
		}
		L = b1-a1;      //updating Value of L/
		
		//printf("Value of a1 = %lf , Value of b1 = %lf , value of xm = %lf , Value of L = %lf\n",a1,b1,xm,L);
		//break;
		fevalu = fevalu + 2;
	}
	//printf("\n");
	//fprintf(fp, "\n");
	//printf("\nMinima point will lie in (%lf,%lf) according to Interval Half Method.\n",a1,b1);        //print interval containing maxima
	//fprintf(fp1, "\nMinima point will lie in (%lf,%lf) according to Interval Half Method.\n",a1,b1); 
	//printf("\nobjective function value f(alpha*) = %lf at minima point alpha* = %lf\n",objectiveFunction(n,a1,x,x_new,s),a1);
	//fprintf(fp1, "\nobjective function value f(alpha*) = %lf at minima point alpha* = %lf\n",objectiveFunction(n,a1,x,x_new,s),a1);
	return a1;
	
//}
//fclose(fp1);

}


double objectiveFunction(int n, double alp, double x[], double x_new[],double s[], double r){       //n = number of variables, alp =value of alpha ,x=vector(initial point), x_new=vector(to store new point),s=vector = SERACH DIIRECTION	//printf("%lf\n",x);
	
	for (int i =1;i<=n;i++){            //loop reading the file and storing initial point
		//fscanf(fp,"%lf",&x[i-1]);           
		//printf("%lf\n",x[i-1]);
		x_new[i-1] = x[i-1] + alp*s[i-1];                
	}
	

	//double f =  hB(n,x_new);
	double f =  Px(n,x_new, r);
	return f;
	
}


double randNum(){                //will return random float between [0,1]
	return ( (double) rand() / (double) RAND_MAX);            //RAND_MAX is max random number generated by the machine
}

double randInterval(double a1, double b1){                 //will return random float between [a1,b1]
	double rand = randNum()*(b1-a1) + a1;
	return rand;

}



int main(){
fp = fopen("input.txt","r");
fp1 = fopen("sum10.txt","w");
int n;
fscanf(fp,"%d",&n);

double x[n];
double grad[n];
double Hessian[n][n];
double upper[n][n];
double matrixInv[n][n];
double copyMatrix[n][n];
double re[n];
double x_new[n];        //use to store new point

double h = 0.00001;              //use to find derivatives

////////////////////// INITIAL INTERVAL for ALPHA//////////////////////////////
	
	
	double a1=10;
	double b1= 10000.0;                                                      //initial interval for function variables as well as ALPHA
	double r=0.01;
	
	

for(int i=0;i<n;i++){
	fscanf(fp,"%lf",&x[i]);
	fprintf(fp1,"x[%d] = %lf\n",i,x[i]);
	printf("x[%d] = %lf\n",i,x[i]);

}

fprintf(fp1,"\n");

//Fx(n,x);
//Gx1(n,x);
//Gx2(n,x);
//omega(n,x,0.1);
//Px(n,x,0.1);
//Pgradient(n,x,grad,h,0.1);
//PfindHessian(n,x,Hessian,h,0.1);
//UpperTriangular(n,Hessian,upper);
newton(n,x,grad,Hessian,upper,matrixInv,copyMatrix,re,x_new,a1,b1,h,r);             //a=vector(initial guess) , g =vector(gradient vector),H=Hessian matrix,u=upper matrix to store UPPER TRIANGULAR for of Hessian, (h=0.0001=small value to find grad,hessian) 
//matrixInverse(n,Hessian,matrixInv,copyMatrix);

return 0;
}

