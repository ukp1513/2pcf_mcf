#include<stdio.h>
#include<math.h>
#include<stdlib.h>

#define D_R (M_PI/180.0)
#define c 300000
#define Omega_M 0.3
#define Omega_L 0.7
#define H0 70.0
#define Omega_k (1-Omega_M-Omega_L)
#define D_H (c/H0)	// in Mpc
#define abs_mag_sun 4.83

double h=H0/100.0;

double one_by_E_z(double z)		// reciprocal of E(z) in the equation for comoving distance
{
	return 1/sqrt(Omega_M*pow(1+z,3)+Omega_L);
}

double D_C(double z)					// comoving distance eq 15 of Hogg 2000
{
	int i,n=10;
	double h1=z/n,f[n],I,s1=0;
	double f0=one_by_E_z(0);
	double fn=one_by_E_z(z);
	for(i=1;i<n;i++)
		s1=s1+one_by_E_z(i*h1);
	I=(h1/2)*(f0+fn+2*s1);		
	return D_H*I*h;					// in h^-1 Mpc
}

double D_M(double z)					// transverse comoving distance eq 16 of Hogg 2000
{
	if(Omega_k == 0)
		return D_C(z);
	else if(Omega_k < 0)
		return D_H*(1/sqrt(fabs(Omega_k)))*sin(sqrt(fabs(Omega_k))*D_C(z)/D_H);
	else if(Omega_k > 0)
		return D_H*(1/sqrt(Omega_k))*sinh(sqrt(Omega_k)*D_C(z)/D_H);
}

double D_L(double z)					// luminosity distance eq 21 of Hogg 2000
{
	return (1+z)*D_M(z);
}

double Dist_Modulus(double z)		// distance modulus eq 25 of Hogg 2000
{
	return 5*log10(D_L(z)*1E6/10);
}

double K_corr(double z)				// k correction yet to be implemented
{
	return 0;
}

double abs_mag(double app_mag,double z)	// absolute magnitude eq 1 of Zehavi et al. 2002
{
	return app_mag-Dist_Modulus(z)-K_corr(z)+(5*log10(h));	// mag corrected for hubble constant
}

double luminosity(double abs_mag)			// luminosity eq 3.7 of Skibba's PhD thesis
{
	return pow(10,-0.4*(abs_mag - abs_mag_sun)); 	// in units of L_sun 
}

double array_max(double array[], int size)	// maximum value in the array
{
	double M=array[0];
	int k;
	for(k=1;k<size;k++)
		if(array[k]>M)
			M=array[k];
	return M;
}

double array_min(double array[], int size)	// minimum value in the array
{
	double m=array[0];
	int k;
	for(k=1;k<size;k++)
		if(array[k]<m)
			m=array[k];
	return m;
}

double array_avg(double array[], int size)	// average value of elements in the array
{
	double sum=0;
	int k;
	for(k=0;k<size;k++)
		sum+=array[k];
	return (sum/size); 
}

void array_sort(double *array , int n)		// sort an array
{ 
	int i=0 , j=0;
	double temp=0;
	for(i=0;i<n;i++)
	{
		for(j=0;j<n-1;j++)
		{
			if(array[j]>array[j+1])
			{
				temp = array[j];
				array[j] = array[j+1];
				array[j+1] = temp;
			}
		}
	}
}

double array_median(double array[] , int n)		// median of an array
{
	array_sort(array,n);
	double median=0;
	    
	if(n%2 == 0)			// if number of elements are even
		median = (array[(n-1)/2] + array[n/2])/2.0;
    
	else						// if number of elements are odd
		median = array[n/2];
    
	return median;
}

double array_var(double array[], int size)	// variance of array
{
	int j;
	double average=array_avg(array,size),sum=0,variance;
	for(j=0;j<size;j++)	
		sum += pow(array[j]-average,2);
	variance = sum / (float) (size-1);
	return variance;
}

double array_std_dev(double array[], int n)		// standard deviation 
{
	double mean=array_avg(array,n);
	int i;
	double sd=0;

	for(i=0;i<n;i++)
		sd += pow((array[i] - mean), 2);

	return sqrt(sd/n);
}

double rand_normal (double mu, double sigma)		// generate random number from normal distribution
{																// Marsaglia polar method
  double U1, U2, W, mult;
  static double X1, X2;
  static int call = 0;
 
  if (call == 1)
    {
      call = !call;
      return (mu + sigma * (double) X2);
    }
 
  do
    {
      U1 = -1 + ((double) rand () / RAND_MAX) * 2;
      U2 = -1 + ((double) rand () / RAND_MAX) * 2;
      W = pow (U1, 2) + pow (U2, 2);
    }
  while (W >= 1 || W == 0);
 
  mult = sqrt ((-2 * log (W)) / W);
  X1 = U1 * mult;
  X2 = U2 * mult;
 
  call = !call;
 
  return (mu + sigma * (double) X1);
}

void rankify(double* A , int n,double* R)	 	// Function to rank an array 
{
	int i,j;
	for (i = 0; i < n; i++)	R[i]=0;
	for (i = 0; i < n; i++) { 
		int r = 1, s = 1; 
		for (j = 0; j < n; j++)	{ 
			if (j != i && A[j] < A[i]) 
				r += 1; 
			if (j != i && A[j] == A[i]) 
				s += 1;      
		} 
		R[i] = r + (double)(s - 1) / (double) 2;   
	}  
} 

void equateArray(double *A, double *B, int n)	// equate array B to A (equivalent to A=B)
{
	for(int i=0;i<n;i++)	A[i]=B[i];
}	

void equateMatrix(double **A, double **B, int n1, int n2)	// equate array B to A (equivalent to A=B)
{
	for(int i=0;i<n1;i++)	{
		for(int j=0;j<n2;j++)	A[i][j]=B[i][j];
	}
}

void displayArray(double *A, int n)	// to display an array on screen
{
	for(int i=0; i<n;i++)	printf("%lf\n",A[i]);
}

void displayMatrix(double **A, int nrow, int ncol)	// to display a 2d array on screen
{
	for(int i=0; i<nrow; i++)	{
		for(int j=0; j<ncol; j++)
			printf("%lf\t",A[i][j]);
		printf("\n");
	}
}

void writeMatrixToFile(int nrow, int ncol, double A[nrow][ncol], FILE *f) // to write matrix A to file f
{
	for(int i=0; i<nrow; i++)	{
		for(int j=0; j<ncol; j++)
			fprintf(f,"%lf\t",A[i][j]);
		fprintf(f,"\n");
	}
}

void readMatrixFromFile(FILE *f, double **A, int nrow, int ncol)	// to read matrix A from file f
{
	for(int i=0; i<nrow; i++)	{
		for(int j=0; j<ncol; j++)
			fscanf(f,"%lf",&A[i][j]);
	}
}


double* alloc_1d(size_t n)	// dynamically allocates a 1d array[n] of type double
{
	double *array = malloc(n*sizeof(double));
	return array;
}

double** alloc_2d(size_t n1, size_t n2)	// dynamically allocates a 2d array[n1][n2] of type double
{
	int i;
	double **array=malloc(n1 * sizeof(double*));
	for(i=0;i<n1;i++)	
		array[i]=malloc(n2 * sizeof(double));
	return array;
}

double*** alloc_3d(size_t n1, size_t n2, size_t n3)	// dynamically allocates a 3d array[n1][n2][n3] of type double
{
	int i,j;
	double ***array=malloc(n1 * sizeof(double**));
	for(i=0;i<n1;i++)	{	
		array[i]=malloc(n2 * sizeof(double*));
		for(j=0;j<n2;j++)	
			array[i][j]=malloc(n3 * sizeof(double));
	}
	return array;
}

void free_1d(double *array)	// to free the 1-d array
{
	free(array);
}

void free_2d(double **array, size_t n1)	// to free the 2-d array of n1 number of rows
{
	int i;
	for(i=0;i<n1;i++)	free(array[i]);
	free(array);
}

void free_3d(double ***array, size_t n1, size_t n2)	// to free the 3-d array of n1 number of rows and n2 number of columns
{
	int i,j;
	for(i=0;i<n1;i++)	{
		for(j=0;j<n2;j++)	free(array[i][j]);
		free(array[i]);
	}
	free(array);
}

int count_lines(char *file)	// to count the number of lines in a file 
{
    FILE *fic;
    int n;
    char com[BUFSIZ];
    
    fic=fopen(file,"r");
    if (fic) {
	n=0;
	while(!feof(fic)) {
	    fscanf(fic,"%[^\n]\n",com);
	    n++;
	}
	fclose(fic);
	return n;
    } else {
        fprintf(stderr,"Can't read %s file !\n",file);
	return 0;
    }
    //TODO note: if there are only \n in the file it hangs
}


int count_cols(char *file)	// to count the number of columns in a file
{
    FILE *fic;
    int n;
    char com[BUFSIZ];
    
    fic=fopen(file,"r");
    if (fic) {
	n=0;
	while(!feof(fic)) {
	    fscanf(fic,"%s",com);
	    n++;
	}
	fclose(fic);
	return n/count_lines(file);
    } else {
        fprintf(stderr,"Can't read %s file !\n",file);
	return 0;
    }
}

int fact(int num)	// to find factorial of a number
{
    int k = 1, i;
    // factorial of 0 is 1
    if (num == 0)
    {
        return(k);
    }
    else
    {
        for (i = 1; i <= num; i++)
    {
            k = k * i;
	}
    }
    return(k);
}


