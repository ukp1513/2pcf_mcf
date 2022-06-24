// one property
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/stat.h>
#include "my_functions.h"
#include<dirent.h>

#define AB 12

int is_bs=0, is_jk=0;
int n_bins, n_bins_mp, n_bins_mp_shuffle, n_bins_tofit, n_copies, n_marks, n_shuffles;




double chi_square(double rp[n_bins_tofit],double CF[n_bins_tofit],double mat[n_bins_tofit][n_bins_tofit],double r0,double gamma);
double chi_square_svd(double rp[n_bins_tofit],double CF[n_bins_tofit],double errCF[n_bins_tofit],double mat[n_bins_tofit][n_bins_tofit],double r0,double gamma);
double determinant(double A[n_bins_tofit][n_bins_tofit], int n);
void inverse(double A[n_bins_tofit][n_bins_tofit], double inverse[n_bins_tofit][n_bins_tofit]);
void adjoint(double A[n_bins_tofit][n_bins_tofit],double adj[n_bins_tofit][n_bins_tofit]);
void getCofactor(double A[n_bins_tofit][n_bins_tofit], double temp[n_bins_tofit][n_bins_tofit], int p, int q, int n);
double gamma_func(double z);
double wp_model(double rp, double r0, double gamma);
void least_square_fit(double rp[n_bins], double wp[n_bins]);

int main(void)
{

	int i,j,k,o;
	
	mkdir("finals",0700);
	
	DIR* dir_bs=opendir("results/bootstraps");
	DIR* dir_jk=opendir("results/jackknifes");
	
	if(dir_bs && dir_jk)	{
		printf("\nBoth jackknife and bootstrap results exist! \n");
		exit(1);
	}
	else if(dir_bs)
	{
		printf("\nBootstrap copies exist\n");
		is_bs=1;
	}
	else if(dir_jk)	{
		printf("\nJackknife copies exist\n");
		is_jk=1;
	}
	else
	{
		printf("\nNo wp(rp) copies exist\n");
	}
	
	// READING wp DATA
	
	FILE *fEff=fopen("biproducts/effective_bins.txt","r");
	int n_bins_eff;
	fscanf(fEff,"%d",&n_bins_eff);
	fclose(fEff);
	
	n_bins = count_lines("results/wpRealAll_filtered.txt");
	n_bins_mp = count_lines("results/mpRealAll_mk1_filtered.txt");
	n_bins_mp_shuffle = count_lines("results/mpRealAllShuffle_mk1_filtered.txt");

	n_copies = count_cols("results/wpRealAll_filtered.txt")-2;	// first two cols are rp and real wp
	n_shuffles = count_cols("results/mpRealAllShuffle_mk1_filtered.txt")-2;
	n_marks = count_lines("propertylist.txt");
	
	n_bins_tofit=n_bins;
	
	printf("\nNumber of bins for wp: %d",n_bins);
	printf("\nNumber of bins for Mp: %d",n_bins_mp);
	printf("\nNumber of bins effectively used for fitting : %d",n_bins_eff);
	printf("\nNumber of copies for wp : %d",n_copies);
	printf("\nNumber of shuffles for mp : %d",n_shuffles);
	printf("\nNumber of properties : %d",n_marks);
	
	
	double **all_wp_rp = alloc_2d(n_bins,n_copies+2);
	FILE *all_wp=fopen("results/wpRealAll_filtered.txt","r");
	if(all_wp == NULL)	{
		fprintf(stderr,"Error opening results/wpRealAll.txt !\n\n");
		exit(1);
	}
	readMatrixFromFile(all_wp, all_wp_rp, n_bins, n_copies+2);
	
	
	
	double rp[n_bins],wp[n_bins],err_wp[n_bins];
	double rp_tofit[n_bins_tofit],wp_tofit[n_bins_tofit];
	
	for(i=0;i<n_bins;i++)
	{	
		rp[i]=all_wp_rp[i][0];
		wp[i]=all_wp_rp[i][1];
	}
	
	
	
	//******** fitting wp ****************
	
	printf("\n\n-----------FITTING wp(rp)-----------\n");
	
	double sig_gamma_min, sig_gamma_max, sig_r0_min, sig_r0_max,sig_r0,sig_gamma,r0Best,gammaBest;
	
	printf("\nLeast square fits: \n");
	least_square_fit(rp,wp);	// to see the appoximate fit parameters (not actual fits)

	if(n_copies>0)	{
		
		
			
		//COMPUTING COVARIANCE MATRIX OF ALL BINS

		FILE *cov;
		double covar[n_bins][n_bins],sum,mean[n_bins],ravoc;

		for(i=0;i<n_bins;i++)
		{
			sum=0;
			for(j=2;j<n_copies+2;j++)
				sum+=all_wp_rp[i][j];
			mean[i]=sum/n_copies;
		}

		for(i=0;i<n_bins;i++)
		{
			for(j=0;j<n_bins;j++)
			{
				ravoc=0;
				for(k=2;k<n_copies+2;k++)
				{
					ravoc+=(all_wp_rp[i][k]-mean[i])*(all_wp_rp[j][k]-mean[j]);
				}
				if(is_bs == 1)	covar[i][j]=ravoc/(n_copies-1);	
				if(is_jk == 1)	covar[i][j]=ravoc*(n_copies-1)/n_copies;		
			}
		}
		
		FILE *fCov = fopen("biproducts/cov_matrix.txt","w");
		for(i=0;i<n_bins;i++)
		{
			for(j=0;j<n_bins;j++)	fprintf(fCov,"%lf\t",covar[i][j]);
			fprintf(fCov,"\n");
		}
		fclose(fCov);

		// ASSIGNING THE ERROR BARS
		
		for(i=0;i<n_bins;i++)
		{
			err_wp[i] = sqrt(covar[i][i]);	// square root of diagonal elements of covariance matrix
		}
		
		// COMPUTING COVARIANCE MATRIX FOR BINS TO FIT

		
		
		
		
		double **all_wp_rp_tofit = alloc_2d(n_bins_tofit,n_copies+2);
		FILE *all_wp_tofit=fopen("results/wpRealAll_filtered_tofit.txt","r");
		if(all_wp_tofit == NULL)	{
			fprintf(stderr,"Error opening results/wpRealAll_tofit.txt !\n\n");
			exit(1);
		}
		readMatrixFromFile(all_wp_tofit, all_wp_rp_tofit, n_bins_tofit, n_copies+2);
		
		for(i=0;i<n_bins_tofit;i++)
		{	
			rp_tofit[i]=all_wp_rp_tofit[i][0];
			wp_tofit[i]=all_wp_rp_tofit[i][1];
		}

		FILE *cov_tofit;
		double covar_tofit[n_bins_tofit][n_bins_tofit],sum_tofit,mean_tofit[n_bins_tofit],ravoc_tofit;

		for(i=0;i<n_bins_tofit;i++)
		{
			sum_tofit=0;
			for(j=2;j<n_copies+2;j++)
				sum_tofit+=all_wp_rp_tofit[i][j];
			mean_tofit[i]=sum_tofit/n_copies;
		}

		for(i=0;i<n_bins_tofit;i++)
		{
			for(j=0;j<n_bins_tofit;j++)
			{
				ravoc_tofit=0;
				for(k=2;k<n_copies+2;k++)
				{
					ravoc_tofit+=(all_wp_rp_tofit[i][k]-mean_tofit[i])*(all_wp_rp_tofit[j][k]-mean_tofit[j]);
				}
				if(is_bs == 1)	covar_tofit[i][j]=ravoc_tofit/(n_copies-1);	
				if(is_jk == 1)	covar_tofit[i][j]=ravoc_tofit*(n_copies-1)/n_copies;		
			}
		}
		
		FILE *fCov_tofit = fopen("biproducts/cov_matrix_tofit.txt","w");
		for(i=0;i<n_bins_tofit;i++)
		{
			for(j=0;j<n_bins_tofit;j++)	fprintf(fCov_tofit,"%lf\t",covar_tofit[i][j]);
			fprintf(fCov_tofit,"\n");
		}
		fclose(fCov_tofit);
		
		double Cinv_SVD[n_bins_tofit][n_bins_tofit];
		FILE *fCinv_SVD=fopen("biproducts/Cinv_SVD.txt","r");
		for(int i=0;i<n_bins_tofit;i++)	{
			for(int j=0;j<n_bins_tofit;j++)	{
				fscanf(fCinv_SVD,"%lf",&Cinv_SVD[i][j]);
			}
		}
		fclose(fCinv_SVD);
		
		
		// INVERTING COVARIANCE MATRIX
		FILE *invcov;
		double invcovar[n_bins_tofit][n_bins_tofit];
		inverse(covar_tofit,invcovar);
	
		
		FILE *fICov = fopen("biproducts/inv_cov_matrix.txt","w");
		for(i=0;i<n_bins_tofit;i++)
		{
			for(j=0;j<n_bins_tofit;j++)	fprintf(fICov,"%lf\t",invcovar[i][j]);
			fprintf(fICov,"\n");
		}
		fclose(fICov);
	
		// CHI-SQUARE FIT
		
		double gammaMin=1.01,gammaMax=2.5,r0Min=1.0,r0Max=10.0,dGamma=0.01,dR0=0.01;
		int gammaParamCount=(gammaMax-gammaMin)/dGamma,r0ParamCount=(r0Max-r0Min)/dR0;
		double gammaSpace[gammaParamCount],r0Space[r0ParamCount];
		double chi2matrix[gammaParamCount][r0ParamCount];
		
		for(i=0;i<gammaParamCount;i++)
			gammaSpace[i]=gammaMin+(i*dGamma);
		for(i=0;i<r0ParamCount;i++)
			r0Space[i]=r0Min+(i*dR0);
			
			
		double chi2min=1e10;
		int bestR0Index,bestGammaIndex;
		for(i=0;i<gammaParamCount;i++)	{	
			for(j=0;j<r0ParamCount;j++)	{
				double chi2 = chi_square_svd(rp_tofit,wp_tofit,err_wp,Cinv_SVD,r0Space[j],gammaSpace[i]);
				chi2matrix[i][j] = chi2;
				if(chi2 < chi2min)	{
					chi2min=chi2;
					bestR0Index=j;
					bestGammaIndex=i;
					r0Best = r0Space[j];
					gammaBest = gammaSpace[i];
				}
			}
		}
		
		FILE *fr0_space = fopen("biproducts/r0_space.txt","w");
		for(i=0;i<r0ParamCount;i++)			fprintf(fr0_space,"%lf\n",r0Space[i]);
		fclose(fr0_space);
		
		FILE *fgamma_space = fopen("biproducts/gamma_space.txt","w");
		for(i=0;i<gammaParamCount;i++)			fprintf(fgamma_space,"%lf\n",gammaSpace[i]);
		fclose(fgamma_space);		
				
		FILE *fChisq = fopen("biproducts/chi_sq.txt","w");
		for(i=0;i<gammaParamCount;i++)
		{
			for(j=0;j<r0ParamCount;j++)	fprintf(fChisq,"%lf\t",chi2matrix[i][j]);
			fprintf(fChisq,"\n");
		}
		fclose(fChisq);
		
		// ESTIMATING ERRORS IN PARAMETERS USING DELTA_chi2matrix2 = 2.3 (1 SIGMA)
		// Refer Section 15.6.5 of Numerical Recipes (Press et al., 1992)
		
		
		double deltaChi2matrix[gammaParamCount][r0ParamCount];
		for(i=0;i<gammaParamCount;i++)		
			for(j=0;j<r0ParamCount;j++)	
				deltaChi2matrix[i][j] = chi2matrix[i][j]-chi2min;
				
		FILE *fdChisq = fopen("biproducts/delta_chi_sq.txt","w");
		for(i=0;i<gammaParamCount;i++)
		{
			for(j=0;j<r0ParamCount;j++)	fprintf(fdChisq,"%lf\t",deltaChi2matrix[i][j]);
			fprintf(fdChisq,"\n");
		}
		fclose(fdChisq);
		
		int imin = gammaParamCount, imax=0;
		int jmin = r0ParamCount, jmax=0;
		for(i=0;i<gammaParamCount;i++)		{
			for(j=0;j<r0ParamCount;j++)	{
				if(deltaChi2matrix[i][j] > 0.0 && deltaChi2matrix[i][j] <= 2.3)	{
					if(i<imin)	imin=i;
					if(i>imax)	imax=i;
					if(j<jmin)	jmin=j;
					if(j>jmax)	jmax=j;
				}
			}
		}
		
		sig_gamma_min=gammaBest-gammaSpace[imin];
		sig_gamma_max=gammaSpace[imax]-gammaBest;
		
		sig_r0_min=r0Best-r0Space[jmin];
		sig_r0_max=r0Space[jmax]-r0Best;
		
		//sig_r0 = (sig_r0_max + sig_r0_min)/2.0;
		//sig_gamma = (sig_gamma_max + sig_gamma_min)/2.0;
		
		if(sig_r0_min > sig_r0_max)	sig_r0 = sig_r0_min;
		else sig_r0 = sig_r0_max;

		if(sig_gamma_min > sig_gamma_max)	sig_gamma = sig_gamma_min;
		else sig_gamma = sig_gamma_max;
				
		printf("\nMinimum value of chi = %0.2f\n\nr0 = %0.2lf +/- %0.2lf\ngamma = %0.2lf +/- %0.2lf\n",chi2min,r0Best,sig_r0,gammaBest,sig_gamma);
				
	}

	else	{
		for(i=0;i<n_bins;i++)	err_wp[i]=0;
	}
	
	FILE *wp_rp_errors=fopen("finals/final_wp.txt","w");
	for(i=0;i<n_bins;i++)
		fprintf(wp_rp_errors,"%lf\t%lf\t%lf\n",rp[i],wp[i],err_wp[i]);
	fclose(wp_rp_errors);
	
	FILE *wp_params=fopen("finals/wp_params.txt","w");
	fprintf(wp_params,"%lf\n%lf\n%lf\n%lf",r0Best,sig_r0,gammaBest,sig_gamma);
	fclose(wp_params);
	
	
	// COMPUTING MCF ERRORS
	
	printf("\n\n-----------Mp(rp)-----------\n");
	
	FILE *fMarks=fopen("propertylist.txt","r");
	

	for(int prop_nr=1;prop_nr<=n_marks;prop_nr++)	{
	
		char prop[20];
		fscanf(fMarks,"%s",prop);
	
		//************ reading the shuffles ***********

		double **all_Mp_rp_shuffle = alloc_2d(n_bins_mp_shuffle, n_shuffles+2);
		char name_mp_shuffle[50];
		sprintf(name_mp_shuffle,"results/mpRealAllShuffle_mk%d_filtered.txt",prop_nr);
		FILE *all_Mp_shuffle=fopen(name_mp_shuffle,"r");
		if(all_Mp_shuffle == NULL)	{
			continue;
		}	

		readMatrixFromFile(all_Mp_shuffle, all_Mp_rp_shuffle, n_bins_mp_shuffle, n_shuffles+2);
		
		double rp_mp_shuffle[n_bins_mp_shuffle],Mp_shuffle[n_bins_mp_shuffle],err_Mp_shuffle[n_bins_mp_shuffle];
		
		for(i=0;i<n_bins_mp_shuffle;i++)
		{	
			rp_mp_shuffle[i]=all_Mp_rp_shuffle[i][0];
			Mp_shuffle[i]=all_Mp_rp_shuffle[i][1];
		}
	
		// ESTIMATING ERROR BARS OF MCF (variance among all mcfs of shuffles)
		if(n_shuffles > 0)	{
					
			for(i=0;i<n_bins_mp_shuffle;i++)
			{
				double bs_row_Mp_shuffle[n_shuffles],var_sum_shuffle=0,var_avg_shuffle;
				for(j=2;j<(n_shuffles+2);j++)
					bs_row_Mp_shuffle[j-2]=all_Mp_rp_shuffle[i][j];
				err_Mp_shuffle[i] = sqrt(array_var(bs_row_Mp_shuffle,n_shuffles));
			}
		}
		else	{
			for(i=0;i<n_bins_mp;i++)	err_Mp_shuffle[i]=0.0;
		}
		
		//************ reading the JK/BS copies ***********

		double **all_Mp_rp = alloc_2d(n_bins_mp, n_copies+2);
		char name_mp[50];
		sprintf(name_mp,"results/mpRealAll_mk%d_filtered.txt",prop_nr);
		FILE *all_Mp=fopen(name_mp,"r");
		if(all_Mp == NULL)	{
			continue;
		}	

		readMatrixFromFile(all_Mp, all_Mp_rp, n_bins_mp, n_copies+2);
		
		double rp_mp[n_bins_mp],Mp[n_bins_mp],err_Mp[n_bins_mp];
		
		for(i=0;i<n_bins_mp;i++)
		{	
			rp_mp[i]=all_Mp_rp[i][0];
			Mp[i]=all_Mp_rp[i][1];
		}
	
		// ESTIMATING ERROR BARS OF MCF (variance among all mcfs of JK/BS copies)
		if(n_copies > 0)	{
					
			for(i=0;i<n_bins_mp;i++)
			{
				double bs_row_Mp[n_copies],var_sum=0,var_avg;
				for(j=2;j<(n_copies+2);j++)
					bs_row_Mp[j-2]=all_Mp_rp[i][j];
				err_Mp[i] = sqrt(array_var(bs_row_Mp,n_copies));
			}
		}
		else	{
			for(i=0;i<n_bins_mp;i++)	err_Mp[i]=0.0;
		}
		
		
		char result_Mp_name[50];
		sprintf(result_Mp_name,"finals/final_Mp_%s.txt",prop);
		FILE *Mp_rp_errors=fopen(result_Mp_name,"w");
		for(i=0;i<n_bins_mp;i++)
			fprintf(Mp_rp_errors,"%lf\t%lf\t%lf\t%lf\n",rp_mp[i],Mp[i],err_Mp_shuffle[i],err_Mp[i]);
		fclose(Mp_rp_errors);
		
	}
	
	fclose(fMarks);
	printf("\n\n");
	return(0);


	

	
}



double compute_Mp(double rp, double wp, double Wp)
{
	double Mp = (1+(Wp/rp))/(1+(wp/rp));
	return Mp;
}

double error_Mp(double rp, double wp, double Wp, double Mp, double err_wp, double err_Wp)
{
	double err_Mp;
	err_Mp = Mp * sqrt(pow(err_Wp/(Wp+rp),2)+pow(err_wp/(wp+rp),2));
	return err_Mp;
}

int count_nbins(char sample[100])
{
	char buf2[255];
	char name_0[50];
	int n_bins_inst=0;
	sprintf(name_0,"%s/Mp_shuffle/Mp_rp0.txt",sample);
	FILE *act_file=fopen(name_0,"r");
	if(act_file == NULL)	{
		fprintf(stderr,"\nError opening actual result file in %s\n",sample);
		exit(1);
	}
	while(!feof(act_file))	{
		fgets(buf2,255,act_file);
		n_bins_inst++;
	}
	n_bins_inst--;
	fprintf(stdout,"\nNumber of bins in %s: %d\n",sample,n_bins_inst);
	fclose(act_file);
	return n_bins_inst;
}


	



 
double chi_square(double rp[n_bins_tofit],double CF[n_bins_tofit],double mat[n_bins_tofit][n_bins_tofit],double r0,double gamma)
{
	double chi2;
	int i,j;

	chi2=0;

	for(i=0;i<n_bins_tofit;i++)
	{
		for(j=0;j<n_bins_tofit;j++)
			chi2+=(wp_model(rp[i],r0,gamma)-CF[i])*mat[i][j]*(wp_model(rp[j],r0,gamma)-CF[j]);
	}

	return chi2;
}

double chi_square_svd(double rp[n_bins_tofit],double CF[n_bins_tofit],double errCF[n_bins_tofit],double mat[n_bins_tofit][n_bins_tofit],double r0,double gamma)
{
	double chi2;
	int i,j;

	chi2=0;

	for(i=0;i<n_bins_tofit;i++)
	{
		for(j=0;j<n_bins_tofit;j++)
			chi2+=((wp_model(rp[i],r0,gamma)-CF[i])/errCF[i])*mat[i][j]*((wp_model(rp[j],r0,gamma)-CF[j])/errCF[j]);
	}

	return chi2;
}

double wp_model(double rp, double r0, double gamma)
{
	return rp*pow(r0/rp,gamma)*gamma_func(0.5)*gamma_func(0.5*(gamma-1))/gamma_func(0.5*gamma);
}


double determinant(double A[n_bins_tofit][n_bins_tofit], int n)
{
	double D = 0; 
	if (n == 1)
		return A[0][0];
	double temp[n_bins_tofit][n_bins_tofit]; 
	int sign = 1;  

	for (int f = 0; f < n; f++)
	{
		getCofactor(A, temp, 0, f, n);
		D += sign * A[0][f] * determinant(temp, n - 1);
		sign = -sign;
	}
	return D;
}

void inverse(double A[n_bins_tofit][n_bins_tofit], double inverse[n_bins_tofit][n_bins_tofit])
{
	double det = determinant(A, n_bins_tofit);
	if (det == 0)
		fprintf(stdout,"\nSingular matrix, can't find its inverse");

	double adj[n_bins_tofit][n_bins_tofit];
	adjoint(A, adj);
	for (int i=0; i<n_bins_tofit; i++)
		for (int j=0; j<n_bins_tofit; j++)
			inverse[i][j] = adj[i][j]/det;
}
void adjoint(double A[n_bins_tofit][n_bins_tofit],double adj[n_bins_tofit][n_bins_tofit])
{
	if (n_bins == 1)
		adj[0][0] = 1;

	int sign = 1;
	double temp[n_bins_tofit][n_bins_tofit];
 
	for (int i=0; i<n_bins_tofit; i++)
	{
		for (int j=0; j<n_bins_tofit; j++)
		{
			getCofactor(A, temp, i, j, n_bins_tofit);
			sign = ((i+j)%2==0)? 1: -1;
			adj[j][i] = (sign)*(determinant(temp, n_bins_tofit-1));
		}
	}
}
void getCofactor(double A[n_bins_tofit][n_bins_tofit], double temp[n_bins_tofit][n_bins_tofit], int p, int q, int n)
{
	int i = 0, j = 0;
	for (int row = 0; row < n; row++)
	{
		for (int col = 0; col < n; col++)
		{
			if (row != p && col != q)
			{
				temp[i][j++] = A[row][col];
				if (j == n - 1)
				{
					j = 0;
					i++;
				}
			}
		}
	}
}

double gamma_func(double z)
{
  const int a = AB;
  static double c_space[AB];
  static double *ca = NULL;
  int k;
  double accm;
 
  if ( ca == NULL ) {
    double k1_factrl = 1.0; 
    ca = c_space;
    ca[0] = sqrt(2.0*M_PI);
    for(k=1; k < a; k++) {
      ca[k] = exp(a-k) * pow(a-k, k-0.5) / k1_factrl;
	  k1_factrl *= -k;
    }
  }
  accm = ca[0];
  for(k=1; k < a; k++) {
    accm += ca[k] / ( z + k );
  }
  accm *= exp(-(z+a)) * pow(z+a, z+0.5); 
  return accm/z;
}

void least_square_fit(double rp[n_bins], double wp[n_bins])
{
	int i,n=n_bins;
	double x[n_bins],y[n_bins],sumx=0,sumy=0,sumxy=0,sumx2=0,a,b,r0_lr,gamma_lr;
	for(i=0;i<n;i++)	{
		x[i]=log10(rp[i]);
		y[i]=log10(wp[i]);
	}
	for(i=0;i<n;i++)	{
        sumx=sumx +x[i];
        sumx2=sumx2 +x[i]*x[i];
        sumy=sumy +y[i];
        sumxy=sumxy +x[i]*y[i];
    }
    a=((sumx2*sumy -sumx*sumxy)*1.0/(n*sumx2-sumx*sumx)*1.0);
    b=((n*sumxy-sumx*sumy)*1.0/(n*sumx2-sumx*sumx)*1.0);
	
	gamma_lr = 1-b;
	r0_lr = pow(pow(10,a)*gamma_func(gamma_lr*0.5)/(gamma_func(0.5)*gamma_func(0.5*(gamma_lr-1))),1/gamma_lr);

	fprintf(stdout,"\nr0=%0.2f \ngamma=%0.2f\n",r0_lr,gamma_lr);
}
