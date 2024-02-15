// This code computes projected 2pt and marked correlation functions 
// of real data and bootstrap/jackknife samples.
// Details in README.
// LAST UPDATE : 24 June 2022.
// -------------- UNNIKRISHNAN SURESHKUMAR -------------------//
// ----------------- ukp1513@gmail.com -----------------------//

#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include <sys/stat.h>
#include "my_functions.h"
#include<omp.h>
#include<string.h>

/*------------ SET PARAMETERS --------------*/

#define non_prop 4				// Number of columns in real file that are not to be used as marks 
								// These columns should be the first columns				
#define boot_meth 2				// 1: individual; 2: blockwise
#define n_boots 0				// number of bootstrap samples
#define n_jacks 24				// Number of jackknife samples 
#define bin_rp 20				// number of bins
#define bin_pi 40				// here goes the pmax

double rp_init=0.01;			// centre of first (smallest) bin in rp
double pi_init=0.5;				// centre of first (smallest) bin in pi
double dlrp=0.24;				// binwidth in log scale
double dpi=1.0;					// binwidth in pi (linear scale)

/*-------------------------------------------*/

#if n_boots > 0
#define n_copies n_boots
#endif

#if n_jacks > 0
#define n_copies n_jacks
#endif

#if n_jacks == 0 && n_boots == 0
#define n_copies 0
#endif

int n_marks;
double rp_bins[bin_rp],pi_bins[bin_pi];
double **RR;
double rp_wp_all[bin_rp][n_copies+2];
double rp_lowcut,rp_highcut,pi_highcut;
FILE *f_summary;

void shuffle(double *array, size_t n)	// to scramble the array having marks
{								
	
    if (n > 1) 
    {
        size_t i;
        for (i = 0; i < n - 1; i++) 
        {
			srand((unsigned) time(NULL));
          size_t j = i + rand() / (RAND_MAX / (n - i) + 1);
          double t = array[j];
          array[j] = array[i];
          array[i] = t;
        }
    }
}

void sep(double ra1,double dec1,double z1,double ra2,double dec2,double z2,double *RP,double *PI, double *TH)
{
	double th,a1,a2,d1,d2,r,v1,v2,dab,s,pi,rp,S;
	double x1,x2,y1,y2,zz1,zz2;

	// converting degree to radian for trigonometric functions in C
	a1=ra1*D_R;					
	d1=dec1*D_R;				
	a2=ra2*D_R;					
	d2=dec2*D_R;

	// converting RA, DEC, Z into cartesian coords. x,y,z
	x1=D_C(z1)*cos(d1)*cos(a1);		
	x2=D_C(z2)*cos(d2)*cos(a2);

	y1=D_C(z1)*cos(d1)*sin(a1);
	y2=D_C(z2)*cos(d2)*sin(a2);

	zz1=D_C(z1)*sin(d1);
	zz2=D_C(z2)*sin(d2);

	s=(x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(zz1-zz2)*(zz1-zz2);
	pi=x1*x1-x2*x2+y1*y1-y2*y2+zz1*zz1-zz2*zz2;
	pi=pi*pi/((x1+x2)*(x1+x2)+(y1+y2)*(y1+y2)+(zz1+zz2)*(zz1+zz2));

	rp=sqrt(s-pi);
	pi=sqrt(pi);
	//s=sqrt(s);
	th=(acos((sin(d1)*sin(d2))+(cos(d1)*cos(d2)*cos(a1-a2))));

	*TH=th;
	*RP=rp;
	*PI=pi;
}

void DD_count(double *ra1,double *dec1,double *z1, unsigned long int n1,unsigned long int count_array_DD[bin_rp][bin_pi])
{
	double rp,pi,th;
	for(int j=0;j<bin_rp;j++)	{
		for(int k=0;k<bin_pi;k++)	{
			count_array_DD[j][k]=0;
		}
	}
	
	FILE *f_sepDD=fopen("cache/real_real_sep.txt","w");
	for(int l=0;l<n1;l++)	{
		for(int m=0;m<n1;m++)	{
			if(l<m)	{
				sep(ra1[l],dec1[l],z1[l],ra1[m],dec1[m],z1[m],&rp,&pi,&th);
				if(rp < rp_lowcut || rp > rp_highcut || pi > pi_highcut)	continue;
				fprintf(f_sepDD,"%lf\t%lf\t%lf\n",th,rp,pi);
				rp=log10(rp);
				
				int bin_of_rp = (rp-log10(rp_lowcut))/dlrp;
				int bin_of_pi = (pi-0.0)/dpi;
				
				count_array_DD[bin_of_rp][bin_of_pi]+=1;
				
				
			}
		}
	}
	fclose(f_sepDD);
	

}

void DR_count(double *ra1,double *dec1,double *z1, unsigned long int n1,double *ra2,double *dec2,double *z2, unsigned long int n2,unsigned long int count_array_DR[bin_rp][bin_pi])
{
	double rp,pi, th;
	for(int j=0;j<bin_rp;j++)	{
		for(int k=0;k<bin_pi;k++)	{
			count_array_DR[j][k]=0;
		}
	}

	FILE *f_sepDR=fopen("cache/real_rand_sep.txt","w");
	for(int l=0;l<n1;l++)	{
		for(int m=0;m<n2;m++)	{
			sep(ra1[l],dec1[l],z1[l],ra2[m],dec2[m],z2[m],&rp,&pi,&th);
			if(rp < rp_lowcut || rp > rp_highcut || pi > pi_highcut)	continue;
			fprintf(f_sepDR,"%lf\t%lf\t%lf\n",th,rp,pi);
			rp=log10(rp);
			int bin_of_rp = (rp-log10(rp_lowcut))/dlrp;
			int bin_of_pi = (pi-0.0)/dpi;
			
			count_array_DR[bin_of_rp][bin_of_pi]+=1;
		}	
	}
	fclose(f_sepDR);	
	

}

void RR_count(double *ra1,double *dec1,double *z1 ,unsigned long int n1,unsigned long int count_array_RR[bin_rp][bin_pi])
{
	
	double rp,pi,th;
	for(int j=0;j<bin_rp;j++)	
		for(int k=0;k<bin_pi;k++)	
			count_array_RR[j][k]=0;

	FILE *f_sepRR=fopen("cache/rand_rand_sep.txt","w");
	for(int l=0;l<n1;l++)	{
		for(int m=0;m<n1;m++)	{
			if(l<m)	{
				sep(ra1[l],dec1[l],z1[l],ra1[m],dec1[m],z1[m],&rp,&pi, &th);
				if(rp < rp_lowcut || rp > rp_highcut || pi > pi_highcut)	continue;
				fprintf(f_sepRR,"%lf\t%lf\t%lf\n",th,rp,pi);
				rp=log10(rp);
				int bin_of_rp = (rp-log10(rp_lowcut))/dlrp;
				int bin_of_pi = (pi-0.0)/dpi;
				
				count_array_RR[bin_of_rp][bin_of_pi]+=1;	
			}
		}	
	}
	fclose(f_sepRR);
	

}

void compute2pCF(int sample_nr, unsigned long int n_d, double *ra_d, double *dec_d, double *z_d, unsigned long int n_r, double *ra_r, double *dec_r, double *z_r)
{

	unsigned long int DD[bin_rp][bin_pi], DR[bin_rp][bin_pi], RR[bin_rp][bin_pi];
	long double xi_2d[bin_rp][bin_pi], wp[bin_rp];


	unsigned long long int ndd=n_d*(n_d-1)*0.5;
	unsigned long long int nrr=n_r*(n_r-1)*0.5;
	unsigned long long int ndr=n_d*n_r;
	
	// ************* COMPUTING RR ****************************
	
	char RRName[100];
	sprintf(RRName,"cache/RR_rp_pi_%d.txt",sample_nr);
	FILE *f_RR=fopen(RRName,"r");

	if(!f_RR)	{
		fprintf(stdout,"\nSample %d: RR not computed! Computing RR and writing to RR_rp_pi.txt\n", sample_nr);
		fprintf(f_summary,"\nSample %d: RR not computed! Computing RR and writing to RR_rp_pi.txt\n", sample_nr);
		RR_count(ra_r,dec_r,z_r,n_r,RR);
		f_RR=fopen(RRName,"w");
		for(int i=0;i<bin_rp;i++)	{
			for(int j=0;j<bin_pi;j++)
				fprintf(f_RR,"%lu\t",RR[i][j]);
			fprintf(f_RR,"\n");
		}
		fclose(f_RR);
	}
	else	{
		fprintf(stdout,"\nReading RR from cache/RR_rp_pi.txt \n");
		fprintf(f_summary,"\nReading RR from cache/RR_rp_pi.txt \n");
		for(int i=0;i<bin_rp;i++)	
			for(int j=0;j<bin_pi;j++)
				fscanf(f_RR,"%lu",&RR[i][j]);
		fclose(f_RR);
	}	
	
	// ************* COMPUTING DD ****************************

	char DDName[100];
	sprintf(DDName,"cache/DD_rp_pi_%d.txt",sample_nr);
	FILE *f_DD=fopen(DDName,"r");
	
	if(!f_DD)	{
		fprintf(stdout,"\nSample %d: DD not computed! Computing DD and writing to DD_rp_pi.txt\n", sample_nr);
		fprintf(f_summary,"\nSample %d: DD not computed! Computing DD and writing to DD_rp_pi.txt\n", sample_nr);
		DD_count(ra_d,dec_d,z_d,n_d,DD);
		f_DD=fopen(DDName,"w");
		for(int i=0;i<bin_rp;i++)	{
			for(int j=0;j<bin_pi;j++)
				fprintf(f_DD,"%lu\t",DD[i][j]);
			fprintf(f_DD,"\n");
		}
		fclose(f_DD);
	}
	else	{
		fprintf(stdout,"\nReading DD from cache/DD_rp_pi.txt \n");
		fprintf(f_summary,"\nReading DD from cache/DD_rp_pi.txt \n");
		for(int i=0;i<bin_rp;i++)	
			for(int j=0;j<bin_pi;j++)
				fscanf(f_DD,"%lu",&DD[i][j]);
		fclose(f_DD);
	}
	
	// ************* COMPUTING DR ****************************

	char DRName[100];
	sprintf(DRName,"cache/DR_rp_pi_%d.txt",sample_nr);
	FILE *f_DR=fopen(DRName,"r");
	
	if(!f_DR)	{
		fprintf(stdout,"\nSample %d: DR not computed! Computing DR and writing to DR_rp_pi.txt\n", sample_nr);
		fprintf(f_summary,"\nSample %d: DR not computed! Computing DR and writing to DR_rp_pi.txt\n", sample_nr);
		DR_count(ra_d,dec_d,z_d,n_d,ra_r,dec_r,z_r,n_r,DR);
		f_DR=fopen(DRName,"w");
		for(int i=0;i<bin_rp;i++)	{
			for(int j=0;j<bin_pi;j++)
				fprintf(f_DR,"%lu\t",DR[i][j]);
			fprintf(f_DR,"\n");
		}
		fclose(f_DR);
	}
	else	{
		fprintf(stdout,"\nReading DR from cache/DR_rp_pi.txt \n");
		fprintf(f_summary,"\nReading DR from cache/DR_rp_pi.txt \n");
		for(int i=0;i<bin_rp;i++)	
			for(int j=0;j<bin_pi;j++)
				fscanf(f_DR,"%lu",&DR[i][j]);
		fclose(f_DR);
	}
	
	// ******** COMPUTING XI(RP,PI) AND W(RP,PI) **********************************

	for(int j=0;j<bin_rp;j++)	{
		for(int k=0;k<bin_pi;k++)	{
			long double DDn = (long double)DD[j][k]/ndd;
			long double DRn = (long double)DR[j][k]/ndr;
			long double RRn = (long double)RR[j][k]/nrr;
				
			if(RR[j][k]==0)	{
				xi_2d[j][k]=NAN;
			}
			else	{	
				xi_2d[j][k]= (DDn-(2*DRn)+RRn)/RRn;
			}
		}
	}
	
	
	// ******** WRITING 2D XI(Rp,PI) AND W(Rp,PI) TO FILE ********************
	
	char xi2dName[100];
	sprintf(xi2dName,"cache/xi2d_%d.txt",sample_nr);
	
	FILE *fxi2d = fopen(xi2dName,"w");
	for(int j=0;j<bin_rp;j++)	{
		for(int k=0;k<bin_pi;k++)	{
			fprintf(fxi2d,"%Lf\t",xi_2d[j][k]);
		}
		fprintf(fxi2d,"\n");
	}
	fclose(fxi2d);
	
	// ******** COMPUTING PROJECTED CF BY INTEGRATION *************************

	for(int j=0;j<bin_rp;j++)	{
		wp[j]=0;
		for(int k=0;k<bin_pi;k++)	{
			wp[j]+=(dpi*xi_2d[j][k]);
		}		
		wp[j]*=2;
	}
	
	// writing rp, wp to file wpReal.txt and adding rp, wp to the matrix rp_wp_all[][]
	if(sample_nr == 0)	{
		FILE *wp_real=fopen("results/wpReal.txt","w");
		for(int i=0;i<bin_rp;i++)	{
			fprintf(wp_real,"%lf\t%Lf\n",pow(10,rp_bins[i]), wp[i]);
			rp_wp_all[i][0]=pow(10,rp_bins[i]);
			rp_wp_all[i][1]=wp[i];
		}
		fclose(wp_real);
	}
	else	{
		char wp_jk_name[100];
		sprintf(wp_jk_name,"results/jackknifes/wpJackknife_jk%d.txt",sample_nr);
		FILE *wp_jk_file=fopen(wp_jk_name,"w");
		
		// writing results rp, wp of the jackknife sample to wpJackknife<>.txt 
		for(int i=0;i<bin_rp;i++)	fprintf(wp_jk_file,"%lf\t%Lf\n",pow(10,rp_bins[i]),wp[i]);
		fclose(wp_jk_file);
	}
	
}

int main(void)
{
	
	clock_t begin=time(NULL);	// to compute total time taken by code
	
	if(n_boots != 0 && n_jacks != 0)	{
		fprintf(stderr,"\nError! Can't have both bootstrap copies and jackknife copies together..! \n");
		exit(1);
	}
	
	// *********************************
	
	mkdir("biproducts",0700);
	mkdir("results",0700);
	mkdir("cache",0700);
	mkdir("results/jackknifes",0700);
	
	f_summary=fopen("biproducts/summary.txt","w");
	
	// marking the current data and time
	time_t t = time(NULL);
	struct tm tm = *localtime(&t);
	fprintf(f_summary,"The code started on %d-%02d-%02d at %02d:%02d:%02d\n", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);
	
	// ********* SETTING UP BINS IN RP AND PI *************************************
	
	for(int i=0;i<bin_rp;i++)	
		rp_bins[i]=log10(rp_init)+(i*dlrp);
	
	for(int i=0;i<bin_pi;i++)
		pi_bins[i]=pi_init+(i*dpi);
		
	
	FILE *frps=fopen("biproducts/rps_log.txt","w");
	for(int i=0;i<bin_rp;i++)	fprintf(frps,"%lf\n",rp_bins[i]);
	fclose(frps);
	
	FILE *frps_lin=fopen("biproducts/rps_lin.txt","w");
	for(int i=0;i<bin_rp;i++)	fprintf(frps_lin,"%lf\n",pow(10,rp_bins[i]));
	fclose(frps_lin);

	FILE *fpis=fopen("biproducts/pis.txt","w");
	for(int i=0;i<bin_pi;i++)	fprintf(fpis,"%lf\n",pi_bins[i]);
	fclose(fpis);
	
	// ******** SETTING UPPER AND LOWER CUTOFFS FOR Rp AND PI *********************
	
	rp_lowcut=pow(10,rp_bins[0]-dlrp*0.5);
	rp_highcut=pow(10,rp_bins[bin_rp-1]+dlrp*0.5);
	pi_highcut=bin_pi+dpi*0.5;
		
			
	// ******** COUNTING THE TOTAL NUMBER OF REAL GALAXIES *************************
	
	unsigned long int n_d_real;
	
	char buf[255];
	FILE *f_realGal=fopen("real_galaxies","r");
	if(f_realGal == NULL)	{
		fprintf(stderr,"\nError opening real data file!");
		exit(1);
	}
	n_d_real=count_lines("real_galaxies")-1;		// -1 to remove the header
	rewind(f_realGal);	

	// ******* COUNTING THE NUMBER OF MARKS PRESENT IN DATA FILE ********** //

  	int ncols= count_cols("real_galaxies");
  	rewind(f_realGal);
  	
	n_marks = ncols-non_prop;
	
	rewind(f_realGal);
  	
	
	fprintf(stdout,"\nTotal number of real galaxies: %ld\n",n_d_real);
	fprintf(f_summary,"\nTotal number of real galaxies: %ld\n",n_d_real);
	
	fprintf(stdout,"\nTotal number of columns in real galaxies file: %d\n",ncols);
	fprintf(f_summary,"\nTotal number of columns in real galaxies file: %d\n",ncols);
	
	fprintf(stdout,"\nNumber of columns with quantities other than galaxy properties: %d\n",non_prop);
	fprintf(f_summary,"\nNumber of columns with quantities other than galaxy properties: %d\n",non_prop);
		
	fprintf(stdout,"\nNumber of properties used in the work: %d\n",n_marks);
	fprintf(f_summary,"\nNumber of properties used in the work: %d\n",n_marks);

	// *********** READING FILES AND COMPUTING CF PARALLELLY HAPPENING FROM HERE ******************* 
	
	#pragma omp parallel for
	for(int smplNr=0;smplNr<=n_jacks;smplNr++)	{
	
		// ******** READING THE REAL GALAXIES' DATA FILE *************************
		
		unsigned long int n_d;
		
		char fileRealName[100];
			
		if(smplNr==0)
			sprintf(fileRealName,"real_galaxies");
		else	
			sprintf(fileRealName,"jackknife_data/jk%d_real_galaxies",smplNr);
			
		FILE *fReal = fopen(fileRealName, "r");
		
		if(fReal == NULL) {
			fprintf(stderr,"\nReal galaxy %d not found!\n",smplNr);
			exit(1);
		}
			
		n_d = count_lines(fileRealName)-1;	// -1 to avoid counting the header
		
		char ch;											 	
		// ignoring the first header line
		do											
	  		ch = fgetc(fReal);			 
		while (ch != '\n');			
		
		double **real_gal = alloc_2d(n_d, ncols);
		readMatrixFromFile(fReal, real_gal, n_d, ncols);
		
		
		double *id_d = alloc_1d(n_d);
		double *ra_d = alloc_1d(n_d);
		double *dec_d = alloc_1d(n_d);
		double *z_d = alloc_1d(n_d);
		double *w_d = alloc_1d(n_d);
		double **prop_2d = alloc_2d(n_d, n_marks);

		for(int i=0;i<n_d;i++)	{
			id_d[i]=real_gal[i][0];
			ra_d[i]=real_gal[i][1];
			dec_d[i]=real_gal[i][2];
			z_d[i]=real_gal[i][3];
			w_d[i]=1.0;
			for(int j=0;j<n_marks;j++)
				prop_2d[i][j]=real_gal[i][j+non_prop];
			
		}
		fclose(fReal);
		free_2d(real_gal,n_d);
		
		// ************* READING RANDOM DATA ****************************
		
		unsigned long int n_r;
		
		char fileRandName[100];
			
		if(smplNr==0)
			sprintf(fileRandName,"random_galaxies");
		else	
			sprintf(fileRandName,"jackknife_data/jk%d_random_galaxies",smplNr);
		
		n_r = count_lines(fileRandName)-1;	// -1 to avoid counting the header
		
		fprintf(stdout,"\nTotal number of Random galaxies : %lu\n",n_r);
		fprintf(f_summary,"\nTotal number of Random galaxies : %lu\n",n_r);
		
		double *id_r = alloc_1d(n_r);
		double *ra_r = alloc_1d(n_r);
		double *dec_r = alloc_1d(n_r);
		double *z_r = alloc_1d(n_r);
		double *w_r = alloc_1d(n_r);
		
		FILE *fRand = fopen(fileRandName, "r");
		
		if(!fRand)	{	// if random file doesn't exist
			fprintf(stderr,"\nRandom data file not found!\n");
			fprintf(f_summary,"\nRandom data file not found!\n");
			exit(1);	
		}
		else	{	// if random file exits
			char ch;							 	
			do	// ignoring the first header line
				ch = fgetc(fRand);			 
			while (ch != '\n');	

															   
			for(int i=0;i<n_r;i++)
				fscanf(fRand,"%lf%lf%lf%lf%lf",&id_r[i],&ra_r[i],&dec_r[i],&z_r[i],&w_r[i]);
			fclose(fRand);
			
			
		}
		
		// compute wp for the real data and JK data
		compute2pCF(smplNr, n_d, ra_d, dec_d, z_d, n_r, ra_r, dec_r, z_r);
	

		free_1d(id_d);
		free_1d(ra_d);
		free_1d(dec_d);
		free_1d(z_d);
		free_1d(w_d);
		free_2d(prop_2d, n_d);
		free_1d(id_r);
		free_1d(ra_r);
		free_1d(dec_r);
		free_1d(z_r);
		free_1d(w_r);
		
	}
	
	
	if(n_boots != 0)	fprintf(f_summary,"\nNumber of bootstrap samples for wp : %d\n",n_boots);
	if(n_jacks != 0)	fprintf(f_summary,"\nNumber of jackknife samples for wp : %d\n",n_jacks);
	
	time_t t_f = time(NULL);
	struct tm tm_f = *localtime(&t_f);
	fprintf(stdout,"\nThe code finished successfully on %d-%02d-%02d at %02d:%02d:%02d\n", tm_f.tm_year + 1900, tm_f.tm_mon + 1, tm_f.tm_mday, tm_f.tm_hour, tm_f.tm_min, tm_f.tm_sec);
	fprintf(f_summary,"\nThe code finished successfully on %d-%02d-%02d at %02d:%02d:%02d\n", tm_f.tm_year + 1900, tm_f.tm_mon + 1, tm_f.tm_mday, tm_f.tm_hour, tm_f.tm_min, tm_f.tm_sec);
	
	clock_t end=time(NULL);
	int time_spent = (int)(end-begin);
	fprintf(stdout,"\nTotal time spent = %d hours %d minutes!!! \n", time_spent/3600,(time_spent%3600)/60);
	fprintf(f_summary,"\nTotal time spent = %d hours %d minutes!!! \n", time_spent/3600,(time_spent%3600)/60);
	fclose(f_summary);
	return 0;
}	





