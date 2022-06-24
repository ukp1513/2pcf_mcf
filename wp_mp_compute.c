// This code computes projected 2pt and marked correlation functions 
// of real data and bootstrap/jackknife samples.
// This code uses the weights of random galaxies too.
// Details are given in the README file.
// LAST UPDATE : 24 March 2022.
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

#define non_prop 4				// Number of columns in real file that are not to be used as marks 
								// These columns should be the first columns
								
#define boot_meth 2				// 1: individual; 2: blockwise

#define n_boots 0							// number of bootstrap samples
#define n_jacks 9				// Number of jackknife samples 
#define n_shuffles 100					// number of shuffled copies of marks 

#define bin_rp 16							// number of bins
#define bin_pi 40							// here goes the pmax

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

double rp_init=0.01;							// first bin
double pi_init=0.5;
double dlrp=0.27;							// binwidth in log scale
double dpi=1.0;

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

void sep(double ra1,double dec1,double z1,double ra2,double dec2,double z2,double *RP,double *PI)
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
	//th=(acos((sin(d1)*sin(d2))+(cos(d1)*cos(d2)*cos(a1-a2))));

	*RP=rp;
	*PI=pi;
}

void DD_count(double *ra1,double *dec1,double *z1, double *w1,int n1,double **count_array_DD)
{
	double rp,pi;
	for(int j=0;j<bin_rp;j++)	{
		for(int k=0;k<bin_pi;k++)	{
			count_array_DD[j][k]=0;
		}
	}
	for(int l=0;l<n1;l++)	{
		for(int m=0;m<n1;m++)	{
			if(l<m)	{
				sep(ra1[l],dec1[l],z1[l],ra1[m],dec1[m],z1[m],&rp,&pi);
				if(rp < rp_lowcut || rp > rp_highcut || pi > pi_highcut)	continue;
				if(rp != 0.0)	{
					rp=log10(rp);
					for(int j=0;j<bin_rp;j++)	{
						if(rp < (rp_bins[j]+dlrp*0.5) && rp > (rp_bins[j]-dlrp*0.5))		{
							for(int k=0;k<bin_pi;k++)	{
								if(pi < (pi_bins[k]+dpi*0.5) && pi > (pi_bins[k]-dpi*0.5))		{
									count_array_DD[j][k]+=w1[l]*w1[m];							
								}
							}
						}
					}
				}
			}
		}
	}
}

void DR_count(double *ra1,double *dec1,double *z1,double *w1,int n1,double *ra2,double *dec2,double *z2, double *w2,int n2,double **count_array_DR)
{
	double rp,pi;
	for(int j=0;j<bin_rp;j++)	{
		for(int k=0;k<bin_pi;k++)	{
			count_array_DR[j][k]=0;
		}
	}

	for(int l=0;l<n1;l++)	{
		for(int m=0;m<n2;m++)	{
			sep(ra1[l],dec1[l],z1[l],ra2[m],dec2[m],z2[m],&rp,&pi);
			if(rp < rp_lowcut || rp > rp_highcut || pi > pi_highcut)	continue;
			rp=log10(rp);
			for(int j=0;j<bin_rp;j++)	{
				if(rp < (rp_bins[j]+dlrp*0.5) && rp > (rp_bins[j]-dlrp*0.5))	{
					for(int k=0;k<bin_pi;k++)	{
						 if(pi < (pi_bins[k]+dpi*0.5) && pi > (pi_bins[k]-dpi*0.5))	{
							count_array_DR[j][k]+=w1[l]*w2[m];					
						}
					}
				}
			}
		}	
	}	
}



void DD_WW_count(double *ra1,double *dec1,double *z1,double *w1,double **m1,int n1,double **count_array_DD,double ***count_array_WW)
{
	double rp,pi;
	for(int j=0;j<bin_rp;j++)	{
		for(int k=0;k<bin_pi;k++)	{
			count_array_DD[j][k]=0;
			for(int l=0;l<n_shuffles+1;l++)
				count_array_WW[j][k][l]=0;
		}
	}
	for(int l=0;l<n1;l++)	{
    for(int m=0;m<n1;m++)	{
			if(l<m)	{
				sep(ra1[l],dec1[l],z1[l],ra1[m],dec1[m],z1[m],&rp,&pi);
				if(rp < rp_lowcut || rp > rp_highcut || pi > pi_highcut)	continue;
				if(rp != 0.0)	{
					rp=log10(rp);
					for(int j=0;j<bin_rp;j++)	{
						if(rp < (rp_bins[j]+dlrp*0.5) && rp > (rp_bins[j]-dlrp*0.5))		{
							for(int k=0;k<bin_pi;k++)	{
								if(pi < (pi_bins[k]+dpi*0.5) && pi > (pi_bins[k]-dpi*0.5))		{
									count_array_DD[j][k]+=w1[l]*w1[m];							
									for(int o=0;o<n_shuffles+1;o++)	
										count_array_WW[j][k][o]+=w1[l]*w1[m]*m1[l][o]*m1[m][o];
								}
							}
						}
					}
				}
			}
		}
	}
}

void DD_WW_count_JK(double *ra1,double *dec1,double *z1,double *w1,double **m1,int n1,double **count_array_DD,double ***count_array_WW)
{
	double rp,pi;
	for(int j=0;j<bin_rp;j++)	{
		for(int k=0;k<bin_pi;k++)	{
			count_array_DD[j][k]=0;
			for(int l=0;l<n_marks;l++)
				count_array_WW[j][k][l]=0;
		}
	}
	for(int l=0;l<n1;l++)	{
		for(int m=0;m<n1;m++)	{
			if(l<m)	{
				sep(ra1[l],dec1[l],z1[l],ra1[m],dec1[m],z1[m],&rp,&pi);
				if(rp < rp_lowcut || rp > rp_highcut || pi > pi_highcut)	continue;
				if(rp != 0.0)	{
					rp=log10(rp);
					for(int j=0;j<bin_rp;j++)	{
						if(rp < (rp_bins[j]+dlrp*0.5) && rp > (rp_bins[j]-dlrp*0.5))		{
							for(int k=0;k<bin_pi;k++)	{
								if(pi < (pi_bins[k]+dpi*0.5) && pi > (pi_bins[k]-dpi*0.5))		{
									count_array_DD[j][k]+=w1[l]*w1[m];							
									for(int o=0;o<n_marks;o++)	
										count_array_WW[j][k][o]+=w1[l]*w1[m]*m1[l][o]*m1[m][o];
								}
							}
						}
					}
				}
			}
		}
	}
}
	

void DR_WR_count(double *ra1,double *dec1,double *z1,double *w1, double **m1,int n1,double *ra2,double *dec2,double *z2, double *w2,int n2,double **count_array_DR,double ***count_array_WR)
{
	double rp,pi;
	for(int j=0;j<bin_rp;j++)	{
		for(int k=0;k<bin_pi;k++)	{
			count_array_DR[j][k]=0;
			for(int l=0;l<n_shuffles+1;l++)
				count_array_WR[j][k][l]=0;
		}
	}

	for(int l=0;l<n1;l++)	{
		for(int m=0;m<n2;m++)	{
			sep(ra1[l],dec1[l],z1[l],ra2[m],dec2[m],z2[m],&rp,&pi);
			if(rp < rp_lowcut || rp > rp_highcut || pi > pi_highcut)	continue;
			rp=log10(rp);
			for(int j=0;j<bin_rp;j++)	{
				if(rp < (rp_bins[j]+dlrp*0.5) && rp > (rp_bins[j]-dlrp*0.5))	{
					for(int k=0;k<bin_pi;k++)	{
						 if(pi < (pi_bins[k]+dpi*0.5) && pi > (pi_bins[k]-dpi*0.5))	{
							count_array_DR[j][k]+=w1[l]*w2[m];					
							for(int o=0;o<n_shuffles+1;o++)				
								count_array_WR[j][k][o]+=(w1[l]*w2[m]*m1[l][o]);
						}
					}
				}
			}
		}	
	}	
}	

void DR_WR_count_JK(double *ra1,double *dec1,double *z1,double *w1, double **m1,int n1,double *ra2,double *dec2,double *z2, double *w2,int n2,double **count_array_DR,double ***count_array_WR)
{
	double rp,pi;
	for(int j=0;j<bin_rp;j++)	{
		for(int k=0;k<bin_pi;k++)	{
			count_array_DR[j][k]=0;
			for(int l=0;l<n_marks;l++)
				count_array_WR[j][k][l]=0;
		}
	}

	for(int l=0;l<n1;l++)	{
		for(int m=0;m<n2;m++)	{
			sep(ra1[l],dec1[l],z1[l],ra2[m],dec2[m],z2[m],&rp,&pi);
			if(rp < rp_lowcut || rp > rp_highcut || pi > pi_highcut)	continue;
			rp=log10(rp);
			for(int j=0;j<bin_rp;j++)	{
				if(rp < (rp_bins[j]+dlrp*0.5) && rp > (rp_bins[j]-dlrp*0.5))	{
					for(int k=0;k<bin_pi;k++)	{
						 if(pi < (pi_bins[k]+dpi*0.5) && pi > (pi_bins[k]-dpi*0.5))	{
							count_array_DR[j][k]+=w1[l]*w2[m];					
							for(int o=0;o<n_marks;o++)				
								count_array_WR[j][k][o]+=(w1[l]*w2[m]*m1[l][o]);
						}
					}
				}
			}
		}	
	}	
}			

void RR_count(double *ra1,double *dec1,double *z1,double *w1,int n1,double **count_array_RR)
{
	
	double rp,pi;
	for(int j=0;j<bin_rp;j++)	
		for(int k=0;k<bin_pi;k++)	
			count_array_RR[j][k]=0;

	for(int l=0;l<n1;l++)	{
		for(int m=0;m<n1;m++)	{
			if(l<m)	{
				sep(ra1[l],dec1[l],z1[l],ra1[m],dec1[m],z1[m],&rp,&pi);
				if(rp < rp_lowcut || rp > rp_highcut || pi > pi_highcut)	continue;
				rp=log10(rp);
				for(int j=0;j<bin_rp;j++)	{
					if(rp < (rp_bins[j]+dlrp*0.5) && rp > (rp_bins[j]-dlrp*0.5))	{
						for(int k=0;k<bin_pi;k++)	{
							if(pi < (pi_bins[k]+dpi*0.5) && pi > (pi_bins[k]-dpi*0.5))	
								count_array_RR[j][k]+=w1[l]*w1[m];							
						}
					}
				}	
			}
		}	
	}
}

// computeJackknife() computes CFs for all jackknife samples
void computeJackknife(int jkNr, long int n_d_jk, double *ra_d_jk, double *dec_d_jk, double *z_d_jk, double *w_d_jk, long int n_r_jk, double *ra_r_jk, double *dec_r_jk, double *z_r_jk, double *w_r_jk)
{
	double **DD_jk=alloc_2d(bin_rp,bin_pi);
	double **DR_jk=alloc_2d(bin_rp,bin_pi);
	double **RR_jk=alloc_2d(bin_rp,bin_pi);
	double **DDn_jk=alloc_2d(bin_rp,bin_pi);
	double **DRn_jk=alloc_2d(bin_rp,bin_pi);
	double **RRn_jk=alloc_2d(bin_rp,bin_pi);
	double **xi_jk=alloc_2d(bin_rp,bin_pi);
	
	double ndd_jk=n_d_jk*(n_d_jk-1)*0.5;
	double nrr_jk=n_r_jk*(n_r_jk-1)*0.5;
	double ndr_jk=n_d_jk*n_r_jk;

	// ************* COMPUTING WW and WR **************************
	
	fprintf(stdout,"\nComputing DD for the jackknife sample %d\n", jkNr+1);
	DD_count(ra_d_jk,dec_d_jk,z_d_jk,w_d_jk,n_d_jk,DD_jk);
	fprintf(stdout,"\nComputing DR for the jackknife sample %d\n", jkNr+1);
	DR_count(ra_d_jk,dec_d_jk,z_d_jk,w_d_jk,n_d_jk,ra_r_jk,dec_r_jk,z_r_jk,w_r_jk,n_r_jk,DR_jk);
	fprintf(stdout,"\nComputing RR for the jackknife sample %d\n", jkNr+1);

	char RR_jk_name[100];
	sprintf(RR_jk_name,"biproducts/RR_jk%d.txt",jkNr+1);

	FILE *f_RR_jk=fopen(RR_jk_name,"r");
	if(!f_RR_jk)	{
		fprintf(stdout,"\nRR_jk%d.txt not found! Computing RR and writing to RR_%d.txt\n", jkNr+1,jkNr+1);
		RR_count(ra_r_jk,dec_r_jk,z_r_jk, w_r_jk,n_r_jk,RR_jk);
		f_RR_jk=fopen(RR_jk_name,"w");
		for(int i=0;i<bin_rp;i++)	{
			for(int j=0;j<bin_pi;j++)
				fprintf(f_RR_jk,"%lf\t",RR_jk[i][j]);
			fprintf(f_RR_jk,"\n");
		}
		fclose(f_RR_jk);
	}
	else	{
		fprintf(stdout,"\nReading RR from RR_jk%d.txt \n",jkNr+1);
		for(int i=0;i<bin_rp;i++)	
			for(int j=0;j<bin_pi;j++)
				fscanf(f_RR_jk,"%lf",&RR_jk[i][j]);
		fclose(f_RR_jk);
	}

	// ******** COMPUTING XI(RP,PI) AND W(RP,PI) **********************************

	for(int j=0;j<bin_rp;j++)	{
		for(int k=0;k<bin_pi;k++)	{
			DDn_jk[j][k] = DD_jk[j][k]/ndd_jk;
			DRn_jk[j][k] = DR_jk[j][k]/ndr_jk;
			RRn_jk[j][k] = RR_jk[j][k]/nrr_jk;
				
			if(RRn_jk==0)	xi_jk[j][k]=0;
			
			else	xi_jk[j][k]= (DDn_jk[j][k]-(2*DRn_jk[j][k])+RRn_jk[j][k])/RRn_jk[j][k];
			
		}
	}
	
	// ******** COMPUTING PROJECTED CF BY INTEGRATION *************************

	double wp_jk[bin_rp];
	for(int j=0;j<bin_rp;j++)	{
		wp_jk[j]=0;
		for(int k=0;k<bin_pi;k++)	wp_jk[j]+=(dpi*xi_jk[j][k]);
		wp_jk[j]*=2;
	}

	// *********** WRITING RESULTS **************************
	
	char wp_jk_name[50];
	sprintf(wp_jk_name,"results/jackknifes/wpJackknife_jk%d.txt",jkNr+1);
	FILE *wp_jk_file=fopen(wp_jk_name,"w");
	
	// writing results rp, wp of the jackknife sample to wpJackknife<>.txt 
	for(int i=0;i<bin_rp;i++)	fprintf(wp_jk_file,"%lf\t%lf\n",pow(10,rp_bins[i]),wp_jk[i]);
	fclose(wp_jk_file);
	
	free_2d(DD_jk,bin_rp);
	free_2d(DR_jk,bin_rp);
	free_2d(RR_jk,bin_rp);
	free_2d(DDn_jk,bin_rp);
	free_2d(DRn_jk,bin_rp);
	free_2d(RRn_jk,bin_rp);
	free_2d(xi_jk,bin_rp);
	
}

// computeBootstraps() computes CFs for all bootstrap samples using method 1
void computeBootstraps_meth1(int bootstrapSampleNr, long int n_d, double **ra_bs_full, double **dec_bs_full, double **z_bs_full, double **w_bs_full, long int n_r, double *ra_r, double *dec_r, double *z_r, double *w_r)
{
	double **DD=alloc_2d(bin_rp,bin_pi);
	double **DR=alloc_2d(bin_rp,bin_pi);
	double **DDn=alloc_2d(bin_rp,bin_pi);
	double **DRn=alloc_2d(bin_rp,bin_pi);
	double **RRn=alloc_2d(bin_rp,bin_pi);
	double **xi=alloc_2d(bin_rp,bin_pi);
	
	double *ra_bs = alloc_1d(n_d);
	double *dec_bs = alloc_1d(n_d);
	double *z_bs = alloc_1d(n_d);
	double *w_bs = alloc_1d(n_d);

	for(int i=0;i<n_d;i++)	{
		ra_bs[i]=ra_bs_full[i][bootstrapSampleNr];
		dec_bs[i]=dec_bs_full[i][bootstrapSampleNr];
		z_bs[i]=z_bs_full[i][bootstrapSampleNr];
		w_bs[i]=w_bs_full[i][bootstrapSampleNr];
	}

	double ndd=n_d*(n_d-1)*0.5;
	double nrr=n_r*(n_r-1)*0.5;
	double ndr=n_d*n_r;

	// ************* COMPUTING WW and WR **************************
	
	fprintf(stdout,"\nComputing DD for the bootstrap sample %d\n", bootstrapSampleNr+1);
	DD_count(ra_bs,dec_bs,z_bs, w_bs,n_d,DD);
	fprintf(stdout,"\nComputing DR for the bootstrap sample %d\n", bootstrapSampleNr+1);
	DR_count(ra_bs,dec_bs,z_bs, w_bs,n_d,ra_r,dec_r,z_r, w_r,n_r,DR);

	free_1d(ra_bs);
	free_1d(dec_bs);
	free_1d(z_bs);
	free_1d(w_bs);
		
	// ******** COMPUTING XI(RP,PI) AND W(RP,PI) **********************************

	

	for(int j=0;j<bin_rp;j++)	{
		for(int k=0;k<bin_pi;k++)	{
			DDn[j][k] = DD[j][k]/ndd;
			DRn[j][k] = DR[j][k]/ndr;
			RRn[j][k] = RR[j][k]/nrr;
				
			if(RRn==0)	xi[j][k]=0;
			
			else	xi[j][k]= (DDn[j][k]-(2*DRn[j][k])+RRn[j][k])/RRn[j][k];
			
		}
	}

	// ******** COMPUTING PROJECTED CF BY INTEGRATION *************************

	double wp[bin_rp];
	for(int j=0;j<bin_rp;j++)	{
		wp[j]=0;
		for(int k=0;k<bin_pi;k++)	wp[j]+=(dpi*xi[j][k]);
		wp[j]*=2;
	}

	// *********** WRITING RESULTS **************************
	
	char wp_bootstrap_name[50];
	sprintf(wp_bootstrap_name,"results/bootstraps/wpBootstrap_Bs%d.txt",bootstrapSampleNr+1);
	FILE *wp_bootstrap=fopen(wp_bootstrap_name,"w");
	
	// writing results rp, wp of the bootstrap sample to wpBootstrap<>.txt 
	for(int i=0;i<bin_rp;i++)	fprintf(wp_bootstrap,"%lf\t%lf\n",pow(10,rp_bins[i]),wp[i]);
	fclose(wp_bootstrap);
	
	free_2d(DD,bin_rp);
	free_2d(DR,bin_rp);
	free_2d(DDn,bin_rp);
	free_2d(DRn,bin_rp);
	free_2d(RRn,bin_rp);
	free_2d(xi,bin_rp);

}

// computeBootstraps() computes CFs for all bootstrap samples using method 2
void computeBootstraps_meth2(int bsNr, long int n_d_bs, double *ra_d_bs, double *dec_d_bs, double *z_d_bs, double *w_d_bs, long int n_r_bs, double *ra_r_bs, double *dec_r_bs, double *z_r_bs,double *w_r_bs)
{

	double **DD_bs=alloc_2d(bin_rp,bin_pi);
	double **DR_bs=alloc_2d(bin_rp,bin_pi);
	double **RR_bs=alloc_2d(bin_rp,bin_pi);
	double **DDn_bs=alloc_2d(bin_rp,bin_pi);
	double **DRn_bs=alloc_2d(bin_rp,bin_pi);
	double **RRn_bs=alloc_2d(bin_rp,bin_pi);
	double **xi_bs=alloc_2d(bin_rp,bin_pi);
	
	double ndd_bs=n_d_bs*(n_d_bs-1)*0.5;
	double nrr_bs=n_r_bs*(n_r_bs-1)*0.5;
	double ndr_bs=n_d_bs*n_r_bs;

	// ************* COMPUTING WW and WR **************************
	
	fprintf(stdout,"\nComputing DD for the bootstrap sample %d\n", bsNr+1);
	DD_count(ra_d_bs,dec_d_bs,z_d_bs,z_d_bs,n_d_bs,DD_bs);
	fprintf(stdout,"\nComputing DR for the bootstrap sample %d\n", bsNr+1);
	DR_count(ra_d_bs,dec_d_bs,z_d_bs,w_d_bs,n_d_bs,ra_r_bs,dec_r_bs,z_r_bs,w_r_bs,n_r_bs,DR_bs);
	fprintf(stdout,"\nComputing RR for the bootstrap sample %d\n", bsNr+1);
	RR_count(ra_r_bs,dec_r_bs,z_r_bs,w_r_bs,n_r_bs,RR_bs);


	// ******** COMPUTING XI(RP,PI) AND W(RP,PI) **********************************

	

	for(int j=0;j<bin_rp;j++)	{
		for(int k=0;k<bin_pi;k++)	{
			DDn_bs[j][k] = DD_bs[j][k]/ndd_bs;
			DRn_bs[j][k] = DR_bs[j][k]/ndr_bs;
			RRn_bs[j][k] = RR_bs[j][k]/nrr_bs;
				
			if(RRn_bs==0)	xi_bs[j][k]=0;
			
			else	xi_bs[j][k]= (DDn_bs[j][k]-(2*DRn_bs[j][k])+RRn_bs[j][k])/RRn_bs[j][k];
			
		}
	}

	// ******** COMPUTING PROJECTED CF BY INTEGRATION *************************

	double wp_bs[bin_rp];
	for(int j=0;j<bin_rp;j++)	{
		wp_bs[j]=0;
		for(int k=0;k<bin_pi;k++)	wp_bs[j]+=(dpi*xi_bs[j][k]);
		wp_bs[j]*=2;
	}

	// *********** WRITING RESULTS **************************
	
	char wp_bootstrap_name[50];
	sprintf(wp_bootstrap_name,"results/bootstraps/wpBootstrap_Bs%d.txt",bsNr+1);
	FILE *wp_bootstrap=fopen(wp_bootstrap_name,"w");
	
	// writing results rp, wp of the bootstrap sample to wpBootstrap<>.txt 
	for(int i=0;i<bin_rp;i++)	fprintf(wp_bootstrap,"%lf\t%lf\n",pow(10,rp_bins[i]),wp_bs[i]);
	fclose(wp_bootstrap);
	
	free_2d(DD_bs,bin_rp);
	free_2d(DR_bs,bin_rp);
	free_2d(RR_bs,bin_rp);
	free_2d(DDn_bs,bin_rp);
	free_2d(DRn_bs,bin_rp);
	free_2d(RRn_bs,bin_rp);
	free_2d(xi_bs,bin_rp);
}

void computeJKwithoutshuffles(int jkNr, long int n_d_jk, double *ra_d_jk, double *dec_d_jk, double *z_d_jk,double *w_d_jk,double **w_mark_jk, long int n_r_jk, double *ra_r_jk, double *dec_r_jk, double *z_r_jk, double *w_r_jk)
{
	double **DD_jk=alloc_2d(bin_rp,bin_pi);
	double **DR_jk=alloc_2d(bin_rp,bin_pi);
	double **RR_jk=alloc_2d(bin_rp,bin_pi);
	double **DDn_jk=alloc_2d(bin_rp,bin_pi);
	double **DRn_jk=alloc_2d(bin_rp,bin_pi);
	double **RRn_jk=alloc_2d(bin_rp,bin_pi);
	double **xi_jk=alloc_2d(bin_rp,bin_pi);
	
	double ***WW_jk=alloc_3d(bin_rp,bin_pi,n_marks);
	double ***WR_jk=alloc_3d(bin_rp,bin_pi,n_marks);
	double ***WWn_jk=alloc_3d(bin_rp,bin_pi,n_marks);
	double ***WRn_jk=alloc_3d(bin_rp,bin_pi,n_marks);
	double ***W_jk=alloc_3d(bin_rp,bin_pi,n_marks);
	
	double ndd_jk=n_d_jk*(n_d_jk-1)*0.5;
	double nrr_jk=n_r_jk*(n_r_jk-1)*0.5;
	double ndr_jk=n_d_jk*n_r_jk;
	
	// ************* COMPUTING WW and WR **************************
	
	fprintf(stdout,"\nComputing DD for the jackknife sample %d\n", jkNr+1);
	DD_WW_count_JK(ra_d_jk,dec_d_jk,z_d_jk,w_d_jk,w_mark_jk,n_d_jk,DD_jk,WW_jk);
	fprintf(stdout,"\nComputing DR for the jackknife sample %d\n", jkNr+1);
	DR_WR_count_JK(ra_d_jk,dec_d_jk,z_d_jk,w_d_jk,w_mark_jk,n_d_jk,ra_r_jk,dec_r_jk,z_r_jk,w_r_jk,n_r_jk,DR_jk,WR_jk);
	
	char RR_jk_name[100];
	sprintf(RR_jk_name,"biproducts/RR_jk%d.txt",jkNr+1);

	FILE *f_RR_jk=fopen(RR_jk_name,"r");
	if(!f_RR_jk)	{
		fprintf(stdout,"\nRR_jk%d.txt not found! Computing RR and writing to RR_%d.txt\n", jkNr+1,jkNr+1);
		RR_count(ra_r_jk,dec_r_jk,z_r_jk, w_r_jk,n_r_jk,RR_jk);
		f_RR_jk=fopen(RR_jk_name,"w");
		for(int i=0;i<bin_rp;i++)	{
			for(int j=0;j<bin_pi;j++)
				fprintf(f_RR_jk,"%lf\t",RR_jk[i][j]);
			fprintf(f_RR_jk,"\n");
		}
		fclose(f_RR_jk);
	}
	else	{
		fprintf(stdout,"\nReading RR from RR_jk%d.txt \n",jkNr+1);
		for(int i=0;i<bin_rp;i++)	
			for(int j=0;j<bin_pi;j++)
				fscanf(f_RR_jk,"%lf",&RR_jk[i][j]);
		fclose(f_RR_jk);
	}


	// ******** COMPUTING XI(RP,PI) AND W(RP,PI) **********************************



	for(int j=0;j<bin_rp;j++)	{
		for(int k=0;k<bin_pi;k++)	{
			DDn_jk[j][k] = DD_jk[j][k]/ndd_jk;
			DRn_jk[j][k] = DR_jk[j][k]/ndr_jk;
			RRn_jk[j][k] = RR_jk[j][k]/nrr_jk;
				
			for(int l=0;l<n_marks;l++)	{
				WWn_jk[j][k][l] = WW_jk[j][k][l]/ndd_jk;
				WRn_jk[j][k][l] = WR_jk[j][k][l]/ndr_jk;
			}

			if(RRn_jk==0)	{
				xi_jk[j][k]=0;
				for(int l=0;l<n_marks;l++)	W_jk[j][k][l]=0;	
			}
			else	{	
				xi_jk[j][k]= (DDn_jk[j][k]-(2*DRn_jk[j][k])+RRn_jk[j][k])/RRn_jk[j][k];
				for(int l=0;l<n_marks;l++)
					W_jk[j][k][l] = (WWn_jk[j][k][l]-(2*WRn_jk[j][k][l])+RRn_jk[j][k])/RRn_jk[j][k];
			}
		}
	}
	
	// ******** COMPUTING PROJECTED CF BY INTEGRATION *************************

	double wp_jk[bin_rp],Wp_jk[bin_rp][n_marks],Mp_jk[bin_rp][n_marks];
	for(int j=0;j<bin_rp;j++)	{
		wp_jk[j]=0;
		for(int l=0;l<n_marks;l++)		Wp_jk[j][l]=0;
		for(int k=0;k<bin_pi;k++)	{
			wp_jk[j]+=(dpi*xi_jk[j][k]);
			for(int l=0;l<n_marks;l++)
				Wp_jk[j][l]+=(dpi*W_jk[j][k][l]);
		}		
		wp_jk[j]*=2;
		for(int l=0;l<n_marks;l++)	Wp_jk[j][l]*=2;
	}
	

	
	for(int j=0;j<bin_rp;j++)	{
 		for(int l=0;l<n_marks;l++)	{
			Mp_jk[j][l] = (1+Wp_jk[j][l]/pow(10,rp_bins[j]))/(1+wp_jk[j]/pow(10,rp_bins[j]));
		}
	}


	// rp_wp_mp_alljk stores the columns rp, wp, mp(mark1), mp(mark2),...
	double rp_wp_mp_alljk[bin_rp][n_marks+2];
	
	for(int i=0;i<bin_rp;i++)	{
		rp_wp_mp_alljk[i][0]=pow(10,rp_bins[i]);
		rp_wp_mp_alljk[i][1]=wp_jk[i];
		for(int l=2;l<n_marks+2;l++)	
			rp_wp_mp_alljk[i][l]=Mp_jk[i][l-2];
	}	
	
	// writing rp, wp, mp(mark1), mp(mark2), ... to file wpJackknife_jk<jkNr>.txt
	
	char wp_jk_name[50];
	sprintf(wp_jk_name,"results/jackknifes/wpMpJackknife_jk%d.txt",jkNr+1);
	FILE *fwp_jk=fopen(wp_jk_name,"w");

	for(int i=0;i<bin_rp;i++)	{
		for(int l=0;l<n_marks+2;l++)	
			fprintf(fwp_jk,"%lf\t",rp_wp_mp_alljk[i][l]);
		fprintf(fwp_jk,"\n");
	}	
	fclose(fwp_jk);
	
	free_2d(DD_jk,bin_rp);
	free_2d(DR_jk,bin_rp);
	free_2d(RR_jk,bin_rp);
	free_2d(DDn_jk,bin_rp);
	free_2d(DRn_jk,bin_rp);
	free_2d(RRn_jk,bin_rp);
	free_2d(xi_jk,bin_rp);
	
	free_3d(WW_jk,bin_rp,bin_pi);
	free_3d(WR_jk,bin_rp,bin_pi);
	free_3d(WWn_jk,bin_rp,bin_pi);
	free_3d(WRn_jk,bin_rp,bin_pi);
	free_3d(W_jk,bin_rp,bin_pi);
	
}

// computeReal() computes CF and MCFs for <n_shuffles> shuffled marks
void computeReal(int mark_nr, long int n_d, double *ra_d, double *dec_d, double *z_d,double *w_d,double **w_mark, long int n_r, double *ra_r, double *dec_r, double *z_r, double *w_r)
{
	double **DD=alloc_2d(bin_rp,bin_pi);
	double **DR=alloc_2d(bin_rp,bin_pi);
	double **DDn=alloc_2d(bin_rp,bin_pi);
	double **DRn=alloc_2d(bin_rp,bin_pi);
	double **RRn=alloc_2d(bin_rp,bin_pi);
	double **xi=alloc_2d(bin_rp,bin_pi);
	
	double ***WW = alloc_3d(bin_rp, bin_pi, n_shuffles+1);
	double ***WR = alloc_3d(bin_rp, bin_pi, n_shuffles+1);
	double ***WWn = alloc_3d(bin_rp, bin_pi, n_shuffles+1);
	double ***WRn = alloc_3d(bin_rp, bin_pi, n_shuffles+1);
	double ***W = alloc_3d(bin_rp, bin_pi, n_shuffles+1);
	
	double ndd=n_d*(n_d-1)*0.5;
	double nrr=n_r*(n_r-1)*0.5;
	double ndr=n_d*n_r;
	
	// ************* COMPUTING WW and WR **************************

	fprintf(stdout,"\nComputing DD and WW for the real sample : Property %d\n",mark_nr+1);
	DD_WW_count(ra_d,dec_d,z_d,w_d,w_mark,n_d,DD,WW);
	fprintf(stdout,"\nComputing DR and WR for the real sample : Property %d\n", mark_nr+1);
	DR_WR_count(ra_d,dec_d,z_d,w_d,w_mark,n_d,ra_r,dec_r,z_r,w_r,n_r,DR,WR);


	// ******** COMPUTING XI(RP,PI) AND W(RP,PI) **********************************



	for(int j=0;j<bin_rp;j++)	{
		for(int k=0;k<bin_pi;k++)	{
			DDn[j][k] = DD[j][k]/ndd;
			DRn[j][k] = DR[j][k]/ndr;
			RRn[j][k] = RR[j][k]/nrr;
				
			for(int l=0;l<n_shuffles+1;l++)	{
				WWn[j][k][l] = WW[j][k][l]/ndd;
				WRn[j][k][l] = WR[j][k][l]/ndr;
			}

			if(RRn==0)	{
				xi[j][k]=0;
				for(int l=0;l<n_shuffles+1;l++)	W[j][k][l]=0;	
			}
			else	{	
				xi[j][k]= (DDn[j][k]-(2*DRn[j][k])+RRn[j][k])/RRn[j][k];
				for(int l=0;l<n_shuffles+1;l++)
					W[j][k][l] = (WWn[j][k][l]-(2*WRn[j][k][l])+RRn[j][k])/RRn[j][k];
			}
		}
	}
	
	// ******** WRITING 2D XI(Rp,PI) AND W(Rp,PI) TO FILE ********************
	
	FILE *fxi2d = fopen("biproducts/xi2d.txt","w");
	FILE *fW2d = fopen("biproducts/W2d.txt","w");
	for(int j=0;j<bin_rp;j++)	{
		for(int k=0;k<bin_pi;k++)	{
			fprintf(fxi2d,"%lf\t",xi[j][k]);
			fprintf(fW2d,"%lf\t",W[j][k][0]);
		}
		fprintf(fxi2d,"\n");
		fprintf(fW2d,"\n");
	}
	fclose(fxi2d);
	fclose(fW2d);

	// ******** COMPUTING PROJECTED CF BY INTEGRATION *************************

	double wp[bin_rp],Wp[bin_rp][n_shuffles+1],Mp[bin_rp][n_shuffles+1];
	for(int j=0;j<bin_rp;j++)	{
		wp[j]=0;
		for(int l=0;l<n_shuffles+1;l++)		Wp[j][l]=0;
		for(int k=0;k<bin_pi;k++)	{
			wp[j]+=(dpi*xi[j][k]);
			for(int l=0;l<n_shuffles+1;l++)
				Wp[j][l]+=(dpi*W[j][k][l]);
		}		
		wp[j]*=2;
		for(int l=0;l<n_shuffles+1;l++)	Wp[j][l]*=2;
	}
	
	for(int j=0;j<bin_rp;j++)	{
 		for(int l=0;l<n_shuffles+1;l++)	{
			Mp[j][l] = (1+Wp[j][l]/pow(10,rp_bins[j]))/(1+wp[j]/pow(10,rp_bins[j]));
		}
	}
	
	// rp_Wp_all stores the columns rp, Wp(real), Wp(shuffle1), Wp(shuffle2),...
	double rp_Wp_all[bin_rp][n_shuffles+2];
	for(int i=0;i<bin_rp;i++)	{
		rp_Wp_all[i][0]=pow(10,rp_bins[i]);
		for(int l=1;l<n_shuffles+2;l++)	
			rp_Wp_all[i][l]=Wp[i][l-1];
	}	
	
	double rp_Mp_all[bin_rp][n_shuffles+2];
	for(int i=0;i<bin_rp;i++)	{
		rp_Mp_all[i][0]=pow(10,rp_bins[i]);
		for(int l=1;l<n_shuffles+2;l++)	
			rp_Mp_all[i][l]=Mp[i][l-1];
	}
	
	if(mark_nr == 0)	{
		// writing rp, wp to file wpReal.txt and adding rp, wp to the matrix rp_wp_all[][]
		FILE *wp_real=fopen("results/wpReal.txt","w");
		for(int i=0;i<bin_rp;i++)	{
			fprintf(wp_real,"%lf\t%lf\n",pow(10,rp_bins[i]), wp[i]);
			rp_wp_all[i][0]=pow(10,rp_bins[i]);
			rp_wp_all[i][1]=wp[i];
		}
		fclose(wp_real);
	}
	
	char Wp_name[50], Mp_name[50];
	// writing rp, Wp(real), Wp(shuffle1), Wp(shuffle2), ... to file WpRealAll.txt
	// and rp, Mp(real), Mp(shuffle1), Mp(shuffle2), ... to file mpRealAll.txt
	sprintf(Wp_name,"results/WpRealAllShuffle_mk%d.txt",mark_nr+1);
	sprintf(Mp_name,"results/mpRealAllShuffle_mk%d.txt",mark_nr+1);
	FILE *Wp_real=fopen(Wp_name,"w");
	FILE *Mp_real=fopen(Mp_name,"w");
	for(int i=0;i<bin_rp;i++)	{
		for(int l=0;l<n_shuffles+2;l++)	{
			fprintf(Wp_real,"%lf\t",rp_Wp_all[i][l]);
			fprintf(Mp_real,"%lf\t",rp_Mp_all[i][l]);
		}
		fprintf(Wp_real,"\n");
		fprintf(Mp_real,"\n");
	}	
	fclose(Wp_real);
	fclose(Mp_real);
	
	free_2d(DD,bin_rp);
	free_2d(DR,bin_rp);
	free_2d(DDn,bin_rp);
	free_2d(DRn,bin_rp);
	free_2d(RRn,bin_rp);
	free_2d(xi,bin_rp);
	
	free_3d(WW,bin_rp,bin_pi);
	free_3d(WR,bin_rp,bin_pi);
	free_3d(WWn,bin_rp,bin_pi);
	free_3d(WRn,bin_rp,bin_pi);
	free_3d(W,bin_rp,bin_pi);
	
}

int main(void)
{
	
	clock_t begin=time(NULL);	// to compute total time taken by code
	
	if(n_boots != 0 && n_jacks != 0)	{
		fprintf(stderr,"\nError! Can't have both bootstrap copies and jackknife copies together..! \n");
		exit(1);
	}
	
	// ***************************************************************
	
	long int n_d;
	int randoms;
	
	mkdir("biproducts",0700);
	mkdir("results",0700);
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
	
	char buf[255];
	FILE *f_realGal=fopen("real_galaxies","r");
	if(f_realGal == NULL)	{
		fprintf(stderr,"\nError opening real data file!");
		exit(1);
	}
	n_d=count_lines("real_galaxies")-1;		// -1 to remove the header
	rewind(f_realGal);	

	// ******* COUNTING THE NUMBER OF MARKS PRESENT IN DATA FILE ********** //

  	int ncols= count_cols("real_galaxies");
  	rewind(f_realGal);
  	
	n_marks = ncols-non_prop;
	
	rewind(f_realGal);
  	
	
	fprintf(stdout,"\nTotal number of real galaxies: %ld\n",n_d);
	fprintf(f_summary,"\nTotal number of real galaxies: %ld\n",n_d);
	
	fprintf(stdout,"\nTotal number of columns in real galaxies file: %d\n",ncols);
	fprintf(f_summary,"\nTotal number of columns in real galaxies file: %d\n",ncols);
	
	fprintf(stdout,"\nNumber of columns with quantities other than galaxy properties: %d\n",non_prop);
	fprintf(f_summary,"\nNumber of columns with quantities other than galaxy properties: %d\n",non_prop);
		
	fprintf(stdout,"\nNumber of properties used in the work: %d\n",n_marks);
	fprintf(f_summary,"\nNumber of properties used in the work: %d\n",n_marks);


	
	// ******** READING THE REAL GALAXIES' DATA FILE *************************
	
	char ch;											 	
	// ignoring the first header line
	do											
  		ch = fgetc(f_realGal);			 
	while (ch != '\n');			
	
	double **real_gal = alloc_2d(n_d, ncols);
	readMatrixFromFile(f_realGal, real_gal, n_d, ncols);
	
	
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
	fclose(f_realGal);
	free_2d(real_gal,n_d);
	
	// ************* READING RANDOM DATA ****************************
	
	RR=alloc_2d(bin_rp,bin_pi);
	
	long int n_r;
	
	
	n_r = count_lines("random_galaxies")-1;		// -1 to remove the header
	

	randoms=n_r/n_d;

	fprintf(stdout,"\nTotal number of Random galaxies : %ld\n",n_r);
	fprintf(f_summary,"\nTotal number of Random galaxies : %ld\n",n_r);
	
	double *id_r = alloc_1d(n_r);
	double *ra_r = alloc_1d(n_r);
	double *dec_r = alloc_1d(n_r);
	double *z_r = alloc_1d(n_r);
	double *w_r = alloc_1d(n_r);
	
	FILE *f_randGal=fopen("random_galaxies","r");
	if(!f_randGal)	{	// if random file doesn't exist
		fprintf(stderr,"\nRandom data file not found!\n");
		fprintf(f_summary,"\nRandom data file not found!\n");
		exit(1);	
	}
	else	{	// if random file exits
		char ch;							 	
		do	// ignoring the first header line
			ch = fgetc(f_randGal);			 
		while (ch != '\n');	

														   
		for(int i=0;i<n_r;i++)
			fscanf(f_randGal,"%lf%lf%lf%lf%lf",&id_r[i],&ra_r[i],&dec_r[i],&z_r[i],&w_r[i]);
		fclose(f_randGal);
		
		double rand_weight_mean = array_avg(w_r,n_r); 
		
		// normalising all random weights to mean
		for(int i=0;i<n_r;i++)	{
				w_r[i]/=rand_weight_mean;
		}
		
	}

	// ************* COMPUTING RR ****************************

	FILE *f_RR=fopen("biproducts/RR.txt","r");
	if(!f_RR)	{
		fprintf(stdout,"\nbiproducts/RR.txt not found! Computing RR and writing to RR.txt\n");
		fprintf(f_summary,"\nbiproducts/RR.txt not found! Computing RR and writing to RR.txt\n");
		RR_count(ra_r,dec_r,z_r,w_r,n_r,RR);
		f_RR=fopen("biproducts/RR.txt","w");
		for(int i=0;i<bin_rp;i++)	{
			for(int j=0;j<bin_pi;j++)
				fprintf(f_RR,"%lf\t",RR[i][j]);
			fprintf(f_RR,"\n");
		}
		fclose(f_RR);
	}
	else	{
		fprintf(stdout,"\nReading RR from biproducts/RR.txt \n");
		fprintf(f_summary,"\nReading RR from biproducts/RR.txt \n");
		for(int i=0;i<bin_rp;i++)	
			for(int j=0;j<bin_pi;j++)
				fscanf(f_RR,"%lf",&RR[i][j]);
		fclose(f_RR);
	}	
			
	// *************  ASSIGNING MARKS AND COMPUTING FOR REAL SAMPLE *********************************

	#pragma omp parallel for
	for(int mark_nr=0;mark_nr<n_marks;mark_nr++)	{
		
		double **prop_2d_mark = alloc_2d(n_d,n_marks);
		equateMatrix(prop_2d_mark, prop_2d,n_d,n_marks);
		
		double *mark_d = alloc_1d(n_d);
		double *prop_d = alloc_1d(n_d);
		
		// extracting marks from prop_2d array
		
		for(int i=0;i<n_d;i++)	prop_d[i]=prop_2d_mark[i][mark_nr];
		
		// making all marks positive if negative (in case of absolute magnitude)
		if(prop_d[0] < 0)	{
			for(int i=0;i<n_d;i++)	prop_d[i]=-1*prop_d[i]; 
		}

		// ranking the properties
		rankify(prop_d,n_d,mark_d); 
		
		//equateArray(mark_d,prop_d,n_d);
		
		// shuffledMark: marks shuffled instantaneously
		double *shuffledMark = alloc_1d(n_d);
	
		double **weight_mark=alloc_2d(n_d,n_shuffles+1);
		
		// first column of weight[][] be the original marks
		for(int i=0;i<n_d;i++)	weight_mark[i][0]=mark_d[i];
		
		// subsequent columns be the shuffled marks each time
		for(int i=0;i<n_d;i++)	shuffledMark[i]=mark_d[i];
		for(int j=0;j<n_shuffles;j++)	{
			shuffle(shuffledMark,n_d);
			for(int i=0;i<n_d;i++)	weight_mark[i][j+1]=shuffledMark[i];
		}
		
		// mean value of the marks in whole sample to normalize marks
		double weightScheme = array_avg(mark_d,n_d); 
		
		// normalising all marks (including shuffled) to mean
		for(int i=0;i<n_d;i++)	{
			for(int j=0;j<n_shuffles+1;j++)
				weight_mark[i][j]/=weightScheme;
		}
		
		// compute wp and MCF (with <nshuffle> repetitions) for the real data
		computeReal(mark_nr, n_d, ra_d, dec_d, z_d, w_d, weight_mark, n_r, ra_r, dec_r, z_r, w_r);

		free_2d(prop_2d_mark,n_d);
		free_2d(weight_mark,n_d);
		free_1d(prop_d);
		free_1d(mark_d);
		free_1d(shuffledMark);
	}
	
	// ************* COMPUTING FOR JACKKNIFE SAMPLES ***********
	
	if(n_jacks > 0)	{
		
		fprintf(stdout,"\nJackknife computation starts here parallely!\n");
		mkdir("results/jackknifes",0700);
		
		// executing jackknife computations parallely
		
		#pragma omp parallel for
		for(int jkNr=0;jkNr<n_jacks;jkNr++)	{
			long int n_d_jk, n_r_jk;
			double *id_d_jk,*ra_d_jk, *dec_d_jk, *z_d_jk, *w_d_jk;
			double *id_r_jk, *ra_r_jk, *dec_r_jk, *z_r_jk, *w_r_jk;
			double **prop_2d_jk;
			char fileJKReal[100],fileJKRan[100];
			
			sprintf(fileJKReal,"jackknife_data/jk%d_real_galaxies",jkNr+1);
			
			FILE *fJKReal = fopen(fileJKReal, "r");
			
			if(fJKReal == NULL) {
				fprintf(stderr,"\nJackknife real galaxy %d not found!\n",jkNr+1);
				exit(1);
			}
			
			n_d_jk = count_lines(fileJKReal)-1;	// -1 to avoid counting the header

			char ch;							 	
			do	// ignoring the first header line
				ch = fgetc(fJKReal);			 
			while (ch != '\n');

			double **real_jk_gal = alloc_2d(n_d_jk, ncols);
			
			if(!real_jk_gal) {
				fprintf(stderr,"\nJackknife real galaxy %d not found!\n",jkNr+1);
				exit(1);
			}
			
			readMatrixFromFile(fJKReal, real_jk_gal, n_d_jk, ncols);
			fclose(fJKReal);

			id_d_jk = alloc_1d(n_d_jk);
			ra_d_jk = alloc_1d(n_d_jk);
			dec_d_jk = alloc_1d(n_d_jk);
			z_d_jk = alloc_1d(n_d_jk);
			w_d_jk = alloc_1d(n_d_jk);
			prop_2d_jk = alloc_2d(n_d_jk, n_marks);

			for(int i=0;i<n_d_jk;i++)	{
				id_d_jk[i]=real_jk_gal[i][0];
				ra_d_jk[i]=real_jk_gal[i][1];
				dec_d_jk[i]=real_jk_gal[i][2];
				z_d_jk[i]=real_jk_gal[i][3];
				w_d_jk[i]=1.0;
				for(int j=0;j<n_marks;j++)
					prop_2d_jk[i][j]=real_jk_gal[i][j+non_prop];
	
			}
			
			sprintf(fileJKRan,"jackknife_data/jk%d_random_galaxies",jkNr+1);
			FILE *fJKRan = fopen(fileJKRan, "r");
			if(fJKRan)	{
				n_r_jk = count_lines(fileJKRan)-1;	// -1 to avoid counting the header
				id_r_jk = alloc_1d(n_r_jk);
				ra_r_jk = alloc_1d(n_r_jk);
				dec_r_jk = alloc_1d(n_r_jk);
				z_r_jk = alloc_1d(n_r_jk);
				w_r_jk=alloc_1d(n_r_jk);
				
				char ch;							 	
				do	// ignoring the first header line
					ch = fgetc(fJKRan);			 
				while (ch != '\n');
				
				for(int i=0;i<n_r_jk;i++)
					fscanf(fJKRan,"%lf%lf%lf%lf%lf",&id_r_jk[i],&ra_r_jk[i],&dec_r_jk[i],&z_r_jk[i],&w_r_jk[i]);
				fclose(fJKRan);
				
			}
			else {
				fprintf(stderr,"\nJackknife random galaxy %d not found!\n",jkNr+1);
				exit(1);
			}
			
			fprintf(f_summary,"\nJK %d: N_real: %ld, N_random: %ld\n",jkNr+1,n_d_jk,n_r_jk);
			
			double rand_weight_mean_jk = array_avg(w_r_jk,n_r_jk); 
		
			// normalising all marks to mean
			for(int i=0;i<n_r_jk;i++)	{
					w_r_jk[i]/=rand_weight_mean_jk;
			}
			
			double **prop_2d_jk_mark = alloc_2d(n_d_jk,n_marks);
			
			for(int mark_nr=0;mark_nr<n_marks;mark_nr++)	{
				
				double *mark_d_jk = alloc_1d(n_d_jk);
				double *prop_d_jk = alloc_1d(n_d_jk);
				
				// extracting marks from prop_2d array
				for(int j=0;j<n_d_jk;j++)	prop_d_jk[j]=prop_2d_jk[j][mark_nr];
				
				// making all marks positive if negative (in case of absolute magnitude)
				if(prop_d_jk[0] < 0)	{
					for(int i=0;i<n_d_jk;i++)	prop_d_jk[i]=-1*prop_d_jk[i];
				}
				
				// ranking the properties
				rankify(prop_d_jk,n_d_jk,mark_d_jk);
				
				double weightScheme = array_avg(mark_d_jk,n_d_jk);
				
				for(int j=0;j<n_d_jk;j++)	mark_d_jk[j]/=weightScheme;
				
				for(int j=0;j<n_d_jk;j++)	prop_2d_jk_mark[j][mark_nr]=mark_d_jk[j];
			
				free_1d(mark_d_jk);
				free_1d(prop_d_jk);
			
			}
			
			computeJKwithoutshuffles(jkNr,n_d_jk,ra_d_jk,dec_d_jk,z_d_jk,w_d_jk,prop_2d_jk_mark,n_r_jk, ra_r_jk,dec_r_jk,z_r_jk, w_r_jk);
			
			free_1d(id_d_jk);
			free_1d(ra_d_jk);
			free_1d(dec_d_jk);
			free_1d(z_d_jk);
			free_1d(w_d_jk);
			free_1d(id_r_jk);
			free_1d(ra_r_jk);
			free_1d(dec_r_jk);
			free_1d(z_r_jk);		
			free_1d(w_r_jk);
			
			
			free_2d(real_jk_gal, n_d_jk);
			free_2d(prop_2d_jk, n_d_jk);
			free_2d(prop_2d_jk_mark, n_d_jk);
		}
	}

	// *************COMPUTING FOR BOOTSTRAP SAMPLES*******************
	
	if(boot_meth == 1)	{
	
		if(n_boots > 0)	{
			fprintf(stdout,"\nBootstrap (individual) computation starts here parallely!\n");
			mkdir("results/bootstraps",0700);
		
			
			double **ra_bs_full = alloc_2d(n_d, n_boots);
			double **dec_bs_full = alloc_2d(n_d, n_boots);
			double **z_bs_full = alloc_2d(n_d, n_boots);
			double **w_bs_full = alloc_2d(n_d, n_boots);
			
			// creating bootstrap samples
			srand((unsigned) time(NULL));
			long int bs_index;
			for(int bootstrapSampleNr = 0;bootstrapSampleNr < n_boots;bootstrapSampleNr++)	{
				for(int i=0;i<n_d;i++)	{
					bs_index=rand()%(n_d);
					ra_bs_full[i][bootstrapSampleNr]=ra_d[bs_index];
					dec_bs_full[i][bootstrapSampleNr]=dec_d[bs_index];
					z_bs_full[i][bootstrapSampleNr]=z_d[bs_index];
					w_bs_full[i][bootstrapSampleNr]=w_d[bs_index];
				}
			}
			
			// executing bootstrap computations parallely
			#pragma omp parallel for
			for(int bootstrapSampleNr=0;bootstrapSampleNr<n_boots;bootstrapSampleNr++)	{
				if(bootstrapSampleNr==0)	fprintf(f_summary,"\nTotal number of parallel threads: %d\n",omp_get_num_threads());
				computeBootstraps_meth1(bootstrapSampleNr,n_d,ra_bs_full,dec_bs_full,z_bs_full,w_bs_full,n_r, ra_r,dec_r,z_r,w_r);		
			}
			
			free_2d(ra_bs_full,n_d);
			free_2d(dec_bs_full,n_d);
			free_2d(z_bs_full,n_d);
			free_2d(w_bs_full,n_d);
		}
	}
	
	if(boot_meth == 2)	{
		if(n_boots > 0)	{
			
			fprintf(stdout,"\nBootstrap (blockwise) computation starts here parallely!\n");
			mkdir("results/bootstraps",0700);
			
			#pragma omp parallel for
			for(int bsNr=0;bsNr<n_boots;bsNr++)	{
				long int n_d_bs, n_r_bs;
				double *id_d_bs, *ra_d_bs, *dec_d_bs, *z_d_bs,*w_d_bs, *id_r_bs, *ra_r_bs, *dec_r_bs, *z_r_bs,*w_r_bs;
				char fileBSReal[100],fileBSRan[100];
				sprintf(fileBSReal,"bootstrap_data/bs%d_real_galaxies.txt",bsNr+1);
				FILE *fBSReal = fopen(fileBSReal, "r");
				if(fBSReal)	{
					n_d_bs = count_lines(fileBSReal);
					id_d_bs = alloc_1d(n_d_bs);
					ra_d_bs = alloc_1d(n_d_bs);
					dec_d_bs = alloc_1d(n_d_bs);
					z_d_bs = alloc_1d(n_d_bs);
					w_d_bs = alloc_1d(n_d_bs);
					
					char ch;							 	
					do	// ignoring the first header line
						ch = fgetc(fBSReal);			 
					while (ch != '\n');
					
					for(int i=0;i<n_d_bs;i++)	{
						fscanf(fBSReal,"%lf%lf%lf%lf",&id_d_bs[i],&ra_d_bs[i],&dec_d_bs[i],&z_d_bs[i]);
						w_d_bs[i]=1.0;
					}
					fclose(fBSReal);
				}
				else {
					fprintf(stderr,"\nBootstrap real galaxy %d not found!\n",bsNr+1);
					exit(1);
				}
				
				sprintf(fileBSRan,"bootstrap_data/bs%d_random_galaxies.txt",bsNr+1);
				FILE *fBSRan = fopen(fileBSRan, "r");
				if(fBSRan)	{
					n_r_bs = count_lines(fileBSRan);
					id_r_bs = alloc_1d(n_r_bs);
					ra_r_bs = alloc_1d(n_r_bs);
					dec_r_bs = alloc_1d(n_r_bs);
					z_r_bs = alloc_1d(n_r_bs);
					w_r_bs = alloc_1d(n_r_bs);
					
					char ch;							 	
					do	// ignoring the first header line
						ch = fgetc(fBSRan);			 
					while (ch != '\n');
					
					for(int i=0;i<n_r_bs;i++)
						fscanf(fBSRan,"%lf%lf%lf%lf%lf",&id_r_bs[i],&ra_r_bs[i],&dec_r_bs[i],&z_r_bs[i],&w_r_bs[i]);
					fclose(fBSRan);
				}
				else {
					fprintf(stderr,"\nBootstrap random galaxy %d not found!\n",bsNr+1);
					exit(1);
				}
				
				fprintf(f_summary,"\nBS %d: N_real: %ld, N_random: %ld\n",bsNr+1,n_d_bs,n_r_bs);
				
				computeBootstraps_meth2(bsNr,n_d_bs,ra_d_bs,dec_d_bs,z_d_bs, w_d_bs,n_r_bs, ra_r_bs,dec_r_bs,z_r_bs,w_r_bs);

				free_1d(id_d_bs);
				free_1d(ra_d_bs);
				free_1d(dec_d_bs);
				free_1d(z_d_bs);
				free_1d(w_d_bs);
				free_1d(id_r_bs);
				free_1d(ra_r_bs);
				free_1d(dec_r_bs);
				free_1d(z_r_bs);		
				free_1d(w_r_bs);
			}
		}
	}
	
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
	
	free_2d(RR,bin_rp);
	
	if(n_boots != 0)	fprintf(f_summary,"\nNumber of bootstrap samples for wp : %d\n",n_boots);
	if(n_jacks != 0)	fprintf(f_summary,"\nNumber of jackknife samples for wp : %d\n",n_jacks);
	fprintf(f_summary,"\nNumber of suffles for Mp : %d\n",n_shuffles);
	
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






