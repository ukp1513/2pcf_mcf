# LAST UPDATED : 24 JUNE 2022

#just checking

import os
import numpy as np
from numpy.linalg import inv
from scipy.special import gamma
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import bokeh.palettes as bp
from matplotlib.ticker import ScalarFormatter
from matplotlib.ticker import FormatStrFormatter
import seaborn as sns
from matplotlib.patches import Ellipse
import scipy.optimize as opt

def wp_model(rp,r0,g):
    return rp*pow(r0/rp,g)*gamma(0.5)*gamma((g-1)*0.5)/gamma(g*0.5)
    
def covmat_to_corrmat(covariance):
    v = np.sqrt(np.diag(covariance))
    outer_v = np.outer(v, v)
    correlation = covariance / outer_v
    correlation[covariance == 0] = 0
    return correlation

# do you want to implement the SVD method?    
to_filter=1	# 1 : yes , 0: no

#max: 150 for dc2, 70 for gama
rpmin_wp,rpmax_wp=0.1,40.0
rpmin_mp,rpmax_mp=0.1,40.0
    
# CHECKING BOOTSTRAP/JACKKNIFE

cwd_fullpath = os.getcwd()
cwd_split=cwd_fullpath.split('/')
sample_name=cwd_split[-1].upper()

print("\nFitting sample %s...." %(sample_name))

is_bs=0
is_jk=0

if(os.path.isdir("results/bootstraps")):
	is_bs=1
if(os.path.isdir("results/jackknifes")):
	is_jk=1
	
# COUNTING THE NUMBER OF MCFs IN RESULT FOLDER

fProp = open("propertylist.txt","r")
lines = fProp.readlines()
mark_count=len(lines)
mark_name=[]
for line in lines: 
	mark_name.append(line.strip())
	
print('Nr. of marks: ',mark_count)
print('Marks: ',mark_name)

# CREATING wpRealAll.txt for JK/BS copies

rp=np.loadtxt('results/wpReal.txt')[:,0]
wpReal=np.loadtxt('results/wpReal.txt')[:,1]

total_nbins = len(rp)
ncopies = len(os.listdir('results/jackknifes'))

print('Nr. of total bins: ', total_nbins)
print('Nr. of JK/BS copies: ', ncopies)

wpRealAll_tofile = np.ndarray(shape=(ncopies+2,total_nbins), dtype=float)

wpRealAll_tofile[:][0]=rp
wpRealAll_tofile[:][1]=wpReal

for copy in range(ncopies):
	wpJK=np.loadtxt('results/jackknifes/wpMpJackknife_jk%d.txt' %(copy+1))[:,1]
	wpRealAll_tofile[:][copy+2]=wpJK

np.savetxt("results/wpRealAll.txt",np.transpose(wpRealAll_tofile),delimiter="\t",fmt='%f')

# CREATING mpRealAll.txt for JK/BS copies 

for mark_nr in range(mark_count):

	mpRealAll_tofile = np.ndarray(shape=(ncopies+2,total_nbins), dtype=float)

	mpRealAll_tofile[:][0]=np.loadtxt('results/mpRealAllShuffle_mk%d.txt' %(mark_nr+1))[:,0]
	mpRealAll_tofile[:][1]=np.loadtxt('results/mpRealAllShuffle_mk%d.txt' %(mark_nr+1))[:,1]

	for copy in range(ncopies):
		mpJK=np.loadtxt('results/jackknifes/wpMpJackknife_jk%d.txt' %(copy+1))[:,mark_nr+2]
		mpRealAll_tofile[:][copy+2]=mpJK

	np.savetxt("results/mpRealAll_mk%d.txt" %(mark_nr+1),np.transpose(mpRealAll_tofile),delimiter="\t",fmt='%f')
	
# FILTERING NAN AND INF VALUES

wpRealAll = np.loadtxt('results/wpRealAll.txt')

nrows,ncols = wpRealAll.shape

total_nbins = nrows
ncopies = ncols-2

nrows=total_nbins
ncols=ncopies

filter_index_wp = []
filter_index_mp = []
for i in range(0,nrows):
	if(np.isnan(wpRealAll[i,1]).any() == True or np.isinf(wpRealAll[i,1]).any() == True or wpRealAll[i,1]<0.):
		filter_index_wp.append(i)
		filter_index_mp.append(i)
	else: 
		for j in range(ncols):
			if(np.isnan(wpRealAll[i,j]).any() == True or np.isinf(wpRealAll[i,j]).any() == True):
				filter_index_wp.append(i)
	if(wpRealAll[i,0] < rpmin_wp or wpRealAll[i,0] > rpmax_wp):
		filter_index_wp.append(i)
	if(wpRealAll[i,0] < rpmin_mp or wpRealAll[i,0] > rpmax_mp):
		filter_index_mp.append(i)

filter_index_wp=list(set(filter_index_wp))	
filter_index_mp=list(set(filter_index_mp))			



#REMOVING NAN BINS FROM WP FILE
wpRealAll=np.delete(wpRealAll, filter_index_wp, axis=0)

nbins_wp=total_nbins-len(filter_index_wp)
nbins_mp=total_nbins-len(filter_index_mp)

np.savetxt('results/wpRealAll_filtered.txt',wpRealAll,delimiter='\t',fmt='%f')
np.savetxt('results/wpRealAll_filtered_tofit.txt',wpRealAll,delimiter='\t',fmt='%f')

#REMOVING THOSE BINS FROM MP FILES ALSO..
for i in range(mark_count):
	mpRealAll = np.loadtxt('results/mpRealAll_mk%d.txt' %(i+1))
	mpRealAllShuffle = np.loadtxt('results/mpRealAllShuffle_mk%d.txt' %(i+1))
		
	mpRealAll = np.delete(mpRealAll, filter_index_mp, axis=0)
	mpRealAllShuffle = np.delete(mpRealAllShuffle, filter_index_mp, axis=0)
	np.savetxt('results/mpRealAll_mk%d_filtered.txt' %(i+1),mpRealAll,delimiter='\t',fmt='%f')
	np.savetxt('results/mpRealAllShuffle_mk%d_filtered.txt' %(i+1),mpRealAllShuffle,delimiter='\t',fmt='%f')


	
# COMPUTING COVARIANCE MATRIX 

if(ncopies > 0):

	allCopiesWps = wpRealAll[:,2:ncopies+2]

	if(is_jk==1):
		cov_mat=np.cov(allCopiesWps, bias=True)	# C = (Njk-1)/Njk x SUM
		cov_mat = (ncopies-1)*cov_mat

	if(is_bs==1):
		cov_mat=np.cov(allCopiesWps, bias=False)	# C = (Njk-1)/Njk x SUM
		
	corr_mat = covmat_to_corrmat(cov_mat)

	# FITTING USING SVD 

	U, Dvector, UT = np.linalg.svd(corr_mat)	# C = U D UT

	Dinv_vec = []

	neff=nbins_wp

	for i in Dvector:
		if(to_filter==1):
			if(i < np.sqrt(2./ncopies)):
				neff-=1
				Dinv_vec.append(0.0)
			else:
				Dinv_vec.append(1./i)	
		else:
				Dinv_vec.append(1./i)
				
	Dinv = np.diag(Dinv_vec)
	Cinv = np.matmul(U,np.matmul(Dinv,UT))
	
	fEff=open("biproducts/effective_bins.txt","w")
	fEff.write(str(neff))
	fEff.close()
	
	np.savetxt("biproducts/Cinv_SVD.txt",np.transpose(Cinv),delimiter="\t",fmt='%f')
