print 'Importing...'
import numpy as np
import pyfits
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from sklearn import mixture
import random
import cosmolopy.distance as cd
import cosmolopy.density 
import esutil
from scipy import spatial,integrate,optimize
from scipy.interpolate import interp1d,interp2d
from scipy.special import erf
from collections import Counter


###OUTPUT DIRECTORY###
dir='afterburner_outputs'
print 'Output Directory:',dir

#read in data
print 'Getting Data'
c=pyfits.getdata('redmapper_v6.4.11_full.fits')
g_1=pyfits.getdata('y1a1_gold_short0.fits',ignore_missing_end=True)
g_2=pyfits.getdata('y1a1_gold_desdmphotoz_short.fits',ignore_missing_end=True)
g_3=pyfits.getdata('y1a1_gold_desdmphotoz_Mr.fits')

rac=c.field('RA')
decc=c.field('DEC')
z=c.field('Z_LAMBDA')
rag1=g_1.field('RA')
decg1=g_1.field('DEC')
zg1=g_2.field('ZP')
ngals=c.field('LAMBDA_CHISQ')
AMAG=g_3.field('Mr')
modest_class=g_1.field('MODEST_CLASS')
mult_niter=g_1.field('MULT_NITER_MODEL')
flags_gold=g_1.field('FLAGS_GOLD')
magg=g_2.field('MAG_AUTO_G')
magr=g_2.field('MAG_AUTO_R')
magi=g_2.field('MAG_AUTO_I')
magz=g_2.field('MAG_AUTO_Z')

#make cuts
zmin=0.1
zmax=0.9

ra1=90
ra2=100
dec1=-50
dec2=-40

w, = np.where((z>zmin) & (z<zmax) & (rac>ra1) & (rac<ra2) & (decc>dec1) & (decc<dec2) & (ngals!=0))
# & (rac>ra1) & (rac<ra2) & (decc>dec1) & (decc<dec2)
c1=c[w]

zmin=0.05
zmax=1.1

ra1=89
ra2=101
dec1=-51
dec2=-39

#'crazy color' cut - don't use galaxies with colors less than -1 or greater than 4
crazy1=-1
crazy2=4
gr=magg-magr
ri=magr-magi
iz=magi-magz

w, = np.where((zg1>zmin) & (zg1<zmax) & (rag1>ra1) & (rag1<ra2) & (decg1>dec1) & (decg1<dec2) & (modest_class==1) & (mult_niter>0) & (flags_gold==0) & (AMAG<=-19.) & (gr>crazy1) & (gr<crazy2) & (ri>crazy1) & (ri<crazy2) & (iz>crazy1) & (iz<crazy2))
# & (rag1>ra1) & (rag1<ra2) & (decg1>dec1) & (decg1<dec2)
g1a=g_1[w]
g1b=g_2[w]
g1c=g_3[w]

print 'total clusters: ',len(c1)
print 'total galaxies: ',len(g1a)

rac=c1.field('RA')
decc=c1.field('DEC')
z=c1.field('Z_LAMBDA')
NGALS=c1.field('LAMBDA_CHISQ')
lambda_r=c1.field('R_LAMBDA')
cid=c1.field('MEM_MATCH_ID')
rag=g1a.field('RA')
decg=g1a.field('DEC')
zg=g1b.field('ZP')
zgerr=g1b.field('ZPE')
galid=g1a.field('COADD_OBJECTS_ID')


###########BACKGROUND GALAXY DENSITY CALCULATION###############
print 
print 'calculating background densities'

unique_galid,inds=np.unique(galid,return_index=True)
galid=list(galid)
zgtemp=zg[inds]
cosmo={'omega_M_0' : 0.23, 'omega_lambda_0' : 0.77, 'h' : 1.}
cosmo=cd.set_omega_k_0(cosmo)

#set redshift bins
zbin1=[]
zbin2=[]
zbin3=[]
zbin4=[]
zbin5=[]
zbin6=[]
zbin7=[]
zbin8=[]
zbin9=[]
zbin10=[]
zmin=min(zg)
zmax=max(zg)
step=(zmax-zmin)/10


zbins=[zbin1,zbin2,zbin3,zbin4,zbin5,zbin6,zbin7,zbin8,zbin9,zbin10]


area=(max(rag)-min(rag))*(max(decg)-min(decg))

#calculate background density in each bin
sum_gals=[]
volume=[]
globaldensity=[]
for i in range(0,len(zbins)):
    unique_id=unique_galid[np.where((zgtemp>zmin+(i*step))&(zgtemp<=zmin+(i+1)*step))]
    sumg=unique_id.size
    sum_gals.append(sumg)
    vol=cd.comoving_volume(zmin+((i+1)*step),**cosmo)-cd.comoving_volume(zmin+(i*step),**cosmo)
    vol=vol*(4*np.pi/41253)*area
    volume.append(vol)
    den=sumg/vol
    globaldensity.append(den)

'''
print 'plotting background densities'
vertlines=[zmin,zmin+step,zmin+2*step,zmin+3*step,zmin+4*step,zmin+5*step,zmin+6*step,zmin+7*step,zmin+8*step,zmin+9*step,zmax]
fig=plt.figure()

a1=fig.add_subplot(111)
#plt.plot(zcenter,rho_b,'r.',label='All Points')
#plt.plot(ave_z,ave_rhob,'b*',ms=10.,label='Binned Averages')
plt.plot(ave_z,globaldensity,'go',ms=10.,label='$\\Sigma (gals) / V$')
plt.xlabel('Redshift')
plt.ylabel('Background Galaxy Density ($gals/(\\frac{Mpc}{h})^3$)')
#plt.ylim(0.,0.3)
ymin,ymax=plt.ylim()
plt.ylim( (ymin,ymax) )
plt.vlines(vertlines,ymin,ymax,colors='k',linestyles='dashed',alpha=0.4,label='Bin Edges')
plt.legend(loc='upper right',prop={'size':10})

plt.tight_layout()
plt.savefig(dir+'/test_htm.png')

exit()'''

###########CALCULATING NGALS PER CLUSTER#############
print 
print 'matching galaxies within test radii'

#get cluster centers etc.
central_ra=c1.field('RA_CENT')
central_ra=[i[0] for i in central_ra]
central_ra=np.array(central_ra)
central_dec=c1.field('DEC_CENT')
central_dec=[i[0] for i in central_dec]
central_dec=np.array(central_dec)
central_z=c1.field('Z_LAMBDA')
central_z_err=c1.field('Z_LAMBDA_E')
cluster_id=c1.field('MEM_MATCH_ID')
cluster_ngals=c1.field('LAMBDA_CHISQ_CENT')
cluster_ngals=[i[0] for i in cluster_ngals]
cluster_ngals=np.array(cluster_ngals)
maskfrac=c1.field('MASKFRAC')
galid=g1a.field('COADD_OBJECTS_ID')

ang_diam_dist=cd.angular_diameter_distance(central_z,z0=0,**cosmo)

rmax=3.0     #maximum test radius in mpc. rmax will always be included as a test radius regardless of rmin,step
rmin=0.1   #minimum test radius in mpc. rmin always included as test radius
step=0.1  #step size, stepping from rmin to rmax, in mpc
radii=np.r_[rmin:rmax:step,rmax]

def gaussian(x,mu,sigma):
    return 1/(sigma*np.sqrt(2*np.pi))*np.exp(-(x-mu)**2/(2*sigma**2))

#old. used for exact integration, but thats hella slow
#def gaussint(xmin,xmax,mu,sigma):
#    return integrate.quad(gaussian,xmin,xmax,args=(mu,sigma))

totalgals=[]
backgroundgals=[]
ngals=[]
density=[]
background_dense=[]
vols=[]
zmatch=[]
depth=10
h=esutil.htm.HTM(depth)
for j in radii:
    degrees=(360/(2*np.pi))*(float(j)/ang_diam_dist) #convert radii to angular sizes
    #match cluster centers to surrounding galaxies
    m1,m2,adist=h.match(central_ra,central_dec,rag,decg,radius=degrees,maxmatch=50000)
    m1uniq=np.unique(m1)#find unique clusters (m1 has 1 cluster index for each matched galaxy)
    clusterid=cluster_id[m1uniq]
    truez=central_z[m1uniq]
    truezerr=central_z_err[m1uniq]
    zmin=truez-0.1
    zmax=truez+0.1
    total=[]
    totalerr=[]
    probz=[]
    for x in range(len(m1uniq)): #get the total number of matched galaxies for each cluster
        w0=np.where(m1==m1uniq[x]) #find points where gals matched to specific cluster
        w=m2[w0]
        membz=zg[w]
        membzerr=zgerr[w]
        intmin=zmin[x] #Integration window minimum
        intmax=zmax[x] #Integration window maximum
    #        proberr=[]
        zpts,zstep=np.linspace(intmin,intmax,20,retstep=True) #split redshift window for approximation
        area=[]
        for i in range(len(zpts)-1): #approximate integral using trapezoidal riemann sum
            zpt1=zpts[i] #zpt1,zpt2 are left,right points for trapezoid, respectively
            zpt2=zpts[i+1]
            gauss1=gaussian(zpt1,membz,membzerr) #gauss1/2 are gaussian values for left,right points, respectively
            gauss2=gaussian(zpt2,membz,membzerr)
            area1=((gauss1+gauss2)/2.)*zstep
            area.append(area1)
        area=np.array(area)
        arflip=area.swapaxes(0,1)
        prob=np.sum(arflip,axis=1)
        #    if x<20:
        #        plt.plot(membz,prob,'.')
        #        ymin,ymax=plt.ylim()
        #        plt.vlines(np.array([intmin,intmax]),ymin,ymax)
        #        plt.show()
        probz.append(prob)
        total.append(np.sum(prob))
    #        totalerr.append(np.sqrt(np.sum(proberr**2)))
    maskfrac1=maskfrac[m1uniq] #find fraction of the cluster that is masked
    mf=1-maskfrac1
    volume=(4./3.)*np.pi*(j**3)*(mf) #calculate total volume of cluster (assuming spherical with radius j)
    total_density=total/volume #calculate total galaxy density
    ztemp=central_z[m1uniq]
    #find background densities for each redshift bin (I know it looks messy, but its quick and it works. Deal with it.)
    background_density=np.where(ztemp<=zmin+step,globaldensity[0],
                                np.where((ztemp>zmin+step)&(ztemp<=zmin+2*step),globaldensity[1],
                                         np.where((ztemp>zmin+2*step)&(ztemp<=zmin+3*step),globaldensity[2],
                                                  np.where((ztemp>zmin+3*step)&(ztemp<=zmin+4*step),globaldensity[3],
                                                           np.where((ztemp>zmin+4*step)&(ztemp<=zmin+5*step),globaldensity[4],
                                                                    np.where((ztemp>zmin+5*step)&(ztemp<=zmin+6*step),globaldensity[5],
                                                                             np.where((ztemp>zmin+6*step)&(ztemp<=zmin+7*step),globaldensity[6],
                                                                                      np.where((ztemp>zmin+7*step)&(ztemp<=zmin+8*step),globaldensity[7],
                                                                                               np.where((ztemp>zmin+8*step)&(ztemp<=zmin+9*step),globaldensity[8],
                                                                                                        np.where(ztemp>zmin+9*step,globaldensity[9],-999))))))))))
    background=background_density*volume
    n=total-background #calculate background subtracted richness
    n_density=total_density-background_density
    #save it all
    totalgals.append(total)
    backgroundgals.append(background)
    background_dense.append(background_density)
    ngals.append(n)
    density.append(n_density)
    vols.append(volume)


print
print 'Making R200/M200 measurements'

#Msat redshift dependence
def logMsat(z,M0=12.33,a=-0.27):
    return M0 + a*z

#Mass-richness relation functions (see Tinker et al 2011)
def ncen(M,log_Mmin,sigma):
    #takes logM_min and logSigma_logM from paper. returns Ncentral from paper
    sigma=10**sigma
    return (1./2.)*(1+erf((np.log10(M)-log_Mmin)/sigma))

def ntot(M,Msat,log_Mmin,sigma,alpha_sat,Mcut):
    #takes logMmin, logSigma_logM, logMsat, logMcut from paper. Returns Ntotal=Ncentral+Nsatellite from paper
    #Msat=logMsat(z,M0,a)
    Msat=10**Msat
    Mcut=10**Mcut
    #Msat_inverse=1./Msat
    #MdivMsat=np.outer(Msat_inverse,M)
    return ncen(M,log_Mmin,sigma)*(1+((M/Msat**alpha_sat)*np.exp(-Mcut/M)))

def hod_mass(N,params):
    #params: logMmin,logMsat,alphasat,logMcut,logsigmalogM directly from table 4 of paper
    Mmin=params[0]
    #M0=params[1]
    #a=params[2]
    Msat=params[1]
    alpha=params[2]
    Mcut=params[3]
    sigma=params[4]
    mass=np.linspace(1*10**7,2*10**17,1*10**6)
    #ztest=np.linspace(0.05,1.1,1*10**4)
    m=interp1d(ntot(mass,Msat,Mmin,sigma,alpha,Mcut),mass)
    #xx,yy=np.meshgrid(mass,Msat)
    #n_tot=ntot(mass,ztest,M0,a,Mmin,sigma,alpha,Mcut)
    #n_a=[]
    #for i in range(len(Msat)):
    #    m=interp1d(ntot(mass,Msat[i],Mmin,sigma,alpha,Mcut),mass)
    #    ni=m(N)
    #    n_a.append(ni)
    #m=interp2d(n_tot,ztest,mass)
    n=m(N)
    #n_a=np.array(n_a)
    #n=np.diagonal(n_a)
    return n




'''
def bins(data,num):
    d=np.sort(data)
    inds=np.argsort(data)
    bins=np.array_split(d,num)
    bin_indices=np.array_split(inds,num)
    return bins,bin_indices

num=10
zbins,zinds=bins(z,num)
'''
bins=np.linspace(0.1,0.9,10)
bin_inds=np.digitize(z,bins)

ngals=np.array(ngals)
density=np.array(density)
volume=np.array(vols)

full_bin_inds=np.array([bin_inds,]*30)

param_0=[11.6,12.,0.85,12.25,-0.69] #parameters for mass conversion - see table 4 in Tinker paper
params=np.array([param_0,]*len(bins))
for i in range(1,len(bins)):
    zrange=np.array([bins[i-1],bins[i]])
    zmid=(np.median(zrange))
    msat=logMsat(zmid)
    params[i][1]=msat



#zmatch=np.array([z,]*30)
n_mass=np.concatenate(ngals)
#z_mass=np.concatenate(zmatch)
t=[]
for i in n_mass:  #set minimum value for mass conversions to prevent code from failing
    if i>=0.9:
        t.append(i)
    elif i<0.9:
        t.append(0.9)

n_mass=np.array(t)
n_mass.shape=ngals.shape
mass=np.zeros_like(n_mass)

for i in range(1,len(bins)):
    ntemp=n_mass[np.where(full_bin_inds==i)]
    mass_i=hod_mass(ntemp,params[i-1])
    mass[np.where(full_bin_inds==i)]=mass_i

#nlow=n_mass[np.where(zmatch<=0.35)]
#nhi=n_mass[np.where(zmatch>0.35)]
#mass=hod_mass(n_mass,params) #calculate mass given ngals (see above functions)
#mass[np.where(zmatch<=0.35)]=hod_mass(nlow,lowz_params)
#mass.shape=ngals.shape
mass_density=mass/volume
test=ngals[0]
#print test.size
ngals=ngals.swapaxes(0,1)
#lambda_ngals=lambda_ngals.swapaxes(0,1)
density=density.swapaxes(0,1)
mass_density=mass_density.swapaxes(0,1)
mass=mass.swapaxes(0,1)
background_dense=np.array(background_dense)
background_dense=background_dense.swapaxes(0,1)
subdir='/smaller_backsub'
subdir2='/sphere_mass'
#print 'subdirs:',subdir,subdir2
rho_crit,rho_0=cosmolopy.density.cosmo_densities(**cosmo) #current critical density
print 'Critical Density:',rho_crit
pc=200*np.ones_like(radii)

#change critical density with redshift
def crit_density(p_c,z,Omega_m,Omega_lambda):
    return p_c*(Omega_m*(1+z)**3 + Omega_lambda)#/(Omega_m*(1+z)**3)*Omega_m

R200_measure=[]
M200_measure=[]
redshift=[]
lambda_ngals=[]
R_lambda=[]
cidsave=[]
rasave=[]
decsave=[]
N_back=[]
X=200 #desired excess over critical density, ex. if X=200, calculates R/M200
dX=1 #acceptance window around X
interpradii=np.append([0],radii)
j=0
for i in range(0,test.size):
    cluster=ngals[i]
    dense=density[i]
    clustermass=mass[i]
    massdense=mass_density[i]
    x=clusterid[i]
    cidsave.append(x)
    c_ngals=NGALS[np.where(cid==x)]
    if c_ngals.size == 0:
        c_ngals=np.array([0])
    c_ngals=float(c_ngals)
    lambda_ngals.append(c_ngals)
    c_r=lambda_r[np.where(cid==x)]
    if c_r.size==0:
        c_r=np.array([0])
    c_r=float(c_r)
    R_lambda.append(c_r)
    c_z=truez[i]
    if c_z.size==0:
        c_z=o
    c_z=float(c_z)
    redshift.append(c_z)
    c_ra=rac[np.where(cid==x)]
    if c_ra.size==0:
        c_ra=np.array([0])
    c_ra=float(c_ra)
    rasave.append(c_ra)
    c_dec=decc[np.where(cid==x)]
    if c_dec.size==0:
        c_dec=np.array([0])
    c_dec=float(c_dec)
    decsave.append(c_dec)
#
    backdense=background_dense[i]
    backdense1=float(backdense[0]) #background density for this cluster
#
    critdense1=crit_density(rho_crit,c_z,0.23,0.77)
    critdense=critdense1*np.ones_like(radii)
    ratio=massdense/critdense
#
    f=interp1d(radii,ratio)
    radii_new=np.linspace(rmin,rmax,10000)
    ratio_new=f(radii_new)
    r200m=radii_new[np.where((ratio_new>=X-dX)&(ratio_new<=X+dX))] #find possible r200s within acceptance range
    if r200m.size > 0:
        r200m=np.mean(r200m) #mean of all possible r200s is measured r200
    else:
        r200m=0. #bogus r200=0 if nothing within acceptance range
        print 'bad cluster:',x
        j=j+1
#
    R200_measure.append(r200m)
#
    vol=(4*np.pi/3.)*(r200m**3) #cluster volume inside R200
    nback=backdense1*vol #expected number of background galaxies within R200
    N_back.append(nback)
    interpmass=np.append([0],clustermass)
    interp=interp1d(interpradii,interpmass)
    m200m=interp(r200m) #measured M200 (mass within R200)
    M200_measure.append(m200m)
'''#
    if c_ngals>=85:
        fig=plt.figure()
#
        ax=fig.add_axes([0.1,0.1,0.65,0.85])
        ln1=ax.plot(radii,clustermass,'b-',label='Mass')
        ax.set_title('Cluster ID:'+str(x))
        ax.set_xlabel('Radius (mpc)')
        ax.set_ylabel('Mass ($M_{\odot}$)')
        textfit='Cluster data: \n' \
            'True z= %.2f \n' \
            '$\\rho_{crit}(z) = $ %.2E \n' \
            'Meas. R200: %.1f \n' \
            'Meas. M200: %.2E \n' \
            %(c_z,critdense1,r200m,m200m)
        plt.figtext(0.81,0.65,textfit,fontsize='9')
#
        a2=ax.twinx()
        ln2=a2.plot(radii,ratio,'r-',label='$\\rho / \\rho_c$')
        ln3=a2.plot(radii,pc,'k--',label='200')
        a2.set_ylabel('$\\rho / \\rho_{crit}$')
        lns=ln1+ln2+ln3
        labels=[l.get_label() for l in lns]
        ax.legend(lns,labels,loc='center right')
        plt.show()
        plt.savefig(dir+'/cluster_'+str(x)+'.png')

#
    fig=plt.figure()
#
    ax=fig.add_axes([0.1,0.1,0.65,0.85])
    ln1=ax.plot(radii,clustermass,'b-',label='Mass')
    ax.set_title('Cluster ID:'+str(x))
    ax.set_xlabel('Radius (mpc)')
    ax.set_ylabel('Mass ($M_{\odot}$)')
    textfit='Cluster data: \n' \
        'True z= %.2f \n' \
        '$\\rho_{crit}(z) = $ %.2E \n' \
        'Meas. R200: %.1f \n' \
        'Meas. M200: %.2E \n' \
        %(c_z,critdense1,r200m,m200m)
    plt.figtext(0.81,0.65,textfit,fontsize='9')
#
    a2=ax.twinx()
    ln2=a2.plot(radii,ratio,'r-',label='$\\rho / \\rho_c$')
    ln3=a2.plot(radii,pc,'k--',label='200')
    a2.set_ylabel('$\\rho / \\rho_{crit}$')
    lns=ln1+ln2+ln3
    labels=[l.get_label() for l in lns]
    ax.legend(lns,labels,loc='center right')
    plt.savefig('largest_massprofile.png')
#
#    plt.savefig(dir+subdir2+'/cluster_'+str(i)+'.png')
'''
print 'Total bad clusters:',j


R200_measure=np.array(R200_measure)
M200_measure=np.array(M200_measure)
lambda_ngals=np.array(lambda_ngals)
R_lambda=np.array(R_lambda)
cidsave=np.array(cidsave)
redshift=np.array(redshift)
rasave=np.array(rasave)
decsave=np.array(decsave)
N_back=np.array(N_back)

col1=pyfits.Column(name='MATCH_MEMB_ID',format='J',array=cidsave)
col2=pyfits.Column(name='RA',format='D',array=rasave)
col3=pyfits.Column(name='DEC',format='D',array=decsave)
col4=pyfits.Column(name='Z',format='E',array=redshift)
col5=pyfits.Column(name='R200',format='E',array=R200_measure)
col6=pyfits.Column(name='M200',format='E',array=M200_measure)
col9=pyfits.Column(name="LAMBDA_CHISQ",format='E',array=lambda_ngals)
cols=pyfits.ColDefs([col1,col2,col3,col4,col5,col6,col9])
tbhdu=pyfits.new_table(cols)
tbhdu.writeto('hodz_test.fit',clobber=True)


'''
#test. ignore this
def lambmass(lamb,z,lamb0=20,logM0=14.06,Flamb=1.3,Gz=-2,z0=0):
    M0=10**logM0
    frac1=lamb/lamb0
    frac2=(1.+z)/(1.+z0)
    return M0*(frac1**Flamb)*(frac2**Gz)

lambdamass=lambmass(lambda_ngals,redshift)
massresid=(lambdamass-M200_measure)
absresid=np.abs(massresid)
minresid=min(absresid)
goodid=cidsave[np.where(absresid==minresid)]
print 'good cluster:',goodid

'''


######SELECTING MEMBER GALAXIES WITHIN R200######
print 'selecting member galaxies'
#degrees=(360/(2*np.pi))*(rmax/ang_diam_dist)
#matchind,adist,zdist=h.cylmatch(central_ra,central_dec,central_z,rag,decg,zg,radius=degrees,dz=0.025,maxmatch=50000)
#m1,m2,adist=h.match(central_ra,central_dec,rag,decg,radius=degrees,maxmatch=50000)


magg=g1b.field('MAG_AUTO_G')
magr=g1b.field('MAG_AUTO_R')
magi=g1b.field('MAG_AUTO_I')
magz=g1b.field('MAG_AUTO_Z')
amagr=g1c.field('Mr')
spread=g1a.field('SPREAD_MODEL_I')
mult=g1a.field('MULT_NITER_MODEL')
zgerr=g1b.field('ZPE')
kcorr=g1c.field('Kir')
gr0=g1c.field('gr0')

print 'calculating probabilities'

def sigma(R,R200,c=10):
    #Radial NFW profile implementation. Takes array of radii, value of R200,
    #and NFW concentration parameter (set to 10 by default)
    if R200>0:
        Rs=float(R200)/float(c)
        r=R/Rs
        try:
            pre=1./((r**2)-1)
            arctan_coeff=2./(np.sqrt(r**2-1))
            arctan_arg=np.sqrt((r-1)/(r+1))
            sigma=np.where(r>1,pre*(1-arctan_coeff*np.arctan(arctan_arg)),1./3.)
        except ZeroDivisionError:
            sigma=1./3.
        return sigma
    else:
        bogusval=-99.*np.ones_like(R)
        return bogusval


#integral method for finding constant of normalization
def integrand(R,sigma,R200,c=10):
    #sets up integrand to use to calculate probability
    return 2*np.pi*R*sigma(R,R200,c)

def prob(r,R200,sigma=sigma,c=10):
    #computes the probability of membership for each galaxy
    #takes array of radii, cluster R200 value, returns probability of membership
    #(max prob ~1/4)
    Rcore=float(R200)/float(c)
    dr=Rcore
    integral=integrate.quad(integrand,Rcore,R200,args=(sigma,R200,c))
    p=np.where(r>Rcore,(2*np.pi*r*dr*sigma(r,R200,c))/(integral[0]),(2*np.pi*Rcore*dr*sigma(Rcore,R200,c))/(integral[0]))
    return p


#member file stuff
galaxyID=[]
hostID=[]
P_radial=[]
angular_dist=[]
P_redshift=[]
hostR200=[]
hostM200=[]
hostN200=[]
galaxyRA=[]
galaxyDEC=[]
galaxyZ=[]
galZerr=[]
galmagG=[]
galmagR=[]
galmagI=[]
galmagZ=[]
galamagR=[]
galmult=[]
galspread=[]
galkcorr=[]
galgr0=[]
#cluster file stuff
R200=[]
M200=[]
N200=[]
N200p=[]
N_background=[]
cluster_ID=[]
cluster_RA=[]
cluster_DEC=[]
cluster_Z=[]
for i in range(0,len(m1uniq)):
    #first do all the cluster stuff
#investigate your indexing here. Try using ind=m1uniq[i], thing=thing0[ind]
#it may not change anything, but it should at least make the code more robust
    cind=clusterid[i]
    cluster_ID.append(cind)
    hostr=R200_measure[i]
    R200.append(hostr)
    hostm=M200_measure[i]
    M200.append(hostm)
    hostback=N_back[i]
    N_background.append(hostback)
    cluster_RA.append(central_ra[i])
    cluster_DEC.append(central_dec[i])
    cluster_Z.append(central_z[i])
    #now do all the galaxy stuff
    match0=np.where(m1==m1uniq[i])
    match=m2[match0]
    add=ang_diam_dist[i]
    angular=adist[match0]
    R=(2.*np.pi/360.)*angular*add
    p_rad=prob(R,hostr)
    #add in redshift probabilities
    p_z=probz[i]
    pdist=p_rad*p_z
    n_p=np.sum(pdist)
#    n=len(match)
#    N200.append(n)
    N200p.append(n_p)
    P_radial.append(p_rad)
    angular_dist.append(angular)
    P_redshift.append(p_z)
    glxid=galid[match]
    glxra=rag[match]
    glxdec=decg[match]
    glxz=zg[match]
    glxmagg=magg[match]
    glxmagr=magr[match]
    glxmagi=magi[match]
    glxmagz=magz[match]
    glxamagr=amagr[match]
    glxspread=spread[match]
    glxmult=mult[match]
    glxzerr=zgerr[match]
    glxkcorr=kcorr[match]
    glxgr0=gr0[match]
    galaxyID.append(glxid)
    tempID=cind*np.ones_like(glxid)
    hostID.append(tempID)
    tempR=hostr*np.ones_like(glxid)
    hostR200.append(tempR)
    tempM=hostm*np.ones_like(glxid)
    hostM200.append(tempM)
    tempn=n_p*np.ones_like(glxid)
    hostN200.append(tempn)
    galaxyRA.append(glxra)
    galaxyDEC.append(glxdec)
    galaxyZ.append(glxz)
    galZerr.append(glxzerr)
    galmagG.append(glxmagg)
    galmagR.append(glxmagr)
    galmagI.append(glxmagi)
    galmagZ.append(glxmagz)
    galamagR.append(glxamagr)
    galmult.append(glxmult)
    galspread.append(glxspread)
    galkcorr.append(glxkcorr)
    galgr0.append(glxgr0)

galaxyID=np.array([item for sublist in galaxyID for item in sublist])
galaxyRA=np.array([item for sublist in galaxyRA for item in sublist])
galaxyDEC=np.array([item for sublist in galaxyDEC for item in sublist])
galaxyZ=np.array([item for sublist in galaxyZ for item in sublist])
hostID=np.array([item for sublist in hostID for item in sublist])
P_radial=np.array([item for sublist in P_radial for item in sublist])
P_redshift=np.array([item for sublist in P_redshift for item in sublist])
angular_dist=np.array([item for sublist in angular_dist for item in sublist])
hostR200=np.array([item for sublist in hostR200 for item in sublist])
hostM200=np.array([item for sublist in hostM200 for item in sublist])
hostN200=np.array([item for sublist in hostN200 for item in sublist])
galZerr=np.array([item for sublist in galZerr for item in sublist])
galmagG=np.array([item for sublist in galmagG for item in sublist])
galmagR=np.array([item for sublist in galmagR for item in sublist])
galmagI=np.array([item for sublist in galmagI for item in sublist])
galmagZ=np.array([item for sublist in galmagZ for item in sublist])
galamagR=np.array([item for sublist in galamagR for item in sublist])
galmult=np.array([item for sublist in galmult for item in sublist])
galspread=np.array([item for sublist in galspread for item in sublist])
galkcorr=np.array([item for sublist in galkcorr for item in sublist])
galgr0=np.array([item for sublist in galgr0 for item in sublist])

mxp=max(P_radial)
P_radial=P_radial*(1./mxp)
Pdist=P_radial*P_redshift

#print 'Good matched galaxies:'
#print galaxyID[np.where(hostID==goodid)]
#exit()

print 
print 'calculating GMM probabilities'

def linear(p,x):
    return p[0]+p[1]*x


#pull expected color vs redshift data
annis=np.loadtxt('red_galaxy_El1_COSMOS_DES_filters.txt')
jimz=[i[0] for i in annis]
jimgr=[i[2] for  i in annis]
jimri=[i[3] for i in annis]
jimiz=[i[4] for i in annis]

jimgr=np.array(jimgr)+0.2
jimri=np.array(jimri)+0.10

treedata=zip(jimz,np.zeros_like(jimz))
tree=spatial.KDTree(treedata)

            #background=np.where(alpha==min(alpha))
            #redmu=np.delete(mu,background)
            #red=np.where(mu==max(redmu))
            #blue=np.where((mu!=max(redmu))&(sigma!=max(sigma)))

def gmmfit(band1,band2,expcol,distlim=0.05,n_components=3,tol=0.0000001,galaxyID=galaxyID,
           hostID=hostID,P_radial=P_radial,P_redshift=P_redshift):
#for g-r, band1=galmagG, band2=galmagR
    Pdist=P_radial*P_redshift
    gmm=mixture.GMM(n_components=n_components,tol=tol,n_iter=500)
    slope=[]
    yint=[]
    mu_r=[]
    mu_b=[]
    mu_bg=[]
    sigma_r=[]
    sigma_b=[]
    sigma_bg=[]
    alpha_r=[]
    alpha_b=[]
    alpha_bg=[]
    Pred=[]
    Pblue=[]
    Pbg=[]
    probgalid=[]
    converged=[]
    for x in cluster_ID:
        glxid=galaxyID[np.where((hostID==x)&(Pdist>=distlim))]
        magg1=band1[np.where((hostID==x)&(Pdist>=distlim))]
        magr1=band2[np.where((hostID==x)&(Pdist>=distlim))]
        rprob1=P_radial[np.where((hostID==x)&(Pdist>=distlim))]
        zprob1=P_redshift[np.where((hostID==x)&(Pdist>=distlim))]
        distprob=rprob1*zprob1
        gr1=magg1-magr1
#        find expected RS color
        zcl=cluster_Z[np.where(cluster_ID==x)]
        qpt=zip(zcl,np.zeros_like(zcl))
        d,ind=tree.query(qpt)
        expred=expcol[ind]
        if len(glxid) > 1:
            gr1.shape=(len(gr1),1)
            distprob.shape=(len(distprob),1)
            fit=gmm.fit(gr1,data_weights=distprob)
            conv=gmm.converged_
            converged.append(conv)
            mu=gmm.means_
            mu.shape=(len(mu),)
            alpha=gmm.weights_
            alpha.shape=(len(alpha),)
            covars=gmm.covars_
            sigma=np.sqrt(covars)
            sigma.shape=(len(sigma),)
            mudif=np.abs(mu-expred)
            red=np.where(mudif==min(mudif))
            blualph=np.delete(alpha,red)
            blue=np.where(alpha==max(blualph))
            background=np.where(alpha==min(blualph))
            alphar=alpha[red]
            alphabg=alpha[background]
            if alphar<0.1: #arbitrarily set, maybe play with this at some point
                if alphabg>alphar:
                    temp=red
                    red=background
                    background=temp
            mur=mu[red]
            mub=mu[blue]
            if mur<mub:
                temp=red
                red=blue
                blue=temp
            mur=mu[red]
            mub=mu[blue]
            mubg=mu[background]
            sigr=sigma[red]
            sigb=sigma[blue]
            sigbg=sigma[background]
            alphar=alpha[red]
            alphab=alpha[blue]
            alphabg=alpha[background]
#
            mu_r.append(mur)
            mu_b.append(mub)
            mu_bg.append(mubg)
            sigma_r.append(sigr)
            sigma_b.append(sigb)
            sigma_bg.append(sigbg)
            alpha_r.append(alphar)
            alpha_b.append(alphab)
            alpha_bg.append(alphabg)
            exr=-((gr1-mur)**2)/(2*(sigr**2))
            exb=-((gr1-mub)**2)/(2*(sigb**2))
            exbg=-((gr1-mubg)**2)/(2*(sigbg**2))
            p_red=(1/(sigr*np.sqrt(2*np.pi)))*np.exp(exr)
            p_blue=(1/(sigb*np.sqrt(2*np.pi)))*np.exp(exb)
            p_bg=(1/(sigbg*np.sqrt(2*np.pi)))*np.exp(exbg)
            maxLred=(1/(sigr*np.sqrt(2*np.pi)))
            maxLblue=(1/(sigb*np.sqrt(2*np.pi)))
            Pred.append(p_red/maxLred)
            Pblue.append(p_blue/maxLblue)
            probgalid.append(glxid)
            magr1.shape=(len(magr1),)
            gr1.shape=(len(gr1),)
            weights=distprob*p_red#*(1-p_blue)
            weights.shape=(len(weights),)
            rfit,info=np.polynomial.polynomial.polyfit(magr1,gr1,deg=1,w=weights,full=True)
            rfit.shape=(2,)
            slope.append(rfit[1])
            yint.append(rfit[0])
#            bf=alphab/(alphar+alphab)
#            if zcl>0.3 and zcl<=0.4 and mur<1.3:
#                fig=plt.figure()
#                fig.add_subplot(111)
#                plt.hist(gr1,facecolor='g',alpha=0.3,normed=True,bins=45,range=(-1,4))
#                plt.plot(y,mlab.normpdf(y,mur,sigr)*alphar,'r-')
#                plt.plot(y,mlab.normpdf(y,mub,sigb)*alphab,'b-')
#                plt.plot(y,mlab.normpdf(y,mubg,sigbg)*alphabg,'k-')
#                plt.title('z='+str(zcl))
#                plt.xlabel('r-i')
#                plt.show()
#                plt.close(fig)
#            chisq=info[0]
#            dof=len(redr)
#            redchisq=chisq/dof
    slope=np.array((slope))
    yint=np.array((yint))
    mu_r=np.array((mu_r))
    mu_b=np.array((mu_b))
    mu_bg=np.array((mu_bg))
    sigma_r=np.array((sigma_r))
    sigma_b=np.array((sigma_b))
    sigma_bg=np.array((sigma_bg))
    alpha_r=np.array((alpha_r))
    alpha_b=np.array((alpha_b))
    alpha_bg=np.array((alpha_bg))
    Pred=np.array(Pred)
    Pblue=np.array(Pblue)
    Pbg=np.array(Pbg)
    probgalid=np.array(probgalid)
    converged=np.array(converged)
    return slope,yint,mu_r,mu_b,mu_bg,sigma_r,sigma_b,sigma_bg,alpha_r,alpha_b,alpha_bg,Pred,Pblue,Pbg,probgalid,converged
 
cluster_Z=np.array(cluster_Z)

grinfo=gmmfit(galmagG,galmagR,jimgr)

riinfo=gmmfit(galmagR,galmagI,jimri)

izinfo=gmmfit(galmagI,galmagZ,jimiz)


grslope=grinfo[0]
gryint=grinfo[1]
grmu_r=grinfo[2]
grmu_b=grinfo[3]
grmu_bg=grinfo[4]
grsigma_r=grinfo[5]
grsigma_b=grinfo[6]
grsigma_bg=grinfo[7]
gralpha_r=grinfo[8]
gralpha_b=grinfo[9]
gralpha_bg=grinfo[10]
grPred=grinfo[11]
grPblue=grinfo[12]
grPbg=grinfo[13]
grprobgalid=grinfo[14]
grconverged=grinfo[15]

rislope=riinfo[0]
riyint=riinfo[1]
rimu_r=riinfo[2]
rimu_b=riinfo[3]
rimu_bg=riinfo[4]
risigma_r=riinfo[5]
risigma_b=riinfo[6]
risigma_bg=riinfo[7]
rialpha_r=riinfo[8]
rialpha_b=riinfo[9]
rialpha_bg=riinfo[10]
riPred=riinfo[11]
riPblue=riinfo[12]
riPbg=riinfo[13]
riprobgalid=riinfo[14]
riconverged=riinfo[15]

izslope=izinfo[0]
izyint=izinfo[1]
izmu_r=izinfo[2]
izmu_b=izinfo[3]
izmu_bg=izinfo[4]
izsigma_r=izinfo[5]
izsigma_b=izinfo[6]
izsigma_bg=izinfo[7]
izalpha_r=izinfo[8]
izalpha_b=izinfo[9]
izalpha_bg=izinfo[10]
izPred=izinfo[11]
izPblue=izinfo[12]
izPbg=izinfo[13]
izprobgalid=izinfo[14]
izconverged=izinfo[15]


#three gaussian testing site
'''

            if mur<=1.35*mub:
                fig=plt.figure()
                fig.add_subplot(111)
                plt.hist(gr1,facecolor='g',alpha=0.3,normed=True,bins=45,range=(-1,4))
                plt.plot(y,mlab.normpdf(y,mur,sigr)*alphar,'r-')
                plt.plot(y,mlab.normpdf(y,mub,sigb)*alphab,'b-')
                plt.plot(y,mlab.normpdf(y,mubg,sigbg)*alphabg,'k-')
                plt.title('z='+str(f))
                plt.xlabel('r-i')
                plt.savefig('ri_odd/hist/cluster_'+str(x)+'.png')
                plt.close(fig)
#
                line=np.linspace(15,23,1000)
                fig=plt.figure()
                fig.add_subplot(111)
                plt.hexbin(magr1,gr1,C=weights,cmap='YlOrRd')
                plt.plot(line,linear(rfit,line),'g-')
                plt.xlabel('r-mag')
                plt.ylabel('$g-r$')
                plt.title('z='+str(f))
                plt.savefig('ri_odd/cmag/cluster_'+str(x)+'.png')
                plt.close(fig)
                


gmm=mixture.GMM(n_components=3,tol=0.0000001,n_iter=500)

z2=[]
sz2=[]
slope=[] #add uncertainties to all of these soon
yint=[]
reduced_chisq=[]
mu_r=[]
mu_b=[]
mu_bg=[]
sigma_r=[]
sigma_b=[]
sigma_bg=[]
alpha_r=[]
alpha_b=[]
alpha_bg=[]
blue_frac=[]
Pred=[]
Pblue=[]
probgalid=[]
y=np.arange(-1,4.0,0.0001)
z_temp=cluster_Z
z_temp=np.array(z_temp)
bla=[]
distlim=0.05
cluster_ID=np.array(cluster_ID)
gm=[]
rm=[]
dataweight=[]
converged=[]
for x in misfit:#[np.where(z_temp<=0.4)]:
    glxid=galaxyID[np.where((hostID==x)&(Pdist>=distlim))]
    magg1=galmagR[np.where((hostID==x)&(Pdist>=distlim))]
    gm.append(magg1)
    magr1=galmagI[np.where((hostID==x)&(Pdist>=distlim))]
    rm.append(magr1)
    rprob1=P_radial[np.where((hostID==x)&(Pdist>=distlim))]
    zprob1=P_redshift[np.where((hostID==x)&(Pdist>=distlim))]
    distprob=rprob1*zprob1
    gr1=magg1-magr1
    if len(glxid) > 1:
        gr1.shape=(len(gr1),1)
        distprob.shape=(len(distprob),1)
        fit=gmm.fit(gr1,data_weights=distprob)
        conv=gmm.converged_
        converged.append(conv)
        mu=gmm.means_
        mu.shape=(len(mu),)
        alpha=gmm.weights_
        alpha.shape=(len(alpha),)
        covars=gmm.covars_
        sigma=np.sqrt(covars)
        sigma.shape=(len(sigma),)
        background=np.where(alpha==min(alpha))
        redmu=np.delete(mu,background)
        red=np.where(mu==max(redmu))
        blue=np.where((mu!=max(redmu))&(sigma!=max(sigma)))
        mur=mu[red]
        mu_r.append(mu[red])
        mub=mu[blue]
        mu_b.append(mu[blue])
        mubg=mu[background]
        mu_bg.append(mu[background])
        sigr=sigma[red]
        sigma_r.append(sigma[red])
        sigb=sigma[blue]
        sigma_b.append(sigma[blue])
        sigbg=sigma[background]
        sigma_bg.append(sigma[background])
        alphar=alpha[red]
        alpha_r.append(alpha[red])
        alphab=alpha[blue]
        alpha_b.append(alpha[blue])
        alphabg=alpha[background]
        alpha_bg.append(alpha[background])
        #index=cluster_ID.index(x)
        index=np.where(cluster_ID==x)
        f=z_temp[index]
        z2.append(f)
        if alphabg>alphar:
            fig=plt.figure()
            fig.add_subplot(111)
            plt.hist(gr1,facecolor='g',alpha=0.3,normed=True,bins=45,range=(-1,4))
            plt.plot(y,mlab.normpdf(y,mur,sigr)*alphar,'r-')
            plt.plot(y,mlab.normpdf(y,mub,sigb)*alphab,'b-')
            plt.plot(y,mlab.normpdf(y,mubg,sigbg)*alphabg,'k-')
            plt.title('z='+str(f))
            plt.xlabel('g-r')
            plt.savefig('ri_odd/hist/cluster_'+str(x)+'.png')
            plt.close(fig)

        if f<0.15:
            plt.savefig('zbins_gmm_clusterplots/1_15_hist/cluster_'+str(x))
        elif f>=0.15 and f<0.2:
            plt.savefig('zbins_gmm_clusterplots/15_2_hist/cluster_'+str(x))
        elif f>=0.2 and f<0.25:
            plt.savefig('zbins_gmm_clusterplots/2_25_hist/cluster_'+str(x))
        elif f>=0.25 and f<0.3:
            plt.savefig('zbins_gmm_clusterplots/25_3_hist/cluster_'+str(x))
        elif f>=0.3 and f<0.35:
            plt.savefig('zbins_gmm_clusterplots/3_35_hist/cluster_'+str(x))
        elif f>=0.35 and f<0.4:
            plt.savefig('zbins_gmm_clusterplots/35_4_hist/cluster_'+str(x))
        else:
            plt.savefig('zbins_gmm_clusterplots/4_end_hist/cluster_'+str(x))
        plt.close(fig)
#
        redgr=[]
        redr=[]
        bluegr=[]
        bluer=[]
        exr=-((gr1-mur)**2)/(2*(sigr**2))
        exb=-((gr1-mub)**2)/(2*(sigb**2))
        exbg=-((gr1-mubg)**2)/(2*(sigbg**2))
        p_red=(1/(sigr*np.sqrt(2*np.pi)))*np.exp(exr)
        p_blue=(1/(sigb*np.sqrt(2*np.pi)))*np.exp(exb)
        p_bg=(1/(sigbg*np.sqrt(2*np.pi)))*np.exp(exbg)
        probgalid.append(glxid)
        magr1.shape=gr1.shape
        #w=np.where((p_red>p_blue)&(gr1>mub))
        #notw=np.where((p_red<=p_blue)|(p_red>p_blue)&(gr1<=mub))
        #redgr=gr1[w]
        #bluegr=gr1[notw]
        #redr=magr1[w]
        #bluer=magr1[notw]
        #if len(redgr)>1:
        Pred.append(p_red)
        Pblue.append(p_blue)
        magr1.shape=(len(magr1),)
        gr1.shape=(len(gr1),)
        weights=distprob*p_red#*(1-p_blue)*(1-p_bg)
        weights.shape=(len(weights),)
        dataweight.append(weights)
        rfit,info=np.polynomial.polynomial.polyfit(magr1,gr1,deg=1,w=weights,full=True)
        rfit.shape=(2,)
        slope.append(rfit[1])
        yint.append(rfit[0])
        chisq=info[0]
        dof=len(redr)
        redchisq=chisq/dof
        reduced_chisq.append(redchisq)
        sz2.append(f)
#
        #total=float(len(redgr)+len(bluegr))
        #blue_frac.append(len(bluegr)/total)
        if alphabg>alphar:
            line=np.linspace(15,23,1000)
            fig=plt.figure()
            fig.add_subplot(111)
            #plt.plot(redr,redgr,'r.',bluer,bluegr,'b.',line,rplot(line),'g-')
            plt.hexbin(magr1,gr1,C=weights,cmap='YlOrRd')
            plt.plot(line,linear(rfit,line),'g-')
            plt.xlabel('r-mag')
            plt.ylabel('$g-r$')
            plt.title('z='+str(f))
            plt.savefig('ri_odd/cmag/cluster_'+str(x)+'.png')
            plt.close(fig)

        if f<0.15:
            plt.savefig('zbins_gmm_clusterplots/1_15_cmag/cluster_'+str(x))
        elif f>=0.15 and f<0.2:
            plt.savefig('zbins_gmm_clusterplots/15_2_cmag/cluster_'+str(x))
        elif f>=0.2 and f<0.25:
            plt.savefig('zbins_gmm_clusterplots/2_25_cmag/cluster_'+str(x))
        elif f>=0.25 and f<0.3:
            plt.savefig('zbins_gmm_clusterplots/25_3_cmag/cluster_'+str(x))
        elif f>=0.3 and f<0.35:
            plt.savefig('zbins_gmm_clusterplots/3_35_cmag/cluster_'+str(x))
        elif f>=0.35 and f<0.4:
            plt.savefig('zbins_gmm_clusterplots/35_4_cmag/cluster_'+str(x))
        else:
            plt.savefig('zbins_gmm_clusterplots/4_end_cmag/cluster_'+str(x))
        plt.close(fig)

    elif len(glxid)==1:
        print '1 richness cluster:',x
        Pred.append(-99)
        Pblue.append(-99)
        probgalid.append(glxid)
        z2.append(-99)
        sz2.append(-99)
        slope.append(-99)
        yint.append(-99)
        mu_r.append(-99)
        mu_b.append(-99)
        sigma_r.append(-99)
        sigma_b.append(-99)
        #blue_frac.append(-99)
        bla.append(-99)
    else:
        print '0 richness cluster:',x
        z2.append(-99)
        sz2.append(-99)
        slope.append(-99)
        yint.append(-99)
        mu_r.append(-99)
        mu_b.append(-99)
        sigma_r.append(-99)
        sigma_b.append(-99)
        blue_frac.append(-99)
        Pred.append(-99)
        Pblue.append(-99)
        #blue_frac.append(-99)
        bla.append(-99)


z2=np.array(z2)
sz2=np.array(sz2)
slope=np.array(slope)
converged=np.array(converged)
yint=np.array(yint)
mu_r=np.array(mu_r)
mu_b=np.array(mu_b)
sigma_r=np.array(sigma_r)
sigma_b=np.array(sigma_b)
alpha_r=np.array(alpha_r)
alpha_b=np.array(alpha_b)
blue_frac=alpha_b/alpha_r
Pred=np.array(Pred)
Pred=np.array([item for sublist in Pred for item in sublist])
Pblue=np.array(Pblue)
Pblue=np.array([item for sublist in Pblue for item in sublist])
probgalid=np.array(probgalid) #tentatively same IDs as existing galaxyID. If things look weird later, look here to fix it
probgalid=np.array([item for sublist in probgalid for item in sublist])
Pred.shape=probgalid.shape
Pblue.shape=probgalid.shape
#gm=np.array(gm)
#rm=np.array(rm)
#dataweight=np.array(dataweight)
'''


#check for one galaxy belonging to multiple clusters
w=np.where(Pdist>0.05)
galid2=galaxyID[w]
hostid2=hostID[w]
galra2=galaxyRA[w]
galdec2=galaxyDEC[w]
galz2=galaxyZ[w]
galzerr2=galZerr[w]
galmagG2=galmagG[w]
galmagR2=galmagR[w]
galmagI2=galmagI[w]
galmagZ2=galmagZ[w]
prad2=P_radial[w]
pz2=P_redshift[w]

#izgalid2=np.array([item for sublist in izprobgalid for item in sublist])

gralpha_r.shape=(gralpha_r.size,)
gralpha_b.shape=(gralpha_b.size,)
grPcolor=(gralpha_r*grPred) + (gralpha_b*grPblue)
grpc2=np.array([item for sublist in grPcolor for item in sublist])
grpc2.shape=(grpc2.size,)
grPred=np.array([item for sublist in grPred for item in sublist])
grPblue=np.array([item for sublist in grPblue for item in sublist])
grPbg=np.array([item for sublist in grPbg for item in sublist])

rialpha_r.shape=(rialpha_r.size,)
rialpha_b.shape=(rialpha_b.size,)
riPcolor=(rialpha_r*riPred)+(rialpha_b*riPblue)
ripc2=np.array([item for sublist in riPcolor for item in sublist])
ripc2.shape=(ripc2.size,)
riPred=np.array([item for sublist in riPred for item in sublist])
riPblue=np.array([item for sublist in riPblue for item in sublist])
riPbg=np.array([item for sublist in riPbg for item in sublist])

izalpha_r.shape=(izalpha_r.size,)
izalpha_b.shape=(izalpha_b.size,)
izPcolor=(izalpha_r*izPred)+(izalpha_b*izPblue)
izpc2=np.array([item for sublist in izPcolor for item in sublist])
izpc2.shape=(izpc2.size,)
izPred=np.array([item for sublist in izPred for item in sublist])
izPblue=np.array([item for sublist in izPblue for item in sublist])
izPbg=np.array([item for sublist in izPbg for item in sublist])


duplicates=[item for item, count in Counter(galid2).iteritems() if count>1]
'''
#this needs fixing
#if len(duplicates)!=0: #add stuff to this mess to include Pred/Pblue
new_hostid=[]
new_Pmemb=[]
#    new_r=[]
#    new_m=[]
#    new_n=[]
new_grPr=[]
new_grPb=[]
indices=[]
for i in duplicates:
    ind=np.where(galid2==i)
    indices.append(ind)
    new_hostid.append(hostid2[ind])
    new_Pmemb.append(prad2[ind])
        #new_r.append(hostR200[ind])
        #new_m.append(hostM200[ind])
        #new_n.append(hostN200[ind])
    new_grPr.append(grpred2[ind])
    new_grPb.append(grpblue2[ind])
#
ind1=[]
ind2=[]
for i in indices:
    ind1.append(i[0][0])
    u=i[0][1:]
    u=u.tolist()
    ind2.append(u)
#
ind2=[item for sublist in ind2 for item in sublist]
#
#    lengths=[len(x) for x in new_hostid]
#    for i in range(0,len(new_hostid)):
#        if len(new_hostid[i])<max(lengths):
#            new_Pmemb[i]=np.append(new_Pmemb[i],np.zeros(max(lengths)-len(new_hostid[i])))
#            new_r[i]=np.append(new_r[i],np.zeros(max(lengths)-len(new_hostid[i])))
#            new_m[i]=np.append(new_m[i],np.zeros(max(lengths)-len(new_hostid[i])))
#            new_n[i]=np.append(new_n[i],np.zeros(max(lengths)-len(new_hostid[i])))
#            new_hostid[i]=np.append(new_hostid[i],np.zeros(max(lengths)-len(new_hostid[i])))
#            new_Pr[i]=np.append(new_Pr[i],np.zeros(max(lengths)-len(new_hostid[i])))
#            new_Pb[i]=np.append(new_Pb[i],np.zeros(max(lengths)-len(new_hostid[i])))
#
    hostID=hostID.tolist()
    P_radial=P_radial.tolist()
    hostR200=hostR200.tolist()
    hostM200=hostM200.tolist()
    hostN200=hostN200.tolist()
    Pred=Pred.tolist()
    Pblue=Pblue.tolist()
    for i in range(0,len(hostID)):
        hostID[i]=np.append(hostID[i],np.zeros(max(lengths)-1))
        P_radial[i]=np.append(P_radial[i],np.zeros(max(lengths)-1))
        hostR200[i]=np.append(hostR200[i],np.zeros(max(lengths)-1))
        hostM200[i]=np.append(hostM200[i],np.zeros(max(lengths)-1))
        hostN200[i]=np.append(hostN200[i],np.zeros(max(lengths)-1))
        Pred[i]=np.append(Pred[i],np.zeros(max(lengths)-1))
        Pblue[i]=np.append(Pblue[i],np.zeros(max(lengths)-1))
#        
    hostID=np.array(hostID)
    new_hostid=np.array(new_hostid)
    P_radial=np.array(P_radial)
    new_Pmemb=np.array(new_Pmemb)
    hostR200=np.array(hostR200)
    new_r=np.array(new_r)
    hostM200=np.array(hostM200)
    new_m=np.array(new_m)
    hostN200=np.array(hostN200)
    new_n=np.array(new_n)
    Pred=np.array(Pred)
    new_Pr=np.array(new_Pr)
    Pblue=np.array(Pblue)
    new_Pb=np.array(new_Pb)
#
    for i in range(0,len(new_hostid)):
        hostID[ind1[i]]=new_hostid[i]
        P_radial[ind1[i]]=new_Pmemb[i]
        hostR200[ind1[i]]=new_r[i]
        hostM200[ind1[i]]=new_m[i]
        hostN200[ind1[i]]=new_n[i]
        Pred[ind1[i]]=new_Pr[i]
        Pblue[ind1[i]]=new_Pb[i]
#        
    galaxyID=np.delete(galaxyID,ind2,0)
    hostID=np.delete(hostID,ind2,0)
    P_radial=np.delete(P_radial,ind2,0)
    angular_dist=np.delete(angular_dist,ind2,0)
#redshift_dist=np.delete(redshift_dist,ind2,0)
    hostR200=np.delete(hostR200,ind2,0)
    hostM200=np.delete(hostM200,ind2,0)
    hostN200=np.delete(hostN200,ind2,0)
    galaxyRA=np.delete(galaxyRA,ind2,0)
    galaxyDEC=np.delete(galaxyDEC,ind2,0)
    galaxyZ=np.delete(galaxyZ,ind2,0)
    galZerr=np.delete(galZerr,ind2,0)
    galmagG=np.delete(galmagG,ind2,0)
    galmagR=np.delete(galmagR,ind2,0)
    galmagI=np.delete(galmagI,ind2,0)
    galmagZ=np.delete(galmagZ,ind2,0)
    galamagR=np.delete(galamagR,ind2,0)
    galmult=np.delete(galmult,ind2,0)
    galspread=np.delete(galspread,ind2,0)
    galkcorr=np.delete(galkcorr,ind2,0)
    galgr0=np.delete(galgr0,ind2,0)
    Pred=np.delete(Pred,ind2,0)
    Pblue=np.delete(Pblue,ind2,0)
#
    l=hostID.shape
    l=l[1]
    #The column format piece needs some work, right now its hard coded to only accept
    #galaxies matched to two halos
    col2=pyfits.Column(name='HOST_HALOID',format='2J',array=hostID)
    col3=pyfits.Column(name='P_RADIAL',format='2E',array=P_radial)
    col17=pyfits.Column(name='HOST_R200',format='2E',array=hostR200)
    col18=pyfits.Column(name='HOST_M200',format='2E',array=hostM200)
    col19=pyfits.Column(name='HOST_N200',format='2E',array=hostN200)
    col21=pyfits.Column(name='P_RED',format='2E',array=Pred)
    col22=pyfits.Column(name='P_BLUE',format='2E',array=Pblue)
else:
    col2=pyfits.Column(name='HOST_HALOID',format='J',array=hostID)
    col3=pyfits.Column(name='P_RADIAL',format='E',array=P_radial)
    col17=pyfits.Column(name='HOST_R200',format='E',array=hostR200)
    col18=pyfits.Column(name='HOST_M200',format='E',array=hostM200)
    col19=pyfits.Column(name='HOST_N200',format='E',array=hostN200)
    col21=pyfits.Column(name='P_RED',format='E',array=Pred)
    col22=pyfits.Column(name='P_BLUE',format='E',array=Pblue)
'''

#col2=pyfits.Column(name='HOST_HALOID',format='J',array=hostID)
#col3=pyfits.Column(name='P_RADIAL',format='E',array=P_radial)
#col17=pyfits.Column(name='HOST_R200',format='E',array=hostR200)
#col18=pyfits.Column(name='HOST_M200',format='E',array=hostM200)
#col19=pyfits.Column(name='HOST_N200',format='E',array=hostN200)
#col21=pyfits.Column(name='GRP_RED',format='E',array=grPred)
#col22=pyfits.Column(name='GRP_BLUE',format='E',array=grPblue)
#col23=pyfits.Column(name='GRP_BACKGROUND',format='E',array=grPbg)
#col24=pyfits.Column(name='P_REDSHIFT',format='E',array=P_redshift)
#col25=pyfits.Column(name='RIP_RED',format='E',array=riPred)
#col26=pyfits.Column(name='RIP_BLUE',format='E',array=riPblue)
#col27=pyfits.Column(name='RIP_BG',format='E',array=riPbg)
#col28=pyfits.Column(name='IZP_RED',format='E',array=izPred)
#col29=pyfits.Column(name='IZP_BLUE',format='E',array=izPblue)
#col30=pyfits.Column(name='IZP_BG',format='E',array=izPbg)


######WRITING OUTPUT FILES######
grPmemb=grpc2*prad2*pz2
riPmemb=ripc2*prad2*pz2
izPmemb=izpc2*prad2*pz2

print 'writing member file'

col1=pyfits.Column(name='COADD_OBJECTS_ID',format='K',array=galid2)
col2=pyfits.Column(name='HOST_HALOID',format='J',array=hostid2)
col3=pyfits.Column(name='RA',format='D',array=galra2)
col4=pyfits.Column(name='DEC',format='D',array=galdec2)
col5=pyfits.Column(name='ZP',format='E',array=galz2)
col6=pyfits.Column(name='ZPE',format='E',array=galzerr2)
col7=pyfits.Column(name='MAG_AUTO_G',format='E',array=galmagG2)
col8=pyfits.Column(name='MAG_AUTO_R',format='E',array=galmagR2)
col9=pyfits.Column(name='MAG_AUTO_I',format='E',array=galmagI2)
col10=pyfits.Column(name='MAG_AUTO_Z',format='E',array=galmagZ2)
col11=pyfits.Column(name='P_RADIAL',format='E',array=prad2)
col12=pyfits.Column(name='P_REDSHIFT',format='E',array=pz2)
col13=pyfits.Column(name='GR_P_COLOR',format='E',array=grpc2)
col14=pyfits.Column(name='RI_P_COLOR',format='E',array=ripc2)
col15=pyfits.Column(name='IZ_P_COLOR',format='E',array=izpc2)
col16=pyfits.Column(name='GR_P_MEMBER',format='E',array=grPmemb)
col17=pyfits.Column(name='RI_P_MEMBER',format='E',array=riPmemb)
col18=pyfits.Column(name='IZ_P_MEMBER',format='E',array=izPmemb)
col19=pyfits.Column(name='AMAG_R',format='E',array=galamagR)
col20=pyfits.Column(name='DIST_TO_CENTER',format='E',array=angular_dist)
col21=pyfits.Column(name='GRP_RED',format='E',array=grPred)
col22=pyfits.Column(name='GRP_BLUE',format='E',array=grPblue)
col23=pyfits.Column(name='GRP_BG',format='E',array=grPbg)
col24=pyfits.Column(name='RIP_RED',format='E',array=riPred)
col25=pyfits.Column(name='RIP_BLUE',format='E',array=riPblue)
col26=pyfits.Column(name='RIP_BG',format='E',array=riPbg)
col27=pyfits.Column(name='IZP_RED',format='E',array=izPred)
col28=pyfits.Column(name='IZP_BLUE',format='E',array=izPblue)
col29=pyfits.Column(name='IZP_BG',format='E',array=izPbg)


cols=pyfits.ColDefs([col1,col2,col3,col4,col5,col6,col7,col8,col9,col10,col11,col12,col13,col14,col15,col16,col17,col18,col19,col20,col21,col22,col23,col24,col25,col26,col27,col28,col29])
tbhdu=pyfits.BinTableHDU.from_columns(cols)
tbhdu.writeto('gold_redmapper_v6.4.11_full_newmembers_v2.fit',clobber=True)

print 'writing cluster file'

R200=np.array(R200)
M200=np.array(M200)
N200=np.array(N200)
N_background=np.array(N_background)
cluster_ID=np.array(cluster_ID)
cluster_RA=np.array(cluster_RA)
cluster_DEC=np.array(cluster_DEC)
cluster_Z=np.array(cluster_Z)

grsep=(grmu_r-grmu_b)/grsigma_r
risep=(rimu_r-rimu_b)/risigma_r
izsep=(izmu_r-izmu_b)/izsigma_r

grw=np.where(grsep>=1.2)
riw=np.where(risep>=1.2)
izw=np.where(izsep>=1.2)

grflag=np.zeros_like(grmu_r)
grflag[grw]=1
riflag=np.zeros_like(rimu_r)
riflag[riw]=1
izflag=np.zeros_like(izmu_r)
izflag[izw]=1


col1=pyfits.Column(name='MEM_MATCH_ID',format='J',array=cluster_ID)
col2=pyfits.Column(name='RA',format='D',array=cluster_RA)
col3=pyfits.Column(name='DEC',format='D',array=cluster_DEC)
col4=pyfits.Column(name='Z',format='E',array=cluster_Z)
col5=pyfits.Column(name='R200',format='E',array=R200)
col6=pyfits.Column(name='M200',format='E',array=M200)
col7=pyfits.Column(name='N200',format='E',array=N200)
col8=pyfits.Column(name='Nback',format='E',array=N_background)
col9=pyfits.Column(name="LAMBDA_CHISQ",format='E',array=lambda_ngals)
col10=pyfits.Column(name="LAMBDA_R",format='E',array=R_lambda)
col11=pyfits.Column(name='GR_SLOPE',format='E',array=grslope)
col12=pyfits.Column(name='GR_INTERCEPT',format='E',array=gryint)
col13=pyfits.Column(name='GRMU_R',format='E',array=grmu_r)
col14=pyfits.Column(name='GRMU_B',format='E',array=grmu_b)
col15=pyfits.Column(name='GRSIGMA_R',format='E',array=grsigma_r)
col16=pyfits.Column(name='GRSIGMA_B',format='E',array=grsigma_b)
col17=pyfits.Column(name='GRW_R',format='E',array=gralpha_r)
col18=pyfits.Column(name='GRW_B',format='E',array=gralpha_b)
col19=pyfits.Column(name='RI_SLOPE',format='E',array=rislope)
col20=pyfits.Column(name='RI_INTERCEPT',format='E',array=riyint)
col21=pyfits.Column(name='RIMU_R',format='E',array=rimu_r)
col22=pyfits.Column(name='RIMU_B',format='E',array=rimu_b)
col23=pyfits.Column(name='RISIGMA_R',format='E',array=risigma_r)
col24=pyfits.Column(name='RISIGMA_B',format='E',array=risigma_b)
col25=pyfits.Column(name='RIW_R',format='E',array=rialpha_r)
col26=pyfits.Column(name='RIW_B',format='E',array=rialpha_b)
col27=pyfits.Column(name='GRMU_BG',format='E',array=grmu_bg)
col28=pyfits.Column(name='GRSIGMA_BG',format='E',array=grsigma_bg)
col29=pyfits.Column(name='GRW_BG',format='E',array=gralpha_bg)
col30=pyfits.Column(name='GR_CONVERGED',format='E',array=grconverged)
col31=pyfits.Column(name='RIMU_BG',format='E',array=rimu_bg)
col32=pyfits.Column(name='RISIGMA_BG',format='E',array=risigma_bg)
col33=pyfits.Column(name='RIW_BG',format='E',array=rialpha_bg)
col34=pyfits.Column(name='RI_CONVERGED',format='E',array=riconverged)
col35=pyfits.Column(name='IZ_SLOPE',format='E',array=izslope)
col36=pyfits.Column(name='IZ_INTERCEPT',format='E',array=izyint)
col37=pyfits.Column(name='IZMU_R',format='E',array=izmu_r)
col38=pyfits.Column(name='IZMU_B',format='E',array=izmu_b)
col39=pyfits.Column(name='IZSIGMA_R',format='E',array=izsigma_r)
col40=pyfits.Column(name='IZSIGMA_B',format='E',array=izsigma_b)
col41=pyfits.Column(name='IZW_R',format='E',array=izalpha_r)
col42=pyfits.Column(name='IZW_B',format='E',array=izalpha_b)
col43=pyfits.Column(name='IZMU_BG',format='E',array=izmu_bg)
col44=pyfits.Column(name='IZSIGMA_BG',format='E',array=izsigma_bg)
col45=pyfits.Column(name='IZW_BG',format='E',array=izalpha_bg)
col46=pyfits.Column(name='IZ_CONVERGED',format='E',array=izconverged)
col47=pyfits.Column(name='GR_SEP_FLAG',format='L',array=grflag)
col48=pyfits.Column(name='RI_SEP_FLAG',format='L',array=riflag)
col49=pyfits.Column(name='IZ_SEP_FLAG',format='L',array=izflag)


cols=pyfits.ColDefs([col1,col2,col3,col4,col5,col6,col7,col8,col9,col10,col11,col12,col13,col14,col15,col16,col17,col18,col19,col20,col21,col22,col23,col24,col25,col26,col27,col28,col29,col30,col31,col32,col33,col34,col35,col36,col37,col38,col39,col40,col41,col42,col43,col44,col45,col46,col47,col48,col49])
tbhdu=pyfits.BinTableHDU.from_columns(cols)
tbhdu.writeto('gold_redmapper_v6.4.11_full_clusterR200_v2.fit',clobber=True)



print 'success!!'
exit()

w=np.where(rialpha_r<0.1)
w=w[0]

misfit=cluster_ID[w]
zmis=cluster_Z[w]
mrmis=rimu_r[w]
mbmis=rimu_b[w]
mbgmis=rimu_bg[w]
srmis=risigma_r[w]
sbmis=risigma_b[w]
sbgmis=risigma_bg[w]
armis=rialpha_r[w]
abmis=rialpha_b[w]
abgmis=rialpha_bg[w]

y=np.linspace(-1,4,1000)
for i in range(len(misfit)):
    rimis=ri[np.where(hostID==misfit[i])]
    f=zmis[i]
    mr=mrmis[i]
    mb=mbmis[i]
    mbg=mbgmis[i]
    sr=srmis[i]
    sb=sbmis[i]
    sbg=sbgmis[i]
    ar=armis[i]
    ab=abmis[i]
    abg=abgmis[i]
    fig=plt.figure()
    fig.add_subplot(111)
    plt.hist(rimis,facecolor='g',alpha=0.3,normed=True,bins=45,range=(-1,4))
    plt.plot(y,mlab.normpdf(y,mr,sr)*ar,'r-')
    plt.plot(y,mlab.normpdf(y,mb,sb)*ab,'b-')
    plt.plot(y,mlab.normpdf(y,mbg,sbg)*abg,'k-')
    plt.xlabel('r-i')
    plt.title(str(f))
    plt.savefig('ri_odd/cluster_'+str(misfit[i])+'.png')
    plt.close(fig)
