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
from scipy.interpolate import interp1d
from scipy.special import erf
from collections import Counter


###FILE INFO###
#inputs
cluster_indir='/data/des30.a/data/bwelch/xray_clusters/'
gal_indir='/data/des30.a/data/bwelch/redmapper_y1a1/'
color_indir='/data/des30.a/data/bwelch/redmapper_y1a1/'
clusterfile=cluster_indir+'Chandra_Y1A1-6.4.16-Feb-26-2017-prelim.fits'
galfile=gal_indir+'y1a1_gold_BPZ_mof_gals_cut_allcolumns.fits'
colorfile=color_indir+'red_galaxy_El1_COSMOS_DES_filters.txt'
#outputs
outdir='/data/des30.a/data/bwelch/xray_clusters'
cluster_outfile=outdir+'bpz_mof_chandra_prelim_clusters.fit'
member_outfile=outdir+'bpz_mof_chandra_prelim_members.fit'
print 'Output Directory:',outdir

#read in data
print 'Getting Data'
c=pyfits.open(clusterfile)
c=c[1].data
g=pyfits.open(galfile)
g=g[1].data
#g_1=pyfits.getdata('y1a1_gold_short0.fits',ignore_missing_end=True)
#g_2=pyfits.getdata('y1a1_gold_bpz_short.fits',ignore_missing_end=True)
#desdmphotoz->bpz for new photoz's
#g_3=pyfits.getdata('y1a1_gold_bpz_Mr.fits')

rac=c.field('r500_ra')
decc=c.field('r500_dec')
z=c.field('redshift')
#ngals=c.field('LAMBDA_CHISQ')
rag1=g.field('RA')
decg1=g.field('DEC')
zg1=g.field('MEDIAN_Z')
AMAG=g.field('Mr')
#modest_class=g.field('MODEST_CLASS')
mult_niter=g.field('MULT_NITER_MODEL')
flags_gold=g.field('FLAGS_GOLD')
magg=g.field('MAG_AUTO_G')
magr=g.field('MAG_AUTO_R')
magi=g.field('MAG_AUTO_I')
magz=g.field('MAG_AUTO_Z')

#make cuts
zmin=0.1
zmax=1.0

ra1=0
ra2=360
dec1=-60
dec2=3

w, = np.where((z>zmin) & (z<zmax) & (rac>ra1) & (rac<ra2) & (decc>dec1) & (decc<dec2) )#& (ngals!=0))
c1=c[w]

zmin=0.05
zmax=1.1

ra1=0
ra2=361
dec1=-61
dec2=4

#'crazy color' cut - don't use galaxies with colors less than -1 or greater than 4
crazy1=-1
crazy2=4
gr=magg-magr
ri=magr-magi
iz=magi-magz

w, = np.where((zg1>zmin) & (zg1<zmax) & (rag1>ra1) & (rag1<ra2) & (decg1>dec1) & (decg1<dec2) & (mult_niter>0) & (flags_gold==0) & (AMAG<=-19.) & (gr>crazy1) & (gr<crazy2) & (ri>crazy1) & (ri<crazy2) & (iz>crazy1) & (iz<crazy2))
#& (modest_class==1)
g1=g[w]
#g1a=g_1[w]
#g1b=g_2[w]
#g1c=g_3[w]

print 'total clusters: ',len(c1)
print 'total galaxies: ',len(g1)

rac=c1.field('r500_ra')
decc=c1.field('r500_dec')
z=c1.field('redshift')
#NGALS=c1.field('LAMBDA_CHISQ')
#lambda_r=c1.field('R_LAMBDA')
#maskfrac=c1.field('MASKFRAC')
cid=c1.field('MEM_MATCH_ID')
rag=g1.field('RA')
decg=g1.field('DEC')
zg=g1.field('MEDIAN_Z')
zgerr=g1.field('Z_SIGMA')
galid=g1.field('COADD_OBJECTS_ID')

magg=g1.field('MAG_AUTO_G')
magr=g1.field('MAG_AUTO_R')
magi=g1.field('MAG_AUTO_I')
magz=g1.field('MAG_AUTO_Z')
amagr=g1.field('Mr')
kcorr=g1.field('Kir')
gr0=g1.field('gr0')



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


###########CALCULATING NGALS PER CLUSTER#############
print 
print 'matching galaxies within test radii'

#get cluster centers etc.
central_ra=rac
central_dec=decc
central_z=z
cluster_id=cid
galid=np.array(galid)

ang_diam_dist=cd.angular_diameter_distance(central_z,z0=0,**cosmo)

rmax=3.0     #maximum test radius in mpc. rmax will always be included as a test radius regardless of rmin,step
rmin=0.1   #minimum test radius in mpc. rmin always included as test radius
step=0.1  #step size, stepping from rmin to rmax, in mpc
radii=np.r_[rmin:rmax:step,rmax]

def gaussian(x,mu,sigma):
    return 1/(sigma*np.sqrt(2*np.pi))*np.exp(-(x-mu)**2/(2*sigma**2))


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
    if len(m1uniq)!=len(central_ra):
        diff=np.diff(m1uniq)
        ind=np.where(diff!=1)
        ind=np.array(ind[0])
        for i in ind:
            p=0
            while p<diff[i]-1:
                addval=m1uniq[i+p]+1
                insert_ind=i+p+1
                m1uniq=np.insert(m1uniq,insert_ind,addval)
                p=p+1
    if len(m1uniq)!=len(central_ra):
        if m1uniq[0]!=0:
            addval=np.arange(0,m1uniq[0],1)
            m1uniq=np.insert(m1uniq,0,addval)
        if m1uniq[-1]!=len(central_ra)-1:
            addval=np.arange(m1uniq[-1],len(central_ra)-1,1)
            m1uniq=np.insert(m1uniq,len(m1uniq),addval)
    clusterid=cluster_id[m1uniq]
    truez=central_z[m1uniq]
    zmin=truez-(0.025*(1+truez)**2)
    zmax=truez+(0.025*(1+truez)**2)
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
        probz.append(prob)
        total.append(np.sum(prob))
    if 'maskfrac' in globals():
        maskfrac1=maskfrac[m1uniq] #find fraction of the cluster that is masked
    else:
        maskfrac1=np.zeros_like(m1uniq) #set zero maskfraction if not provided
    mf=1-maskfrac1
    volume=(4./3.)*np.pi*(j**3)*(mf) #calculate total volume of cluster (assuming spherical with radius j)
    total_density=total/volume #calculate total galaxy density
    ztemp=central_z[m1uniq]
    #find background densities for each redshift bin (I know its messy, but its quick and it works.)
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
    Msat=10**Msat
    Mcut=10**Mcut
    return ncen(M,log_Mmin,sigma)*(1+(((M/Msat)**alpha_sat)*np.exp(-Mcut/M)))

def hod_mass(N,params):
    #params: logMmin,logMsat,alphasat,logMcut,logsigmalogM directly from table 4 of paper
    Mmin=params[0]
    Msat=params[1]
    alpha=params[2]
    Mcut=params[3]
    sigma=params[4]
    mass=np.linspace(1*10**7,2*10**17,1*10**6)
    m=interp1d(ntot(mass,Msat,Mmin,sigma,alpha,Mcut),mass)
    n=m(N)
    return n




####
#commented sections here used for redshift-varying HOD parameter
#We aren't actually using this, but I'll leave it for posterity
#For redshift-varying HOD parameters, un-comment lines below thru for loop
#comment out mass=hod_mass() line
####

#bins=np.linspace(0.1,0.9,10)
#bin_inds=np.digitize(z,bins)

ngals=np.array(ngals)
density=np.array(density)
volume=np.array(vols)

#full_bin_inds=np.array([bin_inds,]*30)

params=[11.6,12.45,1.0,12.25,-0.69]#parameters for mass conversion - see table 4 in Tinker paper
#param_0=[11.6,12.45,0.85,12.25,-0.69] 
#params=np.array([param_0,]*len(bins))
#for i in range(1,len(bins)):
#    zrange=np.array([bins[i-1],bins[i]])
#    zmid=(np.median(zrange))
#    msat=logMsat(zmid)
#    params[i][1]=msat



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
#n_mass.shape=ngals.shape
#mass=np.zeros_like(n_mass)

#for i in range(1,len(bins)):
#    ntemp=n_mass[np.where(full_bin_inds==i)]
#    mass_i=hod_mass(ntemp,params[i-1])
#    mass[np.where(full_bin_inds==i)]=mass_i

mass=hod_mass(n_mass,params) #calculate mass given ngals (see above functions)
#comment out above line for redshift-varying HOD
mass.shape=ngals.shape
mass_density=mass/volume
test=ngals[0]
ngals=ngals.swapaxes(0,1)
density=density.swapaxes(0,1)
mass_density=mass_density.swapaxes(0,1)
mass=mass.swapaxes(0,1)
background_dense=np.array(background_dense)
background_dense=background_dense.swapaxes(0,1)
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
    if 'NGALS' in globals():
        c_ngals=NGALS[np.where(cid==x)]
        if c_ngals.size == 0:
            c_ngals=np.array([0])
        c_ngals=float(c_ngals)
        lambda_ngals.append(c_ngals)
    if 'lambda_r' in globals():
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


print 'Total bad clusters:',j


R200_measure=np.array(R200_measure)
M200_measure=np.array(M200_measure)
if 'lambda_ngals' in globals():
    lambda_ngals=np.array(lambda_ngals)
if 'R_lambda' in globals():
    R_lambda=np.array(R_lambda)
cidsave=np.array(cidsave)
redshift=np.array(redshift)
rasave=np.array(rasave)
decsave=np.array(decsave)
N_back=np.array(N_back)



######SELECTING MEMBER GALAXIES WITHIN R200######
print 'calculating radial probabilities'

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
    angular=adist[match0] #angular distance to cluster center, degrees
    R=(2.*np.pi/360.)*angular*add #distance to cluster center in Mpc
    p_rad=prob(R,hostr)
    #add in redshift probabilities
    p_z=probz[i]
    pdist=p_rad*p_z
    n_p=np.sum(pdist)
    N200p.append(n_p)
    P_radial.append(p_rad)
    angular_dist.append(R) #save distance to cluster center in Mpc
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
#    glxspread=spread[match]
#    glxmult=mult[match]
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
#    galmult.append(glxmult)
#    galspread.append(glxspread)
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
#galmult=np.array([item for sublist in galmult for item in sublist])
#galspread=np.array([item for sublist in galspread for item in sublist])
galkcorr=np.array([item for sublist in galkcorr for item in sublist])
galgr0=np.array([item for sublist in galgr0 for item in sublist])

mxp=max(P_radial)
P_radial=P_radial*(1./mxp)
Pdist=P_radial*P_redshift


print 
print 'calculating GMM probabilities'

def linear(p,x):
    return p[0]+p[1]*x


#pull expected color vs redshift data
annis=np.loadtxt(colorfile)
jimz=[i[0] for i in annis]
jimgr=[i[2] for  i in annis]
jimri=[i[3] for i in annis]
jimiz=[i[4] for i in annis]

jimgr=np.array(jimgr)+0.2
jimri=np.array(jimri)+0.10

treedata=zip(jimz,np.zeros_like(jimz))
tree=spatial.KDTree(treedata)


def clip(data,data2):
    #does sigma clipping step
    m=data.mean();s=data.std();ix=(data>(m-2*s))&(data<(m+2*s))
    dat1=data[ix]
    dat1.shape=(len(dat1),1)
    dat2=data2[ix]
    dat2.shape=(len(dat2),1)
    return dat1,dat2

def pick_colors(mu,sigma,alpha,expred):
    #choose which gaussian fit corresponds to red/blue/background
    mudif=np.abs(mu-expred)
    red=np.where(mudif==min(mudif))
    blualph=np.delete(alpha,red)
    blue=np.where(alpha==max(blualph))
    background=np.where(alpha==min(blualph))
    alphar=alpha[red]
    alphabg=alpha[background]
    if alphar<0.1:
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
    return mur,mub,mubg,sigr,sigb,sigbg,alphar,alphab,alphabg

def sigma_clip(col,distprob,expred):
    #compute GMM fit for data after sigma clipping
    gmm=mixture.GMM(n_components=3,tol=0.0000001,n_iter=500)
    gr3,distprob3=clip(col,distprob)
    gmm.fit(gr3,data_weights=distprob3)
    mu=gmm.means_
    mu.shape=(len(mu),)
    alpha=gmm.weights_
    alpha.shape=(len(alpha),)
    covars=gmm.covars_
    sigma=np.sqrt(covars)
    sigma.shape=(len(sigma),)
    mur,mub,mubg,sigr,sigb,sigbg,alphar,alphab,alphabg=pick_colors(mu,sigma,alpha,expred)
    return gr3,distprob3,mur,mub,mubg,sigr,sigb,sigbg,alphar,alphab,alphabg




def gmmfit(band1,band2,expcol,distlim=0.05,n_components=3,tol=0.0000001,galaxyID=galaxyID,hostID=hostID,P_radial=P_radial,P_redshift=P_redshift):
    #Full GMM calculation for single observed color
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
    p0=0
    p1=0
    p2=0
    p3=0
    p4=0
    p5=0
    q0=0
    q1=0
    q2=0
    q3=0
    q4=0
    q5=0
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
        if len(glxid) >= 3:
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
            mur,mub,mubg,sigr,sigb,sigbg,alphar,alphab,alphabg=pick_colors(mu,sigma,alpha,expred)
            if np.abs(mur-mub)>=1.*sigr:
                p0=p0+1
            if np.abs(mur-mub)<1.*sigr:
                gr2,distprob2,mur,mub,mubg,sigr,sigb,sigbg,alphar,alphab,alphabg=sigma_clip(gr1,distprob,expred)
                p1=p1+1
                if np.abs(mur-mub)<1.*sigr:
                    gr3,distprob3,mur,mub,mubg,sigr,sigb,sigbg,alphar,alphab,alphabg=sigma_clip(gr2,distprob2,expred)
                    p2=p2+1
                    if np.abs(mur-mub)<1.*sigr:
                        gr4,distprob4,mur,mub,mubg,sigr,sigb,sigbg,alphar,alphab,alphabg=sigma_clip(gr3,distprob3,expred)
                        p3=p3+1
                        if np.abs(mur-mub)<1.*sigr:
                            gr5,distprob5,mur,mub,mubg,sigr,sigb,sigbg,alphar,alphab,alphabg=sigma_clip(gr4,distprob4,expred)
                            p4=p4+1
                            if np.abs(mur-mub)<1.*sigr:
                                gr6,distprob6,mur,mub,mubg,sigr,sigb,sigbg,alphar,alphab,alphabg=sigma_clip(gr5,distprob5,expred)
                                p5=p5+1
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
            maxLbg=(1/(sigbg*np.sqrt(2*np.pi)))
            Pred.append(p_red/maxLred)
            Pblue.append(p_blue/maxLblue)
            Pbg.append(p_bg/maxLbg)
            probgalid.append(glxid)
            magr1.shape=(len(magr1),)
            gr1.shape=(len(gr1),)
            weights=distprob*p_red#*(1-p_blue)
            weights.shape=(len(weights),)
            rfit,info=np.polynomial.polynomial.polyfit(magr1,gr1,deg=1,w=weights,full=True)
            rfit.shape=(2,)
            slope.append(rfit[1])
            yint.append(rfit[0])
        elif len(glxid)==0:
            slope.append(-999.)
            yint.append(-999.)
            mu_r.append(-999.)
            mu_b.append(-999.)
            mu_bg.append(-999.)
            sigma_r.append(-999.)
            sigma_b.append(-999.)
            sigma_bg.append(-999.)
            alpha_r.append(-999.)
            alpha_b.append(-999.)
            alpha_bg.append(-999.)
            Pred.append(np.array([]))
            Pblue.append(np.array([]))
            Pbg.append(np.array([]))
            probgalid.append(-999.)
            converged.append(-999.)
        elif len(glxid)==1:
            slope.append(-999.)
            yint.append(-999.)
            mu_r.append(-999.)
            mu_b.append(-999.)
            mu_bg.append(-999.)
            sigma_r.append(-999.)
            sigma_b.append(-999.)
            sigma_bg.append(-999.)
            alpha_r.append(-999.)
            alpha_b.append(-999.)
            alpha_bg.append(-999.)
            Pred.append(np.array([-999.]))
            Pblue.append(np.array([-999.]))
            Pbg.append(np.array([-999.]))
            probgalid.append(-999.)
            converged.append(-999.)
        elif len(glxid)==2:
            slope.append(-999.)
            yint.append(-999.)
            mu_r.append(-999.)
            mu_b.append(-999.)
            mu_bg.append(-999.)
            sigma_r.append(-999.)
            sigma_b.append(-999.)
            sigma_bg.append(-999.)
            alpha_r.append(-999.)
            alpha_b.append(-999.)
            alpha_bg.append(-999.)
            Pred.append(np.array([-999.,-999.]))
            Pblue.append(np.array([-999.,-999.]))
            Pbg.append(np.array([-999.,-999.]))
            probgalid.append(-999.)
            converged.append(-999.)
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
    p=np.array([p0,p1,p2,p3,p4,p5])
    q=np.array([q0,q1,q2,q3,q4,q5])
    return slope,yint,mu_r,mu_b,mu_bg,sigma_r,sigma_b,sigma_bg,alpha_r,alpha_b,alpha_bg,Pred,Pblue,Pbg,probgalid,converged,p,q



def gmm_restframe(color,band2,expcol,distlim=0.05,n_components=3,tol=0.0000001,galaxyID=galaxyID,hostID=hostID,P_radial=P_radial,P_redshift=P_redshift):
    #full GMM fitting for rest frame color
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
    p0=0
    p1=0
    p2=0
    p3=0
    p4=0
    p5=0
    for x in cluster_ID:
        glxid=galaxyID[np.where((hostID==x)&(Pdist>=distlim))]
        magr1=band2[np.where((hostID==x)&(Pdist>=distlim))]
        rprob1=P_radial[np.where((hostID==x)&(Pdist>=distlim))]
        zprob1=P_redshift[np.where((hostID==x)&(Pdist>=distlim))]
        distprob=rprob1*zprob1
        gr1=color[np.where((hostID==x)&(Pdist>=distlim))]
#        find expected RS color
        zcl=cluster_Z[np.where(cluster_ID==x)]
        ind=0
        expred=expcol[ind]
        if len(glxid) >= 3:
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
            mur,mub,mubg,sigr,sigb,sigbg,alphar,alphab,alphabg=pick_colors(mu,sigma,alpha,expred)
            if np.abs(mur-mub)>=1.*sigr:
                p0=p0+1
            if np.abs(mur-mub)<1.*sigr:
                gr2,distprob2,mur,mub,mubg,sigr,sigb,sigbg,alphar,alphab,alphabg=sigma_clip(gr1,distprob,expred)
                p1=p1+1
                if np.abs(mur-mub)<1.*sigr:
                    gr3,distprob3,mur,mub,mubg,sigr,sigb,sigbg,alphar,alphab,alphabg=sigma_clip(gr2,distprob2,expred)
                    p2=p2+1
                    if np.abs(mur-mub)<1.*sigr:
                        gr4,distprob4,mur,mub,mubg,sigr,sigb,sigbg,alphar,alphab,alphabg=sigma_clip(gr3,distprob3,expred)
                        p3=p3+1
                        if np.abs(mur-mub)<1.*sigr:
                            gr5,distprob5,mur,mub,mubg,sigr,sigb,sigbg,alphar,alphab,alphabg=sigma_clip(gr4,distprob4,expred)
                            p4=p4+1
                            if np.abs(mur-mub)<1.*sigr:
                                gr6,distprob6,mur,mub,mubg,sigr,sigb,sigbg,alphar,alphab,alphabg=sigma_clip(gr5,distprob5,expred)
                                p5=p5+1
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
            weights=distprob*p_red
            weights.shape=(len(weights),)
            rfit,info=np.polynomial.polynomial.polyfit(magr1,gr1,deg=1,w=weights,full=True)
            rfit.shape=(2,)
            slope.append(rfit[1])
            yint.append(rfit[0])
        elif len(glxid)==0:
            slope.append(-999.)
            yint.append(-999.)
            mu_r.append(-999.)
            mu_b.append(-999.)
            mu_bg.append(-999.)
            sigma_r.append(-999.)
            sigma_b.append(-999.)
            sigma_bg.append(-999.)
            alpha_r.append(-999.)
            alpha_b.append(-999.)
            alpha_bg.append(-999.)
            Pred.append(np.array([]))
            Pblue.append(np.array([]))
            Pbg.append(np.array([]))
            probgalid.append(-999.)
            converged.append(-999.)
        elif len(glxid)==1:
            slope.append(-999.)
            yint.append(-999.)
            mu_r.append(-999.)
            mu_b.append(-999.)
            mu_bg.append(-999.)
            sigma_r.append(-999.)
            sigma_b.append(-999.)
            sigma_bg.append(-999.)
            alpha_r.append(-999.)
            alpha_b.append(-999.)
            alpha_bg.append(-999.)
            Pred.append(np.array([-999.]))
            Pblue.append(np.array([-999.]))
            Pbg.append(np.array([-999.]))
            probgalid.append(-999.)
            converged.append(-999.)
        elif len(glxid)==2:
            slope.append(-999.)
            yint.append(-999.)
            mu_r.append(-999.)
            mu_b.append(-999.)
            mu_bg.append(-999.)
            sigma_r.append(-999.)
            sigma_b.append(-999.)
            sigma_bg.append(-999.)
            alpha_r.append(-999.)
            alpha_b.append(-999.)
            alpha_bg.append(-999.)
            Pred.append(np.array([-999.,-999.]))
            Pblue.append(np.array([-999.,-999.]))
            Pbg.append(np.array([-999.,-999.]))
            probgalid.append(-999.)
            converged.append(-999.)
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
    p=np.array([p0,p1,p2,p3,p4,p5])
    return slope,yint,mu_r,mu_b,mu_bg,sigma_r,sigma_b,sigma_bg,alpha_r,alpha_b,alpha_bg,Pred,Pblue,Pbg,probgalid,converged,p



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
grp=grinfo[16]
grq=grinfo[17]
print 'g-r p012345'
print grp

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
rip=riinfo[16]
riq=riinfo[17]
print 'r-i p012345'
print rip

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
izp=izinfo[16]
izq=izinfo[17]
print 'i-z p012345'
print izp


restgr=gmm_restframe(galgr0,galmagR,jimgr)

restslope=restgr[0]
restyint=restgr[1]
restmu_r=restgr[2]
restmu_b=restgr[3]
restmu_bg=restgr[4]
restsigma_r=restgr[5]
restsigma_b=restgr[6]
restsigma_bg=restgr[7]
restalpha_r=restgr[8]
restalpha_b=restgr[9]
restalpha_bg=restgr[10]
restPred=restgr[11]
restPblue=restgr[12]
restPbg=restgr[13]
restprobgalid=restgr[14]
restconverged=restgr[15]
restp=restgr[16]
print 'rest p012345'
print restp



#check for one galaxy belonging to multiple clusters
w=np.where(Pdist>=0.05)
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
galgr02=galgr0[w]


gralpha_r.shape=(gralpha_r.size,)
gralpha_b.shape=(gralpha_b.size,)
grPcolor=(gralpha_r*grPred) + (gralpha_b*grPblue)
grpc2=np.array([item for sublist in grPcolor for item in sublist])
grpc2.shape=(grpc2.size,)
grPred=np.array([item for sublist in grPred for item in sublist])
grPred.shape=(grPred.size,)
grPblue=np.array([item for sublist in grPblue for item in sublist])
grPblue.shape=(grPblue.size,)
grPbg=np.array([item for sublist in grPbg for item in sublist])
grPbg.shape=(grPbg.size,)

rialpha_r.shape=(rialpha_r.size,)
rialpha_b.shape=(rialpha_b.size,)
riPcolor=(rialpha_r*riPred)+(rialpha_b*riPblue)
ripc2=np.array([item for sublist in riPcolor for item in sublist])
ripc2.shape=(ripc2.size,)
riPred=np.array([item for sublist in riPred for item in sublist])
riPred.shape=(riPred.size,)
riPblue=np.array([item for sublist in riPblue for item in sublist])
riPblue.shape=(riPblue.size,)
riPbg=np.array([item for sublist in riPbg for item in sublist])
riPbg.shape=(riPbg.size,)

izalpha_r.shape=(izalpha_r.size,)
izalpha_b.shape=(izalpha_b.size,)
izPcolor=(izalpha_r*izPred)+(izalpha_b*izPblue)
izpc2=np.array([item for sublist in izPcolor for item in sublist])
izpc2.shape=(izpc2.size,)
izPred=np.array([item for sublist in izPred for item in sublist])
izPred.shape=(izPred.size,)
izPblue=np.array([item for sublist in izPblue for item in sublist])
izPblue.shape=(izPblue.size,)
izPbg=np.array([item for sublist in izPbg for item in sublist])
izPbg.shape=(izPbg.size,)

restalpha_r.shape=(restalpha_r.size,)
restalpha_b.shape=(restalpha_b.size,)
restPcolor=(restalpha_r*restPred)+(restalpha_b*restPblue)
restpc2=np.array([item for sublist in restPcolor for item in sublist])
restpc2.shape=(restpc2.size,)
restPred=np.array([item for sublist in restPred for item in sublist])
restPred.shape=(restPred.size,)
restPblue=np.array([item for sublist in restPblue for item in sublist])
restPblue.shape=(restPblue.size,)
restPbg=np.array([item for sublist in restPbg for item in sublist])
restPbg.shape=(restPbg.size,)



######WRITING OUTPUT FILES######
grPmemb=grpc2*prad2*pz2
riPmemb=ripc2*prad2*pz2
izPmemb=izpc2*prad2*pz2
restPmemb=restpc2*prad2*pz2

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
col30=pyfits.Column(name='RESTP_RED',format='E',array=restPred)
col31=pyfits.Column(name='RESTP_BLUE',format='E',array=restPblue)
col32=pyfits.Column(name='RESTP_BG',format='E',array=restPbg)
col33=pyfits.Column(name='REST_P_COLOR',format='E',array=restpc2)
col34=pyfits.Column(name='REST_P_MEMBER',format='E',array=restPmemb)
col35=pyfits.Column(name='GR0',format='E',array=galgr02)



cols=pyfits.ColDefs([col1,col2,col3,col4,col5,col6,col7,col8,col9,col10,col11,col12,col13,col14,col15,col16,col17,col18,col19,col20,col21,col22,col23,col24,col25,col26,col27,col28,col29,col30,col31,col32,col33,col34,col35])
tbhdu=pyfits.BinTableHDU.from_columns(cols)
tbhdu.writeto(member_outfile,clobber=True)

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
restsep=(restmu_r-restmu_b)/restsigma_r

grw=np.where(grsep>=1.)
riw=np.where(risep>=1.)
izw=np.where(izsep>=1.)
restw=np.where(restsep>=1.)

grflag=np.zeros_like(grmu_r)
grflag[grw]=1
riflag=np.zeros_like(rimu_r)
riflag[riw]=1
izflag=np.zeros_like(izmu_r)
izflag[izw]=1
restflag=np.zeros_like(restmu_r)
restflag[restw]=1

zeros=np.zeros_like(cluster_ID)


col1=pyfits.Column(name='MEM_MATCH_ID',format='J',array=cluster_ID)
col2=pyfits.Column(name='RA',format='D',array=cluster_RA)
col3=pyfits.Column(name='DEC',format='D',array=cluster_DEC)
col4=pyfits.Column(name='Z',format='E',array=cluster_Z)
col5=pyfits.Column(name='R200',format='E',array=R200)
col6=pyfits.Column(name='M200',format='E',array=M200)
col7=pyfits.Column(name='N200',format='E',array=N200)
if 'lambda_ngals' in globals():
    col8=pyfits.Column(name="LAMBDA_CHISQ",format='E',array=lambda_ngals)
else:
    col8=pyfits.Column(name="LAMBDA_CHISQ",format='E',array=zeros)
col9=pyfits.Column(name='GR_SLOPE',format='E',array=grslope)
col10=pyfits.Column(name='GR_INTERCEPT',format='E',array=gryint)
col11=pyfits.Column(name='GRMU_R',format='E',array=grmu_r)
col12=pyfits.Column(name='GRMU_B',format='E',array=grmu_b)
col13=pyfits.Column(name='GRSIGMA_R',format='E',array=grsigma_r)
col14=pyfits.Column(name='GRSIGMA_B',format='E',array=grsigma_b)
col15=pyfits.Column(name='GRW_R',format='E',array=gralpha_r)
col16=pyfits.Column(name='GRW_B',format='E',array=gralpha_b)
col17=pyfits.Column(name='RI_SLOPE',format='E',array=rislope)
col18=pyfits.Column(name='RI_INTERCEPT',format='E',array=riyint)
col19=pyfits.Column(name='RIMU_R',format='E',array=rimu_r)
col20=pyfits.Column(name='RIMU_B',format='E',array=rimu_b)
col21=pyfits.Column(name='RISIGMA_R',format='E',array=risigma_r)
col22=pyfits.Column(name='RISIGMA_B',format='E',array=risigma_b)
col23=pyfits.Column(name='RIW_R',format='E',array=rialpha_r)
col24=pyfits.Column(name='RIW_B',format='E',array=rialpha_b)
col25=pyfits.Column(name='GRMU_BG',format='E',array=grmu_bg)
col26=pyfits.Column(name='GRSIGMA_BG',format='E',array=grsigma_bg)
col27=pyfits.Column(name='GRW_BG',format='E',array=gralpha_bg)
col28=pyfits.Column(name='RIMU_BG',format='E',array=rimu_bg)
col29=pyfits.Column(name='RISIGMA_BG',format='E',array=risigma_bg)
col30=pyfits.Column(name='RIW_BG',format='E',array=rialpha_bg)
col31=pyfits.Column(name='IZ_SLOPE',format='E',array=izslope)
col32=pyfits.Column(name='IZ_INTERCEPT',format='E',array=izyint)
col33=pyfits.Column(name='IZMU_R',format='E',array=izmu_r)
col34=pyfits.Column(name='IZMU_B',format='E',array=izmu_b)
col35=pyfits.Column(name='IZSIGMA_R',format='E',array=izsigma_r)
col36=pyfits.Column(name='IZSIGMA_B',format='E',array=izsigma_b)
col37=pyfits.Column(name='IZW_R',format='E',array=izalpha_r)
col38=pyfits.Column(name='IZW_B',format='E',array=izalpha_b)
col39=pyfits.Column(name='IZMU_BG',format='E',array=izmu_bg)
col40=pyfits.Column(name='IZSIGMA_BG',format='E',array=izsigma_bg)
col41=pyfits.Column(name='IZW_BG',format='E',array=izalpha_bg)
col42=pyfits.Column(name='GR_SEP_FLAG',format='L',array=grflag)
col43=pyfits.Column(name='RI_SEP_FLAG',format='L',array=riflag)
col44=pyfits.Column(name='IZ_SEP_FLAG',format='L',array=izflag)
col45=pyfits.Column(name='REST_SLOPE',format='E',array=restslope)
col46=pyfits.Column(name='REST_INTERCEPT',format='E',array=restyint)
col47=pyfits.Column(name='RESTMU_R',format='E',array=restmu_r)
col48=pyfits.Column(name='RESTMU_B',format='E',array=restmu_b)
col49=pyfits.Column(name='RESTMU_BG',format='E',array=restmu_bg)
col50=pyfits.Column(name='RESTSIGMA_R',format='E',array=restsigma_r)
col51=pyfits.Column(name='RESTSIGMA_B',format='E',array=restsigma_b)
col52=pyfits.Column(name='RESTSIGMA_BG',format='E',array=restsigma_bg)
col53=pyfits.Column(name='RESTW_R',format='E',array=restalpha_r)
col54=pyfits.Column(name='RESTW_B',format='E',array=restalpha_b)
col55=pyfits.Column(name='RESTW_BG',format='E',array=restalpha_bg)
col56=pyfits.Column(name='REST_SEP_FLAG',format='L',array=restflag)


cols=pyfits.ColDefs([col1,col2,col3,col4,col5,col6,col7,col8,col9,col10,col11,col12,col13,col14,col15,col16,col17,col18,col19,col20,col21,col22,col23,col24,col25,col26,col27,col28,col29,col30,col31,col32,col33,col34,col35,col36,col37,col38,col39,col40,col41,col42,col43,col44,col45,col46,col47,col48,col49,col50,col51,col52,col53,col54,col55,col56])
tbhdu=pyfits.BinTableHDU.from_columns(cols)
tbhdu.writeto(cluster_outfile,clobber=True)



print 'success!!'
exit()
