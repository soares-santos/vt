
# coding: utf-8

# In[23]:

get_ipython().system(u'jupyter nbconvert --to script py.ipynb')


# In[4]:

import matplotlib
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits, ascii
from scipy.special import erf 
from scipy.interpolate import interp1d
import os
from configparser import ConfigParser


# In[21]:

parser = ConfigParser()
parser.read('pipecheck.ini')

dataparser.get('bug_tracker', 'url')


# In[17]:

# configure matplotlib to make figures look nicer
matplotlib.rc("figure", facecolor="white")
matplotlib.rc("legend", frameon=False)
#matplotlib.rcParams.update({'font.size': 24})
#matplotlib.rcParams['axes.linewidth'] = 2
#matplotlib.rcParams.update({'ytick.major.size': 24})
#matplotlib.rcParams.update({'xtick.major.size': 24})
#matplotlib.rcParams['xtick.major.width'] = 2
#matplotlib.rcParams['ytick.major.width'] = 2
#matplotlib.rcParams.update({'xtick.minor.size': 12})
#matplotlib.rcParams.update({'ytick.minor.size': 12})
#matplotlib.rcParams['xtick.minor.width'] = 2
#matplotlib.rcParams['ytick.minor.width'] = 2


# In[ ]:

class DataSetColumnNames:
    def __init__(self):
        self.CLUSTER_ID='MEM_MATCH_ID'
        self.Z='Z'


# In[4]:

class MembAssign:

    def __init__(self,mode='both',datadir='.',simsdir='.',auxfilesdir='.',outputsdir='.',
                 cluster_data_file_0=None,cluster_data_file_1=None,cluster_data_file_2=None,
                 cluster_sigma_phot_file=None,external_data_files=[]):
        
        if mode not in ['real','sims','both']:
            sys.exit("Error: Unknonw mode. Valid options are: real , sims , both ")

        if inifile is not None:    
            
        c=DataSetColumnNames()
            
        if (mode=='real') or (mode=='both'):
            cluster_data_file_0=os.path.join(datadir,cluster_data_file_0)
            cluster_hdulist_0 = fits.open(cluster_data_file_0)
            data0 = cluster_hdulist_0[1].data
            datalabel = cluster_data_file_0.split('/')[-1].split('_clusters.fit')[0]
        
        if (mode=='sims' or (mode=='both')):
            cluster_data_file_1=os.path.join(simsdir,cluster_data_file_1)
            cluster_data_file_2=os.path.join(simsdir,cluster_data_file_2)
            cluster_hdulist_1 = fits.open(cluster_data_file_1)
            cluster_hdulist_2 = fits.open(cluster_data_file_2)
            data1 = cluster_hdulist_1[1].data
            data2 = cluster_hdulist_2[1].data
            simlabel = cluster_data_file_1.split('/')[-1].split('_clusters.fit')[0]
            # match sim results to truth table
            x=data1[c.CLUSTER_ID]
            y=data2[c.CLUSTER_ID]
            ix1 = np.where(np.in1d(x.ravel(), y).reshape(x.shape))  ## data1[ix1]
            ix2 = np.where(np.in1d(y.ravel(), x).reshape(y.shape))  ## data2[ix2]

        if cluster_sigma_phot_file is not None:
            cluster_sigma_phot_file=os.path.join(auxfilesdir,cluster_sigma_phot_file)
            files.append(cluster_sigma_phot_file)
            
        for edf in external_data_files:  
            external_data_files=os.path.join(auxfilesdir,external_data_files)
            files.append(external_data_files)
            
        for f in files:
                if not os.path.exists(f): sys.exit("Error: File not found: "+f)
                        
        if not os.path.isdir(outputsdir): os.makedirs(outputsdir)
        self.odir = outputsdir
            

# set labels for plots and output files

# read cluster data files


# read phot err file and match to real clusters
sigma2phot = ascii.read(cluster_sigma_phot_file) 
x=data0['MEM_MATCH_ID']
y=sigma2phot['haloid']
ix = np.where(np.in1d(x.ravel(), y).reshape(x.shape))  ## data0[ix]
jx = np.where(np.in1d(y.ravel(), x).reshape(y.shape))  ## sigma2phot[jx]

# read results from the literature and sort by redshift 
r1color = ascii.read(sdss_rm_color_file) 
r1slope = ascii.read(sdss_rm_slope_file) 
r1sigma = ascii.read(sdss_rm_sigma_file) 
r1color.sort('z')
r1slope.sort('z')
r1sigma.sort('z')


# In[5]:

# interpolate phot err data as a function of cluster redshift
grphotsigma2 = interp1d(data0['Z'][ix],sigma2phot['sigma^2(g-r)'][jx],bounds_error=False,fill_value='extrapolate')
riphotsigma2 = interp1d(data0['Z'][ix],sigma2phot['sigma^2(r-i)'][jx],bounds_error=False,fill_value='extrapolate')
izphotsigma2 = interp1d(data0['Z'][ix],sigma2phot['sigma^2(i-z)'][jx],bounds_error=False,fill_value='extrapolate')
# compute red sequence instrinsic width
gr_int_red1 = np.sqrt(data1['GRSIGMA_R']**2-grphotsigma2(data1['Z']))
ri_int_red1 = np.sqrt(data1['RISIGMA_R']**2-riphotsigma2(data1['Z']))
iz_int_red1 = np.sqrt(data1['RISIGMA_R']**2-izphotsigma2(data1['Z']))
gr_int_red2 = np.sqrt(data2['GRSIGMA_R']**2-grphotsigma2(data2['Z']))
ri_int_red2 = np.sqrt(data2['RISIGMA_R']**2-riphotsigma2(data2['Z']))
iz_int_red2 = np.sqrt(data2['RISIGMA_R']**2-izphotsigma2(data2['Z']))
# compute blue cloud intrinsic width
gr_int_blue1 = np.sqrt(data1['GRSIGMA_B']**2-grphotsigma2(data1['Z']))
ri_int_blue1 = np.sqrt(data1['RISIGMA_B']**2-riphotsigma2(data1['Z']))
iz_int_blue1 = np.sqrt(data1['RISIGMA_B']**2-izphotsigma2(data1['Z']))
gr_int_blue2 = np.sqrt(data2['GRSIGMA_B']**2-grphotsigma2(data2['Z']))
ri_int_blue2 = np.sqrt(data2['RISIGMA_B']**2-riphotsigma2(data2['Z']))
iz_int_blue2 = np.sqrt(data2['RISIGMA_B']**2-izphotsigma2(data2['Z']))


# In[55]:

# plot phot err vs z
z=np.linspace(0,0.9,20)
plt.figure(figsize=(8,6))
plt.plot(z,np.sqrt(grphotsigma2(z)),label='$g-r$',lw=4,c='g')
plt.plot(z,np.sqrt(riphotsigma2(z)),label='$r-i$',lw=4,c='r')
plt.plot(z,np.sqrt(izphotsigma2(z)),label='$i-z$',lw=4,c='indigo')
plt.xlim(-0.1,1.1)
plt.ylim(0.01,0.7)
plt.yscale('log')
plt.xlabel('redshift')
plt.ylabel('member galaxy color uncertainty')
plt.legend(loc=2)
plt.title(datalabel)
plt.savefig(pdir+'datalabel'+'_colorerr.png')
plt.grid(True,which='both')


# In[ ]:

# define plotting function
def plot_rs(ax,x,y,z,xl='',yl='',tl='',cmap='viridis',lw=6,vmin=None,vmax=None,pfunc='scatter',fontsize=12,alpha=1):
    if pfunc=='hexbin':
        im=ax.hexbin(x,y,C=z,cmap=cmap,lw=lw,vmin=vmin,vmax=vmax,alpha=alpha)
    if pfunc=='scatter':
        im=ax.scatter(x,y,c=z,cmap=cmap,s=lw,alpha=alpha,vmin=vmin,vmax=vmax)
    ax.set_ylabel(yl,fontsize=fontsize)
    ax.set_xlabel(xl,fontsize=fontsize)
    ax.set_title(tl,fontsize=fontsize)
    #ax.tick_params(axis='both', which='major', labelsize=fontsize)
    return im


# In[ ]:

#### old function 
#def plot_rs(ax,x,y,z,xl='',yl='',tl='',cmap='viridis',lw=6,vmin=None,vmax=None):
#    fontsize = 12
#    #hb=ax.hexbin(x,y,C=z, cmap=cmap,lw=lw)
#    hb=ax.scatter(x,y,c=z, cmap=cmap,s=lw,alpha=1,vmin=vmin,vmax=vmax)
#    ax.set_ylabel(yl,fontsize=fontsize)
#    ax.set_xlabel(xl,fontsize=fontsize)
#    ax.set_title(tl,fontsize=fontsize)
#    ax.tick_params(axis='both', which='major', labelsize=fontsize)
#    return hb


# In[ ]:

#### old function 
#def scatter(sigmaobs,sigma2phot):
#    x=sigmaobs**2-sigma2phot
#    return np.sqrt(x)


# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:

a='''
fig = plt.figure(figsize=(15,10))
gs = matplotlib.gridspec.GridSpec(2, 4, width_ratios=[1,1,1,0.1])
ax1 = plt.subplot(gs[0])#(241)
ax2 = plt.subplot(gs[1])#(242)
ax3 = plt.subplot(gs[2])#(243)
ax4 = plt.subplot(gs[4])#(245)
ax5 = plt.subplot(gs[5])#(246)
ax6 = plt.subplot(gs[6])#(247)
ax7 = plt.subplot(gs[3])#(244)
ax8 = plt.subplot(gs[7])#(248)
z1=data1['Z']
z2=data2['Z']
m=np.log10(data1['M200'])
m1=data1['GRW_B']
m2=data2['GRW_B']
tl=['obs','truth']
lw=100
im=plot_rs(ax1,z1,data1['GRMU_R'],m1,tl=tl[0],xl='redshift',yl='peak color (g-r)',lw=lw)
im=plot_rs(ax2,z1,data1['RIMU_R'],m1,tl=tl[0],xl='redshift',yl='peak color (r-i)',lw=lw)
im=plot_rs(ax3,z1,data1['IZMU_R'],m1,tl=tl[0],xl='redshift',yl='peak color (i-z)',lw=lw)
#plot_rs(ax2,z,scatter(data1['GRSIGMA_R'][ix],sigma2phot['sigma^2(g-r)'][jx]),m1,tl=tl[0],xl='redshift',yl='width',lw=lw)
#plot_rs(ax3,z,data1['GR_SLOPE'][ix],m1,tl=tl[0],xl='redshift',yl='slope')
plot_rs(ax4,z2,data2['GRMU_R'],m2,tl=tl[1],xl='redshift',yl='peak color (g-r)',lw=lw)
plot_rs(ax5,z2,data2['RIMU_R'],m2,tl=tl[1],xl='redshift',yl='peak color (r-i)',lw=lw)
plot_rs(ax6,z2,data2['IZMU_R'],m2,tl=tl[1],xl='redshift',yl='peak color (i-z)',lw=lw)
#plot_rs(ax5,z,scatter(data2['GRSIGMA_R'][ix],sigma2phot['sigma^2(g-r)'][jx]),m2,tl=tl[1],xl='redshift',yl='width',lw=lw)
#plot_rs(ax6,z,data2['GR_SLOPE'][ix],m2,tl=tl[1],xl='redshift',yl='slope')
ax1.plot(r1color['z'].data,r1color['color'].data,c='m',lw=2)
#ax2.plot(r1sigma['z'].data,r1sigma['sigma'].data,c='m',lw=2)
#ax3.plot(r1slope['z'].data,r1slope['slope'].data,c='m',lw=2)
ax4.plot(r1color['z'].data,r1color['color'].data,c='m',lw=2)
#ax5.plot(r1sigma['z'].data,r1sigma['sigma'].data,c='m',lw=2)
#ax6.plot(r1slope['z'].data,r1slope['slope'].data,c='m',lw=2)
ax1.set_ylim(0.6,2.2)
ax2.set_ylim(0.2,1.4)
ax3.set_ylim(0.1,0.9)
#ax4.set_ylim(0.6,4)
#ax5.set_ylim(0,1)
#ax6.set_ylim(-0.3,0.3)
fig.colorbar(im,cax=ax7,label='blue fraction')
fig.colorbar(im,cax=ax8,label='blue fraction')
fig.tight_layout()
fig.savefig('rs_params_c10c03.png')
fig.savefig('rs_params_c10c03.pdf')
'''


# In[ ]:




# In[ ]:

#Mass-richness relation functions (see Tinker et al 2011)
def ncen(M,log_Mmin,sigma):
    #takes logM_min and logSigma_logM from paper. returns Ncentral from paper
    sigma=10**sigma
    return (1./2.)*(1+erf((np.log10(M)-log_Mmin)/sigma))

def ntot(M,z=0.,params=[11.6,12.45,1.0,12.25,-0.69]):
    #takes logMmin, logSigma_logM, logMsat, logMcut from paper. Returns Ntotal=Ncentral+Nsatellite from paper
    log_Mmin=params[0]
    log_Msat=logMsat(z,params[1],a=0.)
    alpha=alpha_sat(z,params[2],a=-0.5,z0=0.4)
    log_Mcut=params[3]
    sigma=params[4]
    Msat=10**log_Msat
    Mcut=10**log_Mcut
    return ncen(M,log_Mmin,sigma)*(1+(((M/Msat)**alpha)*np.exp(-Mcut/M)))

#Msat redshift dependence - no longer used, didn't work
def logMsat(z,M0=12.33,a=-0.27):
    return M0 + a*z

#alpha_sat redshift dependence - currently used for redshift varying HOD model
def alpha_sat(z,alpha0=1.0,a=0.0,z0=0.5):
    if z> z0: return alpha0 + a*(z-z0) #+ (alpha0+a*z)
    else: return alpha0

def hod_mass_z(N,z,params=[11.6,12.45,1.0,12.25,-0.69]):
    #params: logMmin,logMsat,alphasat,logMcut,logsigmalogM directly from table 4 of Tinker et al. paper
    Mmin=params[0]
    Msat=params[1]
    alpha=params[2]
    Mcut=params[3]
    sigma=params[4]
    mass=np.logspace(10,16,num=60,dtype=float)
    m200c=np.zeros_like(N)
    for i in range(len(z)):
        Msat=logMsat(z[i],params[1],0.) # make msat to be a function of z
        alpha=alpha_sat(z[i],params[2],-0.5,0.4) # make alpha a functin of z
        m=interp1d(ntot(mass,Msat,Mmin,sigma,alpha,Mcut),mass,
                   bounds_error=False,fill_value='extrapolate')
        for j in range(len(N[i])): N[i][j]=max(N[i][j],0.1) #set minimum value for mass conversions to prevent code from failing
        m200c[i]=m(N[i]/100.)  ### use this with hlin's sim
        #m200c[i]=m(N[i])  ## use this with real data
        for j in range(len(N[i])): print N[i][j],m200c[i][j]
    return m200c


# In[ ]:

print data1['GRW_B'][ix1]-data2['GRW_B'][ix2]


# In[ ]:

#plt.figure(figsize=(15,8))
npts=80000
mass=np.logspace(11,15,num=npts,dtype=float)
reds=np.random.random(npts)
ngals,mhalo,redsh = ([] for i in range(3))
for m,z in zip(mass,reds): 
    ngals.append(ntot(m,z))
    mhalo.append(m)
    redsh.append(z)
modellabel='Tinker13 w/ a(z>0.4)=a(0)+0.5z'    
plt.hexbin(mhalo,ngals,C=redsh,xscale='log',yscale='log',gridsize=500,alpha=1,lw=None,cmap='YlOrRd',vmin=0,vmax=1)
plt.colorbar(label='redshift')
mhalo_out=data1['M200'][ix1]
ngals_sim=(data2['N_B']+data2['N_R'])[ix2]
plt.scatter(mhalo_out,ngals_sim,s=50,lw=0,c='k',label=simlabel) 
#plt.xscale('log')
#plt.yscale('log')
plt.xlim(1e11,1e17)
plt.ylim(5e-1,5e3)
plt.xlabel('M200c')
plt.ylabel('Ngals')
plt.grid(True)
plt.legend(loc=2,frameon=False,borderpad=1) #fontsize=24
plt.title(modellabel, y=1.02)
plt.savefig(pdir+simlabel+'_hod_check.png')


# In[ ]:




# In[ ]:

def quantity2plot(q):
    
    cmap=['RdYlBu','RdYlBu','RdYlBu']

    #rs color
    if q == 'MU_R':
        f1=data1['GRMU_R'][ix1]
        f2=data2['GRMU_R'][ix2]
        f3=data1['RIMU_R'][ix1]
        f4=data2['RIMU_R'][ix2]
        f5=data1['IZMU_R'][ix1]
        f6=data2['IZMU_R'][ix2]
        c1=data1['GRW_B'][ix1]-data2['GRW_B'][ix2]
        c2=data1['RIW_B'][ix1]-data2['RIW_B'][ix2]
        c3=data1['IZW_B'][ix1]-data2['IZW_B'][ix2]
        tl='red sequence color' 
        cl='delta blue fraction (obs-truth)'
        ylim=[-0.1,2.1]
        xlim=[-0.1,1.1]
        clim=[-0.4,0.4]
    # rs sigma
    if q == 'SIGMA_R_OBS':
        f1=data1['GRSIGMA_R'][ix1]
        f2=data2['GRSIGMA_R'][ix2]
        f3=data1['RISIGMA_R'][ix1]
        f4=data2['RISIGMA_R'][ix2]
        f5=data1['IZSIGMA_R'][ix1]
        f6=data2['IZSIGMA_R'][ix2]
        c1=data1['GRW_B'][ix1]-data2['GRW_B'][ix2]
        c2=data1['RIW_B'][ix1]-data2['RIW_B'][ix2]
        c3=data1['IZW_B'][ix1]-data2['IZW_B'][ix2]
        tl='rs width w/ photerr'
        cl='delta blue fraction (obs-truth)'
        ylim=[-0.1,1.1]
        xlim=[-0.1,1.1]
        clim=[-0.4,0.4]
    if q == 'SIGMA_R_INT':
        f1=gr_int_scatter1[ix1] # np.sqrt(data1['GRSIGMA_R']**2-grphotsigma2(data1['Z']))
        f3=ri_int_scatter1[ix1] 
        f5=iz_int_scatter1[ix1] 
        f2=gr_int_scatter2[ix2] 
        f4=ri_int_scatter2[ix2] 
        f6=iz_int_scatter2[ix2]         
        c1=data1['GRW_B'][ix1]-data2['GRW_B'][ix2]
        c2=data1['RIW_B'][ix1]-data2['RIW_B'][ix2]
        c3=data1['IZW_B'][ix1]-data2['IZW_B'][ix2]
        tl='red sequence width'
        cl='delta blue fraction (obs-truth)'
        ylim=[-0.1,1.1]
        xlim=[-0.1,1.1]
        clim=[-0.4,0.4]
        # rs slope 
    if q == 'SLOPE':
        f1=data1['GR_SLOPE'][ix1]
        f2=data2['GR_SLOPE'][ix2]
        f3=data1['RI_SLOPE'][ix1]
        f4=data2['RI_SLOPE'][ix2]
        f5=data1['IZ_SLOPE'][ix1]
        f6=data2['IZ_SLOPE'][ix2]
        c1=data1['GRW_B'][ix1]-data2['GRW_B'][ix2]
        c2=data1['RIW_B'][ix1]-data2['RIW_B'][ix2]
        c3=data1['IZW_B'][ix1]-data2['IZW_B'][ix2]
        tl='red sequence slope'
        cl='delta blue fraction (obs-truth)'
        ylim=[-1.1,1.1]
        xlim=[-0.1,1.1]
        clim=[-0.4,0.4]
    # blue cloud color
    if q == 'MU_B':
        f1=data1['GRMU_B'][ix1]
        f2=data2['GRMU_B'][ix2]
        f3=data1['RIMU_B'][ix1]
        f4=data2['RIMU_B'][ix2]
        f5=data1['IZMU_B'][ix1]
        f6=data2['IZMU_B'][ix2]
        c1=data1['GRW_B'][ix1]-data2['GRW_B'][ix2]
        c2=data1['RIW_B'][ix1]-data2['RIW_B'][ix2]
        c3=data1['IZW_B'][ix1]-data2['IZW_B'][ix2]
        tl='blue cloud color'
        cl='delta blue fraction (obs-truth)'
        ylim=[-0.6,1.6]
        xlim=[-0.1,1.1]
        clim=[-0.4,0.4]
    # blue cloud sigma
    if q == 'SIGMA_B_OBS':
        f1=data1['GRSIGMA_B'][ix1]
        f2=data2['GRSIGMA_B'][ix2]
        f3=data1['RISIGMA_B'][ix1]
        f4=data2['RISIGMA_B'][ix2]
        f5=data1['IZSIGMA_B'][ix1]
        f6=data2['IZSIGMA_B'][ix2]
        c1=data1['GRW_B'][ix1]-data2['GRW_B'][ix2]
        c2=data1['RIW_B'][ix1]-data2['RIW_B'][ix2]
        c3=data1['IZW_B'][ix1]-data2['IZW_B'][ix2]
        tl='bc width w/ photerr'
        cl='delta blue fraction (obs-truth)'
        ylim=[0.01,0.59]
        xlim=[-0.1,1.1]
        clim=[-0.4,0.4]
    if q == 'SIGMA_B_INT':
        f1=gr_int_blue1[ix1] # np.sqrt(data1['GRSIGMA_B']**2-grphotsigma2(data1['Z']))
        f3=ri_int_blue1[ix1] 
        f5=iz_int_blue1[ix1] 
        f2=gr_int_blue2[ix2] 
        f4=ri_int_blue2[ix2] 
        f6=iz_int_blue2[ix2]
        c1=data1['GRW_B'][ix1]-data2['GRW_B'][ix2]
        c2=data1['RIW_B'][ix1]-data2['RIW_B'][ix2]
        c3=data1['IZW_B'][ix1]-data2['IZW_B'][ix2]
        tl='blue cloud width'
        cl='delta blue fraction (obs-truth)'
        ylim=[0.01,0.59]
        xlim=[-0.1,1.1]
        clim=[-0.4,0.4]
    # blue fraction
    if q == 'W_B':
        f1=data1['GRW_B'][ix1]
        f2=data2['GRW_B'][ix2]
        f3=data1['RIW_B'][ix1]
        f4=data2['RIW_B'][ix2]
        f5=data1['IZW_B'][ix1]
        f6=data2['IZW_B'][ix2]
        c1=data1['GRMU_B'][ix1]-data2['GRMU_B'][ix2]
        c2=data1['RIMU_B'][ix1]-data2['RIMU_B'][ix2]
        c3=data1['IZMU_B'][ix1]-data2['IZMU_B'][ix2]
        tl='blue fraction (obs)'
        cl='delta rs color (obs-truth)'
        ylim=[-0.1,1.1]
        xlim=[-0.1,1.1]
        clim=[-0.4,0.4]
        cmap=['PiYG','PiYG','PiYG']
    # red fraction
    if q == 'W_R':
        f1=data1['GRW_R'][ix1]
        f2=data2['GRW_R'][ix2]
        f3=data1['RIW_R'][ix1]
        f4=data2['RIW_R'][ix2]
        f5=data1['IZW_R'][ix1]
        f6=data2['IZW_R'][ix2]
        c1=data1['GRMU_R'][ix1]-data2['GRMU_R'][ix2]
        c2=data1['RIMU_R'][ix1]-data2['RIMU_R'][ix2]
        c3=data1['IZMU_R'][ix1]-data2['IZMU_R'][ix2]
        tl='red fraction (obs)'
        cl='delta rs color (obs-truth)'
        ylim=[-0.1,1.1]
        xlim=[-0.1,1.1]
        clim=[-0.85,0.85]
        cmap=['coolwarm','coolwarm','coolwarm']
    
    # (obs-truth)
    d1=(f1-f2)/f2 #f1-f2
    d2=(f3-f4)/f4
    d3=(f5-f6)/f4

    # plot subtitles
    stl=['(obs)','(truth)','(obs-truth)/(truth)']
    
    
    #adjust vmin,vmax so that they are centered in zero
    vmin=[c1.min(),c2.min(),c3.min()]
    vmax=[c1.max(),c2.max(),c3.max()]
    for i in range(len(vmin)):
        mv=(vmin[i]+vmax[i])/2
        if mv < 0: 
            vmin[i] = round(vmin[i] - 2 * mv,1)
        else: 
            vmax[i] = round(vmax[i] + 2 * mv,1)
                
    return f1,f2,f3,f4,f5,f6,d1,d2,d3,c1,c2,c3,stl,cl,xlim,ylim,clim,vmin,vmax,cmap,tl


# In[ ]:




# In[ ]:

fig = plt.figure(figsize=(15,10))
gs = matplotlib.gridspec.GridSpec(3, 4, width_ratios=[1,1,1,0.1],hspace=0,wspace=0)
ax1 = plt.subplot(gs[0])
ax2 = plt.subplot(gs[1],sharex=ax1,sharey=ax1)
ax3 = plt.subplot(gs[2],sharex=ax1)
ax4 = plt.subplot(gs[3])
ax5 = plt.subplot(gs[4],sharex=ax1,sharey=ax1)
ax6 = plt.subplot(gs[5],sharex=ax1,sharey=ax1)
ax7 = plt.subplot(gs[6],sharex=ax1,sharey=ax3)
ax8 = plt.subplot(gs[7])
ax9 = plt.subplot(gs[8],sharex=ax1,sharey=ax1)
ax10 = plt.subplot(gs[9],sharex=ax1,sharey=ax1)
ax11 = plt.subplot(gs[10],sharex=ax1,sharey=ax3)
ax12 = plt.subplot(gs[11])
lw=50
# x axis
z1=data1['Z'][ix1]
z2=data2['Z'][ix2]
xl='redshift'
# y axis and colorbar
f1,f2,f3,f4,f5,f6,d1,d2,d3,c1,c2,c3,tl,cl,xlim,ylim,clim,vmin,vmax,cmap,figtitle = quantity2plot('W_R')
yl1=figtitle+' (g-r)'
yl2=figtitle+' (r-i)'
yl3=figtitle+' (i-z)'
cl1=cl 
cl2=cl 
cl3=cl 
# gr plots
im1 =plot_rs(ax1,z1,f1,c1,tl=tl[0],xl=xl,cmap=cmap[0],lw=lw,yl=yl1,vmin=vmin[0],vmax=vmax[0])
im2 =plot_rs(ax2,z2,f2,c1,tl=tl[1],xl=xl,cmap=cmap[0],lw=lw)
im3 =plot_rs(ax3,z2,d1,c1,tl=tl[2],xl=xl,cmap=cmap[0],lw=lw)
plt.setp(ax2.get_yticklabels(), visible=False)
ax2.yaxis.set_ticks_position('left')
plt.setp(ax3.get_yticklabels(), visible=False)
ax3.annotate('$+$20%',xy=(0, 0.2),horizontalalignment='center', verticalalignment='center',bbox=dict(facecolor='w', alpha=1,lw=0))
ax3.annotate('0%',xy=(0, 0.0),horizontalalignment='center', verticalalignment='center',bbox=dict(facecolor='w', alpha=1,lw=0))
ax3.annotate('$-$20%',xy=(0,-0.2),horizontalalignment='center', verticalalignment='center',bbox=dict(facecolor='w', alpha=1,lw=0))
ax3.yaxis.set_ticks_position('none')
ax3.xaxis.set_ticks_position('none')
nsteps=5
dv=round((vmax[0]-vmin[0])/nsteps,1)
fig.colorbar(im1,cax=ax4,label=cl1,ticks=np.arange(vmin[0],vmax[0],dv)[1:])
# ri plots
im5 =plot_rs(ax5,z1,f3,c2,xl=xl,cmap=cmap[1],lw=lw,yl=yl2,vmin=vmin[1],vmax=vmax[1])
im6 =plot_rs(ax6,z2,f4,c2,xl=xl,cmap=cmap[1],lw=lw)
im7 =plot_rs(ax7,z2,d2,c2,xl=xl,cmap=cmap[1],lw=lw)
plt.setp(ax6.get_yticklabels(), visible=False)
ax6.yaxis.set_ticks_position('left')
plt.setp(ax7.get_yticklabels(), visible=False)
ax7.annotate('$+$20%',xy=(0, 0.2),horizontalalignment='center', verticalalignment='center',bbox=dict(facecolor='w', alpha=1,lw=0))
ax7.annotate('0%',xy=(0, 0.0),horizontalalignment='center', verticalalignment='center',bbox=dict(facecolor='w', alpha=1,lw=0))
ax7.annotate('$-$20%',xy=(0,-0.2),horizontalalignment='center', verticalalignment='center',bbox=dict(facecolor='w', alpha=1,lw=0))
ax7.yaxis.set_ticks_position('none')
ax7.xaxis.set_ticks_position('none')
##nsteps=6.6
dv=round((vmax[1]-vmin[1])/nsteps,1)
fig.colorbar(im5,cax=ax8,label=cl2,ticks=np.arange(vmin[1],vmax[1],dv)[1:])
# iz plots
im9 =plot_rs(ax9 ,z1,f5,c3,xl=xl,cmap=cmap[2],lw=lw,yl=yl3,vmin=vmin[2],vmax=vmax[2])
im10=plot_rs(ax10,z2,f6,c3,xl=xl,cmap=cmap[2],lw=lw)
im11=plot_rs(ax11,z2,d3,c3,xl=xl,cmap=cmap[2],lw=lw)
plt.setp(ax10.get_yticklabels(), visible=False)
ax10.yaxis.set_ticks_position('left')
plt.setp(ax11.get_yticklabels(), visible=False)
ax11.annotate('$+$20%',xy=(0, 0.2),horizontalalignment='center', verticalalignment='center',bbox=dict(facecolor='w', alpha=1,lw=0))
ax11.annotate('0%',xy=(0, 0.0),horizontalalignment='center', verticalalignment='center',bbox=dict(facecolor='w', alpha=1,lw=0))
ax11.annotate('$+$20%',xy=(0,-0.2),horizontalalignment='center', verticalalignment='center',bbox=dict(facecolor='w', alpha=1,lw=0))
ax11.yaxis.set_ticks_position('none')
ax11.xaxis.set_ticks_position('none')
#nsteps=6
dv=round((vmax[2]-vmin[2])/nsteps,1)
if dv < 0.05: dv = 0.05
fig.colorbar(im9,cax=ax12,label=cl3,ticks=np.arange(vmin[2],vmax[2],dv)[1:])
# adjustments
ax1.set_ylim(ylim)
ax1.set_xlim(xlim)
ax3.set_ylim(clim)
ax3.grid(True)
ax7.grid(True)
ax11.grid(True)
#fig.savefig('bluefrac_c10s3_sim_MU_R.png')
#fig.savefig('bluefrac_c10s3_sim_MU_R.pdf')


# In[ ]:




# In[ ]:

fig = plt.figure(figsize=(15,10))
gs = matplotlib.gridspec.GridSpec(3, 4, width_ratios=[1,1,1,0.2],hspace=0,wspace=0)
ax1 = plt.subplot(gs[0])
ax2 = plt.subplot(gs[1],sharex=ax1,sharey=ax1)
ax3 = plt.subplot(gs[2],sharex=ax1)
ax4 = plt.subplot(gs[3])
ax5 = plt.subplot(gs[4],sharex=ax1,sharey=ax1)
ax6 = plt.subplot(gs[5],sharex=ax1,sharey=ax1)
ax7 = plt.subplot(gs[6],sharex=ax1,sharey=ax3)
ax8 = plt.subplot(gs[7])
ax9 = plt.subplot(gs[8],sharex=ax1,sharey=ax1)
ax10 = plt.subplot(gs[9],sharex=ax1,sharey=ax1)
ax11 = plt.subplot(gs[10],sharex=ax1,sharey=ax3)
ax12 = plt.subplot(gs[11])
lw=100
x=data1['MEM_MATCH_ID']
y=data2['MEM_MATCH_ID']
ix1 = np.where(np.in1d(x.ravel(), y).reshape(x.shape))  ## data1[ix]
ix2 = np.where(np.in1d(y.ravel(), x).reshape(y.shape))  ## sigma2phot[jx]
z1=data1['Z'][ix1]
z2=data2['Z'][ix2]
m1=data1['GRMU_R'][ix1]
m2=data2['GRMU_R'][ix2]
m3=data1['RIMU_R'][ix1]
m4=data2['RIMU_R'][ix2]
m5=data1['IZMU_R'][ix1]
m6=data2['IZMU_R'][ix2]
f1=data1['GRW_B'][ix1]
f2=data2['GRW_B'][ix2]
f3=data1['RIW_B'][ix1]
f4=data2['RIW_B'][ix2]
f5=data1['IZW_B'][ix1]
f6=data2['IZW_B'][ix2]
tl=['blue fraction (obs)','blue fraction (truth)',' delta blue fraction (obs-truth)']
xl='redshift'
yl1='(g-r)'
yl2='(r-i)'
yl3='(i-z)'
cl1='delta red sequence color (g-r)'
cl2='delta red sequence color (r-i)'
cl3='delta red sequence color (i-z)'
cmap='jet'
im1=plot_rs(ax1,z1,f1,m1-m2,tl=tl[0],xl=xl,yl=yl1,cmap=cmap,lw=lw)
im2=plot_rs(ax2,z2,f2,m1-m2,tl=tl[1],xl=xl,cmap=cmap,lw=lw)
im3=plot_rs(ax3,z2,f1-f2,m1-m2,tl=tl[2],xl=xl,cmap=cmap,lw=lw)
fig.colorbar(im1,cax=ax4,label=cl1)#,orientation='horizontal')
im5=plot_rs(ax5,z1,f3,m3-m4,xl=xl,yl=yl2,cmap=cmap,lw=lw)
im6=plot_rs(ax6,z2,f4,m3-m4,xl=xl,cmap=cmap,lw=lw)
im7=plot_rs(ax7,z2,f3-f4,m3-m4,xl=xl,cmap=cmap,lw=lw)
fig.colorbar(im5,cax=ax8,label=cl2)#,orientation='horizontal')
im9=plot_rs(ax9,z1,f5,m5-m6,xl=xl,yl=yl3,cmap=cmap,lw=lw)
im10=plot_rs(ax10,z2,f6,m5-m6,xl=xl,cmap=cmap,lw=lw)
im11=plot_rs(ax11,z2,f5-f6,m5-m6,xl=xl,cmap=cmap,lw=lw)
fig.colorbar(im9,cax=ax12,label=cl3)#,orientation='horizontal')
ax1.set_ylim(-0.1,1.1)
ax1.set_xlim(-0.1,1.1)
ax3.set_ylim(-0.85,0.85)
ax3.grid(True)
ax7.grid(True)
ax11.grid(True)

fig.tight_layout()
fig.savefig('bluefrac_c10s3.png')
fig.savefig('bluefrac_c10s3.pdf')


# In[ ]:

fig = plt.figure(figsize=(12,5))
gs = matplotlib.gridspec.GridSpec(1, 2, width_ratios=[1,0.05])
#ax1 = plt.subplot(gs[0])
#ax2 = plt.subplot(gs[1])
#ax3 = plt.subplot(gs[3])
#ax4 = plt.subplot(gs[4])
ax5 = plt.subplot(gs[0])
ax6 = plt.subplot(gs[1])
z1=data1['Z']
z2=data2['Z']
m1=data1['GRMU_R']
m2=data2['GRMU_R']
f1=data1['GRW_B']
f2=data2['N_B']*1.0/(data2['N_B']+data2['N_R'])
tl=['obs','truth','obs-truth']
idx,=np.where(np.abs(m2-m1)<0.1)
#im3=ax5.scatter(z2[idx],(f2-f1)[idx],c=(m2-m1)[idx],cmap='jet',s=100,vmin=-.1,vmax=.1)
im3=ax5.scatter(z2[idx],(f2-f1)[idx],c=(data2['THETAMAX'])[idx],cmap='jet',s=100)##,vmin=-.1,vmax=.1)

#ax1.set_ylim(0,1)
ax5.set_ylim(-1,1)
ax5.set_ylabel('delta blue fraction')
ax5.set_xlabel('redshift')
ax5.set_title('mock clusters')
#fig.colorbar(im1,cax=ax3,label='red sequence color',orientation='horizontal')
#fig.colorbar(im2,cax=ax4,label='red sequence color',orientation='horizontal')
fig.colorbar(im3,cax=ax6,label='THETAMAX')#,orientation='horizontal')
ax5.grid(True)
fig.tight_layout()
fig.savefig('bluefrac_c10c03.png')
fig.savefig('bluefrac_c10c03.pdf')


# In[ ]:

fig = plt.figure(figsize=(12,5))
gs = matplotlib.gridspec.GridSpec(1, 2, width_ratios=[1,0.05])
ax5 = plt.subplot(gs[0])
ax6 = plt.subplot(gs[1])
z1=data1['Z']
z2=data2['Z']
m1=data1['GRMU_R']
m2=data2['GRMU_R']
f1=data1['GRW_B']
f2=data2['N_B']*1.0/(data2['N_B']+data2['N_R'])
tl=['obs','truth','obs-truth']
idx,=np.where(np.abs(m2-m1)<20)
#im3=ax5.scatter(z2[idx],(f2-f1)[idx],c=(m2-m1)[idx],cmap='jet',s=100,vmin=-.1,vmax=.1)
##im3=ax5.scatter(z2[idx],(f2-f1)[idx],c=(data2['THETAMAX'])[idx],cmap='jet',s=100)##,vmin=-.1,vmax=.1)
im3=ax5.scatter(z2[idx],(m1-m2)[idx],c=(f2-f1)[idx],cmap='jet',s=100)##,vmin=-.1,vmax=.1)

#ax1.set_ylim(0,1)
#ax5.set_ylim(-1,1)
ax5.set_ylabel('delta red sequence color')
ax5.set_xlabel('redshift')
ax5.set_title('mock clusters')
#fig.colorbar(im1,cax=ax3,label='red sequence color',orientation='horizontal')
#fig.colorbar(im2,cax=ax4,label='red sequence color',orientation='horizontal')
fig.colorbar(im3,cax=ax6,label='delta blue fraction')#,orientation='horizontal')
ax5.grid(True)
fig.tight_layout()
fig.savefig('bluefrac_c10c03.png')
fig.savefig('bluefrac_c10c03.pdf')


# In[ ]:

fig = plt.figure(figsize=(5,5))
ax5 = plt.subplot()
z1=data1['Z']
z2=data2['Z']
m1=data1['GRMU_R']
m2=data2['GRMU_R']
f1=data1['GRW_B']
f2=data2['N_B']*1.0/(data2['N_B']+data2['N_R'])
tl=['obs','truth','obs-truth']
idx,=np.where(np.abs(m2-m1)<20)
im3=ax5.scatter(m2[idx],m1[idx],s=100) #,c=(f2-f1)[idx],cmap='jet',s=100)
ax5.set_ylabel('red sequence color (true)')
ax5.set_xlabel('red sequence color (obs)')
ax5.set_title('mock clusters')
#fig.colorbar(im3,cax=ax6,label='delta blue fraction')
#ax5.set_xlim(0,4)
ax5.grid(True)
fig.tight_layout()
fig.savefig('bluefrac_c10c03.png')
fig.savefig('bluefrac_c10c03.pdf')


# In[ ]:

fig = plt.figure(figsize=(12,10))
ax1 = plt.subplot(331)
ax2 = plt.subplot(332)
ax3 = plt.subplot(333)
ax4 = plt.subplot(334)
ax5 = plt.subplot(335)
ax6 = plt.subplot(336)
ax7 = plt.subplot(337)
ax8 = plt.subplot(338)
ax9 = plt.subplot(339)
z=data1['Z'][ix]
m=np.log10(data1['M200'][ix])
##m=m-10.*(data1['RI_SEP_FLAG'][ix])
plot_rs(ax1,z,data1['GRMU_R'][ix],m,tl='g-r',yl='peak color')
plot_rs(ax2,z,data1['RIMU_R'][ix],m,tl='r-i')
plot_rs(ax3,z,data1['IZMU_R'][ix],m,tl='i-z')
plot_rs(ax4,z,scatter(data1['GRSIGMA_R'][ix],sigma2phot['sigma^2(g-r)'][jx]),m,yl='width')
plot_rs(ax5,z,scatter(data1['RISIGMA_R'][ix],sigma2phot['sigma^2(r-i)'][jx]),m)
plot_rs(ax6,z,scatter(data1['IZSIGMA_R'][ix],sigma2phot['sigma^2(i-z)'][jx]),m)
plot_rs(ax7,z,data1['GR_SLOPE'][ix],m,xl='redshift',yl='slope')
plot_rs(ax8,z,data1['RI_SLOPE'][ix],m,xl='redshift')
plot_rs(ax9,z,data1['IZ_SLOPE'][ix],m,xl='redshift')
ax1.plot(r1color['z'].data,r1color['color'].data,c='m',lw=2)
ax4.plot(r1sigma['z'].data,r1sigma['sigma'].data,c='m',lw=2)
ax7.plot(r1slope['z'].data,r1slope['slope'].data,c='m',lw=2)
fig.tight_layout()
fig.savefig('rs_params_c10.png')
fig.savefig('rs_params_c10.pdf')


# In[ ]:

fig = plt.figure(figsize=(15,10))
gs = matplotlib.gridspec.GridSpec(2, 4, width_ratios=[1,1,1,0.1])
ax1 = plt.subplot(gs[0])#(241)
ax2 = plt.subplot(gs[1])#(242)
ax3 = plt.subplot(gs[2])#(243)
ax4 = plt.subplot(gs[4])#(245)
ax5 = plt.subplot(gs[5])#(246)
ax6 = plt.subplot(gs[6])#(247)
ax7 = plt.subplot(gs[3])#(244)
ax8 = plt.subplot(gs[7])#(248)
z=data1['Z'][ix]
m=np.log10(data1['M200'][ix])
m1=data1['RIW_B'][ix]
m2=data2['RIW_B'][ix]
#tl=['c = 3 , sigma clip n = 2','c = 3 , sigma clip n = 3']
plot_rs(ax1,z,data1['RIMU_R'][ix],m1,tl=tl[0],xl='redshift',yl='peak color')
plot_rs(ax2,z,scatter(data1['RISIGMA_R'][ix],sigma2phot['sigma^2(r-i)'][jx]),m1,tl=tl[0],xl='redshift',yl='width')
im=plot_rs(ax3,z,data1['RI_SLOPE'][ix],m1,tl=tl[0],xl='redshift',yl='slope')
plot_rs(ax4,z,data2['RIMU_R'][ix],m2,tl=tl[1],xl='redshift',yl='peak color')
plot_rs(ax5,z,scatter(data2['RISIGMA_R'][ix],sigma2phot['sigma^2(r-i)'][jx]),m2,tl=tl[1],xl='redshift',yl='width')
plot_rs(ax6,z,data2['RI_SLOPE'][ix],m2,tl=tl[1],xl='redshift',yl='slope')
#ax1.plot(r1color['z'].data,r1color['color'].data,c='m',lw=2)
#ax2.plot(r1sigma['z'].data,r1sigma['sigma'].data,c='m',lw=2)
#ax3.plot(r1slope['z'].data,r1slope['slope'].data,c='m',lw=2)
#ax4.plot(r1color['z'].data,r1color['color'].data,c='m',lw=2)
#ax5.plot(r1sigma['z'].data,r1sigma['sigma'].data,c='m',lw=2)
#ax6.plot(r1slope['z'].data,r1slope['slope'].data,c='m',lw=2)
ax1.set_ylim(0.3,1.4)
#ax2.set_ylim(0,1)
ax3.set_ylim(-0.3,0.3)
ax4.set_ylim(0.3,1.4)
#ax5.set_ylim(0,1)
ax6.set_ylim(-0.3,0.3)
fig.colorbar(im,cax=ax7,label='blue fraction')
fig.colorbar(im,cax=ax8,label='blue fraction')
fig.tight_layout()
fig.savefig('rs_params_c10c03.png')
fig.savefig('rs_params_c10c03.pdf')


# In[ ]:

fig = plt.figure(figsize=(10,5))
gs = matplotlib.gridspec.GridSpec(1, 3, width_ratios=[1,1,0.1])
ax1 = plt.subplot(gs[0])
ax2 = plt.subplot(gs[1])
ax3 = plt.subplot(gs[2])
z=data1['Z'][ix]
#m=np.log10(data1['M200'][ix])
m1=data1['RIMU_R'][ix]
m2=data2['RIMU_R'][ix]
ix1, = np.where(m1<2) 
ix2, = np.where(m2<2)
#tl=['c = 10','c = 3']
im=plot_rs(ax1,z[ix1],data1['RIW_B'][ix][ix1],m1[ix1],tl=tl[0],xl='redshift',yl='blue fraction',cmap='jet')
plot_rs(ax2,z[ix2],data2['RIW_B'][ix][ix2],m2[ix2],tl=tl[1],xl='redshift',yl='blue fraction',cmap='jet')
ax1.set_ylim(0,1)
ax2.set_ylim(0,1)
fig.colorbar(im,cax=ax3,label='red sequence color')
fig.tight_layout()
fig.savefig('bluefrac_c10c03.png')
fig.savefig('bluefrac_c10c03.pdf')


# In[ ]:

fig = plt.figure(figsize=(15,10))
gs = matplotlib.gridspec.GridSpec(2, 4, width_ratios=[1,1,1,0.1])
ax1 = plt.subplot(gs[0])#(241)
ax2 = plt.subplot(gs[1])#(242)
ax3 = plt.subplot(gs[2])#(243)
ax4 = plt.subplot(gs[4])#(245)
ax5 = plt.subplot(gs[5])#(246)
ax6 = plt.subplot(gs[6])#(247)
ax7 = plt.subplot(gs[3])#(244)
ax8 = plt.subplot(gs[7])#(248)
z=data1['Z'][ix]
m=np.log10(data1['M200'][ix])
m1=data1['IZW_B'][ix]
m2=data2['IZW_B'][ix]
#tl=['c = 3 , sigma clip n = 2','c = 3 , sigma clip n = 3']
plot_rs(ax1,z,data1['IZMU_R'][ix],m1,tl=tl[0],xl='redshift',yl='peak color')
plot_rs(ax2,z,scatter(data1['IZSIGMA_R'][ix],sigma2phot['sigma^2(i-z)'][jx]),m1,tl=tl[0],xl='redshift',yl='width')
im=plot_rs(ax3,z,data1['IZ_SLOPE'][ix],m1,tl=tl[0],xl='redshift',yl='slope')
plot_rs(ax4,z,data2['IZMU_R'][ix],m2,tl=tl[1],xl='redshift',yl='peak color')
plot_rs(ax5,z,scatter(data2['IZSIGMA_R'][ix],sigma2phot['sigma^2(i-z)'][jx]),m2,tl=tl[1],xl='redshift',yl='width')
plot_rs(ax6,z,data2['IZ_SLOPE'][ix],m2,tl=tl[1],xl='redshift',yl='slope')
#ax1.plot(r1color['z'].data,r1color['color'].data,c='m',lw=2)
#ax2.plot(r1sigma['z'].data,r1sigma['sigma'].data,c='m',lw=2)
#ax3.plot(r1slope['z'].data,r1slope['slope'].data,c='m',lw=2)
#ax4.plot(r1color['z'].data,r1color['color'].data,c='m',lw=2)
#ax5.plot(r1sigma['z'].data,r1sigma['sigma'].data,c='m',lw=2)
#ax6.plot(r1slope['z'].data,r1slope['slope'].data,c='m',lw=2)
#ax1.set_ylim(0.3,1.4)
ax2.set_ylim(0,1)
ax3.set_ylim(-0.3,0.3)
#ax4.set_ylim(0.3,1.4)
ax5.set_ylim(0,1)
ax6.set_ylim(-0.3,0.3)
fig.colorbar(im,cax=ax7,label='blue fraction')
fig.colorbar(im,cax=ax8,label='blue fraction')
fig.tight_layout()
fig.savefig('rs_params_c10c03.png')
fig.savefig('rs_params_c10c03.pdf')


# In[ ]:

fig = plt.figure(figsize=(10,5))
gs = matplotlib.gridspec.GridSpec(1, 3, width_ratios=[1,1,0.1])
ax1 = plt.subplot(gs[0])
ax2 = plt.subplot(gs[1])
ax3 = plt.subplot(gs[2])
z=data1['Z'][ix]
#m=np.log10(data1['M200'][ix])
m1=data1['IZMU_R'][ix]
m2=data2['IZMU_R'][ix]
ix1, = np.where(m1<2) 
ix2, = np.where(m2<2)
#tl=['c = 10','c = 3']
im=plot_rs(ax1,z[ix1],data1['IZW_B'][ix][ix1],m1[ix1],tl=tl[0],xl='redshift',yl='blue fraction',cmap='jet')
plot_rs(ax2,z[ix2],data2['IZW_B'][ix][ix2],m2[ix2],tl=tl[1],xl='redshift',yl='blue fraction',cmap='jet')
ax1.set_ylim(0,1)
ax2.set_ylim(0,1)
fig.colorbar(im,cax=ax3,label='red sequence color')
fig.tight_layout()
fig.savefig('bluefrac_c10c03.png')
fig.savefig('bluefrac_c10c03.pdf')


# In[ ]:

fig = plt.figure(figsize=(15,10))
gs = matplotlib.gridspec.GridSpec(3, 4, width_ratios=[1,1,1,0.1])
#g-r
ax1 = plt.subplot(gs[0])
ax2 = plt.subplot(gs[1],sharex=ax1)
ax3 = plt.subplot(gs[2],sharex=ax1)
ax4 = plt.subplot(gs[3])
#r-i
ax5 = plt.subplot(gs[4],sharex=ax1)
ax6 = plt.subplot(gs[5],sharex=ax1)
ax7 = plt.subplot(gs[6],sharex=ax1)
ax8 = plt.subplot(gs[7])
#i-z
ax9 = plt.subplot(gs[8],sharex=ax1)
ax10 = plt.subplot(gs[9],sharex=ax1)
ax11 = plt.subplot(gs[10],sharex=ax1)
ax12 = plt.subplot(gs[11])
#
d=data1
z=d['Z'][ix]
m1=d['GRW_B'][ix]
m2=d['RIW_B'][ix]
m3=d['IZW_B'][ix]
c1=d['GRMU_R'][ix]
c2=d['RIMU_R'][ix]
c3=d['IZMU_R'][ix]
w1=scatter(d['GRSIGMA_R'][ix],sigma2phot['sigma^2(g-r)'][jx])
w2=scatter(d['RISIGMA_R'][ix],sigma2phot['sigma^2(r-i)'][jx])
w3=scatter(d['IZSIGMA_R'][ix],sigma2phot['sigma^2(i-z)'][jx])
s1=d['GR_SLOPE'][ix]
s2=d['RI_SLOPE'][ix]
s3=d['IZ_SLOPE'][ix]
tl=['g-r','r-i','i-z']
lw=100
cmap='RdBu'
cmap='jet_r'
#
plot_rs(ax1,z,c1,m1,tl=tl[0],xl='redshift',yl='peak color',lw=lw,cmap=cmap)
plot_rs(ax2,z,w1,m1,tl=tl[0],xl='redshift',yl='width',lw=lw,cmap=cmap)
im1=plot_rs(ax3,z,s1,m1,tl=tl[0],xl='redshift',yl='slope',lw=lw,cmap=cmap)
plot_rs(ax5,z,c2,m2,tl=tl[1],xl='redshift',yl='peak color',lw=lw,cmap=cmap)
plot_rs(ax6,z,w2,m2,tl=tl[1],xl='redshift',yl='width',lw=lw,cmap=cmap)
im2=plot_rs(ax7,z,s2,m2,tl=tl[1],xl='redshift',yl='slope',lw=lw,cmap=cmap)
plot_rs(ax9,z,c3,m3,tl=tl[2],xl='redshift',yl='peak color',lw=lw,cmap=cmap)
plot_rs(ax10,z,w3,m3,tl=tl[2],xl='redshift',yl='width',lw=lw,cmap=cmap)
im3=plot_rs(ax11,z,s3,m3,tl=tl[2],xl='redshift',yl='slope',lw=lw,cmap=cmap)
#
#ax1.plot(r1color['z'].data,r1color['color'].data,c='m',lw=2)
#ax2.plot(r1sigma['z'].data,r1sigma['sigma'].data,c='m',lw=2)
#ax3.plot(r1slope['z'].data,r1slope['slope'].data,c='m',lw=2)
#
ax1.set_xlim(0.,1.)
ax1.set_ylim(0.6,3.6)
ax2.set_ylim(0,0.6)
ax3.set_ylim(-0.3,0.3)
ax5.set_ylim(0.2,2.2)
ax6.set_ylim(0,0.6)
ax7.set_ylim(-0.3,0.3)
ax9.set_ylim(0.1,1.1)
ax10.set_ylim(0,0.6)
ax11.set_ylim(-0.3,0.3)
#
fig.colorbar(im1,cax=ax4,label='blue fraction')
fig.colorbar(im2,cax=ax8,label='blue fraction')
fig.colorbar(im3,cax=ax12,label='blue fraction')
fig.tight_layout()
fig.savefig('rs_params_c03.png')
fig.savefig('rs_params_c03.pdf')


# In[ ]:

fig = plt.figure(figsize=(15,5))
gs = matplotlib.gridspec.GridSpec(2, 3, height_ratios=[1,0.1])
ax1 = plt.subplot(gs[0])
ax2 = plt.subplot(gs[1],sharex=ax1,sharey=ax1)
ax3 = plt.subplot(gs[2],sharex=ax1,sharey=ax1)
ax4 = plt.subplot(gs[3])
ax5 = plt.subplot(gs[4])
ax6 = plt.subplot(gs[5])

#
d=data1
z=d['Z'][ix]
m1=d['GRW_B'][ix]
m2=d['RIW_B'][ix]
m3=d['IZW_B'][ix]
c1=d['GRMU_R'][ix]
c2=d['RIMU_R'][ix]
c3=d['IZMU_R'][ix]
tl=['g-r','r-i','i-z']
lw=100
ix1, = np.where(m1<20) 
ix2, = np.where(m2<20)
ix3, = np.where(m3<20)
cmap='jet'
im1=plot_rs(ax1,z[ix1],m1[ix1],c1[ix1],tl=tl[0],xl='redshift',yl='blue fraction',cmap=cmap,lw=lw)
im2=plot_rs(ax2,z[ix2],m2[ix2],c2[ix2],tl=tl[1],xl='redshift',yl='blue fraction',cmap=cmap,lw=lw)
im3=plot_rs(ax3,z[ix3],m3[ix3],c3[ix3],tl=tl[2],xl='redshift',yl='blue fraction',cmap=cmap,lw=lw)
ax1.set_ylim(0,1)
fig.colorbar(im1,cax=ax4,label='red sequence color (g-r)',orientation='horizontal')
fig.colorbar(im2,cax=ax5,label='red sequence color (r-i)',orientation='horizontal')
fig.colorbar(im3,cax=ax6,label='red sequence color (i-z)',orientation='horizontal')
fig.tight_layout()
fig.savefig('bluefrac_c03.png')
fig.savefig('bluefrac_c03.pdf')


# In[ ]:




# In[ ]:

print data1['MEM_MATCH_ID'][0:10]
print data2['MEM_MATCH_ID'][0:10]
#print data1['M200'][0:10]
print sigma2phot


# In[ ]:

cluster_hdulist_2[1].columns


# In[ ]:

cluster_hdulist_1[1].columns


# In[ ]:

#match on id
x=data1['MEM_MATCH_ID']
goodvalues=sigma2phot['haloid']
ix = np.where(np.in1d(x.ravel(), goodvalues).reshape(x.shape))
jx = np.where(np.in1d(goodvalues.ravel(), x).reshape(goodvalues.shape))


# In[ ]:

data1['MEM_MATCH_ID'][ix][:5]


# In[ ]:

sigma2phot['haloid'][jx][:5]


# In[ ]:

def scatter(sigmaobs,sigma2phot):
    x=sigmaobs**2-sigma2phot
    return np.sqrt(x)


# In[ ]:

def intrinsic_scatter(sigmaobs,sigma2phot,redshift):
    x=sigmaobs**2-sigma2phot['sigma^2(g-r)']
    idx=redshift>=0.3
    x[idx]=(sigmaobs**2-sigma2phot['sigma^2(r-i)'])[idx]
    idx=redshift>=0.7
    x[idx]=(sigmaobs**2-sigma2phot['sigma^2(i-z)'])[idx]
    return np.sqrt(x)


# In[ ]:

x=data1['Z'][ix]
y3=intrinsic_scatter(data1['RESTSIGMA_R'][ix],sigma2phot[jx],data1['Z'][ix])
y1=data1['RESTSIGMA_R'][ix]
plt.scatter(x,y1,c='orange',lw=0,s=60,label='original fit, observed')
plt.scatter(x,y3,c='b',lw=0,label='original fit, intrinsic')
plt.xlabel('z')
plt.legend(loc='upper left')
plt.ylabel('red sequence scatter (rest frame colors)')
#plt.ylim(-0.1,0.6)
plt.show()


# In[ ]:

x=data2['Z'][ix]
y3=intrinsic_scatter(data2['RESTSIGMA_R'][ix],sigma2phot[jx],data2['Z'][ix])
y1=data2['RESTSIGMA_R'][ix]
plt.scatter(x,y1,c='orange',lw=0,s=60,label='2nd fit, observed')
plt.scatter(x,y3,c='b',lw=0,label='2nd fit, intrinsic')
plt.xlabel('z')
plt.legend(loc='upper left')
plt.ylabel('red sequence scatter (rest frame colors)')
plt.ylim(0,0.2)
plt.show()


# In[ ]:

x=data1['Z'][ix]
y3=scatter(data1['RISIGMA_R'][ix],sigma2phot['sigma^2(r-i)'][jx])
y2=sigma2phot['sigma^2(r-i)'][jx]
y1=data1['RISIGMA_R'][ix]
plt.scatter(x,y1,c='b',alpha=0.3,s=30)
plt.scatter(x,y3,c='r')
plt.scatter(x,y2, c='g')
plt.xlabel('z')
plt.ylabel('red sequence scatter (r-i)')
plt.show()


# In[ ]:

x=data1['Z'][ix]
y3=scatter(data1['GRSIGMA_R'][ix],sigma2phot['sigma^2(g-r)'][jx])
y2=sigma2phot['sigma^2(g-r)'][jx]
y1=data1['GRSIGMA_R'][ix]
#plt.scatter(x,y1,c='b',alpha=0.3,s=30)
plt.scatter(x,y3,c='r')
#plt.scatter(x,y2, c='g')
plt.xlabel('z')
plt.ylabel('red sequence scatter (g-r)')
plt.show()


# In[ ]:

# x=data1['Z']
y=data1['GRMU_B']
plt.scatter(x, y, c=m1,cmap='Blues',lw=0,s=40)
plt.xlabel('z')
plt.ylabel('blue cloud color (i-z)')
plt.show()


# In[ ]:

idx=data1['GR_SEP_FLAG']==True
x2=data2['Z']
y2=data2['GRMU_R']
x1=data1['Z']
y1=data1['GRMU_R']
m1=np.log10(data1['M200'])
x1t=data1['Z'][idx]
y1t=data1['GRMU_R'][idx]
y3=data1['GRMU_B']
y4=data1['GRMU_BG']
r1x=r1color['z'].data
r1y=r1color['color'].data
plt.scatter(x1, y1, c='w',alpha=0.2,s=40)
#plt.scatter(x1t, y1t, c='b',alpha=0.2)
plt.scatter(x2, y2, c=m1,cmap='Reds',lw=0,s=40)
plt.scatter(x2, y3, c=m1,cmap='Blues',lw=0,s=40)
plt.scatter(x2, y4, c=m1,cmap='Greys',lw=0,s=40)
plt.plot(r1x,r1y,c='m',lw=2)
plt.xlabel('z')
plt.ylabel('color (g-r)')
plt.xlim(0,1)
##plt.ylim(0.7,2.4)
#plt.colorbar(label='log (M/h)')
plt.show()


# In[ ]:

idx=data1['GR_SEP_FLAG']==True
x2=data2['Z']
y2=data2['GRMU_R']
x1=data1['Z']
y1=data1['GRMU_R']
m1=np.log10(data1['M200'])
x1t=data1['Z'][idx]
y1t=data1['GRMU_R'][idx]
r1x=r1color['z'].data
r1y=r1color['color'].data
plt.scatter(x1, y1, c='w',alpha=0.2,s=40)
#plt.scatter(x1t, y1t, c='b',alpha=0.2)
plt.scatter(x2, y2, c=m1,cmap='Reds',lw=0,s=40)
plt.plot(r1x,r1y,c='m',lw=2)
plt.xlabel('z')
plt.ylabel('color (g-r)')
plt.xlim(0,1)
##plt.ylim(0.7,2.4)
plt.colorbar(label='log (M/h)')
plt.show()


# In[ ]:

x2=data2['Z']
y2=data2['GR_SLOPE']
x1=data1['Z']
y1=data1['GR_SLOPE']
m1=np.log10(data1['M200'])
r1x=r1slope['z'].data
r1y=r1slope['slope'].data
plt.scatter(x1, y1, c='w',alpha=0.2,s=40)
plt.scatter(x2, y2, c=m1,cmap='Reds',lw=0,s=40)
plt.plot(r1x,r1y,c='m',lw=2)
plt.xlabel('z')
plt.ylabel('slope (g-r)')
plt.xlim(0,1)
plt.ylim(-0.1,0.08)
plt.show()


# In[ ]:

x2=data2['Z']
y2=data2['GRSIGMA_R']
x1=data1['Z']
y1=data1['GRSIGMA_R']
y3=scatter(data2['GRSIGMA_R'],sigma2phot['sigma^2(g-r)'][ix])
r1x=r1sigma['z'].data
r1y=r1sigma['sigma'].data
##plt.scatter(x1, y2, c='w',alpha=0.2,s=40)
plt.scatter(x2, y3, c=m1,cmap='Reds',lw=0,s=40)
#plt.scatter(x1, y3, c='r',s=10,lw=0)
plt.plot(r1x,r1y,c='m',lw=2)
plt.xlabel('z')
plt.ylabel('sigma (g-r)')
plt.xlim(0,1)
plt.ylim(-0.1,0.6)
plt.show()


# In[ ]:

x2=data2['Z']
y2=data2['RESTMU_R']
x1=data1['Z']
y1=data1['RESTMU_R']
m1=np.log10(data1['M200'])
x1t=data1['Z'][idx]
y1t=data1['RESTMU_R'][idx]
r1x=r1color['z'].data
r1y=r1color['color'].data
plt.scatter(x1, y1, c='w',alpha=0.2,s=20)
#plt.scatter(x1t, y1t, c='b',alpha=0.2)
plt.scatter(x2, y2, c=m1,cmap='Reds',lw=0,s=20)
##plt.plot(r1x,r1y,c='m',lw=2)
plt.xlabel('z')
plt.ylabel('color rest frame (g-r)')
plt.xlim(0,1)
##plt.ylim(0.7,2.4)
plt.colorbar(label='log (M/h)')
plt.show()


# In[ ]:

x2=data2['Z']
y2=data2['RESTMU_R']
x1=data1['Z']
y1=data1['RESTMU_R']
y3=data1['RESTMU_B']
y4=data1['RESTMU_BG']
m1=np.log10(data1['M200'])
x1t=data1['Z'][idx]
y1t=data1['RESTMU_R'][idx]
r1x=r1color['z'].data
r1y=r1color['color'].data
plt.scatter(x1, y1, c='w',alpha=0.2,s=20)
#plt.scatter(x1t, y1t, c='b',alpha=0.2)
plt.scatter(x2, y2, c=m1,cmap='Reds',lw=0,s=20)
plt.scatter(x2, y3, c=m1,cmap='Blues',lw=0,s=20)
plt.scatter(x2, y4, c=m1,cmap='Greys',lw=0,s=20)
##plt.plot(r1x,r1y,c='m',lw=2)
plt.xlabel('z')
plt.ylabel('color rest frame (g-r)')
plt.xlim(0,1)
##plt.ylim(0.7,2.4)
#plt.colorbar(label='log (M/h)')
plt.show()


# In[ ]:

x2=data2['Z']
y2=data2['RIMU_R']
x1=data1['Z']
y1=data1['RIMU_R']
y3=data1['RIMU_B']
y4=data1['RIMU_BG']
m1=np.log10(data1['M200'])
x1t=data1['Z'][idx]
y1t=data1['RIMU_R'][idx]
r1x=r1color['z'].data
r1y=r1color['color'].data
plt.scatter(x1, y1, c='w',alpha=0.2,s=20)
#plt.scatter(x1t, y1t, c='b',alpha=0.2)
plt.scatter(x2, y2, c=m1,cmap='Reds',lw=0,s=20)
plt.scatter(x2, y3, c=m1,cmap='Blues',lw=0,s=20)
plt.scatter(x2, y4, c=m1,cmap='Greys',lw=0,s=20)
##plt.plot(r1x,r1y,c='m',lw=2)
plt.xlabel('z')
plt.ylabel('color (r-i)')
plt.xlim(0,1)
plt.ylim(0.,1)
#plt.colorbar(label='log (M/h)')
plt.show()


# In[ ]:

x2=data2['Z']
y2=data2['IZMU_R']
x1=data1['Z']
y1=data1['IZMU_R']
y3=data1['IZMU_B']
y4=data1['IZMU_BG']
m1=np.log10(data1['M200'])
x1t=data1['Z'][idx]
y1t=data1['IZMU_R'][idx]
r1x=r1color['z'].data
r1y=r1color['color'].data
plt.scatter(x1, y1, c='w',alpha=0.2,s=20)
#plt.scatter(x1t, y1t, c='b',alpha=0.2)
plt.scatter(x2, y2, c=m1,cmap='Reds',lw=0,s=20)
plt.scatter(x2, y3, c=m1,cmap='Blues',lw=0,s=20)
plt.scatter(x2, y4, c=m1,cmap='Greys',lw=0,s=20)
##plt.plot(r1x,r1y,c='m',lw=2)
plt.xlabel('z')
plt.ylabel('color (i-z)')
plt.xlim(0,1)
plt.ylim(0.,1)
#plt.colorbar(label='log (M/h)')
plt.show()


# In[ ]:

x2=data2['Z']
y2=data2['REST_SLOPE']
x1=data1['Z']
y1=data1['REST_SLOPE']
m1=np.log10(data1['M200'])
x1t=data1['Z'][idx]
y1t=data1['REST_SLOPE'][idx]
r1x=r1slope['z'].data
r1y=r1slope['slope'].data
plt.scatter(x1, y1, c='w',alpha=0.2,s=20)
#plt.scatter(x1t, y1t, c='b',alpha=0.2)
plt.scatter(x2, y2, c=m1,cmap='Reds',lw=0,s=20)
##plt.plot(r1x,r1y,c='m',lw=2)
plt.xlabel('z')
plt.ylabel('slope, rest frame (g-r)')
plt.xlim(0,1)
##plt.ylim(0.7,2.4)
plt.colorbar(label='log (M/h)')
plt.show()


# In[ ]:

plt.p


# In[ ]:

y3=intrinsic_scatter(data2['RESTSIGMA_R'],sigma2phot[ix],data2['Z'])
y4=intrinsic_scatter(data1['RESTSIGMA_R'],sigma2phot[ix],data1['Z'])
x2=data2['Z']
y2=data2['RESTSIGMA_R']
x1=data1['Z']
y1=data1['RESTSIGMA_R']
m1=np.log10(data1['M200'])
x1t=data1['Z'][idx]
y1t=data1['RESTSIGMA_R'][idx]
r1x=r1color['z'].data
r1y=r1color['color'].data
#plt.scatter(x1, y1, c='orange',s=60,lw=0,label='original fit, observed')
#plt.scatter(x1t, y1t, c='b',alpha=0.2)
#plt.scatter(x2, y2, c=m1,cmap='Reds',lw=0,s=20)
plt.scatter(x2, y2, c='gray',lw=0,s=60,label='new fit, observed')
plt.scatter(x2, y3, c='yellow',lw=1,s=20,label='new fit, intrinsic')
#plt.scatter(x2, y4, c='b',s=20,lw=0,label='original, intrinsic')
##plt.plot(r1x,r1y,c='m',lw=2)
plt.xlabel('z')
plt.ylabel('sigma, rest frame')
plt.xlim(0,1)
plt.ylim(0.,0.35)
plt.legend(loc='upper left')
#plt.colorbar(label='log (M/h)')
plt.show()


# In[ ]:

#x2=data2['Z']
#y2=data2['GRW_R']
x1=data1['Z']
y1=data1['RESTW_R'] ##/(data1['GRW_R']+data1['GRW_B'])
y2=data1['RESTW_B']
y3=data1['RESTW_BG']
m1=np.log10(data1['M200'])
plt.scatter(x1, y3, c=m1,cmap='Greys',lw=0,s=20)
plt.scatter(x1, y1, c=m1,cmap='Reds',lw=0,s=20)
plt.scatter(x1, y2, c=m1,cmap='Blues',lw=0,s=20)
##plt.plot(r1x,r1y,c='m',lw=2)
plt.xlabel('z')
plt.ylabel('weights (rest frame)')
plt.xlim(0,1)
##plt.ylim(0.7,2.4)
##plt.colorbar(label='log (M/h)')
plt.show()


# In[ ]:

#x2=data2['Z']
#y2=data2['GRW_R']
x1=data1['Z']
y1=data1['GRW_R'] ##/(data1['GRW_R']+data1['GRW_B'])
y2=data1['GRW_B']
y3=data1['GRW_BG']
m1=np.log10(data1['M200'])
m1=m1+(data1['GR_SEP_FLAG'][ix])
plt.scatter(x1, y3, c=m1,cmap='Greys',lw=0,s=20)
plt.scatter(x1, y1, c=m1,cmap='Reds',lw=0,s=20)
plt.scatter(x1, y2, c=m1,cmap='Blues',lw=0,s=20)
##plt.plot(r1x,r1y,c='m',lw=2)
plt.xlabel('z')
plt.ylabel('weights (g-r)')
plt.xlim(0,1)
##plt.ylim(0.7,2.4)
##plt.colorbar(label='log (M/h)')
plt.show()


# In[ ]:

#x2=data2['Z']
#y2=data2['GRW_R']
x1=data1['Z']
y1=data1['RIW_R'] ##/(data1['GRW_R']+data1['GRW_B'])
y2=data1['RIW_B']
y3=data1['RIW_BG']
m1=np.log10(data1['M200'])
m1=m1+(data1['RI_SEP_FLAG'][ix])
plt.scatter(x1, y3, c=m1,cmap='Greys',lw=0,s=20)
plt.scatter(x1, y1, c=m1,cmap='Reds',lw=0,s=20)
plt.scatter(x1, y2, c=m1,cmap='Blues',lw=0,s=20)
##plt.plot(r1x,r1y,c='m',lw=2)
plt.xlabel('z')
plt.ylabel('weights (r-i)')
plt.xlim(0,1)
##plt.ylim(0.7,2.4)
##plt.colorbar(label='log (M/h)')
plt.show()


# In[ ]:

#x2=data2['Z']
#y2=data2['GRW_R']
x1=data1['Z']
y1=data1['IZW_R'] ##/(data1['GRW_R']+data1['GRW_B'])
y2=data1['IZW_B']
y3=data1['IZW_BG']
m1=np.log10(data1['R200'])
m1=m1+(data1['IZ_SEP_FLAG'][ix])
plt.scatter(x1, y3, c=m1,cmap='Greys',lw=0,s=20)
plt.scatter(x1, y1, c=m1,cmap='Reds',lw=0,s=20)
plt.scatter(x1, y2, c=m1,cmap='Blues',lw=0,s=20)
##plt.plot(r1x,r1y,c='m',lw=2)
plt.xlabel('z')
plt.ylabel('weights (i-z)')
plt.xlim(0,1)
##plt.ylim(0.7,2.4)
##plt.colorbar(label='log (M/h)')
plt.show()


# In[ ]:

x1=data1['Z']
y1=data1['GRW_B'] 
y2=(data1['GRW_B']+data1['GRW_R'])
y3=y1/y2
m1=np.log10(data1['M200'])
m1=m1+(data1['GR_SEP_FLAG'][ix])
plt.scatter(x1, y3, c=m1,cmap='Purples',lw=0,s=40)
plt.xlabel('z')
plt.ylabel('$f_b$ (g-r)')
plt.show()


# In[ ]:

x1=data1['Z']
y1=data1['RIW_B'] 
y2=(data1['RIW_B']+data1['RIW_R'])
y3=y1/y2
m1=np.log10(data1['M200'])
m1=m1+(data1['RI_SEP_FLAG'][ix])
plt.scatter(x1, y3, c=m1,cmap='Purples',lw=0,s=40)
plt.xlabel('z')
plt.ylabel('$f_b$ (r-i)')
plt.show()


# In[ ]:

x1=data1['Z']
y1=data1['IZW_B'] ##/(data1['GRW_R']+data1['GRW_B'])
y2=(data1['IZW_B']+data1['IZW_R'])
y3=y1/y2
m1=np.log10(data1['M200'])
m1=m1+(data1['IZ_SEP_FLAG'][ix])
plt.scatter(x1, y3, c=m1,cmap='Purples',lw=0,s=40)
plt.xlabel('z')
plt.ylabel('$f_b$ (i-z)')
plt.show()


# In[ ]:

#x2=data2['Z']
#y2=data2['GRW_R']
x1=data1['Z']
y1=data1['RESTW_R'] ##/(data1['GRW_R']+data1['GRW_B'])
y2=1-data1['RESTW_BG'] ####(data1['RESTW_B']+data1['RESTW_R'])
y3=1-y1/y2
m1=np.log10(data1['M200'])
#x1t=data1['Z'][idx]
#y1t=data1['RESTSIGMA_R'][idx]
#r1x=r1color['z'].data
#r1y=r1color['color'].data
#plt.scatter(x1, y1, c='w',alpha=0.2,s=20)
#plt.scatter(x1t, y1t, c='b',alpha=0.2)
plt.scatter(x1, y3, c=m1,cmap='Purples',lw=0,s=20)
#plt.scatter(x1, y1, c=m1,cmap='Reds',lw=0,s=20)
#plt.scatter(x1, y2, c=m1,cmap='Blues',lw=0,s=20)
##plt.plot(r1x,r1y,c='m',lw=2)
plt.xlabel('z')
plt.ylabel('$f_b$ (rest frame)')
#plt.xlim(0,1)
##plt.ylim(0.7,2.4)
##plt.colorbar(label='log (M/h)')
plt.show()


# In[ ]:

x2=data1['Z']
#y2=data2['GRW_R']
x1=data1['DEC']
y1=data1['RA'] ##/(data1['GRW_R']+data1['GRW_B'])
#y2=(data1['IZW_B']+data1['IZW_R'])
#y3=y1/y2
#m1=np.log10(data1['M200'])
#x1t=data1['Z'][idx]
#y1t=data1['RESTSIGMA_R'][idx]
#r1x=r1color['z'].data
#r1y=r1color['color'].data
#plt.scatter(x1, y1, c='w',alpha=0.2,s=20)
#plt.scatter(x1t, y1t, c='b',alpha=0.2)
plt.scatter(x1, y1, c=x2,lw=0,s=20)
#plt.scatter(x1, y1, c=m1,cmap='Reds',lw=0,s=20)
#plt.scatter(x1, y2, c=m1,cmap='Blues',lw=0,s=20)
##plt.plot(r1x,r1y,c='m',lw=2)
plt.xlabel('DEC')
plt.ylabel('RA')
#plt.xlim(0,1)
##plt.ylim(0.7,2.4)
plt.colorbar(label='redshfit')
plt.show()


# In[ ]:

x2=data1['Z']
x1=np.log10(data1['R200']) #*(1+data1['Z'])**0.5)
y1=np.log10(data1['M200'])###/(1+data1['Z'])**1.5) 
m1=np.log10(data1['M200'])
plt.scatter(10**x1, 10**y1, c=x2,lw=0,s=20)
plt.xlabel('$R_{200}$')
plt.ylabel('$M_{200}$')
#plt.xlim(0.55,2.55)
#plt.ylim(1.4*10**13,1.4*10**15)
plt.xscale('log')
plt.yscale('log')
plt.colorbar(label='redshift')
plt.show()


# In[ ]:

y1=data1['IZW_B'] 
y2=(data1['IZW_B']+data1['IZW_R'])
y3=y1/y2
x2=data1['Z']
x1=data1['M200']  #/(1+data1['Z'])**1.5
m1=(data1['IZ_SEP_FLAG'][ix])
##x1=x1*m1
r1x=np.zeros(2)
r1y=np.zeros(2)
r1x[0]=5.0*10**14
r1x[1]=1.0*10**14
r1y[0]=0.05
r1y[1]=0.1
plt.scatter(x1, y3, c=x2,lw=0,s=20)
plt.plot(r1x,r1y,c='k',lw=2)
plt.xlabel('$M_{200}/(1+z)^{3/2}$')
plt.xlabel('$M_{200}$')
plt.ylabel('$f_{b}$ (i-z)')
plt.ylim(-0.1,1.1)
plt.xlim(0.3*10**14,1.1*10**15)
plt.xscale('log')
plt.colorbar(label='redshift')
plt.show()


# In[ ]:

y1=data1['IZW_B'] 
y2=(data1['IZW_B']+data1['IZW_R'])
y3=y1/y2
x2=data1['Z']
x1=data1['LAMBDA_CHISQ']
m1=(data1['IZ_SEP_FLAG'][ix])
##x1=x1*m1
plt.scatter(x1, y3, c=x2,lw=0,s=20)
plt.ylabel('$f_{b}$ (i-z)')
plt.xlabel('$\lambda$')
plt.ylabel('$f_{b}$')
#plt.ylim(0.,1.)
plt.xscale('log')
#plt.yscale('log')
plt.colorbar(label='redshift')
plt.show()


# In[ ]:

x2=data1['Z']
x1=np.log10(data1['LAMBDA_CHISQ']) 
y1=np.log10(data1['M200'])###/(1+data1['Z'])**1.5) 
m1=np.log10(data1['M200'])
plt.scatter(10**x1, 10**y1, c=x2,lw=0,s=20)
plt.xlabel('$\lambda$')
plt.ylabel('$M_{200}$ (new, fixed z window)')
#plt.xlim(0.55,2.55)
#plt.ylim(1.4*10**13,1.4*10**15)
plt.xscale('log')
plt.yscale('log')
plt.colorbar(label='redshift')
plt.show()


# In[ ]:

x2=data1['Z']
#y2=data2['GRW_R']
x1=data1['R200']
y1=data1['RESTW_B'] ##/(data1['GRW_R']+data1['GRW_B'])
y2=(data1['RESTW_B']+data1['RESTW_R'])
y3=y1/y2
m1=np.log10(data1['M200'])
#x1t=data1['Z'][idx]
#y1t=data1['RESTSIGMA_R'][idx]
#r1x=r1color['z'].data
#r1y=r1color['color'].data
#plt.scatter(x1, y1, c='w',alpha=0.2,s=20)
#plt.scatter(x1t, y1t, c='b',alpha=0.2)
plt.scatter(x1, y3, c=x2,lw=0,s=20)
#plt.scatter(x1, y1, c=m1,cmap='Reds',lw=0,s=20)
#plt.scatter(x1, y2, c=m1,cmap='Blues',lw=0,s=20)
##plt.plot(r1x,r1y,c='m',lw=2)
plt.xlabel('$R_{200}$')
plt.ylabel('$f_{b}$')
plt.ylim(0.,1.)
#plt.xlim(0.8*10**14,1.4*10**15)
plt.xscale('linear')
#plt.yscale('log')
plt.colorbar(label='redshift')
plt.show()


# In[ ]:

x2=data1['Z']
#y2=data2['GRW_R']
x1=data1['M200']/(1+data1['Z'])**1.5
y1=data1['RIW_R'] ##/(data1['GRW_R']+data1['GRW_B'])
y2=1-data1['RIW_BG']
y3=1-y1/y2
m1=np.log10(data1['M200'])
#x1t=data1['Z'][idx]
#y1t=data1['RESTSIGMA_R'][idx]
r1x=np.zeros(2)
r1y=np.zeros(2)
r1x[0]=5.0*10**14
r1x[1]=1.0*10**14
r1y[0]=0.05
r1y[1]=0.1
#plt.scatter(x1, y1, c='w',alpha=0.2,s=20)
#plt.scatter(x1t, y1t, c='b',alpha=0.2)
#plt.scatter(x1[ix], y1[ix], c=x2[ix],lw=0,s=20)
plt.scatter(x1, y3, c=x2,lw=0,s=20)
#plt.scatter(x1, y2, c=m1,cmap='Blues',lw=0,s=20)
plt.plot(r1x,r1y,c='k',lw=2)
plt.xlabel('$M_{200}/(1+z)^{3/2}$')
plt.ylabel('$f_{b} (r-i)$')
plt.ylim(0.01,1.01)
plt.xlim(0.3*10**14,1.4*10**15)
plt.xscale('log')
#plt.yscale('log')
plt.colorbar(label='redshift')
plt.show()


# In[ ]:

x2=data1['Z']
#y2=data2['GRW_R']
x1=data1['M200']/(1+data1['Z'])**1.5
y1=data1['IZW_R'] ##/(data1['GRW_R']+data1['GRW_B'])
y2=1-data1['IZW_BG']
y3=1-y1/y2
m1=np.log10(data1['M200'])
#x1t=data1['Z'][idx]
#y1t=data1['RESTSIGMA_R'][idx]
r1x=np.zeros(2)
r1y=np.zeros(2)
r1x[0]=5.0*10**14
r1x[1]=1.0*10**14
r1y[0]=0.05
r1y[1]=0.1
#plt.scatter(x1, y1, c='w',alpha=0.2,s=20)
#plt.scatter(x1t, y1t, c='b',alpha=0.2)
#plt.scatter(x1[ix], y1[ix], c=x2[ix],lw=0,s=20)
plt.scatter(x1, y3, c=x2,lw=0,s=20)
#plt.scatter(x1, y2, c=m1,cmap='Blues',lw=0,s=20)
plt.plot(r1x,r1y,c='k',lw=2)
plt.xlabel('$M_{200}/(1+z)^{3/2}$')
plt.ylabel('$f_{b} (i-z)$')
plt.ylim(0.01,1.01)
plt.xlim(0.3*10**14,1.4*10**15)
plt.xscale('log')
#plt.yscale('log')
plt.colorbar(label='redshift')
plt.show()


# In[ ]:

x2=data1['Z']
#y2=data2['GRW_R']
x1=data1['M200']/(1+data1['Z'])**1.5
y1=data1['GRW_B'] ##/(data1['GRW_R']+data1['GRW_B'])
y2=1-data1['GRW_BG']
y3=1-y1/y2
m1=np.log10(data1['M200'])
#x1t=data1['Z'][idx]
#y1t=data1['RESTSIGMA_R'][idx]
r1x=np.zeros(2)
r1y=np.zeros(2)
r1x[0]=5.0*10**14
r1x[1]=1.0*10**14
r1y[0]=0.05
r1y[1]=0.1
#plt.scatter(x1, y1, c='w',alpha=0.2,s=20)
#plt.scatter(x1t, y1t, c='b',alpha=0.2)
#plt.scatter(x1[ix], y1[ix], c=x2[ix],lw=0,s=20)
plt.scatter(x1, y3, c=x2,lw=0,s=20)
#plt.scatter(x1, y2, c=m1,cmap='Blues',lw=0,s=20)
plt.plot(r1x,r1y,c='k',lw=2)
plt.xlabel('$M_{200}/(1+z)^{3/2}$')
plt.ylabel('$f_{b} (g-r)$')
plt.ylim(0.01,1.01)
plt.xlim(0.3*10**14,1.4*10**15)
plt.xscale('log')
#plt.yscale('log')
plt.colorbar(label='redshift')
plt.show()


# In[ ]:

x2=data1['Z']
#y2=data2['GRW_R']
x1=data1['M200']/(1+data1['Z'])**1.5
y1=data1['RESTW_R'] ##/(data1['GRW_R']+data1['GRW_B'])
y2=1-data1['RESTW_BG']
y3=1-y1/y2
m1=np.log10(data1['M200'])
#x1t=data1['Z'][idx]
#y1t=data1['RESTSIGMA_R'][idx]
r1x=np.zeros(2)
r1y=np.zeros(2)
r1x[0]=5.0*10**14
r1x[1]=1.0*10**14
r1y[0]=0.05
r1y[1]=0.1
#plt.scatter(x1, y1, c='w',alpha=0.2,s=20)
#plt.scatter(x1t, y1t, c='b',alpha=0.2)
#plt.scatter(x1[ix], y1[ix], c=x2[ix],lw=0,s=20)
plt.scatter(x1, y3, c=x2,lw=0,s=20)
#plt.scatter(x1, y2, c=m1,cmap='Blues',lw=0,s=20)
plt.plot(r1x,r1y,c='k',lw=2)
plt.xlabel('$M_{200}/(1+z)^{3/2}$')
plt.ylabel('$f_{b}$ (rest frame)')
plt.ylim(0.01,1.5)
plt.xlim(0.3*10**14,1.4*10**15)
plt.xscale('log')
#plt.yscale('log')
plt.colorbar(label='redshift')
plt.show()


# In[ ]:

x2=data1['Z']
#y2=data2['GRW_R']
x1=data1['M200']#/(1+data1['Z'])**1.5
y1=data1['RESTW_R'] ##/(data1['GRW_R']+data1['GRW_B'])
y2=1-data1['RESTW_BG']#(data1['IZW_B']+data1['IZW_R'])
y3=1-y1/y2
m1=np.log10(data1['M200'])
#x1t=data1['Z'][idx]
#y1t=data1['RESTSIGMA_R'][idx]
r1x=np.zeros(2)
r1y=np.zeros(2)
r1x[0]=5.0*10**14
r1x[1]=1.0*10**14
r1y[0]=0.05
r1y[1]=0.1
#plt.scatter(x1, y1, c='w',alpha=0.2,s=20)
#plt.scatter(x1t, y1t, c='b',alpha=0.2)
plt.scatter(x1, y3, c=x2,lw=0,s=20)
#plt.scatter(x1, y1, c=m1,cmap='Reds',lw=0,s=20)
#plt.scatter(x1, y2, c=m1,cmap='Blues',lw=0,s=20)
plt.plot(r1x,r1y,c='k',lw=2)
plt.xlabel('$M_{200}/(1+z)^{3/2}$')
plt.ylabel('$f_{b}$ (i-z)')
plt.ylim(-0.1,1.1)
plt.xlim(0.3*10**14,1.1*10**15)
plt.xscale('log')
#plt.yscale('log')
plt.colorbar(label='redshift')
plt.show()


# In[ ]:

import random


# In[ ]:

np.arange(0.,1.,3.)


# In[ ]:




# In[ ]:

x = np.linspace(0.,2.,100)
ngals = 3
redshift_err = random.sample(np.linspace(0.001,0.03, num=100), ngals)
redshift = random.sample(np.linspace(0.,1.0, num=100), ngals)
print redshift
print redshift_err
def pdf(x,x0,sigma): 
    return np.exp(-(x-x0)**2/(2*sigma))/np.sqrt(2*np.pi)

plt.plot(x, pdf(x,x0=redshift[0],sigma=redshift_err[0]))
plt.plot(x, pdf(x,x0=redshift[1],sigma=redshift_err[1]))
plt.plot(x, pdf(x,x0=redshift[2],sigma=redshift_err[2]))

#print pdf(x,x0=redshift[1],sigma=redshift_err[1])

y=x*0.0
for i in np.arange(ngals):
    y=y+pdf(x,x0=redshift[i],sigma=redshift_err[i])
    
y=y/ngals
plt.plot(x,y,ls='-.',lw=4)


# In[ ]:

x = np.linspace(0.,1.,100)
plt.plot(x,0.025*(1+x)**2)


# In[ ]:

x5=data1['Z']
y5=data1['IZ_PBLUE']
plt.scatter(x5, y5,lw=0,s=20)
plt.xlabel('z')
plt.ylabel('i-z')
plt.show()


# In[ ]:

np.log10(1.3)


# In[ ]:

import weightedstats as ws


# In[ ]:

np.cos(np.degrees(10.0))


# In[ ]:



