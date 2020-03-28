import numpy as np
import scipy.linalg as lng
import pylab as pl
from matplotlib import colors, ticker, cm
from scipy.integrate import quad
import GreenAnisotropic2D

pl.rcParams['text.usetex'] = True
pl.rcParams['axes.labelsize'] = 17.
pl.rcParams['legend.fontsize']=18.
#pl.rcParams['legend.frameon'] = 'False'
pl.rcParams['legend.fontsize']=13.
pl.rcParams['xtick.labelsize']=13.
pl.rcParams['ytick.labelsize']=13.
pl.rcParams['legend.numpoints']=1

#fig_format='.eps'
fig_format='.png'

# Select symmetry
isym=0
if (isym==0):
	# Anisotropic case
	p=[.3,3.6,1.,.7,2.1,2.0]
if (isym==1):
	# Orthotropic case 
	p=[.3,0.,.6,.8,1.2,2.5]
if (isym==2):
	# R0-Orthotropic case
	p=[.3,.5,1.2,2.7]
if (isym==3):
	# Square symmetric case
	p=[.3,.8,1.2,2.5]
if (isym==4):
	# Polar isotropic case
	p=[1.2,.3]
if (isym==5):
	# Isotropic case
	p=[1.2,.3]

#
# Define a medium object for the material
mat=GreenAnisotropic2D.medium()
#
# Set up symmetry
mat.set_sym(isym,p)

#
# Verify strain energy is positive and the symmetry is properly identified
def check_pos(mat):
	if (mat.isym==0):
		flag1=True
		flag2=True
		if not (mat.T0-mat.R0>0):
			flag1=False
		if not (mat.T1*(mat.T0**2-mat.R0**2)-2.*mat.R1**2*(mat.T0-mat.R0*np.cos(4.*(mat.P0-mat.P1)))>0): 
			flag1=False
		if not (mat.R0>=0): 
			flag1=False
		if not (mat.R1>=0): 
			flag1=False
		if ((mat.R0==0)|(mat.R1==0)|(np.sin(4.*(mat.P0-mat.P1))==0)): 
			flag2=False				
	elif (mat.isym==1):
		flag1=True
		flag2=True
		if not (mat.T0-mat.R0>0): 
			flag1=False
		if not (mat.T1*(mat.T0+(-1.)**mat.K*mat.R0)-2.*mat.R1**2>0): 
			flag1=False
		if not (mat.R0>=0): 
			flag1=False
		if not (mat.R1>=0): 
			flag1=False
		if ((mat.R0==0)|(mat.R1==0)): 
			flag2=False				
	elif (mat.isym==2): 
		flag1=True
		flag2=True
		if not (mat.T0>0): 
			flag1=False
		if not (mat.T1*mat.T0-2.*mat.R1**2>0): 
			flag1=False
		if not (mat.R1>=0): 
			flag1=False
		if (mat.R1==0): 
			flag2=False	
	elif (mat.isym==3): 
		flag1=True
		flag2=True
		if not (mat.T0-mat.R0>0): 
			flag1=False
		if not (mat.T1*(mat.T0-mat.R0)>0): 
			flag1=False
		if not (mat.R0>=0): 
			flag1=False
		if (mat.R0==0): 
			flag2=False	
	elif (mat.isym==4): 
		return 0
	elif (mat.isym==5):
		if not (mat.k>0): 
			flag1=False
		if not (mat.m>0): 
			flag1=False
	#
	if (flag1&flag2): 
		return 0
	elif (not flag1):
		return 1
	elif (not flag2):
		return 2

status=check_pos(mat)
if (status==0):
	print "Strain energy is positive and symmetry is properly identified."
	#
	# Compute first complete Barnett-Lothe integral
	mat.get_S()
	# Compute second complete Barnett-Lothe integral
	mat.get_H()
elif (status==1):	
	print "Strain energy is negative."
elif (status==2):	
	print "Symmetry not properly identified."


#
# Simulate a random curve
np.random.seed(1307341095)
ns=15
cs=.1-np.random.rand(ns)*.5
#
# Curve
def r_crv(cs,th):
	ns=len(cs)
	crv=1.
	for k in range(ns):
		crv+=cs[k]*np.cos(th*k)/ns	
	return crv
# 
# Not normalized tangent
def d_crv(cs,th):
	ns=len(cs)
	d_crv_x=-(np.sin(th)+1./ns*np.sum([cs[k]*(k*np.sin(k*th)*np.cos(th)+np.cos(k*th)*np.sin(th)) for k in range(ns)]))
	d_crv_y=np.cos(th)+1./ns*np.sum([cs[k]*(np.cos(k*th)*np.cos(th)-k*np.sin(k*th)*np.sin(th)) for k in range(ns)])
	return (d_crv_x,d_crv_y)
#
# Arc length
def arc_length(cs,th):
	_d_crv=d_crv(cs,th)	
	_d_mag=np.sqrt(_d_crv[0]**2+_d_crv[1]**2)
	return _d_mag
#
# Second derivative of the cruve
def dd_crv(cs,th):
	ns=len(cs)
	dd_crv_x=-np.cos(th)+1./ns*np.sum([cs[k]*(2.*k*np.sin(k*th)*np.sin(th)-(1.+k**2)*np.cos(k*th)*np.cos(th)) for k in range(ns)])
	dd_crv_y=-(np.sin(th)+1./ns*np.sum([cs[k]*(2.*k*np.sin(k*th)*np.cos(th)+(1.+k**2)*np.cos(k*th)*np.sin(th)) for k in range(ns)]))
	return (dd_crv_x,dd_crv_y)
#
# Derivative of the arc length
def dmag_d_crv(cs,th):
	ns=len(cs)
	_d_crv=d_crv(cs,th)
	_mag_d_crv=np.sqrt(_d_crv[0]**2+_d_crv[1]**2)
	_dmag_d_crv=-(1.+1./ns*np.sum([(1.-k**2)*cs[k]*np.cos(k*th) for k in range(ns)]))*np.sum([k*cs[k]*np.sin(k*th) for k in range(ns)])
	_dmag_d_crv/=(ns*_mag_d_crv)
	return _dmag_d_crv
#
# Unit outward normal
def n_crv2(cs,th):
	_d_crv=d_crv(cs,th)
	_mag_d_crv=np.sqrt(_d_crv[0]**2+_d_crv[1]**2)
	return (_d_crv[1]/_mag_d_crv,-_d_crv[0]/_mag_d_crv)
	
#
# Stress fields
# sijt:= Stress field sij due to a unit force concentrated at the origin along e_t
def n1(t): return np.cos(t)
def n2(t): return np.sin(t)
def m1(t): return -np.sin(t)
def m2(t): return np.cos(t)
def s111(r,t,mat): 
	return mat.L1111*mat.dnGi(1,[1,1,1],r,t)+mat.L1112*mat.dnGi(1,[1,1,2],r,t)+mat.L1121*mat.dnGi(1,[2,1,1],r,t)+mat.L1122*mat.dnGi(1,[2,1,2],r,t)
def s112(r,t,mat): 
	return mat.L1111*mat.dnGi(1,[1,2,1],r,t)+mat.L1112*mat.dnGi(1,[1,2,2],r,t)+mat.L1121*mat.dnGi(1,[2,2,1],r,t)+mat.L1122*mat.dnGi(1,[2,2,2],r,t)
def s221(r,t,mat): 
	return mat.L2211*mat.dnGi(1,[1,1,1],r,t)+mat.L2212*mat.dnGi(1,[1,1,2],r,t)+mat.L2221*mat.dnGi(1,[2,1,1],r,t)+mat.L2222*mat.dnGi(1,[2,1,2],r,t)
def s222(r,t,mat): 
	return mat.L2211*mat.dnGi(1,[1,2,1],r,t)+mat.L2212*mat.dnGi(1,[1,2,2],r,t)+mat.L2221*mat.dnGi(1,[2,2,1],r,t)+mat.L2222*mat.dnGi(1,[2,2,2],r,t)
def s121(r,t,mat): 
	return mat.L1211*mat.dnGi(1,[1,1,1],r,t)+mat.L1212*mat.dnGi(1,[1,1,2],r,t)+mat.L1221*mat.dnGi(1,[2,1,1],r,t)+mat.L1222*mat.dnGi(1,[2,1,2],r,t)
def s122(r,t,mat): 
	return mat.L1211*mat.dnGi(1,[1,2,1],r,t)+mat.L1212*mat.dnGi(1,[1,2,2],r,t)+mat.L1221*mat.dnGi(1,[2,2,1],r,t)+mat.L1222*mat.dnGi(1,[2,2,2],r,t)
def s211(r,t,mat): 
	return mat.L2111*mat.dnGi(1,[1,1,1],r,t)+mat.L2112*mat.dnGi(1,[1,1,2],r,t)+mat.L2121*mat.dnGi(1,[2,1,1],r,t)+mat.L2122*mat.dnGi(1,[2,1,2],r,t)
def s212(r,t,mat): 
	return mat.L2111*mat.dnGi(1,[1,2,1],r,t)+mat.L2112*mat.dnGi(1,[1,2,2],r,t)+mat.L2121*mat.dnGi(1,[2,2,1],r,t)+mat.L2122*mat.dnGi(1,[2,2,2],r,t)

#
# Taction fields on random curves
# tij_on_crv := Traction component ti due to a unit force concentrated at the origin along e_j	
def tij_on_crv(t,cs,i,j,mat): 
	r=r_crv(cs,t)
	n=n_crv2(cs,t)	
	if (i==1)&(j==1):
		return s111(r,t,mat)*n[0]+s121(r,t,mat)*n[1]
	elif (i==2)&(j==2):
		return s212(r,t,mat)*n[0]+s222(r,t,mat)*n[1]
	elif (i==1)&(j==2):
		return s112(r,t,mat)*n[0]+s122(r,t,mat)*n[1]
	elif (i==2)&(j==1):
		return s211(r,t,mat)*n[0]+s221(r,t,mat)*n[1]

#
# Tractions multiplied by arc length. 
# Used for integration of the traction field on the curve.
def tij_times_arc_length(t,cs,i,j,mat): 
	ds=arc_length(cs,t)		
	if (i==1)&(j==1):
		tij=tij_on_crv(t,cs,1,1,mat)
	elif (i==2)&(j==2):
		tij=tij_on_crv(t,cs,2,2,mat)
	elif (i==1)&(j==2):
		tij=tij_on_crv(t,cs,1,2,mat)
	elif (i==2)&(j==1):
		tij=tij_on_crv(t,cs,2,1,mat)
	return tij*ds



#
# Plot components of the Green's function along some specific directions
def plot_green_gradients(mat,th,name_flag='',rmax=7.,lim_y=10.**9,thresh_y=.0000005):
	Gvals1=[]; Gvals2=[]; Gvals3=[]; Gvals4=[]
	Gvals5=[]; Gvals6=[]; Gvals7=[]; Gvals8=[]	
	rvals=np.linspace(.1,rmax,60)
	ylabels=r'$G^{(n)}_{12,k_1\dots k_n}(r,$'+str(th)+r'$)$'
	dg_labels=[
	r'$G^{(1)}_{12,1}(r,\theta)$', r'$G^{(2)}_{12,12}(r,\theta)$',
	r'$G^{(3)}_{12,121}(r,\theta)$', r'$G^{(4)}_{12,1212}(r,\theta)$',
	r'$G^{(5)}_{12,12121}(r,\theta)$', r'$G^{(6)}_{12,121212}(r,\theta)$',
	r'$G^{(7)}_{12,1212121}(r,\theta)$', r'$G^{(8)}_{12,12121212}(r,\theta)$']
	ith=0
	lwidth=1.5
	Gvals1.append(np.array([mat.dnGi(1,[1,2,1],r,th) for r in rvals]))
	Gvals2.append(np.array([mat.dnGi(2,[1,2,1,2],r,th) for r in rvals]))
	Gvals3.append(np.array([mat.dnGi(3,[1,2,1,2,1],r,th) for r in rvals]))
	Gvals4.append(np.array([mat.dnGi(4,[1,2,1,2,1,2],r,th) for r in rvals]))
	Gvals5.append(np.array([mat.dnGi(5,[1,2,1,2,1,2,1],r,th) for r in rvals]))
	Gvals6.append(np.array([mat.dnGi(6,[1,2,1,2,1,2,1,2],r,th) for r in rvals]))
	Gvals7.append(np.array([mat.dnGi(7,[1,2,1,2,1,2,1,2,1],r,th) for r in rvals]))
	Gvals8.append(np.array([mat.dnGi(8,[1,2,1,2,1,2,1,2,1,2],r,th) for r in rvals]))
	#
	fig=pl.figure()
	dg1,=pl.plot(rvals,Gvals1[ith],lw=lwidth)
	dg2,=pl.plot(rvals,Gvals2[ith],lw=lwidth)
	dg3,=pl.plot(rvals,Gvals3[ith],lw=lwidth)
	dg4,=pl.plot(rvals,Gvals4[ith],lw=lwidth)
	dg5,=pl.plot(rvals,Gvals5[ith],lw=lwidth)
	dg6,=pl.plot(rvals,Gvals6[ith],lw=lwidth)
	dg7,=pl.plot(rvals,Gvals7[ith],lw=lwidth)	
	dg8,=pl.plot(rvals,Gvals8[ith],lw=lwidth)
	#
	pl.yscale('symlog',linthreshy=thresh_y)
	pl.ylabel(ylabels)
	ith+=1
	pl.xlabel(r'$r>0$')
	pl.xlim(.1,rmax)
	pl.ylim(-lim_y,lim_y)
	leg1=pl.legend([dg1,dg2,dg3,dg4],dg_labels[:4],loc=1)
	ax=pl.gca().add_artist(leg1)
	leg2=pl.legend([dg5,dg6,dg7,dg8],dg_labels[4:8],loc=4)
	pl.savefig('figdG'+name_flag+fig_format,bbox_inches='tight')
	pl.close(fig)

plot_green_gradients(mat,np.pi/3.,name_flag='_isym'+str(isym)+'_')


#
# Get a Mandel representation of the stiffness
def get_Lmat(mat):
	Lmat=np.array([[mat.L1111,mat.L1122,np.sqrt(2)*mat.L1112],
				  [mat.L2211,mat.L2222,np.sqrt(2)*mat.L2212],
				  [np.sqrt(2)*mat.L1211,np.sqrt(2)*mat.L1222,2.*mat.L1212]])
	return Lmat
#	
# Generalized Young's modulus
def gen_E(th,mat,Smat):
	m2vec=np.array([np.cos(th)**2,np.sin(th)**2,np.sqrt(2)*np.cos(th)*np.sin(th)])
	return np.dot(np.dot(m2vec,Smat),m2vec)**-1
#
# Generalized shear modulus
def gen_mu(th,mat,Smat):
	mbyp_vec=np.array([np.cos(th)*(-np.sin(th)),np.sin(th)*np.cos(th),np.sqrt(2)*np.cos(th)*np.cos(th)])
	mbyp_vec+=np.array([np.cos(th)*(-np.sin(th)),np.sin(th)*np.cos(th),np.sqrt(2)*np.sin(th)*(-np.sin(th))])
	mbyp_vec/=2.
	return np.dot(4.*np.dot(mbyp_vec,Smat),mbyp_vec)**-1
#
# Generalized Poisson's ratio
def gen_nu(th,mat,Smat):
	m2vec=np.array([np.cos(th)**2,np.sin(th)**2,np.sqrt(2)*np.cos(th)*np.sin(th)])
	p2vec=np.array([(-np.sin(th))**2,np.cos(th)**2,np.sqrt(2)*(-np.sin(th))*np.cos(th)])
	return -np.dot(np.dot(p2vec,Smat),m2vec)/np.dot(np.dot(m2vec,Smat),m2vec)
#
# Generalized absoulte value of Poisson's ratio
def gen_nu(th,mat,Smat):
	m2vec=np.array([np.cos(th)**2,np.sin(th)**2,np.sqrt(2)*np.cos(th)*np.sin(th)])
	p2vec=np.array([(-np.sin(th))**2,np.cos(th)**2,np.sqrt(2)*(-np.sin(th))*np.cos(th)])
	return -np.dot(np.dot(p2vec,Smat),m2vec)/np.dot(np.dot(m2vec,Smat),m2vec)
#
# Plot polar diagram of generalized moduli
def plot_polar(mat,fname,add_line_flag,angle):
	nvals=500
	tvals=np.linspace(0,2.*np.pi,nvals)		
	#
	Lmat=get_Lmat(mat)
	Smat=lng.inv(Lmat)
	E_iso=quad(gen_E,0,2.*np.pi,args=(mat,Smat))[0]/2./np.pi
	E_vals=np.array([gen_E(t,mat,Smat) for t in tvals])/E_iso
	mu_iso=quad(gen_mu,0,2.*np.pi,args=(mat,Smat))[0]/2./np.pi
	mu_vals=np.array([gen_mu(t,mat,Smat) for t in tvals])/mu_iso
	nu_iso=quad(gen_nu,0,2.*np.pi,args=(mat,Smat))[0]/2./np.pi
	nu_vals=np.array([gen_nu(t,mat,Smat) for t in tvals])/nu_iso	
	#
	fig=pl.figure()
	ax = pl.subplot(111, projection='polar')
	ax.plot(tvals, E_vals,lw=1.5,label=r'$E(\theta)/E$')
	ax.plot(tvals, mu_vals,lw=1.5,label=r'$\mu(\theta)/\mu$')
	ax.plot(tvals[nu_vals>=0], nu_vals[nu_vals>=0],lw=1.5,label=r'$\nu(\theta)/\nu$')
	ax.plot(tvals[nu_vals<0], nu_vals[nu_vals<0],lw=1.5)
	max_mag=max(1.1*max(E_vals),1.1*max(mu_vals),1.1*max(nu_vals))
	if add_line_flag:
		ax.plot(2*[angle],[0,max_mag],lw=2.,color='k')
	ax.set_rmax(max_mag)
	pl.legend()
	pl.savefig('figPolar'+fname+fig_format,bbox_inches='tight')
	pl.close(fig)



def plot_traction_on_crv_with_vonMises(cs,f1,f2,mat,name_flag='',normalizer=1.,leg_label=''):
	#
	# Discretize cruve
	nvals=500
	tvals=np.linspace(0,2.*np.pi,nvals)	
	curve=np.array([r_crv(cs,t) for t in tvals])
	xpts=curve*np.cos(tvals)
	ypts=curve*np.sin(tvals)	
	# 
	# Discretize plane
	N=500
	x=np.linspace(-1.25,1.25,N)
	y=np.linspace(-1.25,1.25,N)
	X,Y=np.meshgrid(x,y)
	R=np.sqrt(X**2+Y**2)
	TH=np.arctan2(Y,X)
	#
	# Compute deviatoric stress component due concentrated force
	D11t=np.zeros((len(x),len(y)))
	D12t=np.zeros((len(x),len(y)))	
	for ix in range(len(x)):
		for iy in range(len(y)):
			if (np.sqrt(x[ix]**2+y[iy]**2)<=r_crv(cs,TH[ix][iy])):
				D11t[ix][iy]=f1*.5*(s111(R[ix,iy],TH[ix,iy],mat)-s221(R[ix,iy],TH[ix,iy],mat))\
							+f2*.5*(s112(R[ix,iy],TH[ix,iy],mat)-s222(R[ix,iy],TH[ix,iy],mat))
				D12t[ix][iy]=f1*s121(R[ix,iy],TH[ix,iy],mat)+f2*s122(R[ix,iy],TH[ix,iy],mat)		
	D22t=-D11t
	#
	# Compute von Mises stress
	vonM=np.sqrt(D11t*D11t+D22t*D22t+2.*D12t*D12t)
	vonM=np.ma.masked_where(vonM<=0,vonM)
	#
	#
	fig, ax = pl.subplots(sharex=True,sharey=True)
	ax.set_frame_on(False)
	ax.set_axis_off()
	ax.set_aspect('equal')
	#
	# Create contour plot of the von Mises stress
	cs_plot = ax.contourf(X, Y, vonM, locator=ticker.LogLocator(base=2.), cmap=cm.YlOrRd)
	#
	# Add colorbar
	#cbar = fig.colorbar(cs_plot,shrink=.8)
	#cbar.set_label(leg_label)
	#
	# Set scaling factor
	sc=.2	
	#
	# Plot concentrated force
	pl.plot(0,0,"ko")
	ax.arrow(0, 0,sc*f1,sc*f2,head_width=0.05, head_length=0.1, fc='k', ec='k')
	sc=.7	
	#
	# Plot the curve and its interior
	pl.plot(xpts,ypts,lw=2,color='k')
	#
	# Plot 20 tractions on the curve
	for tv in range(0,nvals,20):
		n=n_crv2(cs,tvals[tv])
		_mag_d_crv=arc_length(cs,tvals[tv])
		ax.arrow(xpts[tv],ypts[tv],\
		sc*(f1*tij_times_arc_length(tvals[tv],cs,1,1,mat)/_mag_d_crv+f2*tij_times_arc_length(tvals[tv],cs,1,2,mat)/_mag_d_crv),\
		sc*(f1*tij_times_arc_length(tvals[tv],cs,2,1,mat)/_mag_d_crv+f2*tij_times_arc_length(tvals[tv],cs,2,2,mat)/_mag_d_crv),\
		head_width=0.025,head_length=0.05,fc='r', ec='r')	
	#
	pl.xlim(-1.2,1.2)
	pl.savefig('figTvM'+name_flag+fig_format,bbox_inches='tight')
	pl.close(fig)


#
# Define concentrated forces
f1=np.sqrt(2.)/2.
f2=np.sqrt(2.)/2.
#
# Compute isotropic shear modulus
Lmat=get_Lmat(mat)
Smat=lng.inv(Lmat)
mu_iso=quad(gen_mu,0,2.*np.pi,args=(mat,Smat))[0]/2./np.pi
#
plot_traction_on_crv_with_vonMises(cs,f1,f2,mat,name_flag='_isym'+str(isym)+'_',normalizer=mu_iso,leg_label=r'$\|\mathbf{s}(\underline{x})\|/\mu$')
