import sys, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
d=np.load(sys.argv[1]); out=sys.argv[2]
q=d['qpg'][:,0]; f=d['fxc'][:,0].real
fin=q>1e-4
f00=f[:,0,0]; a0=-f00[0]*q[0]**2
g=-f00-a0/q**2
c=np.polyfit(q[fin],g[fin],2)
qq=np.linspace(0.04,0.72,200)
fig,ax=plt.subplots(1,3,figsize=(13.5,3.9))
ax[0].plot(q[fin],-f00[fin],'o',label='exported (27 q, IBZ)')
ax[0].plot(qq,a0/qq**2+np.polyval(c,qq),'-',label=r'$\alpha_0/q^2+\gamma(q)$')
ax[0].plot(qq,a0/qq**2,'--',label=r'$\alpha_0/q^2$, $\alpha_0$=%.3f (optical)'%a0)
ax[0].set_yscale('log'); ax[0].set_xlabel(r'$|q|$ (bohr$^{-1}$)'); ax[0].set_ylabel(r'$-f_{xc}^{00}(q,\omega=0)$ (Ha bohr$^3$)')
ax[0].set_title('head'); ax[0].legend()
ax[1].plot(q[fin],g[fin],'o')
ax[1].plot(qq,np.polyval(c,qq),'-',label=r'fit %.1f + %.1f$q$ %+.1f$q^2$'%(c[2],c[1],c[0]))
ax[1].set_xlabel(r'$|q|$ (bohr$^{-1}$)'); ax[1].set_ylabel(r'$\gamma(q)=-f_{xc}^{00}-\alpha_0/q^2$ (Ha bohr$^3$)')
ax[1].set_title('head minus long-range part'); ax[1].legend()
for j in range(1,7): ax[2].plot(q,-f[:,j,j],'o',ms=3,color='C%d'%(j-1))
off=np.array([np.abs(x[1:,1:]-np.diag(np.diag(x[1:,1:]))).max() for x in f])
ax[2].plot(q,off,'k^',ms=4,label=r'max $|f_{xc}^{GG\prime}|$, $G\neq G\prime$')
ax[2].plot([],[],'o',color='gray',label=r'$-f_{xc}^{GG}$, 6 in-plane $G$')
ax[2].set_xlabel(r'$|q|$ (bohr$^{-1}$)'); ax[2].set_ylabel(r'Ha bohr$^3$'); ax[2].set_title('body (first in-plane star)'); ax[2].legend()
fig.tight_layout(); fig.savefig(out,dpi=150)
print('alpha0=%.4f  gamma fit %s'%(a0,c[::-1]))
