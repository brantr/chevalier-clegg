import numpy as np
import matplotlib.pyplot as plt
from array_io import *

fname = "chevalier.txt"

x, M, u_star, P_star, rho_star = read_five_arrays(fname)

xmin = -0.5
xmax =  0.5
ymin = -4.0
ymax =  1.0

fsize = 5.0

bg_size = 0.1
tg_size = 0.02
lg_size = 0.16
rg_size = 0.03

wx = 1.0 - lg_size - rg_size
wy = 1.0 - bg_size - tg_size

plt.figure(figsize=(fsize,fsize*wx/wy))
a0  = plt.axes([lg_size,bg_size,wx,wy])
plt.xlim(xmin,xmax)
plt.ylim(ymin,ymax)
plt.xlabel(r'$\mathrm{log}_{10}\,r/R$')
plt.ylabel(r'$\mathrm{log}_{10}\,y$')

log_x   = np.log10(x)
log_M   = np.log10(M)
log_u   = np.log10(u_star)
log_P   = np.log10(P_star)
log_rho = np.log10(rho_star)
plt.plot(log_x,log_M,'-',color="0.0")
plt.text(0.2,0.26,r'$\mathrm{log}_{10}\,M$',color='0.0')
plt.plot(log_x,log_u,'-',color="red")
plt.text(0.2,-0.2,r'$\mathrm{log}_{10}\,(u_\star)$',color='red')
plt.plot(log_x,log_P,'-',color="blue")
plt.text(0.2,-3.5,r'$\mathrm{log}_{10}\,(P_\star)$',color='blue')
plt.plot(log_x,log_rho,'-',color="green")
plt.text(0.2,-2.2,r'$\mathrm{log}_{10}\,(\rho_\star)$',color='green')

s = 'chevalier.png'
plt.savefig(s,bbox_inches='tight')
#plt.show()



