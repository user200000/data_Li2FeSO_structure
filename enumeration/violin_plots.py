import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid.inset_locator import (inset_axes, InsetPosition,
                                                  mark_inset)
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams.update({
    "pdf.use14corefonts": True
})
from matplotlib import rc
#rc('font', family='Helvetica Neue')
from matplotlib import ticker
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp

import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
fig,axs=plt.subplots(nrows=1, ncols=1,sharey=False)
fig.set_size_inches(9.5, 4)
ece=np.loadtxt('ece.txt')
eew=np.loadtxt('eew.txt')/40
n_poly=np.loadtxt('npoly.txt')
cis_lists=np.loadtxt('ncis.txt')
l4_list=np.loadtxt('l4.txt')
degeneracy=np.loadtxt('deg.txt')
eew/=(4.746907+4.876905+4.726912)/3
fig.set_size_inches(4, 4)

enli=np.array(ece)[np.array(l4_list)/np.array(n_poly)==1]
cli=np.array(cis_lists)[np.array(l4_list)/np.array(n_poly)==1]
pli=np.array(n_poly)[np.array(l4_list)/np.array(n_poly)==1]
data=[]
for nb in np.unique(cli/pli):
    data.append((enli[cli/pli==nb]-np.min(ece))*1000)
parts=axs.violinplot(data, positions=np.unique(cli/pli), showmeans=True, bw_method=0.25,widths=0.1)
for pc in parts['bodies']:
    pc.set_facecolor('#F22525')
    pc.set_edgecolors('#F22525')

for partname in ('cbars','cmins','cmaxes','cmeans'):
    vp = parts[partname]
    vp.set_edgecolor('#F22525')
    vp.set_linewidth(1)
axs.hlines(0,xmin=-0.05,xmax=1, lw=1, linestyles='--',color='black')
#axs.set_xlabel(r'$\operatorname{n}_{\operatorname{cis-OFe}_{2}\operatorname{Li}_{4}}/\operatorname{n}_{\operatorname{OFe}_{2}\operatorname{Li}_{4}}$',size=14)
plt.tight_layout()
labelx = -0.3  # axes coords
#plt.ylabel(r'Energy / eV atom$^{-1}$',size=14)



#plt.sca(axs[1])
#plt.yticks(np.arange(-5.29,-5.269,.01))
#axs[0,count].text(0.7, 0.93, "PC.A", transform=axs[0,count].transAxes,
#  fontsize=14, va='center')
#axs[1,count].text(0.02, 0.93, "PC.B", transform=axs[1,count].transAxes,
#  fontsize=14, va='center')
ecepc=ece
plt.ylim(-5,30)
#axs.set_ylabel(r'Energy / meV atom$^{-1}$',size=14)
#axs.set_title(r'DFT Relaxation',size=14)
plt.tight_layout()
plt.yticks(np.arange(0,65,20),fontsize=14.5 )
plt.xticks(np.arange(0,1.1,.25),fontsize=14.5)
axs.set_xlabel(r'fraction $\operatorname{cis-OFe}_{2}\operatorname{Li}_{4}$',size=14)

axs.set_ylabel(r'$\Delta E / \operatorname{meV} \operatorname{atom}^{-1}$',size=14)


plt.tight_layout()
plt.savefig('cluster_expansion_energies.pdf')

import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
fig,axs=plt.subplots(nrows=1, ncols=1,sharey=False)
fig.set_size_inches(9.5, 4)
ece=np.loadtxt('ece.txt')
eew=np.loadtxt('eew.txt')/40
n_poly=np.loadtxt('npoly.txt')
cis_lists=np.loadtxt('ncis.txt')
l4_list=np.loadtxt('l4.txt')
degeneracy=np.loadtxt('deg.txt')
eew/=(4.746907+4.876905+4.726912)/3


fig.set_size_inches(4, 4)

enli=np.array(eew)[np.array(l4_list)/np.array(n_poly)==1]
cli=np.array(cis_lists)[np.array(l4_list)/np.array(n_poly)==1]
pli=np.array(n_poly)[np.array(l4_list)/np.array(n_poly)==1]
data=[]
for nb in np.unique(cli/pli):
    data.append((enli[cli/pli==nb]-np.min(eew))*1000)
parts=axs.violinplot(data, positions=np.unique(cli/pli), showmeans=True, bw_method=0.25,widths=0.1)
for pc in parts['bodies']:
    pc.set_facecolor('#F22525')
    pc.set_edgecolors('#F22525')

for partname in ('cbars','cmins','cmaxes','cmeans'):
    vp = parts[partname]
    vp.set_edgecolor('#F22525')
    vp.set_linewidth(1)
axs.hlines(0,xmin=-0.05,xmax=1, lw=1, linestyles='--',color='black')
axs.set_xlabel(r'fraction $\operatorname{cis-OFe}_{2}\operatorname{Li}_{4}$',size=14)
plt.tight_layout()
labelx = -0.3  # axes coords
#plt.ylabel(r'Energy / eV atom$^{-1}$',size=14)




ecepc=ece
plt.ylim(-5,65)
axs.set_ylabel(r'$\Delta E / \operatorname{meV} \operatorname{atom}^{-1}$',size=14)
#axs.set_title(r'DFT Relaxation',size=14)
plt.tight_layout()
plt.yticks(np.arange(0,65,20),fontsize=14.5 )
plt.xticks(np.arange(0,1.1,.25),fontsize=14.5)
plt.tight_layout()
plt.savefig('ewald_energies.pdf')
