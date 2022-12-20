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
v_op_list=np.loadtxt('v_op_list')
#rc('font', family='Helvetica Neue')
from matplotlib import ticker
CIS_ORDER_DFT=[]
FOUR_ORDER_DFT=[]
VORDER_DFT=[]
CIS_ORDER_EL=[]
FOUR_ORDER_EL=[]
VORDER_EL=[]
v_op_list=np.array(v_op_list)
eew=np.loadtxt("eew.txt")
eew/=(4.746907+4.876905+4.726912)/3

ece=np.loadtxt("ece.txt")
deg=np.loadtxt("deg.txt")
cis_lists=np.loadtxt('ncis.txt')
l4_list=np.loadtxt('l4.txt')
n_poly=np.loadtxt('npoly.txt')
breadth_of_state_ece=[]
import scipy.constants as pc
kb=pc.physical_constants['Boltzmann constant in eV/K'][0]
FE=[]
EN=[]
POP=[]
for temp in np.arange(0,10000,25):
    weights=deg*np.exp(-np.array(ece-np.min(ece),dtype=np.float128)/(kb*temp)*40,dtype=np.float128)
    populations=weights/np.sum(weights)
    energy=np.sum(populations*ece)
    fe=-kb*temp*np.log(np.sum(weights))/40
    FE.append(fe)
    EN.append(energy)
    POP.append(populations)
    breadth_of_state_ece.append(np.sum(populations>=deg/np.sum(deg)))
for line in POP:
    CIS_ORDER_DFT+=[np.sum(line*cis_lists/n_poly)]
    FOUR_ORDER_DFT+=[np.sum(line*l4_list/n_poly)]
    VORDER_DFT+=[np.sum(line*v_op_list)]
import scipy.constants as pc
FE=[]
EN=[]
POP=[]
breadth_of_state_eew=[]
for temp in np.arange(0,10000,25):
    weights=deg*np.exp(-np.array(eew-np.min(eew),dtype=np.float128)/(kb*temp),dtype=np.float128)
    populations=weights/np.sum(weights)
    energy=np.sum(populations*eew/40)
    fe=-kb*temp*np.log(np.sum(weights))/40
    FE.append(fe)
    EN.append(energy)
    POP.append(populations)
    breadth_of_state_eew.append(np.sum(populations>=deg/np.sum(deg)))
for line in POP:
    CIS_ORDER_EL+=[np.sum((line)*cis_lists/n_poly)]
    FOUR_ORDER_EL+=[np.sum((line)*l4_list/n_poly)]
    VORDER_EL+=[np.sum(line*v_op_list)]
fig, axs = plt.subplots(2, 2,figsize=(8.5,8.5))
# Create first axes, the top-left plot with green plot
import scipy.constants as pc
plt.rcParams["font.size"] = "10"
fig.subplots_adjust(hspace=.35,wspace=.35)
formatter = ticker.ScalarFormatter(useMathText=True)
formatter.set_scientific(True)
formatter.set_powerlimits((-1,1))
from scipy.stats import gaussian_kde
x = np.linspace(-.5,14,10000)
kb=pc.physical_constants['Boltzmann constant in eV/K'][0]
eew=np.loadtxt("eew.txt")
ece=np.loadtxt("ece.txt")
deg=np.loadtxt("deg.txt")
ece*=40
eew/=(4.746907+4.876905+4.726912)/3
eew-=np.min(eew)
ece-=np.min(ece)


axs[0,0].set_ylabel(r'density / states $\operatorname{eV}^{-1}$',size=12)
axs[1,0].set_ylabel(r'density / states $\operatorname{eV}^{-1}$',size=12)
axs[1,1].set_ylabel(r'density / states $\operatorname{eV}^{-1}$',size=12)
weight=np.exp(-1/kb/1025*x)
weight[weight>1]=1

density=gaussian_kde(ece,weights=deg/np.sum(deg),bw_method=.03)
y=density(x)
weight=np.exp(-1/kb/1025*x)
weight[weight>1]=1
axs[1,0].plot(x,y,c='grey',linewidth=2)
axs[1,0].plot(x,y*weight,c='red',linewidth=2)


axs[1,0].set_ylim(0,0.015)
axs[1,0].tick_params(axis="both", labelsize=12)
axs[1,0].fill_between(x[:600*5],(y*weight)[:600*5],y[:600*5],color='grey',alpha=0.5)
axs[1,0].fill_between(x[:600*5],(y*weight)[:600*5],color='red',alpha=0.5)

axs[1,1].ticklabel_format(axis="y", style="scientific", scilimits=(0,0))
axs[1,0].ticklabel_format(axis="y", style="scientific", scilimits=(0,0))

axs[0,1].plot(np.arange(0,10000,25),3*np.array(VORDER_DFT),color='navy',linewidth=2)
axs[0,1].plot(np.arange(0,10000,25),3*np.array(VORDER_EL),color='red',linewidth=2,linestyle='dashed')
axs[0,1].set_xlim(0,2500)

axs[0,1].set_xlabel('Temperature/K',size=12)
axs[0,1].set_ylabel(r'$\langle\Phi_{SR}\rangle$',size=12)

axs[0,0].plot(np.arange(0,10000,25),CIS_ORDER_DFT,color='navy',linewidth=2)
axs[0,0].plot(np.arange(0,10000,25),1-np.array(CIS_ORDER_EL),color='red',linestyle='dashed',linewidth=2)

axs[0,0].set_xlabel('Temperature/K',size=12)
axs[0,0].set_ylabel(r'$\langle\Phi_{LR}\rangle$',size=12)
axs[0,0].set_xlim(0,2500)


density=gaussian_kde(eew,weights=deg/np.sum(deg),bw_method=.03)
y=density(x)
weight=np.exp(-1/kb/1025*x)
weight[weight>1]=1
axs[1,1].plot(x[:600*5],np.clip(y[:600*5],0,0.00035),c='grey',linewidth=2)
axs[1,1].plot(x[:600*5],(y*weight)[:600*5],c='red',linewidth=2)
axs[1,1].set_xlim(-0.1,1)
axs[1,1].set_ylim(0,0.00025)
axs[0,1].set_ylim(-0.05,1.05)
axs[0,0].set_ylim(-0.05,1.05)

axs[1,1].tick_params(axis="both", labelsize=12)
axs[1,1].fill_between(x[:600*5],(y*weight)[:600*5],np.clip(y[:600*5],0,0.00035),color='grey',alpha=0.5)
axs[1,1].fill_between(x[:600*5],(y*weight)[:600*5],color='red',alpha=0.5)
axs[1,0].set_xlim(-0.1,1)
axs[1,0].set_xlabel(r'$\Delta E / \operatorname{eV}$',size=12)
axs[1,1].set_xlabel(r'$\Delta E / \operatorname{eV}$',size=12)

plt.savefig("Order_paramaters_and_cdos.pdf")