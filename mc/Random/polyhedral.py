import mchammer 
from tqdm import tqdm
from pymatgen.io.ase import AseAtomsAdaptor
from polyhedral_analysis.configuration import Configuration
from polyhedral_analysis.polyhedra_recipe import PolyhedraRecipe
from polyhedral_analysis.octahedral_analysis import isomer_is_cis, isomer_is_fac, isomer_is_mer, isomer_is_trans

import numpy as np
import matplotlib.pyplot as plt
from collections import Counter
from mchammer import DataContainer

dc=DataContainer.read('/Users/swc57/Dipoles_Distortions_Analytical_Workflow/mc/Sulfur/Li2FeSO_temp1025_rep=1.mc')
traj = dc.get('trajectory')
import copy
struc=traj[-1]
an=struc.get_atomic_numbers()
cations=np.append(np.where(an==3),np.where(an==26))
random_structures=[]

recipe = PolyhedraRecipe(method='distance cutoff', 
                                coordination_cutoff=2.0, 
                                central_atoms='O',
                                vertex_atoms=['Li','Fe'])

for random_count in range(1000):
    struc.symbols[cations]='Fe'
    indi=np.random.choice(np.array(cations),int(len(cations)/3*2),replace=False)
    struc.symbols[indi]='Li'
    s=copy.deepcopy(struc)
    random_structures.append(s)
all_c=[]
all_n_li2_trans=[]
all_n_li3_mer=[]
all_n_li4_trans=[]
logic=[]
from polyhedral_analysis.octahedral_analysis import opposite_vertex_pairs
  
structures = [AseAtomsAdaptor.get_structure(a) for a in tqdm(random_structures)]

configs = [Configuration(structure=s, recipes=[recipe]) for s in tqdm(structures)]

polyhedra = [p for c in configs for p in c.polyhedra]

for p in polyhedra:
    assert p.coordination_number == 6

pops = [p.vertex_count['Li'] for p in polyhedra]

li2 = [p for p in polyhedra if p.vertex_count['Li'] == 2]
li3 = [p for p in polyhedra if p.vertex_count['Li'] == 3]
li4 = [p for p in polyhedra if p.vertex_count['Li'] == 4]


c = Counter(pops)

n_li2_trans = sum([isomer_is_trans(p,check=False) for p in tqdm(li2)])
n_li3_mer = sum([isomer_is_mer(p,check=False) for p in tqdm(li3)])
n_li4_trans = sum([isomer_is_trans(p,check=False) for p in tqdm(li4)])
all_c.append(c)
all_n_li2_trans.append(n_li2_trans)
all_n_li3_mer.append(n_li3_mer)
all_n_li4_trans.append(n_li4_trans) 
count_store = np.zeros(7)
n4t=0
n2t=0
n3m=0
for count in range(1):
    for countLi in range(7):
        count_store[countLi] += (all_c[count][(countLi)])
    n4t+=all_n_li4_trans[count]
    n2t+=all_n_li2_trans[count]
    n3m+=all_n_li3_mer[count]
colorblind=["#0173B2", "#DE8F05", "#029E73", "#D55E00", "#CC78BC", "#CA9161", "#FBAFE4", "#949494", "#ECE133", "#56B4E9"]
col1 = "#d1cfcf"
col2 = "#00af82"
col3 = "#ef5843"
plt.figure(figsize=(4,4))
plt.bar(0,  count_store[0]/np.sum(count_store), color=col3,label='Single isomer')
for i in [1,2,5,6]:
    plt.bar(i, count_store[i]/np.sum(count_store), color=col3)
for i in range(2,4):
    plt.bar(i, count_store[i]/np.sum(count_store), color=col2)
for i in range(4,5):
    plt.bar(i, count_store[i]/np.sum(count_store), color=col2,label='cis/fac')
plt.bar(2, n2t/np.sum(count_store), color=col1,label='trans/mer')
plt.bar(3, n3m/np.sum(count_store), color=col1)
plt.bar(4, n4t/np.sum(count_store), color=col1)
plt.ylim(0,1.05)
plt.xticks([0,1,2,3,4,5,6])
plt.tick_params(axis="both", labelsize=14)
#plt.legend(fontsize=14,loc=2)
plt.xlabel(r'x in OLi$_{\operatorname{x}}$Fe$_{6-\operatorname{x}}$',size=16)
plt.ylabel(r'Prevelence',size=16)
plt.ylabel(r'Probability',size=16)


plt.tight_layout()

plt.savefig("polyhedra_random.pdf")
