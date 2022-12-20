from pymatgen.core import Structure, Lattice
from polyhedral_analysis.configuration import Configuration
from polyhedral_analysis.polyhedra_recipe import PolyhedraRecipe
from polyhedral_analysis.octahedral_analysis import isomer_is_cis, isomer_is_fac, isomer_is_mer, isomer_is_trans
import numpy as np
from ase import Atoms
from icet.tools import map_structure_to_reference
from icet import (ClusterSpace, StructureContainer, ClusterExpansion)
from tqdm import tqdm
from polyhedral_analysis.octahedral_analysis import opposite_vertex_pairs

lattice=Lattice.from_parameters(3.91395,3.91395,3.91395,90.,90.,90.)
structure= Structure.from_spacegroup(sg='Pm-3m', lattice=lattice,species=['Fe','S','O'],
                                    coords=[[0.,.5,.5],
                                            [0.,0.,0.],
                                            [.5,.5,.5]])
from bsym.interface.pymatgen import unique_structure_substitutions
s222=structure*[2,2,2]
structs=unique_structure_substitutions(structure=s222,
                              to_substitute='Fe',
                              site_distribution={'Fe':8,'Li':16},
                              verbose=True,
                              show_progress=True)
degeneracy=[s.number_of_equivalent_configurations for s in structs]

ordered_structs=[s.get_sorted_structure() for s in structs]

from pymatgen.io.ase import AseAtomsAdaptor
import ase

aseMe=[ AseAtomsAdaptor.get_atoms(s) for s in ordered_structs]
ce=ClusterExpansion.read('../fitting/Sulfur/Li2FeSO.ce')

ece=[]
for structNow in aseMe:
    ece+=[ce.predict(structNow)]
from pymatgen.analysis.ewald import EwaldSummation
from pymatgen.transformations.standard_transformations import OxidationStateDecorationTransformation
ox_states = OxidationStateDecorationTransformation({'Fe':2,'Li':1,'O':-2,'S':-2})
eew=[EwaldSummation(ox_states.apply_transformation(stroNo)).total_energy for stroNo in ordered_structs ]

eew=np.array(eew)
structures = [AseAtomsAdaptor.get_structure(a) for a in tqdm(aseMe)]

cis_lists=[]
pops_lists=[]
n_poly=[]
l4_list=[]
trans_list=[]
v_op_list=[]
# Just processing the first 100 atoms
nstrucList = [a.repeat(2) for a in tqdm(aseMe)]
structures = [AseAtomsAdaptor.get_structure(a) for a in tqdm(nstrucList)]
for s in tqdm(structures):
    recipe = PolyhedraRecipe(method='distance cutoff',
                         coordination_cutoff=2,
                         central_atoms='O',
                         vertex_atoms=['Li','Fe'])
    configs = [Configuration(structure=s, recipes=[recipe])]

    polyhedra = [p for c in configs for p in c.polyhedra]
    for p in polyhedra:
        assert p.coordination_number == 6
    pops = [p.vertex_count['Li'] for p in polyhedra]
    li4 = [p for p in polyhedra if p.vertex_count['Li'] == 4]
    n_li4_cis = sum([isomer_is_cis(p,check=False) for p in li4])
    n_li4_trans= sum([isomer_is_trans(p,check=False) for p in li4])
    logic=[]
    for p in polyhedra:
        for k in range(3):
            logic.append((opposite_vertex_pairs(p,check=False)[k][0].label==opposite_vertex_pairs(p,check=False)[k][1].label)*('Fe'==opposite_vertex_pairs(p,check=False)[k][1].label))
    cis_lists.append(n_li4_cis)
    trans_list.append(n_li4_trans)
    pops_lists.append(pops)
    n_poly.append(len(pops))
    l4_list.append(len(li4))
    v_op_list.append(np.sum(logic)/len(logic))
np.savetxt('ece.txt',ece)
np.savetxt('eew.txt',eew)
np.savetxt('npoly.txt',n_poly)
np.savetxt('ncis.txt',cis_lists)
np.savetxt('l4.txt',l4_list)
np.savetxt('deg.txt',degeneracy)
np.savetxt('v_op_list',v_op_list)

