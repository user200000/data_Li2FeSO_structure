from pymatgen.core import Structure, Lattice
from pymatgen.io.ase import AseAtomsAdaptor
import numpy as np
from ase import Atoms
import ase
from icet.tools import map_structure_to_reference
from icet import (ClusterSpace, StructureContainer, ClusterExpansion)
from trainstation import CrossValidationEstimator
from pymatgen.analysis.ewald import EwaldSummation
from pymatgen.transformations.standard_transformations import OxidationStateDecorationTransformation
primitive_structure = ase.io.read('primitive_structure.cif')
lattice=Lattice.from_parameters(3.91395,3.91395,3.91395,90.,90.,90.)
structure= Structure.from_spacegroup(sg='Pm-3m', lattice=lattice,species=['Fe','S','O'],
                                    coords=[[0.,.5,.5],
                                            [0.,0.,0.],
                                            [.5,.5,.5]])

from bsym.interface.pymatgen import unique_structure_substitutions
continuing_structs=[]
s222=structure*[2,2,2]
structs=unique_structure_substitutions(structure=s222,
                              to_substitute='Fe',
                              site_distribution={'Fe':8,'Li':16},
                              verbose=True,
                              show_progress='notebook')
degeneracy=[s.number_of_equivalent_configurations for s in structs]
ordered_structs=[s.get_sorted_structure() for s in structs]
aseMe=[ AseAtomsAdaptor.get_atoms(s) for s in ordered_structs]
ox_states = OxidationStateDecorationTransformation({'Fe':2,'Li':1,'O':-2,'S':-2})
ew_raise_me_up=[EwaldSummation(ox_states.apply_transformation(stroNo)).total_energy for stroNo in ordered_structs ]
scaled_energies1=np.array(ew_raise_me_up)/(4.746907+4.876905+4.726912)*3/40
continuing_structs.append(structs)
s221=structure*[2,2,1]
structs=unique_structure_substitutions(structure=s221,
                              to_substitute='Fe',
                              site_distribution={'Fe':4,'Li':8},
                              verbose=True,
                              show_progress='notebook')
degeneracy=[s.number_of_equivalent_configurations for s in structs]
ordered_structs=[s.get_sorted_structure() for s in structs]
ox_states = OxidationStateDecorationTransformation({'Fe':2,'Li':1,'O':-2,'S':-2})
ew_raise_me_up=[EwaldSummation(ox_states.apply_transformation(stroNo)).total_energy for stroNo in ordered_structs ]
scaled_energies2=np.array(ew_raise_me_up)/(4.746907+4.876905+4.726912)*3/40
continuing_structs.append(structs)
cs = ClusterSpace(structure=primitive_structure,
                  cutoffs=[14,9,5],
                  chemical_symbols=[['S'],['Li', 'Fe'],['Li', 'Fe'],['Li', 'Fe'],['O']])
sc = StructureContainer(cluster_space=cs)
scale=scaled_energies1
scale=np.append(scale,scaled_energies2)
for k in range(len(aseMe)):
    sc.add_structure(structure=aseMe[k], user_tag=str(k), properties={'potential_energy':scale[k]})
opt = CrossValidationEstimator(fit_data=sc.get_fit_data(key='potential_energy'), fit_method='rfe',estimator='lasso',  seed= np.random.randint(100),n_features=45)
opt.validate()
opt.train()
ce = ClusterExpansion(cluster_space=cs, parameters=opt.parameters)
ce.write('Ewald.ce')
