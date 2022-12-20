from ase.db import connect
from icet import (ClusterSpace, StructureContainer, ClusterExpansion)
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from trainstation import CrossValidationEstimator 

from pymatgen.io.ase import AseAtomsAdaptor
import ase
from pymatgen.analysis.ewald import EwaldSummation
from pymatgen.transformations.standard_transformations import OxidationStateDecorationTransformation
from icet.tools import map_structure_to_reference
from ase import Atoms

primitive_structure=ase.io.read('primitive_structure.cif')

cs = ClusterSpace(structure=primitive_structure,
                  cutoffs=[15,9,5],
                  chemical_symbols=[['S'],['Li', 'Fe'],['Li', 'Fe'],['Li', 'Fe'],['O']])

sc = StructureContainer(cluster_space=cs)

co=0
import json
with open('sacr.json') as json_file:
     dft_data = json.load(json_file)
for calculation in dft_data:
        co+=1
        atoms = Atoms(numbers=calculation['numbers'],positions=calculation['positions'],pbc=calculation['pbc'],cell=calculation['cell'])
        energy = calculation['energy'] 
        ideal_structure, info = map_structure_to_reference(atoms, primitive_structure)
        if len(info['warnings'])==0:
            sc.add_structure(structure=ideal_structure, user_tag=calculation['tag'], properties={'potential_energy':calculation['energy']/np.sum(ideal_structure.get_atomic_numbers()!=0)})

with open('hscr.json') as json_file:
     dft_data = json.load(json_file)
for calculation in dft_data:
        co+=1
        atoms = Atoms(numbers=calculation['numbers'],positions=calculation['positions'],pbc=calculation['pbc'],cell=calculation['cell'])
        energy = calculation['energy'] 
        ideal_structure, info = map_structure_to_reference(atoms, primitive_structure)
        if len(info['warnings'])==0:
            sc.add_structure(structure=ideal_structure, user_tag=calculation['tag'], properties={'potential_energy':calculation['energy']/np.sum(ideal_structure.get_atomic_numbers()!=0)})

with open('racr.json') as json_file:
     dft_data = json.load(json_file)
for calculation in dft_data:
        co+=1
        atoms = Atoms(numbers=calculation['numbers'],positions=calculation['positions'],pbc=calculation['pbc'],cell=calculation['cell'])
        energy = calculation['energy'] 
        ideal_structure, info = map_structure_to_reference(atoms, primitive_structure)
        if len(info['warnings'])==0:
            sc.add_structure(structure=ideal_structure, user_tag=calculation['tag'], properties={'potential_energy':calculation['energy']/np.sum(ideal_structure.get_atomic_numbers()!=0)})


with open('s2cr.json') as json_file:
     dft_data = json.load(json_file)
for calculation in dft_data:
        co+=1
        atoms = Atoms(numbers=calculation['numbers'],positions=calculation['positions'],pbc=calculation['pbc'],cell=calculation['cell'])
        energy = calculation['energy'] 
        ideal_structure, info = map_structure_to_reference(atoms, primitive_structure)
        if len(info['warnings'])==0:
            sc.add_structure(structure=ideal_structure, user_tag=calculation['tag'], properties={'potential_energy':calculation['energy']/np.sum(ideal_structure.get_atomic_numbers()!=0)})
with open('data_last.json') as json_file:
     dft_data = json.load(json_file)
for calculation in dft_data:
        co+=1
        atoms = Atoms(numbers=calculation['numbers'],positions=calculation['positions'],pbc=calculation['pbc'],cell=calculation['cell'])
        energy = calculation['energy']
        ideal_structure, info = map_structure_to_reference(atoms, primitive_structure)
        if len(info['warnings'])==0:
            sc.add_structure(structure=ideal_structure, user_tag=calculation['tag'], properties={'potential_energy':calculation['energy']/np.sum(ideal_structure.get_atomic_numbers()!=0)})

opt = CrossValidationEstimator(fit_data=sc.get_fit_data(key='potential_energy'), fit_method='rfe',estimator='lasso',  seed= 1,n_features=30)
opt.validate()
opt.train()
print(opt)

ce = ClusterExpansion(cluster_space=cs, parameters=opt.parameters)

ce.write("Li2FeSO.ce")
