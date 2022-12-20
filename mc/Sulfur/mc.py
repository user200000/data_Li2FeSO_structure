
from mchammer.ensembles import CanonicalAnnealing, CanonicalEnsemble
from icet import ClusterExpansion
from mchammer.calculators import ClusterExpansionCalculator
import numpy as np
replica_count=0
for seed in [42,66,3,19,12]:
        replica_count+=1
        temp=1025
        ce=ClusterExpansion.read('Li2FeSO.ce')
        primitive_structure=ce.primitive_structure
        primitive_structure[2].symbol="Li"
        primitive_structure[3].symbol="Li"
        print((primitive_structure))
        def generate_randomized_supercell(primitive_structure, size):
                supercell = primitive_structure.repeat(size)
                an=supercell.get_atomic_numbers()
                return supercell

        supercell_structure = generate_randomized_supercell(primitive_structure, 8)
        calc = ClusterExpansionCalculator(supercell_structure, ce)


        mc=CanonicalAnnealing(supercell_structure,calc,20000,temp,500000,random_seed=seed)
        dcmc=mc.data_container

        mc.run()

        from mchammer.ensembles import CanonicalEnsemble
        calc = ClusterExpansionCalculator(dcmc.get_trajectory()[-1], ce)
        mc=CanonicalEnsemble(dcmc.get_trajectory()[-1], calculator=calc,temperature=temp,random_seed=seed)
        mc.run(1000000)
        dcmc=mc.data_container
        dcmc.write('Li2FeSO_'+"temp"+str(temp)+'_rep='+str(replica_count)+'.mc')
