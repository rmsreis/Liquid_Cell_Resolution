from ase import Atoms
from ase.build import bulk, molecule, add_vacuum

def create_structure():
    # Create SiNx layer
    si_nx_layer = bulk('Si', 'diamond', a=5.43)  # Example, should be amorphous
    si_nx_layer *= (5, 5, 1)  # Replicate to get desired thickness
    add_vacuum(si_nx_layer, 20)  # Add vacuum to create layer structure

    # Create water layer
    water_layer = molecule('H2O')
    water_layer *= (10, 10, 1)  # Adjust replication for desired thickness

    # Create Au nanoparticles
    au_nanoparticle = bulk('Au', 'fcc', a=4.08)  # Example, should be a nanoparticle
    au_nanoparticle *= (2, 2, 2)

    # Combine layers
    structure = si_nx_layer + water_layer + si_nx_layer
    structure += au_nanoparticle  # Add nanoparticles to the structure

    # Save the structure
    structure.center(vacuum=10)
    structure.write('tem_liquid_cell_resolution/data/structure.xyz')

def main():
    create_structure()

if __name__ == "__main__":
    main()
