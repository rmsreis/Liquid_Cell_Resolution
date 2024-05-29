import abtem
from abtem import Potential, Wave
import matplotlib.pyplot as plt

def simulate_tem():
    # Load the structure
    atoms = abtem.io.read('tem_liquid_cell_resolution/data/structure.xyz')

    # Create potential from the structure
    potential = Potential(atoms, sampling=0.1)  # Adjust sampling as needed

    # Define the electron wave
    wave = Wave(shape=potential.shape, sampling=potential.sampling, energy=200e3)  # 200 keV

    # Propagate the wave through the potential
    wave.multislice(potential)

    # Get the intensity
    image = wave.intensity

    # Save and plot the image
    plt.imshow(image, cmap='gray')
    plt.colorbar(label='Intensity')
    plt.title('TEM Image')
    plt.savefig('tem_liquid_cell_resolution/data/results/tem_image.png')
    plt.show()

def main():
    simulate_tem()

if __name__ == "__main__":
    main()
