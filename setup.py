from setuptools import setup, find_packages

setup(
    name='tem_liquid_cell_resolution',
    version='0.1.0',
    packages=find_packages(),
    install_requires=[
        'ase',
        'abtem',
        'matplotlib',
        'numpy'
    ],
    entry_points={
        'console_scripts': [
            'create_structure=tem_liquid_cell_resolution.create_structure:main',
            'simulate_tem=tem_liquid_cell_resolution.simulate_tem:main'
        ]
    },
)
