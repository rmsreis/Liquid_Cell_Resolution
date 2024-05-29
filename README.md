# TEM Liquid Cell Resolution Plotter

This project aims to plot the resolution of Transmission Electron Microscopy (TEM) in liquid cells. It considers various imaging modes, solvents, membrane thicknesses, and aberration corrections.

## Features
- Supports TEM and STEM modes
- Includes aberration correction up to the 3rd order
- Plots dose-limited resolution for different solvents and membrane thicknesses

## Installation

```sh
git clone <repository-url>
cd TEM_Liquid_Cell_Resolution
python -m venv venv
source venv/bin/activate  # On Windows use `venv\Scripts\activate`
pip install -r requirements.txt
python setup.py install
