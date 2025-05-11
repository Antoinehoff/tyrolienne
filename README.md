# Tyrolienne Simulation

This repository contains a simulation and analysis of the motion of a rigid bar sliding along an inclined cable, modeled as a zipline. The project includes numerical solutions, energy analysis, and visualizations.

## Contents

- **`tyrolienne.ipynb`**: Jupyter Notebook for Python-based simulation and visualization.
- **`doc/exercice_EN.tex`**: LaTeX document containing the problem statement and detailed derivations in English.
- **`doc/exercice_FR.tex`**: LaTeX document containing the problem statement and detailed derivations in French.
- **`ExoFig/`**: Directory containing figures used in the LaTeX documents.
- **`.gitignore`**: Specifies files to ignore in version control (e.g., generated GIFs).

## Features

- Analytical derivations of equations of motion and energy conservation.
- Numerical integration of the system's equations using Python.
- Energy analysis (potential, kinetic, rotational, and mechanical).
- Visualization of the system's motion and energy evolution.
- Animation of the bar's motion along the zipline.

## Requirements

### Python
- Python 3.8 or later.
- Required libraries: `numpy`, `matplotlib`, `scipy`.

### LaTeX
- A LaTeX distribution (e.g., TeX Live, MiKTeX) with support for `exam` and `tikz` packages.

## Usage

### Python
1. Install the required Python libraries:
   ```bash
   pip install numpy matplotlib scipy
   ```
2. Open `tyrolienne.ipynb` in Jupyter Notebook.
3. Run the cells sequentially to simulate the system, analyze energy, and generate visualizations.

### LaTeX
1. Navigate to the `doc/` directory.
2. Compile `exercice_EN.tex` or `exercice_FR.tex` using a LaTeX editor or command-line tools:
   ```bash
   pdflatex exercice_EN.tex
   ```
3. View the generated PDF for the problem statement and derivations.

## Output

- **Plots**: Time evolution of angle, position, and energy.
- **Animation**: A GIF showing the motion of the bar along the zipline.
- **PDFs**: Detailed problem statements and derivations in English and French.

## Author

code: Antoine C.D. Hoffmann
exercice: Antoine C.D. Hoffmann, Louis Stenger, Paolo Ricci

## License

This project is for educational purposes and is based on an exam question from the 2022 Mechanics I course by Prof. Paolo Ricci.
