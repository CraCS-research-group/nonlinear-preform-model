# nonlinear-preform-model
This repository contains a subroutine and parameter identification scripts for modelling nonlinear material properties in composite process simulations using ABAQUS/Explicit. The subroutine can be used to model an asymmetric constitutive law for biaxial textile materials to represent a non-constant bending stiffness, and the parameter identification scripts can be used to find parameter inputs to ABAQUS' builtin cohesive laws for modelling a nonlinear transverse shear behaviour of a whole preform consisting of multiple layers of fabric.

The code was developed at **Department of Materials and Production, Aalborg University** and at **Bristol Composites Institute, University of Bristol** by Peter H. Broberg, Adam J. Thompson, Jonathan P-H. Belnoue, Stephen R. Hallett, Esben Lindgaard, and Brian L.V. Bak.

The implementation of the subroutine and scrips are described in the publication: \
*Broberg PH, Lindgaard E, Thompson AJ, Belnoue JP-H, Hallett SR, and Bak BLV (2003). "An accurate forming model for capturing the nonlinear material behaviour of multilayered binder-stabilised fabrics and predicting fibre wrinkling" (Manuscript submitted)*

The subroutine is based on the code for the fibre tracking algorithm and hypoelastic constitutive law described in the publication: \
*Thompson AJ, Belnoue JP-H, and Hallett SR (2020). "Modelling defect formation in textiles during the double diaphragm forming process", Composites Part B: Engineering, 202:108357. doi:https://doi.org/10.1016/j.compositesb.2020.108357* \
This subroutine can be found on the repository: https://bristolcompositesinstitute.github.io/HypoDrape/

Please feel free to use and adapt the codes, but remember to give proper attribution.

## Features 

### VUMAT subroutine - vumat_hypodrape_asymmod:
This subroutine can be used to model the behaviour of fabric materials in ABAQUS/Explicit. The subroutine is based on the hypoelastic material law and fibre tracking algorithm from https://bristolcompositesinstitute.github.io/HypoDrape/. The constitutive relation in the fibre direction is described with a multilinear stress-strain relation that can be used to model nonlinear bending behaviour. See the usage below for more information. 

### Python Files for Parameter Identification: 
The two Python files and the module can be used to obtain the parameters needed to run the subroutine and the interface law are included in this repository. The Python file bending.py can be used to obtain the necessary parameters for the VUMAT subroutine based on an inputted moment-curvature relation. The Python file shear.py can be used to obtain the table with damage evolution based on the effective separation needed to model preform transverse shear stiffness with the ABAQUS in-built cohesive interface.

### Examples
Two example files are included in this repository. They are used to exemplify how to use the VUMAT subroutine for modelling a non-constant bending stiffness and the cohesive interface with tabular damage values for modelling a non-constant shear stiffness. The example model_cantilever.inp is a model of the cantilever bending test of the fabric material in the reference. The example model_transverse_shear.inp is a model of the transverse shear test of the preform described in the references. 

## Usage
Below are detailed instructions on how to use the components in the repository.

### VUMAT Subroutine
The subroutine can be used directly in ABAQUS Explicit. 
9 state variables and 21 material constants are needed in the model. A description of the material constants and state variables is given below. 
#### Material properties:
1-3: Initial fibre direction \
4: Tensile stiffness in direction 1 \
5: Tensile stiffness in direction 2 \
6-11: Parameters for the in-plane shear stiffness \
12: Compressive stiffness in direction 1 \
13: Softening stiffness in direction 1 \
14: Softening strain in direction 1 \
15: Locking strain in direction 1 \
16: Compressive stiffness in direction 2 \
17: Softening stiffness in direction 2 \
18: Softening strain in direction 2 \
19: Locking strain in direction 2 \
20: Strain scaling \
21: Viscous damping

#### State variables are stored as:
1-4: Stresses in fibre directions \
5: Shear angle \
6: Strain in fibre direction 1 \
7: Strain in fibre direction 2 \
8: Stiffness in direction 1 at previous increment \
9: Stiffness in direction 2 at previous increment

An example of how to include the user material in your model is given below:

    *Depvar
     9,
    *User Material, constants=21
    1, 0, 0.0, 1.4e9, 1.4e9, 0.0, 0.0, 0.0
    0.0, 0.0, 3.0e9, 1.4e9, 4.9e6, 0.0006, 0.24, 1.4e9
    4.9e6, 0.0006, 0.24, 1000000, 0.0001

### Python Files for Parameter Identification
Two Python files are included for parameter identification. 

The file bending.py can determine the parameters needed to run the VUMAT file from the three variables below:
- t: Thickness of the ply
- E_t: Tensile modulus Pa 
- test_name: Name of the file containing the characterised curve

Furthermore, the following parameters can be adjusted to obtain a better fit to the characterised curve:
- start_fit: Entry to start fitting the moment-curvature relation. The default value is 90
- scale_stiff: Scale on the stiffness in the least squares fit. The default value is 1e6 
- scale_strain: Scale on the strain in the least squares fit. The default value is 1e-5
- x0: Starting guess for the least squares fit. The default value is [1.0, 1.0]

Ensure that the characterised moment-curvature relation is in .csv file and located in the correct folder, and then press run to obtain the desired parameters.

The file shear.py can determine the parameters needed to run the inbuilt cohesive law with the parameters below:
- h: Thickness of an individual ply

- tau_m: Mode-I onset traction
- Gc: Mode-I critical energy release rate
- KI: Mode-I initial stiffness

- test_name: Name of the file containing the characterised shear curve
- slip_end: Mode-II separation at failure - see references for details
- lin_start: Mode-II separation at max load - see references for details
- delta_t: Mode-II separation from max load used to determine slope - see references for details

Ensure that the characterised shear stress-strain relation is in .csv file and located in the correct folder, and then press run to obtain the desired parameters and the table for the damage variables.
#### Packages needed for the parameter identification (version used during the development of the coding)
- python (3.8.8)
- numpy (1.20.1)
- matplotlib (3.3.4)
- pandas (1.2.4)
- cv2 (4.0.1)
- scipy (1.6.2)

### Examples in .inp Format
The examples are a model of a cantilever bending test of the fabric (model_cantilever.inp) and a transverse shear test of the preform (model_transverse_shear.inp). Both models can be run with the subroutine in this repository. 

To run the examples through the command prompt, copy the VUMAT subroutine to the examples folder and use the following commands in the prompt. For the cantilever example use:

    abaqus job=model_cantilever user=vumat_hypodrape_asymmod -interactive

For the transverse shear test example use:

    abaqus job=model_transverse_shear user=vumat_hypodrape_asymmod -interactive

## Documentation
The modules contain docstring and comments describing the functionalities, input variables and outputs.

## Licensing
The software in this repository is distributed under the MIT license. Please see the license file. Additionally, the license for the hypodrape subroutine (https://bristolcompositesinstitute.github.io/HypoDrape/), which the vumat_hypodrape_asymmod subroutine is based on, is located in the src folder.
