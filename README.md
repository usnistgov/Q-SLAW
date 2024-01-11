# qDIC-based Strain Localization Analysis with Wavelets (Q-SLAW)

## Description and Purpose
This dataset is an archival version of an algorithm implementation for strain localization analysis in digital image correlation data. The original application is to detect the occurrence of localization (i.e., collective buckling, hotspots, etc.) in uniaxial stress testing of polymer foams instrumented with 2D digital image correlation. The full-field strain maps are processed with a Marr wavelet and regions of excess strain are identified. These are used to compute an overall localization intensity factor with units of strain, which quantifies the degree of excess strain in localization regions above the nominal background. A Monte Carlo-like (MC) simulation mode is included in the scripts to provide performance validation for the metric under procedurally generated fields with known localization. It is set up for using input data specifically from qDIC (https://github.com/FranckLab/qDIC) but could be straightforwardly adapted to other DIC output data. A run-script with example data from mechanical testing of an open-cell elastomeric foam is included.

## Content and dependancies

The code consists of several scripts, all of which are in the root directory of the project. Subfolders contain example input data and outputs.

### Main runnable function
- `localization_detection_testing` This script that sets up and runs the localization detection in MC mode to test algorithm perfomance.
- *`run_localization_detection_on_foam_data`* This is the main script to used to run the localization detction on experimental data.

### Supporting functions:
- `calculateEij_2d` Used to compute the 2d full-field strain tensors from the deformation gradiant.
- `calculateFij_2d` Used to compute the 2d full-field deformation gradiant tensor from displacement field.
- `compute_localization_residual` Used to compute the least-squares error for the known test cases in the MC mode.
- `gen_synth_localized_strain_maps` Used to procedurally generate synthetic fields with known localization for the MC mode.
- `gradientN` Used to cpmpute the N-dimensional gradient of a field.
- `inc2cum` Used to convert incremental-mode qDIC to cumulative data using a finite-deformation mathematics.
- *`localization_detection`* This is the core function to extract localization information from the displacement/strain field.
- `plotting_qdic_results` Basic visualization tool for qDIC data.

### Directories
- `example_dic_data` contains two example datasets for a nominally 144 kg/m^3 (9 lb/ft^3) impact protection foam.
- `output` contains saved results.

### Dependancies

Core Matlab (developed and tested on R2023B for Win10 x64) with three toolboxes is required. These toolbxes are:
- Image Processing Toolbox: multiple uses for filtering and image handling.
- Wavelet Toolbox: for the wavelet analysis.
- Signal Processing Toolbox: minor usage for peak detection. This could straightforwardly be remove with only minor impact on performance.

## Running the Q-SLAW code
See also the `HowToRun_LocalizationAnalysis.txt` document.

Briefly, the major steps are as follows:
1. Prepare input data from DIC. If you use qDIC the "results" .mat file should have all you need. For other code you'll need to either reformat data into the qDIC format or adjust the run-script's data ingestion. 
2. Update user-changable parameters (data direcotry, threshold, strain filter).
3. Run the script and follow the prompts.
4. Save output, visualize results.

## Other links
For the data used in the development and calibration of the model see: <add data doi here>.

The initial publication this is develop for is: J Tao**, AK Landauer∗∗, Z Yan∗∗, X Li, C Franck, DL Henann "Large-deformation constitutive modeling of viscoelastic foams: Application to a open-cell foam material" (In preperation)

Please cite as:
> [add code citation/doi here].


## Contact and support
For questions, please open a new entry in the "Issues" tab in GitHub. If needed, you can also find authors' contact information on the data description website. 

Alexander Landauer (NIST, Material Measurement Lab, Materials Measurement Science Division) is the developer of the code.

## License

### NIST Software Licensing Statement

NIST-developed software is provided by NIST as a public service. You may use, copy, and distribute copies of the software in any medium, provided that you keep intact this entire notice. You may improve, modify, and create derivative works of the software or any portion of the software, and you may copy and distribute such modifications or works. Modified works should carry a notice stating that you changed the software and should note the date and nature of any such change. Please explicitly acknowledge the National Institute of Standards and Technology as the source of the software.

NIST-developed software is expressly provided "AS IS." NIST MAKES NO WARRANTY OF ANY KIND, EXPRESS, IMPLIED, IN FACT, OR ARISING BY OPERATION OF LAW, INCLUDING, WITHOUT LIMITATION, THE IMPLIED WARRANTY OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NON-INFRINGEMENT, AND DATA ACCURACY. NIST NEITHER REPRESENTS NOR WARRANTS THAT THE OPERATION OF THE SOFTWARE WILL BE UNINTERRUPTED OR ERROR-FREE, OR THAT ANY DEFECTS WILL BE CORRECTED. NIST DOES NOT WARRANT OR MAKE ANY REPRESENTATIONS REGARDING THE USE OF THE SOFTWARE OR THE RESULTS THEREOF, INCLUDING BUT NOT LIMITED TO THE CORRECTNESS, ACCURACY, RELIABILITY, OR USEFULNESS OF THE SOFTWARE.

You are solely responsible for determining the appropriateness of using and distributing the software and you assume all risks associated with its use, including but not limited to the risks and costs of program errors, compliance with applicable laws, damage to or loss of data, programs or equipment, and the unavailability or interruption of operation. This software is not intended to be used in any situation where a failure could cause risk of injury or damage to property. The software developed by NIST employees is not subject to copyright protection within the United States.

### Disclaimers
Certain commercial equipment, instruments, or materials are identified in this paper in order to specify the experimental procedure adequately. Such identification is not intended to imply recommendation or endorsement by NIST, nor is it intended to imply that the materials or equipment identified are necessarily the best available for the purpose.

The opinions, recommendations, findings, and conclusions in the manuscript do not necessarily reflect the views or policies of NIST or the United States Government.