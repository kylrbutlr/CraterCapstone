#   Visualizing Footprints of Space Observations
Senior Capstone project contributing to the Psyche mission. Outputs footprints of an asteroid body from loading and utilizing information from respective SPICE kernel(s).

##  Getting Started
These instructions will get you a copy of the project and allow you to run it on your local machine for development and testing purposes. 

### Prerequisites
Install or update to Python 3.7.1 or better - [Python installation instructions](https://www.python.org/downloads/):
```
install python3
```

Check that you are using Python 3.7.1 or better:
```
python --V
```

Install SpiceyPy - [SpiceyPy installation instructions](https://github.com/AndrewAnnex/SpiceyPy/blob/master/docs/installation.rst):
```
pip install spiceypy
```

Install numpy, scipy packages from SciPy - [SciPy installation instructions](https://scipy.org/install.html)
```
python -m install --user numpy scipy
```

**or run from the executable, cli.exe, which has all the dependencies loaded in already.**

### Usage

To call the program through the script from the command-line of a bash script use this command.

NOTE: The new lines are only included for readability purposes. Do not include in actual call
```
python cli.py --kernel_path=<spice_kernel_path> --lbl_path=<image_lbl_file>
--observing_body=<spacecraft> --target_body=<target_body_name> 
--target_body_frame=<target_reference_frame> --sample_size_scalar=<sample_size_scalar> 
--boundary_alpha_value=<boundary_alpha_value> --file_location =<output_file_path>
--file_name=<output_file_name>
```

#### Inputs

##### kernel_path
This is the file path of the Spice meta-kernel. 

Ex: ..\ROS_OPS.TM

##### lbl_path
This is the file path of the label file associated with the image taken by the spacecraft.

Ex: ..\ros_cam1_20150408t061457.LBL

##### observing_body
This is the spacecraft or instrument name that took the image. This should be the name of the spacecraft that is defined in the spice kernels.

Ex. ROS_NAVCAM-A

##### target_body
This is the body that the instrument is taking pictures of. This should be the name of the body that is defined in the spice kernels.

Ex. 67P/C-G

##### target_body_reference_frame
This is the reference frame of the target body. This should be the name of the frame that is defined in the spice kernels.

Ex. 67P/C-G_CK

##### sample_size_scalar
This is the value that will be multiplied by the size of the image file. The default value is 5.

Ex. 7

##### boundary_alpha_value
This is the alpha value that is used to construct the footprint from the ray intercept points. The lower the value, the more accurate the footprint will be, at the cost of run time efficiency. The default value of this is .15.

Ex. .35

##### file_location
This is the location of the output file. Be sure to include the '\' at the end of this input.
Ex. ..\Output\

##### file_name
This is the name of the output csv file containing the poly-line that describes the footprint.

Ex. output123


#### Output
The output file will be of a csv file which is simply a txt file except that all data members will be comma separated. The numbers represent latitude/longitude coordinates of the image footprint on the asteroid.
Here is the first five points of a sample output:

NOTE: The new lines are only included for readability purposes. It is not include in actual output
```
-1.8327362518869574 -1.459331365020557 -0.23076158687333326,
-1.783247584269684 -1.4870922205781625 -0.1634770855841008,
-1.6810440180945718 -1.4652817790899286 -0.12528752227700846,
-1.6203105367414006 -1.4354717174434024 -0.20363404810525032,
-1.5036997552646678 -1.449186665917284 -0.2539225945162913
```

## Scripts

### System Overview

![System Diagram](https://user-images.githubusercontent.com/21960860/55708490-edbe0600-599a-11e9-8a49-c0adfced004e.jpg)

### spice
This script declares the Spice function wrappers that are used in the project. The functions in this script are renamed in order to improve readability.
### util
This script provides utility functions to assist in miscellaneous tasks in the project. 
### geometry
This script executes the ray intercept functionality. Only instruments with rectangular field of views are implemented.
### boundary
This script constructs the footprint when given a set of ray intercepts. In order to improve accuracy at the cost of run time efficiency, lower the alpha value. 
### cli
This script formats the command line interface of the program.

### Built With
* [SPICE Toolkit](https://naif.jpl.nasa.gov/naif/toolkit.html)
* [SpiceyPy](https://github.com/AndrewAnnex/SpiceyPy)

