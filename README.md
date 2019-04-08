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

Install numpy, scipy, matplotlib packages from SciPy - [SciPy installation instructions](https://scipy.org/install.html)
```
python -m install --user numpy scipy matplotlib
```

Install enum - [PyPi installation instructions](https://pypi.org/project/enum/):
```
pip install enum
```

### Usage

To call the program from the command-line of a bash script use this command
```
python CLI.py --kernel_path=<spice_kernel> --lbl_path=<image_lbl_file> --observing_body=<spacecraft> --target_body=<target_body_name> --target_body_frame=<Target_Reference_Frame> --file_name=<Output_file_name>
```
The output file will be of a csv file which is simply a txt file except that all data members will be comma separated. The numbers represent latitude/longitude coordinates of the image footprint on the asteroid.
Here is a sample output:
```
-2.2340040825904297,-1.2131569590397908,0.20181703414737676,-2.208470929212443,-1.238154617410746,0.1561554069037831
-2.0422742907838907,-0.2732697992565831,1.0860770236333153,-2.41196962714863,-0.5861310008846876,0.12111238545570191
-2.208470929212443,-1.238154617410746,0.1561554069037831,-2.134881035594569,-1.3405743020385485,0.011955008636660054
-0.5156179584999905,-1.2989043049391609,-0.29572687630842526,-0.48695763450715734,-1.2930375403409662,-0.3329612980149817
-0.8760975081578861,-1.4174894911479516,0.06931115224107753,-0.8091802281626755,-1.4144662738042877,0.06133346588940053
0.25716820192068435,0.684360026371913,0.702189624032953,0.20133403360933347,0.6093580423550025,0.8320790029247621
-2.058208459321073,-1.4114558717538093,-0.12342568774306055,-2.021225004659663,-1.4247769592064343,-0.09705012715923522
-1.9300945481345881,-1.455534209940224,-0.1528250354748773,-1.9032778691342371,-1.4664905853951322,-0.1944750281689534
-1.6840089435261587,-1.4735092374141607,-0.0853916140177629,-1.503010616138758,-1.4493355943362132,-0.24259920263331042
0.3408268738689623,-1.0269054578901968,0.2323814196467036,0.2544949429207982,-0.8794033165717003,0.2510301167488057
-0.7512652714782715,-1.3977641953864388,-0.011782006667912459,-0.7223457151215458,-1.389726306732751,-0.04842439574344898
-2.134881035594569,-1.3405743020385485,0.011955008636660054,-2.0845439819289946,-1.3954036973844295,-0.08036521318969041
-2.4285227428895393,-0.7781684100805244,0.316326193982164,-2.430671436849508,-0.8472321363801328,-0.12623347114271358
-2.0845439819289946,-1.3954036973844295,-0.08036521318969041,-2.058208459321073,-1.4114558717538093,-0.12342568774306055
-0.15673788780499937,0.49324990995005624,1.2782535741095418,-0.2431316284219842,0.4890925747841579,1.3967392873833244
```

## Scripts

### Spice
This script declares the Spice function wrappers that are used in the project. The functions in this script are renamed in order to improve readability.
### Util
This script provides utility functions to assist in miscellaneous tasks in the project. 
### Geometry
This script executes the ray intercept functionality. Only instruments with rectangular field of views are implemented.
### Boundary
This script constructs the footprint when given a set of ray intercepts. In order to improve accuracy at the cost of run time efficiency, lower the alpha value. 
### CLI
This script formats the command line interface of the program.

### Built With
* [SPICE Toolkit](https://naif.jpl.nasa.gov/naif/toolkit.html)
* [SpiceyPy](https://github.com/AndrewAnnex/SpiceyPy)

