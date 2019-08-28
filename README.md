# budgetECCO
Code to close tracer budgets offline using output from ECCO Version 4 Release 3 as described in Piecuch (2017), https://dspace.mit.edu/handle/1721.1/111094, written in the MATLAB software environment

budgetECCO

README file last updated by Christopher Piecuch, cpiecuch-at-whoi-dot-edu, Wed Aug 28 2019

Basic description

Citation

This code is for use in closing tracer (e.g., heat, mass, salt, etc.) budgets offline using output from ECCO (Estimating the Circulation and Climate of the Ocean) Version 4 Release 3.  It follows the pseudo-code descriptions given in:

Piecuch, C. G. (2017), “A Note on Practical Evaluation of Budgets in ECCO Version 4 Release 3,” http://hdl.handle.net/1721.1/111094.

Please cite that above reference when using this code. 

This code was generated to produce the main results presented in the main texts of:

Piecuch, C. G., R. M. Ponte, C. M. Little, M. W. Buckley, and I. Fukumori (2017), “Mechanisms underlying recent decadal changes in North Atlantic Ocean heat content,” J. Geophys. Res.-Oceans, 122(9), 7181-7197, https://doi.org/10.1002/2017JC012845.

Ponte, R. M., and C. G. Piecuch (2018), “Mechanisms Controlling Global Mean Sea Surface Temperature Determined From a State Estimate,” Geophys. Res. Lett., 45(7), 3221-3327, https://doi.org/10.1002/2017GL076821.

Piecuch, C. G., P. R. Thompson, R. M. Ponte, M. A. Merrifield, and B. D. Hamlington (under revision), “What Caused Recent Shifts in Tropical Pacific Decadal Sea-Level Trends,” J. Geophys. Res.-Oceans, TBD, https://doi.org/10.1002/2019JC015339 (activated upon publication).


Contents

PDF documents
•	evaluating_budgets_in_eccov4r3_updated.pdf: Outline of pseudo-code, model output needed to do the calculations, and where to get it .

Text files
•	Copyright: copyright statement
•	License: license statement

MATLAB .m files
•	eccov4r3_budget_github.m: this is the main driver code.  Simply execute “eccov4r3_budget_github” in the MATLAB Command Window from the budgetECCO directory, and this code should run “out of the box,” making use of Gael Forget’s “gcmfaces” MATLAB package (https://github.com/gaelforget).  Note that this driver code assumes (1.) you have unzipped the provided version of gcmfaces, and (2.) you have downloaded all the necessary output (see PDF document and more description below).
•	Everything in gcmfaces (see https://github.com/gaelforget/gcmfaces)

Subdirectories 
•	gcmfaces: gcmfaces code package of Gael Forget (https://github.com/gaelforget)

Note that you yourself should create the following two (initially empty) subdirectories within the main base directory
•	nctiles_monthly: Directory for holding monthly averaged nctiles model output
•	nctiles_monthly_snapshots: Directory for holding instantaneous monthly snapshots nctiles model output

Note on usage of code 
I’ve tried to make this code as user-friendly and out-of-the-box useable as possible.  The one thing I haven’t done for you is provide the ECCO Version 4 Release 3 model output necessary for doing the calculations (much too big; hundreds of GB of space would be needed). In Tables 2-4 in the included PDF document, I specify all of the model output files you’ll need to do the calculations, and where to find those outputs on the web.  All you need to do is download them, and put them in the correct “nctiles” subdirectory (see Subdirectories above), which are now empty.  The subdirectory and file naming convention is simple and straightforward, and should be apparent from the MATLAB code.  For example, for monthly averages of the ETAN variable (the elevation of the ocean’s free surface at the air-sea or air-ice interface relative to the reference ellipsoid/geoid [same thing in the model]), you would:
(1.)	Create a subdirectory called “ETAN” within the “nctiles_monthly” subdirectory.
(2.)	You’d place the 13 individual ETAN “nctiles” NetCDF files found here (https://ecco.jpl.nasa.gov/drive/files/Version4/Release3/nctiles_monthly/ETAN) or here (https://web.corral.tacc.utexas.edu/OceanProjects/ECCO/ECCOv4/Release3/nctiles_monthly/ETAN/) within that ETAN subdirectory.

Once you’ve done that for all required outputs, the code should run.  Let me know if you have any problems.




