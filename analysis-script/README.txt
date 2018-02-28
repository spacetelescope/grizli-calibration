This is a script to run the entire grizli-calibration analysis from
start to finish. It will:
    1. Organize data into visit pairs.
    2. Preprocess files.
    2. Make a drizzled mosaic and segmentation map.
    3. Walk you through sources that may have been split with 
    proper motion. 
    4. Do a second one-to-one grism/direct file matching.
    5ish. Optional -- reset DQ flags that may be flagging grisms.
    5. Create beams and models.
    6ish. Optional -- renormalize a given stellar spectra to use 
    as a model input. 
    6. Write out plots that show the grism extraction and how 
    well it fits the model.

This might be a lot more than you need to just do grism data reduction --
use blindly at your own risk!

Also -- this requires some directory filled with unprocessed data (in this
cased FLT files) and will build an "analysis-name" directory in the same
location to to write the processed files and a "plots_beam-id" directory to write
figures to. 

If you run this multiple times, unless you specify a different naming 
convention your processed files and plots will be overwritten.

This script can be run interactively to combines seg-id's. If so it will wait for
a response while you check/edit your segmentation map. 


This can be executed from the command line as such:
::
    python run_full_analysis.py [-p|--path path] [-dn|--data_name data_name] 
    [-s|--spectrum spectrum] [-dq|--dq_reset| dq_reset] [-f|--filter filt] 
    [-i| --interactive_segmap]


