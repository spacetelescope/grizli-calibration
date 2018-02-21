This series of scripts will run the entire analysis 
for the grizli calibration project -- data-preprocessing 
through plotting throughput.

It requres :
1. Up to date grizli installation, config/calib files, and exported paths.
2. Up to date numpy and astropy.
3. WFC3/IR grism data in the form of FLT files.  (Don't do more than one 
field in the same directory please.)

To run : 

>>> python run_full_analysis.py path/to/raw/files

To run WITHOUT correcting DQ flags : 

>>> python run_full_analysis.py path/to/raw/files dq_reset=False

