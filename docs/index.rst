.. grizli-tutorial documentation master file, created by
   sphinx-quickstart on Mon Aug 14 10:55:31 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

#################
`grizli` Cookbook
#################


.. warning::
    ``grizli`` is still in development. We'll try to keep on top of it, but times are changing.
    This extended example was conducted with `v0.3.0`. (Go check out the regular ``grizli`` docs
    `here <http://grizli.readthedocs.io/en/master/>`_ .

``grizli`` Reduction with Minimal Tears
=======================================

This cookbook details a *fairly* general guide to WFC3/IR grism reduction
with ``grizli``. This example was for the purposes of a WFC3 calibration project
examining higher order grism extractions and generally poking against the limits of 
grism analysis with ``grizli``. This analysis was applied to two standard stars (initially 
GD-71 and later GD-153) through the G102 and G141 grism. The analyzed data was non-proprietary 
MAST data. 

PLEASE NOTE -- all of this analysis MAY OR MAY NOT be useful for your project. Follow blindly 
at your own risk. 

To avoid overload -- this cookbook has been organized into 5 sections:


 .. toctree::
    :maxdepth: 2
    
    data_organization
    visit_matching
    drizzling_extracting_etc
    objects_models_outputs
    extraction_and_next_level_analysis

.. toctree::
   :maxdepth: 3
   :hidden:
    index
    data_organization
    visit_matching
    drizzling_extracting_etc
    objects_models_outputs
    extraction_and_next_level_analysis

