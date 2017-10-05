.. grizli-tutorial documentation master file, created by
   sphinx-quickstart on Mon Aug 14 10:55:31 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

#######################
Extraction and Analysis
#######################

After creating ``beams``, the ``optimal_extract`` method provides the best extraction
of the grism data. 

``optimal_extract`` for Your Extraction Needs
===========================================
    
The following extracts the flux, wavelength, and error profile for a given grism data cutout::

    # Continuing on from the beams in place
    wavelength, flux, err = beam.beam.optimal_extract(beam.grism['SCI'] - beam.contam, 
                                                      ivar=beam.ivar)

``optimal_extract`` can extract models as well as data::
    
    # And now, models!
    wave_mod, flux_mod, error_mod = beam.beam.optimal_extract(beam.model, ivar=beam.ivar)

An example of the model (in green), data in black, and a flux ratio of the two with wavelength follows.

 .. raw:: html 

    <img src="./_static/beam_A.png"  width="900"/>
   
.. toctree::
   :maxdepth: 2
   :hidden:
    index.rst
    data_organization
    visit_matching
    drizzling_extracting_etc
    objects_models_outputs
    extraction_and_next_level_analysis
