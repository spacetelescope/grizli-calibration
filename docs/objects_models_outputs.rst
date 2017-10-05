.. grizli-tutorial documentation master file, created by
   sphinx-quickstart on Mon Aug 14 10:55:31 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

####################################
Groups, Multifits, Beams, and Traces
####################################

``grizli`` has a lot of built in classes that store models and data in useful ways. 

``GroupFLT`` Objects with ``MultiFit``
======================================

``GroupFLT`` holds the set of grism observations for one grism/direct filter
set, binding together every `flt` file. The methods include modelling functions 
that can fully model a grism extraction or serve as a starting point upon which
to build -- depending in the complexity of the data. 

Creating an instance of a ``GroupFLT`` class first requires the one-to-one grism to direct exposure 
file matching described earlier. The following performs that last matching::
    
    ## -- IMPORTS
    import os
    import glob
    import grizli
    from grizli import prep
    
    ## -- RUN
    os.chdir('Prep/')
    files = glob.glob('*_flt.fits')
    info = grizli.utils.get_flt_info(files)

    visits, filters = grizli.utils.parse_flt_files(info=info, uniquename=False, 
                                                   get_footprint=True, use_visit=True)
    pairs = grizli.utils.parse_grism_associations(visits, best_direct={'G102':'F105W', 
                                                                       'G141':'F140W'})
    direct_files = []
    grism_files = []
    
    for pair in pairs:
        # Skip G102 so that you'll only get one filter pair's worth of files
        if 'g102' in pair['grism']['product']:
            continue
        
        match = grizli.prep.find_direct_grism_pairs(pair['direct'], pair['grism'])
        direct_files.extend(list(match.keys()))
        for direct in direct_files:

            # This hopefully isn't neccesary but some matching may have gone awry
            try:
                grism_files.extend(match[direct])
                print('It ran through.')
            except (IndexError, KeyError):
                print('It did not run through.')
                print(pair['grism']['product'])

Then for the final group creation::
        
    ## -- IMPORTS
    import drizzlepac

    ## -- RUN
    # Specify seg map and catalog.
    # (Named and created in the drizzling.)
    seg = 'gd71-f105w_seg.fits'
    cat = 'gd71-f105w.cat'

    # Make the group
    grp = grizli.multifit.GroupFLT(grism_files=grism_files, sci_extn=1, 
                                   direct_files=direct_files, pad=200, 
                                   group_name='group', ref_file=None, ref_ext=0, 
                                   seg_file=seg, shrink_segimage=True, verbose=True, 
                                   cpu_count=0, catalog=cat)

From this point ``grizli`` can add models to the ``GroupFLT``. They'll get carried along as attributes 
and will overwrite each other as new models are created::

    # Make a generalized model
    grp.compute_full_model()

``grizli`` computes a model for each file during this, printing something like::
    
    ibwq1asfq_flt.fits: _compute_model Done
    ibwq1asmq_flt.fits: _compute_model Done

as it goes along -- listing each `flt` file it models. 

.. warning::
    It's not very pythonic to import ``drizzlepac`` when we don't use it. 
    Sorry about the lack of Guido-flow, but this will come crashing down should it be neglected. 

.. raw:: html 

    <img src="./_static/drizzlepac_error.png"  width="900"/>

Beams and Initial Testing
=========================

The ``GroupFLT`` class stores the original files, but for the purposes of looking at specific
sources and orders looking at cutouts is significantly more practical. ``GroupFLT`` has a method
``get_beams`` to create cutouts around a source and order.

Presently, grism orders are mapped to letters (A=1, B=0 -- the point source
that makes it to the grism image, C=2, D=3, E=-1 for G141 and G102). 

A regular analysis will likely not use more than first order A, but we're honing
the models for other orders in this example -- for those times when the higher orders come
in handy and for purposes of better contamination modles. 

To make beam object::
    
    # Specify seg id and beam order
    # (Get this from seg map created by grizli)
    id = 16
    beam_id = 'A'

    # Generate the beams in that image
    beams = grp.get_beams(id, size=24, beam_id=beam_id) 
    # A smaller size cutout often goes poorly

``beams`` contains every instance of the extracted grism cutout for the images within
``GroupFLT``, as well as the models of the grisms and contamination. This allows
for some excellent initial testing. 

To view the model and the beam cutout, and see how closely they match::

    ## -- IMPORTS
    import matplotlib.pyplot as plt
    from matplotlib import rc

    ## -- RUN
    # Set up some nice plotting labels
    rc('font', **{'family': 'serif', 'serif':['Computer Modern Roman'],'size':14})
    rc('text', usetex=True)
    
    # Extract the pieces from each beam
    beam = beams[0]
    sci = beam.grism['SCI']
    mod = beam.model
    res = sci - mod - beam.contam
    
    # Plot
    ax = plt.figure(figsize=(15,5))
    
    p1 = ax.add_subplot(1,3,1)
    p2 = ax.add_subplot(1,3,2)
    p3 = ax.add_subplot(1,3,3)

    p1.imshow(mod)
    p1.set_title('Model')
    p2.imshow(sci)
    p2.set_title('Grism Data')
    p3.imshow(res)
    p3.set_title('Residual')

    plt.show()

For the example we've been working with:

.. raw:: html 

    <img src="./_static/ex_mod_res.png"  width="900"/>

.. warning::
    Something looks pretty wrong here. The grism data 
    it a tiny spec, and logically that means subtracting the residual
    looks bonkers. This is due to runaway DQ flagging trying to mark grism
    observations as bad data, and then getting automatically incorporated by 
    ``grizli``.

    Should data look like this -- or potentially be displaced from the 
    grism trace, tweaks may be required for drizzling and matching.


Stopping the Runaway DQ Flag Train
==================================

To reset the DQ flags to exclude (or include) grisms::
    
    ## -- IMPORTS
    from astropy.io import fits
    import numpy as np

    ## -- RUN
    # Continuing on from before where grism_files was defined

    # Write an empty dq array
    with fits.open(grism_files[0]) as hdu:
        dq = np.ones_like(hdu['DQ'].data, dtype=int)*4096
        
    for infile in grism_files:
        with fits.open(infile) as hdu:
            dq &= (hdu['DQ'].data & 4096)
    
    for infile in grism_files:
        with fits.open(infile, mode='update') as hdu:
            cr = (hdu['DQ'].data & 4096) > 0
            hdu['DQ'].data[cr] -= 4096
            hdu['DQ'].data[dq > 0] |= 4096
            hdu.flush()

After resetting DQ flags and remaking group objects things should be back on track.

Creating Next Level Models with ``MultiBeam``
=============================================

The ``MultiBeam`` class facilitates the creation of more complex models ::

    # Using the beams from earlier
    mb = grizli.multifit.MultiBeam(beams, fcontam=0.2, psf=True)

``MultiBeam`` allows for beam specific models based on a user-input spectral profile
for the source. The model spectrum for most standard stars and many more should
be available from STScI through STIS/COS data and MAST retrievals. We used the 
`X-Shooter Spectral Library <http://xsl.u-strasbg.fr/>`_
in the end to download `GD-71 <ftp://ftp.eso.org/pub/stecf/standards/Xshooter/fGD71.dat>`_ .

Using spectra as a model requires renormalizing the spectra with ``pysynphot``::

    ## -- IMPORTS
    import pysynphot as S
    
    ## -- RUN
    # Direct filter
    filt = 'f141w'
    
    # Renorm the object to the given filter
    sp = S.Filespectrum('gd71-xshooter.dat')
    bp = S.ObsBandpass('wfc3,ir,' + filt)
    rn = sp.renorm(1, 'flam', bp)
    
    # Cast the pysynphot data to floats
    rn_wave = np.cast[float](rn.wave)
    rn_flux = np.cast[float](rn.flux)

And then creating the model with ``compute_model``::
    
    # Compute model
    mb.compute_model(spectrum_1d=[rn_wave, rn_flux])

    # Use new beams
    better_beams = mb.beams

``better_beams`` should now store models and grism data through 
the same ``beam.grism['SCI']`` and ``beam.model`` convention as before. 
This should provide a much closer residual. 

.. raw:: html 

    <img src="./_static/better_mod_res.png"  width="900"/>

Now with some good models established, :doc:`extraction_and_next_level_analysis`
can take place. 


.. toctree::
   :maxdepth: 2
   :hidden:
    index.rst
    data_organization
    visit_matching
    drizzling_extracting_etc
    objects_models_outputs
    extraction_and_next_level_analysis
