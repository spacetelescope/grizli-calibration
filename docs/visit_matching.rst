.. grizli-tutorial documentation master file, created by
   sphinx-quickstart on Mon Aug 14 10:55:31 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

#################################
Visit Matching and Pre-Processing
#################################

Some Background
===============

Grism and direct observations are taken in pairs of hopefully at least three exposures. These pairs are easy
enough to spot because they have the same WCS alignment. ``grizli`` requires a one-to-one match between direct and grism exposures. 
However, there may be extra data, direct exposures in multiple filters per grism, etc. Alternatively, there might not be enough data 
to run a few of these grism extractions. ``grizli`` has to make these matches before the extraction and analysis begins.

Visit Matching Tools
====================

Visit matching with ``grizli`` requires the creation of a few tables and dictionaries that are later fed into the preprocessing steps,
but these tables and dictionary structures are also a useful look at the data at hand.
First is the ``info`` table::
    
    ## -- IMPORTS
    import glob
    import os
    import grizli

    ## -- RUN
    os.chdir('Prep/')
    files = glob.glob('../RAW/*_flt.fits')
    info = grizli.utils.get_flt_info(files)
    
    # To create a smaller section with an individual filter
    red_grism = info[info['FILTER'] == 'G141']
    

Just like any table, any column can filter or sort data, but this is also useful to split G102 and G141 in the analysis.

Next the ``visits`` table::
    
    # Picking up from the created info table. 
    visits, filters = grizli.utils.parse_flt_files(info=info, uniquename=True, 
                                                   get_footprint=True)

``visits`` (assuming everything is going smoothly) consists of an ordered dictionary which will have all 
of the files that match in filter and WCS orient angle,
as well as a ``product`` (to be expanded upon soon) and a ``footprint`` which makes it easier to match visits later on. 

Finally, matching up visits provides a dictionary structure much like the last, but with pairs of grism and direct visits::
    
    # Picking up from the created visits dictionary
    matches = grizli.utils.parse_grism_associations(visits)
    
    # For some specific match reference the direct or grism visit
    match = matches[0]
    grism = match['grism']
    direct = match['direct']

Anatomy of the Product
======================

As the visits are organized and matched ``grizli`` will spit information into the terminal for each visit which will 
look something like this:: 
    GD-71-bbt-01-091.0-F105W 9
    GD-71-bll-14-105.0-G102 5

Similarly, each visit will come with a ``product`` that will look something like : ``gd-71-bll-14-105.0-g102``.
Some of the pieces of the product come from the original rootname or visit structure in which the image was taken,
whereas other are important information added by ``grizli``. For example:

- The ``gd-71`` is the object observed -- in this case the white dwarf GD-71.
- The ``105.0`` is the orient WCS angle.
- The ``g102`` is the filter.

.. warning::
    If any of these pieces of your product overlap -- i.e. above the WCS orient angle of 105 is also a filter (F105)
    -- the ``grizli`` matching will fail! We're working on a fix for now, but for the time being please do your analysis without this visit. 

Actually Matching and Preprocessing Visits 
==========================================

The following code matches visits and completes the intial preprocessing.  

.. warning::
    This will break if ``grizli`` cannot find the
    grism specific calibration data it needs in a pre-defined ``iref`` and ``jref`` folder, as described
    `here <http://grizli.readthedocs.io/en/master/grizli/install.html>`_ .

When matching and preprocessing, it can often be a good idea to work with only one direct/grism filter pair at a time, 
and the following creates matches and preprocess ``G141`` and ``F140W``::
    
    ## -- IMPORTS
    import glob
    import os
    import grizli
    from grizli.prep import process_direct_grism_visit

    ## -- RUN
    # Read in raw files and organize into grism visits
    files = glob.glob(os.path.join(raw_path), '*_flt.fits')
    info = grizli.utils.get_flt_info(files)
    
    # Here we'll sort for G141
    grism = 'G141'
    direct = 'F140W'
    
    # Collect in pairs
    color_filter = (info['FILTER'] == direct) | (info['FILTER'] == grism)
    grism_visits, filters = grizli.utils.parse_flt_files(info=info[color_filter], 
                                                          uniquename=True, 
                                                          get_footprint=True)
    grism_matches = grizli.utils.parse_grism_associations(grism_visits)

    # Run the preprocessing
    for pair in grism_matches:
        direct_visit = pair['direct']
        grism_visit = pair['grism']

        # Check that there are enough files
        if (len(direct_visit['files']) > 3) and (len(grism_visit['files']) > 3):
            process_direct_grism_visit(direct=direct_visit, grism=grism_visit, 
                                       align_mag_limits = [14,23])

The ``process_direct_grism_visit_function`` can also take an ``radec`` file (or a ``.dat`` file with RA and DEC coordinates
for the sources in the images) but leaving this blank will simply have ``grizli`` create this for you from ``gaia``, ``pan-STARRS``,
or ``WISE`` -- whichever finds the source first. 

.. warning::
    It's quite easy to set ``direct_visit`` and ``grism_visit`` to the same thing. This results in the error::
        
        UnboundLocalError: local variable 'bg_fixed' referenced before assignment 

This, having completed successfully completed, will write pre-processed (aligned and drizzled images) to
the ``Prep`` folder. Next is :doc:`drizzling_extracting_etc`.



.. toctree::
   :maxdepth: 2
   :hidden:
    index.rst
    data_organization
    visit_matching
    drizzling_extracting_etc
    objects_models_outputs
    extraction_and_next_level_analysis
