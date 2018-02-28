.. grizli-tutorial documentation master file, created by
   sphinx-quickstart on Mon Aug 14 10:55:31 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

#############################################
Making Drizzled Mosaics and Segmentation Maps
#############################################


Making a drizzled mosaic
========================

In order to pick out individual objects to model and extract, we need
to drizzle all of the direct exposures for one filter into an image.

Running the matching and limiting to matches with at least three files in the last
steps should assure that ONLY the correct files for each filter are in the ``Prep/`` 
folder. (Extra files simply require rematching before any additional drizzling.)

With only the correct exposures in ``Prep`` the following creates a drizzled mosaic 
for F140W ::
    
    ## -- IMPORTS
    import glob
    import os
    import grizli
    from grizli import prep
    
    ## -- RUN
    os.chdir('Prep/')
    files = glob.glob('*_flt.fits')
    info = grizli.utils.get_flt_info(files)

    # Make F140W image
    filt = 'f140w'
    group = {'files':list(info['FILE'][info['FILTER'] == filt.upper()]),
            'product':'gd71-{0}'.format(filt)} # Setting a name convention. 
            #In this case we named it gd71 after the star. 
         
    # Make drizzled image.
    grizli.prep.drizzle_overlaps([group], scale=0.06, pixfrac=0.8, skysub=False)


Making a segmentation map and catalog
=====================================

Then, the following will generate a catalog and segmentation map that ``grizli``
will use for the rest of its analysis. ::

    # Make segmentation image and catalog
    cat = grizli.prep.make_drz_catalog(group['product'], threshold=2)


Correcting seg map for proper motion
====================================

.. warning::
    With enough proper motion, the segmentation map may have two object IDs for a single source.
    Looking at the segmentation map grizli creates makes it obvious 
    if there is or isn't a unique ID for the intended source. The seg map is named by the convention
    specified in the call to ``drizzle_overlaps``, for example, the one created with this example
    is ``gd71-f140w_seg.fits``.
    
    If this has happened for a source in wont of extractraction half of the source will need to be changed
    to match the other half.
    
One segmentation map ID absorbing the other is fairly trivial with ``astropy``::
    
    ## -- IMPORTS
    from astropy.io import fits
    import numpy as np

    ## -- RUN
    infile = 'split_seg_map.fits'
    total_id = 16
    partial_id = 15
    
    # Open seg map
    with fits.open(infile) as hdu:
        # Reset all of the partial_id pixels
        seg = np.array(hdu[0].data)
        seg[seg == partial_id] = total_id
        
        # Sanity check to make sure the partial id has been fully eradicated
        print('Is there any trace of the proper-motioned object?')
        print(partial_id in seg)

        # Write the fixed seg map out
        hdu[0].data = seg

        hdu.writeto('single_source_seg.fits')

Checking by eye again to makes sure that the source is now a single object in the segmentation map.

If you'd rather not do this manually, you can also set a pixel area around the
coordinates of the target to your seg_id. (Note that you need to pick a seg_id
already in the seg_map.) ::
    ## -- IMPORTS
    from astropy.io import fits
    from astropy import wcs
    import numpy as np

    ## -- RUN
    seg_map = 'split_seg_map.fits'
       
    # Open the seg map
    with fits.open(seg_map) as hdu:
        hdr = hdu[0].header
        seg_dat = hdu[0].data
            
    # Convert the RA/DEC to pix
    w = wcs.WCS(hdr)
    ra, dec = hdr['RA_TARG'], hdr['DEC_TARG']
    pix = w.wcs_world2pix([[ra,dec]], 1)
    pix_ra, pix_dec = pix[0][0], pix[0][1]
    
    # Pick a seg_id within the region consumed
    # A 100x100 px region seemed to compensate for any error in conversion
    seg_id = np.max(seg_dat[int(pix_dec-100):int(pix_dec+100), int(pix_ra-100):int(pix_ra+100)])
    
    # Make sure there was a seg_id in that region
    if seg_id > 0:
        
        # Update the seg map and save a backup of the old one
        seg_dat[int(pix_dec-100):int(pix_dec+100), int(pix_ra-100):int(pix_ra+100)] = seg_id
    
        hdu.writeto('backup_' + seg_map, overwrite=True)
        hdu[0].data = seg_dat
        hdu.writeto(seg_map, overwrite=True)

        print('A new 200x200 px seg id around the RA/DEC will be used.')
    
    else:
        print('There was no seg_id in that region. Do it yourself.')


Next is a brief exploration of the classes that make up ``grizli`` in :doc:`objects_models_outputs`.

.. toctree::
   :maxdepth: 2
   :hidden:
    index.rst
    data_organization
    visit_matching
    drizzling_extracting_etc
    objects_models_outputs
    extraction_and_next_level_analysis
