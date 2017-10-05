.. grizli-tutorial documentation master file, created by
   sphinx-quickstart on Mon Aug 14 10:55:31 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

#################
Data Organization
#################

Required Data
=============

``grizli`` reduction requires a starting point of IR grism data -- which means grism exposures (in filters G102 and G141) 
and their matching direct exposures (in F105 and F098M or F140 and F160 respectively). These filters span
nearly equivalent wavelength ranges as the grism exposures they match to.

For the purposes of this example we ran two separate analyses of
`GD-71 <http://archive.stsci.edu/hst/search.php?action=Search&sci_instrume=WFC3&sci_targname=GD-71&sci_pep_id=11552,11926,11936,12333,12357,12702,13092,13579,14024,14386&sci_aper_1234=IR,G*>`_ 
and `GD-153 <http://archive.stsci.edu/hst/search.php?action=Search&sci_instrume=WFC3&sci_targname=GD-153&sci_pep_id=11552,11926,11936,12333,12357,12702,13092,13579,14024,14386&sci_aper_1234=IR,G*>`_ 
(MAST data retrieval links are hyperlinked) of which we started by preproccessing the calibrated ``FLT`` files. 

We organized the data by each target in a directory structure with a ``Prep/`` and ``RAW/`` folder. The ``RAW/`` folder
contained all of calibrated the MAST data and the ``Prep/`` folder was where the preproccessing and analysis took place.
In code snippets coming forward there will be commands directing to and reading files from these directories, and some
``grizli`` functions will depend on files being in your working directory.

Next comes :doc:`visit_matching`.



.. toctree::
   :maxdepth: 2
   :hidden:
    index.rst

