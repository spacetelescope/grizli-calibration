"""
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

Authors
-------
    Jules Fowler, Gabe Brammer, 2018
   
Use
---
    
    This can be executed from the command line as such:
    ::
        python run_full_analysis.py [-p|--path path] [-dn|--data_name data_name] 
        [-s|--spectrum spectrum] [-dq|--dq_reset| dq_reset] [-f|--filter filt] 
        [-i| --interactive_segmap]

Notes
-----
    This might be a lot more than you need to just do grism data reduction --
    use blindly at your own risk!

    Also -- this requires some directory filled with unprocessed data (in this
    cased FLT files) and will build an "analysis" directory in the same
    location to to write the processed files and a "plots" directory to write
    figures to. 

    If you run this multiple times, unless you specify a different naming 
    convention your processed files and plots will be overwritten.

    This script can be interactive where it combines seg-id's. If so it will wait for
    a response while you check/edit your segmentation map. 

    Grizli provides A LOT of output. Right now this gets written to a log file
    so that what's printed to the screen is a little neater. 

"""

## -- IMPORTS
import argparse
import glob
import os
import sys
import time

from astropy import wcs
from astropy.io import fits
import drizzlepac
import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np
import pysynphot as S

import grizli
from grizli.prep import process_direct_grism_visit

## -- FUNCTIONS

def main(path, data_name, spectrum='skip', dq_reset=True, interactive=False, filt='best'):
    """ Main wrapper function to run all the grizli
    analysis.

    Parameters
    ----------
    path: str
        Path to the "raw" data files.
    data_name : str
        The name (target, nickname, whatever you please)
        of the data to write files to. 
    spectrum: str, optional
        The data file of a stellar spectra to read in.
    dq_reset: bool, optional
        Whether or not to reset DQ flags. Set to true.
    filt: str, optional 
        The grism/direct filter match to use. Set to best
        match, but otherwise needs list in form of
        "direct1:grism1,direct2:grism2". 
    """
    
    # First some testing to make sure everything is in shape.
    filepath = os.path.join(path, '*flt.fits')
    if len(glob.glob(filepath)) == 0:
        print('STOP THE PRESSES! You have no flt files to work with.')
        return
    
    # Assign filter pairs
    if filt == 'best':
        direct = ['F140W', 'F105W']
        grism = ['G141', 'G102']

    else:
        splits = filt.split(',')
        direct, grism = [], []
        direct = [split.split(':')[0] for split in splits]
        grism = [split.split(':')[1] for split in splits]
    
    for index, direct_filt in enumerate(direct):
        grism_filt = grism[index]
        filt_in = direct_filt.lower()
        grism_in = grism_filt.lower()

        # Navigate to proper directory
        os.chdir(path)
        os.chdir('..')
        if not os.path.exists('analysis-{}'.format(grism_in)):
            os.makedirs('analysis-{}'.format(grism_in))
        os.chdir('analysis-{}'.format(grism_in))
        
        # Preprocess
        preprocess(filepath, direct_filt, grism_filt)
        
        # Make seg map and drizzle mosaic
        make_mosaic_and_seg(filt_in, data_name)

        # Correct seg map
        name = '{0}-{1}'.format(data_name, filt_in)
        
        source_id = fix_seg_map(name, interactive)

        # Initial one-to-one file matching
        direct_files, grism_files = direct_grism_matching(grism_filt)
        
        # DQ flagging
        if dq_reset:
            reset_dq_flags(grism_files)
    
        # Making groups/beams
        grp = create_group(direct_files, grism_files, name)

        # Modelling -- needs to loop through each fringe/beam_id
        beam_ids = ['A', 'B', 'C', 'D']
        for beam_id in beam_ids:
            try:
                beams = add_models_and_beams(grp, source_id, beam_id, grism[index].lower(), spectrum=spectrum)
            
                # Plotting
                if not os.path.exists('plots'):
                    os.makedirs('plots')
                os.chdir('plots')
                plot_traces_and_beams(beams, name, beam_id)
            
            except ValueError: #IndexError) as e:
                print('Looks like your observation is not going to support fringe order {}.'.format(beam_id))
        

def plot_traces_and_beams(beams, name, beam_id):
    """Makes a 3-part subplot for each model/grism/residual
    and a plot of each ratio.

    Paramters
    ---------
    beams : list of beam object
        List of grizli object holding cutout data.
    name : str
        The name convention to use for plots.
    beam_id : str
        The beam_id we're plotting.
    """

    rc('font', **{'family': 'serif', 'serif':['Computer Modern Roman'],'size':14})
    rc('text', usetex=True)
    dgreen = '#117733'

    # Plot the residual-model-science cutouts.
    for index, beam in enumerate(beams):
        sci = beam.grism['SCI']
        mod = beam.model
        res = sci - mod - beam.contam

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
        
        plt.savefig('{0}-beam{1}_{2}.png'.format(name, beam_id, index))
        plt.clf()

    # Plot the ratios
    fig = plt.figure(figsize=[8,8])
    
    ax_flux = fig.add_subplot(211)
    ax_ratio = fig.add_subplot(212)
    
    ymax = 0
     
    for beam in beams:
        ww, ff, ee = beam.beam.optimal_extract(beam.grism['SCI'] - beam.contam, ivar=beam.ivar)
        ax_flux.errorbar(ww/1.e4, ff, ee, color='k', alpha=0.1, marker='.', linestyle='None')
    
        wm, fm, em = beam.beam.optimal_extract(beam.model, ivar=beam.ivar)
        ax_flux.plot(wm/1.e4, fm, color=dgreen, alpha=0.2, zorder=10)
        ymax = np.maximum(ymax, fm[np.isfinite(fm)].max())
        
        # Ratio
        ax_ratio.plot(wm/1.e4, ff/fm, color=dgreen, alpha=0.2)
    
    ax_ratio.set_ylim(0,2)
    ax_flux.set_xticklabels([])
    ax_flux.set_ylim(-0.2*ymax, 1.2*ymax)
    
    for ax in fig.axes:
        ax.grid()
      #  ax.set_xlim(0.95, 1.75)
    
    fig.tight_layout(pad=0.5)
    plt.savefig('{0}-beam{1}-ratio.png'.format(name, beam_id))
    print('Plot plotted.')
    plt.clf()

def add_models_and_beams(grp, source_id, beam_id, grism_filt, spectrum):
    """ Add models and create beam cutouts.

    Parameters
    ----------
    grp : GroupFLT object
        The GroupFLT grizli object. A collection of files,
        models, extensions, etc.
    source_id : int
        The the seg id of the source that needs extracting.
    beam_id : str
        The fringe/beam_id to model.
    grism_filt : str
        The grism filter. 
    spectrum : str
        Either 'skip' to not do a more complex model, or the
        path to a data file with a source spectra.

    Returns
    -------
    beams : list of beam objects
        List of beam objects, which is like a GroupFLT but 
        centered around these beam cutouts tracing the 
        grism.
    """

    # Start with general model
    grp.compute_full_model()

    # Build beam objects 
    beams = grp.get_beams(source_id, size=24, beam_id=beam_id)
    print('-----------------------------------------------------------')
    print(beams)
    print(len(beams))
    print('-----------------------------------------------------------')
    if spectrum=='skip':
        return beams

    else:
        
        mb = grizli.multifit.MultiBeam(beams, fcontam=0.2, psf=True)
        
        # Read in alternate data
        sp = S.FileSpectrum(spectrum)
        bp = S.ObsBandpass('wfc3,ir,' + grism_filt)
        rn = sp.renorm(1, 'flam', bp)
        
        # Cast to floats
        rn_wave = np.cast[float](rn.wave)
        rn_flux = np.cast[float](rn.flux)

        mb.compute_model(spectrum_1d=[rn_wave, rn_flux])
        beams = mb.beams

        return beams

    
def create_group(direct_files, grism_files, name):
    """ Build the group object that will allow for
    extraction/modelling/plotting.

    Parameters
    ----------
    direct_files: list of str
        List of direct files to include in the group.
    grism_files: list of str
        List of grism files to include in the group.
    name: str
        the name convention for the seg map/cat.
    
    Returns
    -------
    grp : GroupFLT object
        The GroupFLT grizli object. A collection of files,
        models, extensions, etc.
    """

    seg = '{}_seg.fits'.format(name)
    cat = '{}.cat'.format(name)

    grp = grizli.multifit.GroupFLT(grism_files=grism_files, sci_extn=1,
            direct_files=direct_files, pad=200, group_name='group',
            ref_file=None, ref_ext=0, seg_file=seg, shrink_segimage=True,
            verbose=True, cpu_count=0, catalog=cat)

    return grp


def direct_grism_matching(grism_filt):
    """Do the one-to-one direct-grism file matching.
    (This may look a lot like the preprocessing but 
    the files need to be re-sorted as they are 
    newly written after they're preprocessed.)
    
    Parameters
    ----------
    grism_filt : str
        The grism filter.
    
    Returns
    -------
    direct_files : list of str
        List of file names with the direct filter.
    grism_files : list of str
        List of file names with the grism filter.
    """

    files = glob.glob('./*_flt.fits')
    info = grizli.utils.get_flt_info(files)
    filt_info = info[info['FILTER'] == grism_filt]
    
    visits, filters = grizli.utils.parse_flt_files(info=info, uniquename=False, 
                                                   get_footprint=True, use_visit=True)
    pairs = grizli.utils.parse_grism_associations(visits,
                                                  best_direct={'G102':'F105W', 'G141':'F140W'})

    direct_files, grism_files = [], []
    for pair in pairs:

        match = grizli.prep.find_direct_grism_pairs(pair['direct'], pair['grism'])
        direct_files.extend(list(match.keys()))
    
        for direct in direct_files:
            try:
                grism_files.extend(match[direct])
            except (IndexError, KeyError):
                print('This grism: {} got a little turned around in the matching.'.format(pair['grism']['product']))

    return direct_files, grism_files


def reset_dq_flags(grism_files):
    """ Resets DQ flags that may have accidentally included 
    grisms.

    Paramters
    ---------
    grism_files : list of str
        The list of filesnames to run through and correct.
    """

    # Open files to make empty of the right size.
    with fits.open(grism_files[0]) as hdu:
        dq = np.ones_like(hdu['DQ'].data, dtype=int)*4096
    
    # Find all the points set to 4096 -- CR flag
    for infile in grism_files:
        with fits.open(infile) as hdu:
            dq &= (hdu['DQ'].data & 4096)
    
    # Reset those flags to zero.
    for infile in grism_files:
        with fits.open(infile, mode='update') as hdu:
            cr = (hdu['DQ'].data & 4096) > 0
            hdu['DQ'].data[cr] -= 4096
            hdu['DQ'].data[dq > 0] |= 4096
            hdu.flush()


def make_mosaic_and_seg(filt, name):
    """ Make drizzled mosaic and segmentation map.

    Parameters
    ----------
    filt : str
        Direct filter.
    name : str
        Naming convention.
    """

    files = glob.glob('*_flt.fits')
    info = grizli.utils.get_flt_info(files)

    group = {'files': list(info['FILE'][info['FILTER'] == filt.upper()]), 'product': name + '-' + filt}
    grizli.prep.drizzle_overlaps([group], scale=0.06, pixfrac=0.8, skysub=False)
    cat = grizli.prep.make_drz_catalog(group['product'], threshold=2)
    
    print('Drizzled mosaic and seg map written under the name: ' + name)

def fix_seg_map(name, interactive):
    """ Runs an interactive back and forth to check and
    correct the segmentation map if it's split sources due
    to proper motion of the star.

    Parameters
    ----------
    name : str
        The full (data + filter) name the seg is written
        to.

    interactive : bool/str
        Whether or not run it automatically or interactively. 

    Returns
    -------
    seg_id : int
        The seg id of the new combined source.
    """

    seg_map = name + '_seg.fits'
    
    if interactive == "True" or interactive == True:
        # Check if we need to fix up the seg map.
        print('Please take a moment to take a peep at the segmentation map written to: ' + seg_map)
        single_source = str(input('Does the seg map require reformatting for the source you intend to extract? [y/n] \n'))
    
        if 'n' in single_source: 
            seg_id = int(input('Great! In that case, we just need the seg id of the source you hope to extract. \n'))

        elif 'y' in single_source:
            print('Why do bad things happen to good seg maps?')
            doubles = str(input('Give the two seg ids, comma separated, and they will be combined into the first. \n'))
            seg_consume = int(doubles.split(',')[0])
            seg_consumed = int(doubles.split(',')[1])
            correct_seg_map(seg_map, seg_consume, seg_consumed)
            seg_id = seg_consume

        else:
            print('Looks like the input is off. Try again? \n')
            fix_seg_map(name, interactive)
            return 
    
    else: 
        
        # Open the seg map
        with fits.open(seg_map) as hdu:
            hdr = hdu[0].header
            seg_dat = hdu[0].data
        
        # Convert the RA/DEC to pix
        w = wcs.WCS(hdr)
        ra, dec = hdr['RA_TARG'], hdr['DEC_TARG']
        pix = w.wcs_world2pix([[ra,dec]], 1)
        pix_ra, pix_dec = pix[0][0], pix[0][1]
        
        seg_id = np.max(seg_dat)*2
        new_id = np.ones((200,200))*seg_id
        seg_dat[int(pix_dec-100):int(pix_dec+100), int(pix_ra-100):int(pix_ra+100)] = new_id
         
        #hdu.writeto('backup_' + seg_map, overwrite=True)
        hdu[0].data = seg_dat
        hdu.writeto(seg_map, overwrite=True)

        print('A new 200x200 px seg id around the RA/DEC will be used.')
        print('The old seg map was saved to: backup_' + seg_map)
    
    return seg_id


def correct_seg_map(seg_map, seg_consume, seg_consumed):
    """Corrects the segmentation map.

    Parameters
    ----------
    seg_map : str
        The name of the seg map.
    seg_consume: int
        The source id to absorb the other.
    seg_consumed: int
        The source id that gets absorbed.
    """

    with fits.open(seg_map) as hdu:
        seg = np.array(hdu[0].data)
        seg[seg == seg_consumed] = seg_consume
        
        test = (seg_consumed in seg)
        if ~test:
        
            hdu.writeto('backup_' + seg_map, overwrite=True)
            hdu[0].data = seg
            hdu.writeto(seg_map, overwrite=True)
            
            print('The seg map was corrected.')
            print('The old seg map was saved to: backup_' + seg_map)

        else:
            print('Something went a little funny. You will need to correct it manually.')
            wait = input('Hit any key once you have corrected the seg map. \n')


def parse_args():
    """ This arg does a parse. """

    # Add arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-p' '--path',
                        dest='path',
                        action='store',
                        required=True,
                        help='The path to the "raw" data files.')
    parser.add_argument('-dn' '--data_name',
                        dest='data_name',
                        action='store',
                        required=True,
                        help='The naming convention for plots/data.')
    parser.add_argument('-s' '--spectrum',
                        dest='spectrum',
                        action='store',
                        required=False,
                        default='skip',
                        help='Path to stellar model if you would rather not rely on CDBS.')
    parser.add_argument('-dq' '--dq_reset',
                        dest='dq_reset',
                        action='store',
                        required=False,
                        default=True,
                        help='Whether or not to reset DQ flags.')
    parser.add_argument('-f' '--filter',
                        dest='filt',
                        action='store',
                        required=False,
                        default='best',
                        help='The grism/direct filter matches to use if not G102-F105 and G141-F140. Of the form "direct1:grism1,direct2:grism2".')
    parser.add_argument('-i' '--interactive_segmap',
                        dest='interactive',
                        action='store',
                        required=False,
                        default=False,
                        help="Whether or not to interactively check the segmap ids.")
    args = parser.parse_args()

    return args


def preprocess(path, direct, grism):
    """ Preprocesses data for a given set of "raw" grism
    data.

    Parameters
    ----------
    path: str
        Path to the "raw" data files.
    direct : str
        The direct filter to check.
    grism : str
        The grism filter to check.
    """

    # Collect files and info
    files = glob.glob(path)
    info = grizli.utils.get_flt_info(files)
    
    # Match visits
    color_filter = (info['FILTER'] == direct) | (info['FILTER'] == grism)
     
    grism_visits, filters = grizli.utils.parse_flt_files(info=info[color_filter], uniquename=True, get_footprint=True)
    safe_visits = remove_conflicting_visits(grism_visits, direct, grism)
    grism_matches = grizli.utils.parse_grism_associations(safe_visits)

    for pair in grism_matches:
        direct_visit = pair['direct']
        grism_visit = pair['grism']
        
        # Run the preprocessing
        if (len(direct_visit['files']) > 3) and (len(grism_visit['files']) > 3):
            process_direct_grism_visit(direct=direct_visit, grism=grism_visit, align_mag_limits=[14,23])
    
    print('Files for: ' + grism + '/' + direct + ' have been processed.')


def remove_conflicting_visits(visits, direct_filt, grism_filt):
    """ Remove anything that might trip up the visit matching.

    Parameters
    ----------
    visits: list of OrderedDicts 
        The visits as created by grizli.
    direct_filt: str
        The direct filter.
    grism_filt : str
        The grism filter.

    Returns
    -------
    safe_visits: list of OrderedDicts
        A replacement list with any confusing visits
        removed.
    """
    direct_n = direct_filt[1:-1]
    grism_n = grism_filt[1:]
    safe_visits = []

    for visit in visits:
        product = visit['product']
        if direct_n in product:
            test = product.split(direct_n)[0]
        else:
            test = product.split(grism_n)[0]
        if (direct_n in test) or (grism_n in test):
            print(product + ' has been removed. It will confuse the poor bear.')
        else:
            safe_visits.append(visit)

    return safe_visits



## -- RUN
if __name__ == "__main__":

#    args = sys.argv
#    dq_reset = True
#    path = sys.argv[1]
#    if args > 2:
#        dq_reset = eval(sys.argv[2])
#    path, data_name, spectrum, dq_reset, filt = args[1:] 
    path ='/grp/hst/wfc3t/jfowler/bright_objects/gd-153-analysis/RAW/'
    data_name = 'gd-153'

    args = parse_args()
    main(args.path, args.data_name, spectrum=args.spectrum, dq_reset=args.dq_reset, interactive=args.interactive, filt=args.filt)

