'''
To image spectral cube.
'''
import os
import sys
from casatasks import uvcontsub, imstat, exportfits, uvsub, importfits
from casatools import table
import argparse
import subprocess
from configparser import ConfigParser
from .sofia_mask import convert_sofia_mask_to_casa, find_sofia_files


def make_config_file():

    parser = argparse.ArgumentParser(description='Make a configuration file for SLIP')
    parser.add_argument('--config', '-c', default='slip.ini', help='Configuration file name')
    parser.add_argument('msfile', type=str, help='Measurement set file')
    args = parser.parse_args()

    config = args.config
    msfile = args.msfile

    config = os.path.join(os.getcwd(), config)

    module_dir = os.path.dirname(os.path.abspath(__file__))
    template = os.path.join(module_dir, 'slip_parameters.ini')
    sofia_template = os.path.join(module_dir, 'slip_sofia_template.par')

    sofia_config = os.path.join(os.getcwd(), config.split('.ini')[0] + '_sofia.par')

    os.system('cp ' + template + ' ' + config)
    os.system('cp ' + sofia_template + ' ' + sofia_config)

    config_parser = ConfigParser()
    config_parser.read(config)

    msfile = os.path.abspath(msfile)
    name = msfile.split('.ms')[0]
    name = os.path.basename(name)

    config_parser.set('DEFAULT', 'name', name)

    config_parser.set('DEFAULT', 'msfile', msfile)

    with open(config, 'w') as configfile:
        config_parser.write(configfile)

def read_parameters(file_path):
    
    config_parser = ConfigParser()
    config_parser.read(file_path)
    
    parameters = {}

    for section in config_parser.sections():
            
        for key in config_parser[section]:

            value = config_parser[section][key]
                
            # Store the key-value pair in the dictionary
            parameters[key] = value
    
    return parameters
# ------------------------------

def run_tclean(out_tclean, ms_uvsub, field, datacolumn, spw_cube, outframe, veltype, restfreq, tcell, tuvrange, tuvtaper, imsize, weighting, niter, cycleniter, nsigma, robust, width, deconvolver, usemask, mask, nproc, nmajor):

    cmd = 'rm -rf ' + out_tclean + '*'
    os.system(cmd)

    try:
        from casampi import MPIEnvironment, MPICommandClient
        parallel = True
    except ImportError:
        parallel = False

    print(f"Running tclean with parameters:")
    print(f"  MS UVSUB: {ms_uvsub}")
    print(f"  Output Tclean: {out_tclean}")
    print(f"  Field: {field}")
    print(f"  Data Column: {datacolumn}")
    print(f"  SPW Cube: {spw_cube}")
    print(f"  Outframe: {outframe}")
    print(f"  Veltype: {veltype}")
    print(f"  Restfreq: {restfreq}")
    print(f"  Cell: {tcell}")
    print(f"  Uvrange: {tuvrange}")
    print(f"  Uvtaper: {tuvtaper}")
    print(f"  Imsize: {imsize}")
    print(f"  Weighting: {weighting}")
    print(f"  Niter: {niter}")
    print(f"  Cycleniter: {cycleniter}")
    print(f"  Nsigma: {nsigma}")
    print(f"  Robust: {robust}")
    print(f"  Width: {width}")
    print(f"  Deconvolver: {deconvolver}")
    print(f"  Usemask: {usemask}")
    print(f"  Mask: {mask}")
    print(f"  Nmajor: {nmajor}")
    print(f"  Nproc: {nproc}")

    if parallel:

        os.environ['OMP_NUM_THREADS'] = '1'

        command = ['mpirun', '-nq', f'{nproc}', 'python', '-m', 'slip.imager', ms_uvsub, out_tclean, field, datacolumn, spw_cube, outframe, veltype, restfreq, tcell, tuvrange, tuvtaper, imsize, weighting, niter, cycleniter, nsigma, robust, width, deconvolver, usemask, mask, str(nmajor)]
    
    else:

        command = ['python', '-m', 'slip.imager', ms_uvsub, out_tclean, field, datacolumn, spw_cube, outframe, veltype, restfreq, tcell, tuvrange, tuvtaper, imsize, weighting, niter, cycleniter, nsigma, robust, width, deconvolver, usemask, mask, str(nmajor)]

    subprocess.run(command)

    # Estimating max SNR 
    tname = out_tclean + '.image'
    noise_stat = imstat(imagename=tname, chans='20')
    rms = noise_stat['rms'][0]
    # noise_stats = imstat(imagename=tname)
    peak = noise_stat['max'][0]
    snr = peak/rms
    print('+++++++++++++++++++++++++')
    print('+++++++++++++++++++++++++')
    print('+++++++++++++++++++++++++')
    print(f'For cube {out_tclean}')
    print('RMS = %4.1f mJy'%(rms*1e3))
    print('Max SNR is %4.1f'%snr)
    print('+++++++++++++++++++++++++')
    print('+++++++++++++++++++++++++')
    print('+++++++++++++++++++++++++')

    # Exporting as a fits file for SoFiA
    exportfits(imagename=out_tclean+'.image',fitsimage = out_tclean + '_image.fits', velocity=False, bitpix=-32, minpix=0,maxpix=-1, overwrite=True, dropstokes=False, stokeslast=True, history=True, dropdeg=False)

    return out_tclean + '_image.fits'

def run_sofia(sofia_infile, sofia_thresh, name, i, nproc, central_mask):

    fitsimage = name + '_cube_' + str(i) + '_image.fits'
    filename = name + '_cube_' + str(i)

    os.environ['OMP_NUM_THREADS'] = f'{nproc}'

    with open(sofia_infile, "r") as f:
        par_lines = f.readlines()

    new_lines = [ f"input.data = {fitsimage}\n" if line.startswith("input.data") else line for line in par_lines]

    new_lines = [f"scfind.threshold = {sofia_thresh}\n" if line.startswith("scfind.threshold") else line for line in new_lines]

    outdir = fitsimage.split('_image.fits')[0] + '_sofia_' + str(i)

    if os.path.exists(outdir):

        os.system('rm -rf ' + outdir)
    else:
        os.system('mkdir ' + outdir)

    new_lines = [f"output.directory = {outdir}\n" if line.startswith("output.directory") else line for line in new_lines]
    new_lines = [f"output.filename = {filename}\n" if line.startswith("output.filename") else line for line in new_lines]

    # Write to a temporary par file
    with open(sofia_infile, "w") as temp_f:
        temp_f.writelines(new_lines)
    # Run SOFIA (replace 'sofia' with the actual command if needed)
    subprocess.run(["sofia", sofia_infile])

    sofia_files = find_sofia_files(outdir)

    mask_file = sofia_files['mask']
    mom0_file = sofia_files['mom0']

    binary_mask = name + '_cube_sofia_' + str(i+1) + '_mask.fits'

    # Convert mask
    result = convert_sofia_mask_to_casa(
        mask_file=mask_file,
        mom0_file=mom0_file,
        output_file=binary_mask,
        central_mask=central_mask)

    casa_mask = name + '_cube_sofia_' + str(i+1) + '_mask.image'

    importfits(fitsimage=binary_mask, imagename=casa_mask, overwrite=True, defaultaxes= True, defaultaxesvalues = ['', '', '', 'I'])

    return casa_mask


def main():

    # Parse the command line arguments

    parser = argparse.ArgumentParser(description='Spectral Line Imaging Pipeline')
    parser.add_argument('infile', type=str, help='Input parameter file')
    args = parser.parse_args()

    infile = args.infile
    sofia_infile = infile.split('.ini')[0] + '_sofia.par'

    # Read the parameters from the file
    par = read_parameters(infile)
    name = par['name']
    bchan = par['bchan']
    echan = par['echan']
    bchan_line = par['bchan_line']
    echan_line = par['echan_line']
    spw = par['spw']
    fitorder = int(par['fit_order'])
    uvsub_flag = par['do_uvsub']
    uvcontsub_flag = par['do_uvcontsub']

    if uvsub_flag.lower() == 'true':
        uvsub_flag = True
    else:
        uvsub_flag = False

    if uvcontsub_flag.lower() == 'true':
        uvcontsub_flag = True
    else:
        uvcontsub_flag = False

    # Performing uvsub
    msfile = par['msfile']
    datacolumn = 'data'
    ms_uvsub = msfile

    if uvsub_flag:
        t = table(msfile, readonly=False)
        uvsub_done = t.getkeywords().get("UVSUB_DONE", False)

        if uvsub_done:
            print('UVSUB already done. Not doing uvsub again.')
        else:
            print('Running UVSUB')
            uvsub(vis=msfile)
            t.putkeyword('UVSUB_DONE', True)

        t.close()
        datacolumn = 'corrected_data'

    if uvcontsub_flag:
        ms_uvsub = msfile.split('.ms')[0] + '_uvcontsub.ms'
        contchans = str(spw) + ':' + str(bchan) + '~' + str(bchan_line) + ';' + str(echan_line) + '~' + str(echan)

        if not(os.path.exists(ms_uvsub)):
            print("Running UVCONTSUB")
            
            uvcontsub(vis=msfile, outputvis=ms_uvsub, spw = str(spw), datacolumn = datacolumn, fitspec=contchans, fitorder=fitorder, field='0')
            datacolumn = 'data'
        else:
            print('UVCONTSUB file exist. Not doing uvcontsub.')
        # ------------------------------

    # Performing line imaging
    uvrange = par['uvrange']
    uvtaper = par['uvtaper']
    cell = par['cell']
    imsize = par['imsize']
    niter = par['niter']
    cycleniter = par['cycleniter']
    restfreq = par['restfreq']
    weighting = par['weighting']
    robust = par['robust']
    nsigma = par['nsigma']
    field = par['field']
    width = par['width']
    spw_cube = str(spw) + ':' + str(bchan) + '~' + str(echan)
    tcell = str(cell) + 'arcsec'
    deconvolver = par['deconvolver']
    if int(uvrange) > 0:
        tuvrange = '0~' + str(uvrange) + 'klambda'
    else:
        tuvrange = ''
    if int(uvtaper) > 0:
        tuvtaper = str(uvtaper) + 'klambda'
    else:
        tuvtaper = ''

    outframe = par['outframe']
    veltype = par['veltype']
    nproc = par['nproc']

    out_tclean = name + '_cube_0'

    run_tclean(out_tclean = out_tclean, ms_uvsub = ms_uvsub, field = field, datacolumn = datacolumn, spw_cube = spw_cube, outframe = outframe, veltype = veltype, restfreq = restfreq, tcell = tcell, tuvrange = tuvrange, tuvtaper = tuvtaper, imsize = imsize, weighting = weighting, niter = '1000', cycleniter = cycleniter, nsigma = nsigma, robust = robust, width = width, deconvolver = 'hogbom', usemask = 'pb', mask = '', nproc = nproc, nmajor = 0)

    sofia_niter = int(par['sofia_niter'])
    sofia_thresh = par['thresh']
    central_mask = par['centre_mask']

    if central_mask.lower() == 'true':
        central_mask = True
    else:
        central_mask = False

    for i in range(sofia_niter):

        casa_mask = run_sofia(sofia_infile, sofia_thresh, name, i, nproc, central_mask)

        out_tclean = name + '_cube_' + str(i+1)

        run_tclean(out_tclean = out_tclean, ms_uvsub = ms_uvsub, field = field, datacolumn = datacolumn, spw_cube = spw_cube, outframe = outframe, veltype = veltype, restfreq = restfreq, tcell = tcell, tuvrange = tuvrange, tuvtaper = tuvtaper, imsize = imsize, weighting = weighting, niter = niter, cycleniter = cycleniter, nsigma = nsigma, robust = robust, width = width, deconvolver = deconvolver, usemask = 'user', mask = casa_mask, nproc = nproc, nmajor= -1)

    # os.system('rm -rf test_sofia/')
    # os.system('mkdir test_sofia')
    # os.system('sofia slip_sofia_template.par')

    # out_fits = out_tclean + '.fits'
    # cmd = 'mv test_image.fits ' + out_fits
    # os.system(cmd)

    # out_sofia_dir = name + '_sofia'
    # cmd = 'rm -rf ' + out_sofia_dir
    # os.system(cmd)
    # cmd = 'mv test_sofia ' + out_sofia_dir
    # os.system(cmd)
    # # ------------------------------

if __name__ == '__main__':

    main()