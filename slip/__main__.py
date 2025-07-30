'''
To image spectral cube.
'''
import os
from casatasks import tclean, uvcontsub, imstat, exportfits, uvsub
from casatools import table
import argparse
import subprocess
from configparser import ConfigParser

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

    os.system('cp ' + template + ' ' + config)
    os.system('cp ' + sofia_template + ' .')

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

def main():

    # Parse the command line arguments

    parser = argparse.ArgumentParser(description='Spectral Line Imaging Pipeline')
    parser.add_argument('infile', type=str, help='Input parameter file')
    args = parser.parse_args()

    infile = args.infile

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
        datacolumn = 'data'

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
    maskcx = par['maskcx']
    maskcy = par['maskcx']
    maskrad = par['maskcx']
    width = par['width']
    spw_cube = str(spw) + ':' + str(bchan) + '~' + str(echan)
    tcell = str(cell) + 'arcsec'
    if int(uvrange) > 0:
        tuvrange = '0~' + str(uvrange) + 'klambda'
    else:
        tuvrange = ''
    if int(uvtaper) > 0:
        tuvtaper = str(uvtaper) + 'klambda'
    else:
        tuvtaper = ''
    # timsize = '[' + str(imsize) + ',' + str(imsize) + ']'
    tmask = 'circle[' + str(maskcx) + 'pix,' + str(maskcy) + 'pix],' + str(maskrad) + 'pix]'
    twidth = str(width) + 'km/s'
    outframe = 'BARY'
    veltype = 'optical'
    #tweighting = 'natural'
    #robust = -2.0 # -2 uniform, +2 natural


    # We will image it with nsigma criteria
    # Imaging with 5k lambda

    out_tclean = name + '_cube'
    cmd = 'rm -rf ' + out_tclean + '*'
    os.system(cmd)

    try:
        from casampi import MPIEnvironment, MPICommandClient
        parallel = True
    except ImportError:
        parallel = False

    if parallel:

        nproc = par['nproc']

        os.environ['OMP_NUM_THREADS'] = '1'

        command = ['mpirun', '-nq', f'{nproc}', 'python', '-m', 'slip.imager', ms_uvsub, out_tclean, field, datacolumn, spw_cube, outframe, veltype, restfreq, tcell, tuvrange, tuvtaper, imsize, weighting, niter, cycleniter, nsigma, robust, width]
    
    else:

        command = ['python', '-m', 'slip.imager', ms_uvsub, out_tclean, field, datacolumn, spw_cube, outframe, veltype, restfreq, tcell, tuvrange, tuvtaper, imsize, weighting, niter, cycleniter, nsigma, robust, width]

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
    print('RMS = %4.1f mJy'%(rms*1e3))
    print('Max SNR is %4.1f'%snr)
    print('+++++++++++++++++++++++++')
    print('+++++++++++++++++++++++++')
    print('+++++++++++++++++++++++++')

    # Exporting as a fits file for SoFiA
    exportfits(imagename=out_tclean+'.image',fitsimage='test_image.fits', velocity=True, optical=True, bitpix=-32, minpix=0,maxpix=-1, overwrite=True, dropstokes=True, stokeslast=True, history=True, dropdeg=False)

    os.environ['OMP_NUM_THREADS'] = f'{nproc}'

    os.system('rm -rf test_sofia/')
    os.system('mkdir test_sofia')
    os.system('sofia slip_sofia_template.par')

    out_fits = out_tclean + '.fits'
    cmd = 'mv test_image.fits ' + out_fits
    os.system(cmd)

    out_sofia_dir = name + '_sofia'
    cmd = 'rm -rf ' + out_sofia_dir
    os.system(cmd)
    cmd = 'mv test_sofia ' + out_sofia_dir
    os.system(cmd)
    # ------------------------------

if __name__ == '__main__':

    main()