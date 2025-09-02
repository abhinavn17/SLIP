#!/usr/bin/env python3
"""
SLIP - Spectral Line Imaging Pipeline
Enhanced with iterative tclean-SoFiA masking

Usage:
    python slip.py config.ini
    python slip.py --make-config data.ms
"""

import os
import sys
import shutil
import logging
import argparse
import subprocess
from pathlib import Path
from configparser import ConfigParser
from typing import Dict, Any, Optional

# CASA tasks and tools
from casatasks import tclean, uvcontsub, imstat, exportfits, uvsub, importfits, makemask
from casatools import table

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class SLIPPipeline:
    """Spectral Line Imaging Pipeline with SoFiA integration."""
    
    def __init__(self, config_file: str):
        """
        Initialize pipeline with configuration file.
        
        Parameters:
        -----------
        config_file : str
            Path to configuration file
        """
        self.config_file = config_file
        self.config = self._read_config()
        self.setup_paths()
        
    def _read_config(self) -> Dict[str, Any]:
        """Read and parse configuration file."""
        if not os.path.exists(self.config_file):
            raise FileNotFoundError(f"Configuration file not found: {self.config_file}")
            
        config_parser = ConfigParser()
        config_parser.read(self.config_file)
        
        # Convert to flat dictionary
        parameters = {}
        for section_name in config_parser.sections():
            section = config_parser[section_name]
            for key, value in section.items():
                parameters[key] = value
                
        # Also check DEFAULT section
        for key, value in config_parser.defaults().items():
            if key not in parameters:
                parameters[key] = value
        
        logger.info(f"Read configuration from: {self.config_file}")
        return parameters
    
    def setup_paths(self):
        """Setup output paths and directories."""
        self.name = self.config['name']
        self.msfile = self.config['msfile']
        self.output_dir = Path(self.config.get('output_dir', './'))
        self.output_dir.mkdir(exist_ok=True)
        
        # Setup working directories
        self.sofia_dir = self.output_dir / f"{self.name}_sofia"
        self.masks_dir = self.output_dir / "masks"
        self.masks_dir.mkdir(exist_ok=True)
        
        logger.info(f"Output directory: {self.output_dir}")
        logger.info(f"Working with MS: {self.msfile}")
    
    def _str_to_bool(self, value: str) -> bool:
        """Convert string to boolean."""
        return value.lower() in ('true', '1', 'yes', 'on')
    
    def _get_config_value(self, key: str, default=None, value_type=str):
        """Get configuration value with type conversion."""
        value = self.config.get(key, default)
        if value is None:
            return None
        
        if value_type == bool:
            return self._str_to_bool(value)
        elif value_type == int:
            return int(value)
        elif value_type == float:
            return float(value)
        else:
            return str(value)
    
    def perform_uv_subtraction(self) -> str:
        """
        Perform UV plane continuum subtraction if requested.
        
        Returns:
        --------
        str : Path to processed MS file
        """
        uvsub_flag = self._get_config_value('do_uvsub', 'false', bool)
        uvcontsub_flag = self._get_config_value('do_uvcontsub', 'false', bool)
        
        processed_ms = self.msfile
        datacolumn = 'data'
        
        # Perform UVSUB
        if uvsub_flag:
            logger.info("Checking for UVSUB...")
            
            tb = table()
            tb.open(self.msfile, nomodify=False)
            
            try:
                uvsub_done = tb.getkeywords().get("UVSUB_DONE", False)
                
                if uvsub_done:
                    logger.info("UVSUB already done. Skipping.")
                else:
                    logger.info("Running UVSUB...")
                    uvsub(vis=self.msfile)
                    tb.putkeyword('UVSUB_DONE', True)
                    logger.info("UVSUB completed.")
                
                datacolumn = 'corrected_data'
                
            finally:
                tb.close()
        
        # Perform UVCONTSUB
        if uvcontsub_flag:
            logger.info("Performing UVCONTSUB...")
            
            processed_ms = self.msfile.replace('.ms', '_uvcontsub.ms')
            
            if os.path.exists(processed_ms):
                logger.info(f"UVCONTSUB file exists: {processed_ms}")
            else:
                # Setup continuum channels
                spw = self._get_config_value('spw', '0')
                bchan = self._get_config_value('bchan')
                echan = self._get_config_value('echan') 
                bchan_line = self._get_config_value('bchan_line')
                echan_line = self._get_config_value('echan_line')
                fitorder = self._get_config_value('fit_order', 1, int)
                
                contchans = f"{spw}:{bchan}~{bchan_line};{echan_line}~{echan}"
                
                logger.info(f"Continuum channels: {contchans}")
                
                uvcontsub(
                    vis=self.msfile,
                    outputvis=processed_ms,
                    spw=spw,
                    datacolumn=datacolumn,
                    fitspec=contchans,
                    fitorder=fitorder,
                    field='0'
                )
                
                logger.info(f"UVCONTSUB completed: {processed_ms}")
            
            datacolumn = 'data'
        
        return processed_ms, datacolumn
    
    def run_initial_tclean(self, ms_file: str, datacolumn: str) -> str:
        """
        Run initial shallow tclean for SoFiA detection.
        
        Parameters:
        -----------
        ms_file : str
            Path to MS file
        datacolumn : str
            Data column to use
            
        Returns:
        --------
        str : Path to initial cleaned image
        """
        logger.info("=== Running Initial TCLEAN ===")
        
        # Get imaging parameters
        cell = self._get_config_value('cell', '4')
        imsize = self._get_config_value('imsize', '512', int)
        weighting = self._get_config_value('weighting', 'briggs')
        robust = self._get_config_value('robust', '0.5', float)
        niter_initial = self._get_config_value('niter_initial', '10000', int)
        threshold_initial = self._get_config_value('threshold_initial', '3sigma')
        
        # Spectral setup
        spw = self._get_config_value('spw', '0')
        bchan = self._get_config_value('bchan')
        echan = self._get_config_value('echan')
        restfreq = self._get_config_value('restfreq', '1.420405752GHz')
        width = self._get_config_value('width', '1')
        
        # UV constraints
        uvrange = self._get_config_value('uvrange', '0')
        uvtaper = self._get_config_value('uvtaper', '0')
        
        # Setup parameters
        spw_cube = f"{spw}:{bchan}~{echan}"
        tcell = f"{cell}arcsec"
        
        tuvrange = f"0~{uvrange}klambda" if int(uvrange) > 0 else ''
        tuvtaper = f"{uvtaper}klambda" if int(uvtaper) > 0 else ''
        
        # Output name
        initial_image = str(self.output_dir / f"{self.name}_initial")
        
        # Clean up existing files
        for ext in ['.image', '.residual', '.model', '.psf', '.pb', '.sumwt']:
            if os.path.exists(f"{initial_image}{ext}"):
                shutil.rmtree(f"{initial_image}{ext}")
        
        # Run tclean
        logger.info(f"Initial imaging: {initial_image}")
        
        tclean(
            vis=ms_file,
            imagename=initial_image,
            field=self._get_config_value('field', '0'),
            spw=spw_cube,
            specmode='cube',
            outframe='BARY',
            veltype='optical',
            restfreq=restfreq,
            gridder='standard',
            mosweight=False,
            cfcache='',
            computepastep=360.0,
            rotatepastep=360.0,
            pblimit=0.2,
            normtype='flatnoise',
            deconvolver='hogbom',
            scales=[],
            nterms=1,
            smallscalebias=0.6,
            imsize=imsize,
            cell=tcell,
            phasecenter='',
            stokes='I',
            projection='SIN',
            startmodel='',
            weighting=weighting,
            robust=robust,
            uvrange=tuvrange,
            uvtaper=tuvtaper if tuvtaper else [],
            niter=niter_initial,
            gain=0.1,
            threshold=threshold_initial,
            nsigma=0.0,
            cycleniter=-1,
            cyclefactor=1.0,
            minpsffraction=0.05,
            maxpsffraction=0.8,
            interactive=False,
            usemask='auto-multithresh',
            sidelobethreshold=2.0,
            noisethreshold=4.25,
            lownoisethreshold=1.5,
            negativethreshold=15.0,
            smoothfactor=1.0,
            minbeamfrac=0.3,
            cutthreshold=0.01,
            growiterations=75,
            dogrowprune=True,
            minpercentchange=-1.0,
            verbose=False,
            fastnoise=True,
            restart=True,
            savemodel='none',
            calcres=True,
            calcpsf=True,
            parallel=False
        )
        
        # Get statistics
        noise_stats = imstat(imagename=f"{initial_image}.image", chans='0~19')
        rms = noise_stats['rms'][0]
        peak_stats = imstat(imagename=f"{initial_image}.image")
        peak = peak_stats['max'][0]
        snr = peak / rms
        
        logger.info(f"Initial clean statistics:")
        logger.info(f"  RMS: {rms*1e3:.1f} mJy")
        logger.info(f"  Peak: {peak*1e3:.1f} mJy") 
        logger.info(f"  Max SNR: {snr:.1f}")
        
        return initial_image
    
    def run_sofia_detection(self, image_cube: str) -> str:
        """
        Run SoFiA source detection on image cube.
        
        Parameters:
        -----------
        image_cube : str
            Path to CASA image cube
            
        Returns:
        --------
        str : Path to SoFiA output directory
        """
        logger.info("=== Running SoFiA Source Detection ===")
        
        # Export to FITS
        fits_cube = str(self.output_dir / f"{self.name}_for_sofia.fits")
        
        exportfits(
            imagename=f"{image_cube}.image",
            fitsimage=fits_cube,
            velocity=True,
            optical=True,
            bitpix=-32,
            minpix=0,
            maxpix=-1,
            overwrite=True,
            dropstokes=True,
            stokeslast=True,
            history=True,
            dropdeg=False
        )
        
        logger.info(f"Exported FITS cube: {fits_cube}")
        
        # Create SoFiA parameter file
        sofia_par = self._create_sofia_parameter_file(fits_cube)
        
        # Run SoFiA
        logger.info("Running SoFiA...")
        
        # Clean up previous SoFiA output
        if self.sofia_dir.exists():
            shutil.rmtree(self.sofia_dir)
        
        # Run SoFiA
        result = subprocess.run(
            ['sofia', sofia_par], 
            capture_output=True, 
            text=True,
            cwd=str(self.output_dir)
        )
        
        if result.returncode != 0:
            logger.error(f"SoFiA failed: {result.stderr}")
            raise RuntimeError("SoFiA source detection failed")
        
        logger.info("SoFiA completed successfully")
        
        return str(self.sofia_dir)
    
    def _create_sofia_parameter_file(self, fits_cube: str) -> str:
        """Create SoFiA parameter file."""
        
        sofia_par = str(self.output_dir / "sofia_config.par")
        
        # Get SoFiA parameters from config
        threshold = self._get_config_value('sofia_threshold', '3.5', float)
        kernels = self._get_config_value('sofia_kernels', 
            "[[0,0,0,'b'],[0,0,3,'b'],[0,0,7,'b'],[3,3,0,'b'],[3,3,3,'b']]")
        
        config_content = f"""
# SoFiA Configuration File
# Generated by SLIP Pipeline

# Input
input.data = {fits_cube}
input.region = 
input.gain = 
input.noise = 

# Flagging
flag.regions = 
flag.auto = false
flag.threshold = 5.0
flag.log = false

# Continuum subtraction  
contsub.enable = false

# Noise measurement
reliability.enable = true
reliability.threshold = 0.9
reliability.scaleKernel = 0.4

# Source finding
scfind.enable = true
scfind.kernelsFile = 
scfind.kernels = {kernels}
scfind.threshold = {threshold}
scfind.replacement = 2.0
scfind.statistic = std
scfind.fluxRange = positive

# Merging
merge.enable = true
merge.radiusX = 3
merge.radiusY = 3  
merge.radiusZ = 3
merge.minSizeX = 5
merge.minSizeY = 5
merge.minSizeZ = 5

# Reliability
reliability.enable = true
reliability.threshold = 0.9
reliability.scaleKernel = 0.4

# Parameterisation
parameter.enable = true
parameter.wcs = true
parameter.physical = false
parameter.prefix = SoFiA
parameter.offset = false

# Output
output.directory = {self.sofia_dir}/
output.filename = {self.name}
output.writeCatASCII = true
output.writeCatXML = false
output.writeCatSQL = false
output.writeNoise = false
output.writeFiltered = false  
output.writeMask = true
output.writeMask2d = false
output.writeRawMask = false
output.writeMoments = true
output.writeCubelets = false
output.marginCubelets = 5
output.thresholdMom12 = 0.0
output.overwrite = true
"""
        
        with open(sofia_par, 'w') as f:
            f.write(config_content)
        
        logger.info(f"Created SoFiA parameter file: {sofia_par}")
        return sofia_par
    
    def convert_sofia_mask_to_casa(self, iteration: int = 0) -> Optional[str]:
        """
        Convert SoFiA mask to CASA format.
        
        Parameters:
        -----------
        iteration : int
            Iteration number for filename
            
        Returns:
        --------
        str : Path to CASA mask or None if no mask found
        """
        logger.info("=== Converting SoFiA mask to CASA format ===")
        
        # Find SoFiA mask file
        mask_patterns = [
            f"{self.name}_mask.fits",
            f"{self.name}*mask*.fits",
            "*mask.fits"
        ]
        
        sofia_mask = None
        for pattern in mask_patterns:
            mask_files = list(self.sofia_dir.glob(pattern))
            if mask_files:
                sofia_mask = mask_files[0]
                break
        
        if sofia_mask is None:
            logger.warning("No SoFiA mask found")
            return None
        
        logger.info(f"Found SoFiA mask: {sofia_mask}")
        
        # Convert to CASA format
        casa_mask_name = f"{self.name}_sofia_mask_iter{iteration:02d}"
        casa_mask_path = str(self.masks_dir / f"{casa_mask_name}.mask")
        
        # Import FITS mask to CASA
        temp_mask_image = str(self.masks_dir / f"temp_mask_iter{iteration:02d}.image")
        
        importfits(
            fitsimage=str(sofia_mask),
            imagename=temp_mask_image,
            overwrite=True
        )
        
        # Get reference image for coordinate system
        initial_image = str(self.output_dir / f"{self.name}_initial.image")
        
        # Create proper CASA mask
        makemask(
            mode='copy',
            inpimage=initial_image,
            inpmask=temp_mask_image,
            output=casa_mask_path,
            overwrite=True
        )
        
        # Clean up temporary file
        if os.path.exists(temp_mask_image):
            shutil.rmtree(temp_mask_image)
        
        logger.info(f"Created CASA mask: {casa_mask_path}")
        return casa_mask_path
    
    def run_deep_tclean_with_mask(self, ms_file: str, datacolumn: str, 
                                 mask_file: str, iteration: int = 1) -> str:
        """
        Run deep tclean with SoFiA mask.
        
        Parameters:
        -----------
        ms_file : str
            Path to MS file
        datacolumn : str
            Data column to use
        mask_file : str
            Path to CASA mask
        iteration : int
            Iteration number
            
        Returns:
        --------
        str : Path to deep cleaned image
        """
        logger.info(f"=== Running Deep TCLEAN (Iteration {iteration}) ===")
        
        # Get deep cleaning parameters
        niter_deep = self._get_config_value('niter_deep', '100000', int)
        threshold_deep = self._get_config_value('threshold_deep', '0.5sigma')
        
        # Setup other parameters (reuse from initial clean)
        cell = self._get_config_value('cell', '4')
        imsize = self._get_config_value('imsize', '512', int)
        weighting = self._get_config_value('weighting', 'briggs')
        robust = self._get_config_value('robust', '0.5', float)
        
        spw = self._get_config_value('spw', '0')
        bchan = self._get_config_value('bchan')
        echan = self._get_config_value('echan')
        restfreq = self._get_config_value('restfreq', '1.420405752GHz')
        
        uvrange = self._get_config_value('uvrange', '0')
        uvtaper = self._get_config_value('uvtaper', '0')
        
        # Setup parameters
        spw_cube = f"{spw}:{bchan}~{echan}"
        tcell = f"{cell}arcsec"
        tuvrange = f"0~{uvrange}klambda" if int(uvrange) > 0 else ''
        tuvtaper = f"{uvtaper}klambda" if int(uvtaper) > 0 else ''
        
        # Output name
        deep_image = str(self.output_dir / f"{self.name}_deep_iter{iteration:02d}")
        
        # Clean up existing files
        for ext in ['.image', '.residual', '.model', '.psf', '.pb', '.sumwt']:
            if os.path.exists(f"{deep_image}{ext}"):
                shutil.rmtree(f"{deep_image}{ext}")
        
        logger.info(f"Deep cleaning with mask: {mask_file}")

        if parallel:

            nproc = par['nproc']

            os.environ['OMP_NUM_THREADS'] = '1'

            command = ['mpirun', '-nq', f'{nproc}', 'python', '-m', 'slip.imager', ms_uvsub, out_tclean, field, datacolumn, spw_cube, outframe, veltype, restfreq, tcell, tuvrange, tuvtaper, imsize, weighting, niter, cycleniter, nsigma, robust, width]
        
        else:

            command = ['python', '-m', 'slip.imager', ms_uvsub, out_tclean, field, datacolumn, spw_cube, outframe, veltype, restfreq, tcell, tuvrange, tuvtaper, imsize, weighting, niter, cycleniter, nsigma, robust, width]
        # Run tclean with mask
        tclean(
            vis=ms_file,
            imagename=deep_image,
            field=self._get_config_value('field', '0'),
            spw=spw_cube,
            specmode='cube',
            outframe='BARY',
            veltype='optical',
            restfreq=restfreq,
            gridder='standard',
            deconvolver='hogbom',
            imsize=imsize,
            cell=tcell,
            weighting=weighting,
            robust=robust,
            uvrange=tuvrange,
            uvtaper=tuvtaper if tuvtaper else [],
            niter=niter_deep,
            gain=0.1,
            threshold=threshold_deep,
            cycleniter=-1,
            cyclefactor=1.0,
            interactive=False,
            usemask='user',
            mask=mask_file,
            restart=True,
            savemodel='none',
            calcres=True,
            calcpsf=True,
            parallel=False
        )
        
        # Get statistics
        noise_stats = imstat(imagename=f"{deep_image}.image", chans='19~39')
        rms = noise_stats['rms'][0]
        peak_stats = imstat(imagename=f"{deep_image}.image")
        peak = peak_stats['max'][0]
        snr = peak / rms
        
        logger.info(f"Deep clean iteration {iteration} statistics:")
        logger.info(f"  RMS: {rms*1e3:.1f} mJy")
        logger.info(f"  Peak: {peak*1e3:.1f} mJy")
        logger.info(f"  Max SNR: {snr:.1f}")
        
        return deep_image
    
    def run_iterative_clean_sofia_loop(self) -> str:
        """
        Run the complete iterative tclean-SoFiA loop.
        
        Returns:
        --------
        str : Path to final cleaned image
        """
        logger.info("=== Starting Iterative TCLEAN-SoFiA Loop ===")
        
        # Step 1: UV plane processing
        processed_ms, datacolumn = self.perform_uv_subtraction()
        
        # Step 2: Initial shallow clean
        initial_image = self.run_initial_tclean(processed_ms, datacolumn)
        
        # Step 3: Run SoFiA on initial clean
        sofia_dir = self.run_sofia_detection(initial_image)
        
        # Step 4: Convert SoFiA mask and run deep clean
        max_iterations = self._get_config_value('max_iterations', '2', int)
        
        current_image = initial_image
        
        for iteration in range(1, max_iterations + 1):
            logger.info(f"=== Iteration {iteration}/{max_iterations} ===")
            
            # Convert SoFiA mask to CASA format
            casa_mask = self.convert_sofia_mask_to_casa(iteration)
            
            if casa_mask is None:
                logger.warning(f"No mask available for iteration {iteration}")
                break
            
            # Run deep clean with mask
            current_image = self.run_deep_tclean_with_mask(
                processed_ms, datacolumn, casa_mask, iteration
            )
            
            # Optionally run SoFiA again on deeper image for next iteration
            if iteration < max_iterations:
                logger.info("Running SoFiA on deeper image for next iteration...")
                sofia_dir = self.run_sofia_detection(current_image)
        
        # Step 5: Export final products
        self.export_final_products(current_image)
        
        logger.info("=== Pipeline Complete ===")
        return current_image
    
    def export_final_products(self, final_image: str):
        """Export final products."""
        logger.info("=== Exporting Final Products ===")
        
        # Export final FITS cube
        final_fits = str(self.output_dir / f"{self.name}_final.fits")
        
        exportfits(
            imagename=f"{final_image}.image",
            fitsimage=final_fits,
            velocity=True,
            optical=True,
            bitpix=-32,
            minpix=0,
            maxpix=-1,
            overwrite=True,
            dropstokes=True,
            stokeslast=True,
            history=True,
            dropdeg=False
        )
        
        logger.info(f"Final FITS cube: {final_fits}")
        
        # Run final SoFiA analysis
        logger.info("Running final SoFiA analysis...")
        final_sofia_dir = self.run_sofia_detection(final_image.replace('_deep_iter', '_final'))
        
        logger.info(f"Final SoFiA products in: {final_sofia_dir}")
    
    def run_pipeline(self) -> str:
        """Run the complete pipeline."""
        try:
            return self.run_iterative_clean_sofia_loop()
        except Exception as e:
            logger.error(f"Pipeline failed: {e}")
            raise


def make_config_template():
    """Create a template configuration file."""
    parser = argparse.ArgumentParser(description='Create SLIP configuration template')
    parser.add_argument('--config', '-c', default='slip.ini', 
                       help='Configuration file name')
    parser.add_argument('msfile', type=str, help='Measurement set file')
    args = parser.parse_args()
    
    msfile = os.path.abspath(args.msfile)
    name = os.path.basename(msfile).replace('.ms', '')
    
    template_content = f"""
[DEFAULT]
# Basic setup
name = {name}
msfile = {msfile}
output_dir = ./

# Data processing
do_uvsub = false
do_uvcontsub = false
spw = 0
field = 0

# Continuum subtraction (if uvcontsub = true)
bchan = 0
echan = 50
bchan_line = 200
echan_line = 800
fit_order = 1

# Imaging parameters
cell = 4
imsize = 512
weighting = briggs
robust = 0.5
uvrange = 0
uvtaper = 0
restfreq = 1.420405752GHz
width = 1

# Initial clean (shallow)
niter_initial = 10000
threshold_initial = 3sigma

# Deep clean (with mask)
niter_deep = 100000
threshold_deep = 0.5sigma

# Iterative loop
max_iterations = 2

# SoFiA parameters
sofia_threshold = 3.5
sofia_kernels = [[0,0,0,'b'],[0,0,3,'b'],[0,0,7,'b'],[3,3,0,'b'],[3,3,3,'b']]

# Parallel processing
nproc = 4
"""
    
    with open(args.config, 'w') as f:
        f.write(template_content)
    
    print(f"Created configuration template: {args.config}")
    print(f"Edit the file and run: python slip.py {args.config}")


def main():
    """Main function."""
    if len(sys.argv) > 1 and sys.argv[1] == '--make-config':
        make_config_template()
        return
    
    parser = argparse.ArgumentParser(description='Spectral Line Imaging Pipeline')
    parser.add_argument('config', help='Configuration file')
    parser.add_argument('--verbose', '-v', action='store_true', 
                       help='Enable verbose logging')
    
    args = parser.parse_args()
    
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    # Run pipeline
    pipeline = SLIPPipeline(args.config)
    final_image = pipeline.run_pipeline()
    
    print(f"\n{'='*50}")
    print("SLIP Pipeline completed successfully!")
    print(f"Final image: {final_image}")
    print(f"Output directory: {pipeline.output_dir}")
    print(f"{'='*50}\n")


if __name__ == '__main__':
    main()
