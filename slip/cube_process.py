#!/usr/bin/env python3

import sys
import os
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy import units as u
from astropy.constants import c
import argparse
import warnings
from spectral_cube import SpectralCube
from radio_beam import Beam
from astropy.convolution import Gaussian2DKernel

# Default HI rest frequency (Hz)
DEFAULT_RESTFREQ = 1420.40575178e6

# Default pattern for parsing through during batch conversion
DEFAULT_PATTERN = ['*.fits', '*.image']


SPECSYS_MAP = {
    'TOPO': 'TOPOCENT',
    'LSRK': 'LSRK',
    'LSRD': 'LSRD',
    'BARY': 'BARYCENT',
    'GALACTO': 'GALACTOC',
    'HELIOCENT': 'HELIOCEN',
    'GEO': 'GEOCENTR',
    'CMB': 'CMBDIPOL',
    'LGROUP': 'LOCALGRP',
    'REST': 'SOURCE',
    'TOPOCENT': 'TOPOCENT',
    'BARYCENT': 'BARYCENT',
    'GEOCENTR': 'GEOCENTR',
    'HELIOCEN': 'HELIOCEN',
    'GALACTOC': 'GALACTOC',
    'LOCALGRP': 'LOCALGRP',
    'CMBDIPOL': 'CMBDIPOL',
    'SOURCE': 'SOURCE'
}
    

def specsys_correction(specsys):
    """
    Map CASA outframe specsys names to standard FITS-compatible names
    """
    if specsys.upper() in SPECSYS_MAP:
        return SPECSYS_MAP[specsys.upper()]
    else:
        warnings.warn(f"Unknown SPECSYS '{specsys}', using as-is.")
        return specsys.upper()

def parse_freq_string(freq_str):
    freq_str = str(freq_str).lower().strip()
    try:
        if freq_str.endswith(('hz', 'khz', 'mhz', 'ghz')):
            if freq_str.endswith('ghz'):
                freq_value = float(freq_str[:-3]) * 1e9
            elif freq_str.endswith('mhz'):
                freq_value = float(freq_str[:-3]) * 1e6
            elif freq_str.endswith('khz'):
                freq_value = float(freq_str[:-3]) * 1e3
            elif freq_str.endswith('hz') and ('k' not in freq_str and 'm' not in freq_str and 'g' not in freq_str):
                freq_value = float(freq_str[:-2])
        else:
            # If no valid suffix, attempt to directly convert to float
            freq_value = float(freq_str)
        return freq_value
    except ValueError:
        raise ValueError(f"Invalid frequency string: '{freq_str}'")

def detect_spectral_axis_type(cube):
    """
    Detect whether spectral axis is frequency or velocity
    
    Returns: 'frequency', 'velocity', or 'unknown'
    """
    spectral_axis = cube.spectral_axis
    
    # Check the unit of the spectral axis
    if spectral_axis.unit.is_equivalent(u.Hz):
        return 'frequency'
    elif spectral_axis.unit.is_equivalent(u.km/u.s) or spectral_axis.unit.is_equivalent(u.m/u.s):
        return 'velocity'
    else:
        return 'unknown'

def freq_to_velocity(freq, rest_freq, convention='optical'):
    """
    Convert frequency to velocity using different conventions
    
    Parameters:
    freq: frequency array (Hz)
    rest_freq: rest frequency (Hz)
    convention: 'optical', 'radio', or 'relativistic'
    
    Returns:
    velocity in km/s
    """
    c_kms = c.to(u.km/u.s).value  # speed of light in km/s
    rest_freq = parse_freq_string(rest_freq)
    
    if convention == 'optical':
        # Optical convention: v = c * (λ - λ0) / λ0 = c * (f0 - f) / f
        velocity = c_kms * (rest_freq - freq) / freq
    elif convention == 'radio':
        # Radio convention: v = c * (f0 - f) / f0
        velocity = c_kms * (rest_freq - freq) / rest_freq
    elif convention == 'relativistic':
        # Relativistic convention: v = c * ((f0/f)^2 - 1) / ((f0/f)^2 + 1)
        ratio = rest_freq / freq
        velocity = c_kms * (ratio**2 - 1) / (ratio**2 + 1)
    else:
        raise ValueError("Convention must be 'optical', 'radio', or 'relativistic'")
    
    return velocity

def convert_freq_to_vel_cube(input_file,
                              rest_freq=1420.40575178e6,
                              convention='radio',
                              specsys='TOPOCENT',
                              smooth=None):
    """
    Convert frequency cube to velocity cube using spectral-cube
    
    Parameters:
    input_file: Input FITS file with frequency axis
    rest_freq: Rest frequency in Hz (default: HI 21cm line)
    convention: Velocity convention ('optical', 'radio', 'relativistic')
    specsys: Spectral reference frame
    smooth: Tuple of (bmaj, bmin, bpa) in arcsec and degrees for smoothing
    """
    print(f"Converting frequency cube: {input_file}")
    rest_freq = parse_freq_string(rest_freq)
    print(f"Rest frequency: {rest_freq/1e6:.6f} MHz")
    print(f"Velocity convention: {convention}")
    
    # Load cube with spectral-cube
    cube = SpectralCube.read(input_file)
    cube.allow_huge_operations=True
    print(f"Input cube shape: {cube.shape}")
    print(f"Spectral axis: {cube.spectral_axis.min()} to {cube.spectral_axis.max()}")
    
    # Apply spatial smoothing if requested
    if smooth:
        bmaj, bmin, bpa = smooth
        print(f"Smoothing to beam: {bmaj}\" x {bmin}\" @ {bpa}°")
        
        target_beam = Beam(major=bmaj*u.arcsec, 
                          minor=bmin*u.arcsec, 
                          pa=bpa*u.deg)
        
        # Convolve to target resolution
        print("Applying spatial smoothing...")
        cube = cube.convolve_to(target_beam)
        print(f"Smoothed beam: {cube.beam}")
    
    # Convert spectral axis to velocity
    print(f"Converting spectral axis to velocity ({convention})...")
    
    if convention == 'optical':
        # Optical velocity convention
        cube_velo = cube.with_spectral_unit(u.km/u.s, 
                                             velocity_convention='optical',
                                             rest_value=rest_freq*u.Hz)
    elif convention == 'radio':
        # Radio velocity convention
        cube_velo = cube.with_spectral_unit(u.km/u.s,
                                             velocity_convention='radio',
                                             rest_value=rest_freq*u.Hz)
    elif convention == 'relativistic':
        # Relativistic velocity convention
        cube_velo = cube.with_spectral_unit(u.km/u.s,
                                             velocity_convention='relativistic',
                                             rest_value=rest_freq*u.Hz)
    else:
        raise ValueError("Convention must be 'optical', 'radio', or 'relativistic'")
    
    print(f"Velocity range: {cube_velo.spectral_axis.min():.1f} to {cube_velo.spectral_axis.max():.1f}")
    
    # Update header with spectral reference frame
    header = cube_velo.header.copy()
    header['SPECSYS'] = (specsys_correction(specsys), 'Spectral reference frame')
    header['RESTFRQ'] = (rest_freq, 'Rest Frequency (Hz)')
    
    # Generate output filename
    base = (input_file.replace('.fits', '') if input_file.endswith('.fits') else input_file)
    if smooth:
        base += '_smooth'
    
    if convention == 'optical':
        output_file = base + '_vopt.fits'
    else:
        output_file = base + '_vrad.fits'
    
    # Save velocity cube
    print(f"Saving velocity cube: {output_file}")
    cube_velo.write(output_file, overwrite=True)
    
    print(f"Successfully created: {output_file}")
    return output_file

def convert_vel_to_freq_cube(input_file,
                              rest_freq=1420.40575178e6,
                              specsys='TOPOCENT',
                              smooth=None):
    """
    Convert velocity cube back to frequency cube using spectral-cube
    
    Parameters:
    input_file: Input FITS file with velocity axis
    rest_freq: Rest frequency in Hz (default: HI 21cm line)
    specsys: Spectral reference frame
    smooth: Tuple of (bmaj, bmin, bpa) in arcsec and degrees for smoothing
    """
    print(f"Converting velocity cube: {input_file}")
    rest_freq = parse_freq_string(rest_freq)
    print(f"Rest frequency: {rest_freq/1e6:.6f} MHz")
    
    # Load cube with spectral-cube
    cube = SpectralCube.read(input_file)
    print(f"Input cube shape: {cube.shape}")
    print(f"Spectral axis: {cube.spectral_axis.min()} to {cube.spectral_axis.max()}")
    
    # Apply spatial smoothing if requested
    if smooth:
        bmaj, bmin, bpa = smooth
        print(f"Smoothing to beam: {bmaj}\" x {bmin}\" @ {bpa}°")
        
        target_beam = Beam(major=bmaj*u.arcsec, 
                          minor=bmin*u.arcsec, 
                          pa=bpa*u.deg)
        
        # Convolve to target resolution
        print("Applying spatial smoothing...")
        cube = cube.convolve_to(target_beam)
        print(f"Smoothed beam: {cube.beam}")
    
    # Convert spectral axis to frequency
    print(f"Converting spectral axis to frequency...")
    
    # Convert to frequency (Hz)
    cube_freq = cube.with_spectral_unit(u.Hz, rest_value=rest_freq*u.Hz)
    
    print(f"Frequency range: {cube_freq.spectral_axis.min().to(u.MHz):.3f} to {cube_freq.spectral_axis.max().to(u.MHz):.3f}")
    
    # Update header with spectral reference frame
    header = cube_freq.header.copy()
    header['SPECSYS'] = (specsys_correction(specsys), 'Spectral reference frame')
    header['RESTFRQ'] = (rest_freq, 'Rest Frequency (Hz)')
    
    # Generate output filename
    base = (input_file.replace('.fits', '') if input_file.endswith('.fits') else input_file).replace('_vopt', '').replace('_vrad', '')
    if smooth:
        base += '_smooth'
    
    output_file = base + '_freq.fits'
    
    # Save frequency cube
    print(f"Saving frequency cube: {output_file}")
    cube_freq.write(output_file, overwrite=True)
    
    print(f"Successfully created: {output_file}")
    return output_file

def convert_cube(input_file,
                 rest_freq=1420.40575178e6,
                 convention='radio',
                 specsys='TOPOCENT',
                 smooth=None,
                 direction='auto'):
    """
    Convert cube between frequency and velocity with automatic detection
    
    Parameters:
    direction: 'auto', 'freq2vel', or 'vel2freq'
    """
    # Load cube to detect spectral axis type
    cube = SpectralCube.read(input_file)
    axis_type = detect_spectral_axis_type(cube)
    
    print(f"Detected spectral axis type: {axis_type}")
    
    # Determine conversion direction
    if direction == 'auto':
        if axis_type == 'frequency':
            direction = 'freq2vel'
        elif axis_type == 'velocity':
            direction = 'vel2freq'
        else:
            raise ValueError(f"Cannot auto-detect spectral axis type. Please specify --direction")
    
    # Perform conversion
    if direction == 'freq2vel':
        if axis_type == 'velocity':
            print("Warning: Input appears to already be in velocity. Converting anyway...")
        return convert_freq_to_vel_cube(input_file, rest_freq, convention, specsys, smooth)
    
    elif direction == 'vel2freq':
        if axis_type == 'frequency':
            print("Warning: Input appears to already be in frequency. Converting anyway...")
        return convert_vel_to_freq_cube(input_file, rest_freq, specsys, smooth)
    
    else:
        raise ValueError("Direction must be 'auto', 'freq2vel', or 'vel2freq'")

def batch_convert_directory(input_dir, rest_freq=1420.40575178e6,
                           convention='radio', specsys='TOPOCENT',
                           pattern=['*.fits', '*.image'], smooth=None, direction='auto'):
    """
    Convert all FITS files in a directory between frequency and velocity
    """
    import glob
    
    in_files = []
    for p in pattern:
        in_files.extend(glob.glob(os.path.join(input_dir, p)))

    if not in_files:
        print(f"No relevant image (FITS or CASA) files found in {input_dir}")
        return
    
    print(f"Found {len(in_files)} image files to convert")
    
    converted_files = []
    failed_files = []
    
    for i, input_file in enumerate(in_files, 1):
        print(f"\n[{i}/{len(in_files)}] Processing: {os.path.basename(input_file)}")
        try:
            output_file = convert_cube(input_file,
                                       rest_freq=rest_freq,
                                       convention=convention,
                                       specsys=specsys,
                                       smooth=smooth,
                                       direction=direction)
            converted_files.append(output_file)
        except Exception as e:
            print(f"Error converting {input_file}: {e}")
            failed_files.append(input_file)
    
    print(f"\nBatch conversion complete!")
    print(f"Successfully converted: {len(converted_files)} files")
    print(f"Failed conversions: {len(failed_files)} files")
    
    if failed_files:
        print("Failed files:")
        for f in failed_files:
            print(f"  {f}")

def main():
    parser = argparse.ArgumentParser(
        description='Convert FITS/CASA image cubes between frequency and velocity and Export as FITS files',
        epilog='Examples:\n'
               '  python cube_process.py cube.fits                                       # Auto-detect and convert\n'
               '  python cube_process.py cube.image --direction freq2vel                 # Force freq to vel\n'
               '  python cube_process.py cube.fits --convention optical                  # Use optical definition of velocity\n'
               '  python cube_process.py cube.image --smooth 15.0 15.0 0.0               # Convert with smoothing\n'
               '  python cube_process.py --batch --fitsonly data_/data/image_dir         # Batch convert only .FITS files in directory\n',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )


    parser.add_argument('input', help='Input FITS/CASA image file or directory (for batch mode)')
    parser.add_argument('-d', '--direction', choices=['auto', 'freq2vel', 'vel2freq'],
                       default='auto', 
                       help='Conversion direction (default: auto-detect)')
    parser.add_argument('-c', '--convention', choices=['optical', 'radio', 'relativistic'],
                       default='radio', 
                       help='Velocity convention for freq2vel (default: radio)')
    parser.add_argument('-s', '--specsys', default='BARYCENT', choices=np.sort(list(SPECSYS_MAP.keys())),
                       help='Spectral reference frame (default: BARYCENT)')
    parser.add_argument('--smooth', nargs=3, type=float, default=None, 
                       metavar=('BMAJ', 'BMIN', 'BPA'),
                       help='Smoothing kernel: bmaj bmin (arcsec) and PA (deg) (default: no smoothing)')
    
    parser.add_argument('--batch', action='store_true',
                       help='Process all FITS/CASA image files in input directory')
    
    # Batch file pattern shortcuts
    group_pattern = parser.add_mutually_exclusive_group()
    group_pattern.add_argument('--all', action='store_const', dest='pattern',
                       const=['*.fits', '*.image'], help='Use Both .FITS and .image files (For batch mode).')
    group_pattern.add_argument('--fitsonly', action='store_const', dest='pattern',
                       const=['*.fits'], help='Use Only .FITS files (For batch mode).')
    group_pattern.add_argument('--imageonly', action='store_const', dest='pattern',
                       const=['*.image'], help='Use Only .image files (For batch mode).')
    
    # Frequency shortcuts
    group_freq = parser.add_mutually_exclusive_group()
    group_freq.add_argument('--hi', action='store_const', dest='rest_freq',
                       const=1420.40575178e6, help='Use HI 21cm rest frequency (1420.40575178 MHz). <-- default')
    group_freq.add_argument('--oh1665', action='store_const', dest='rest_freq',
                       const=1665.4018e6, help='Use OH 1665 MHz rest frequency.')
    group_freq.add_argument('--oh1667', action='store_const', dest='rest_freq',
                       const=1667.3590e6, help='Use OH 1667 MHz rest frequency.')
    group_freq.add_argument('--rest-freq', type=float, help='Custom rest frequency in Hz.')


    args = parser.parse_args()
    
    rest_freq = args.rest_freq if args.rest_freq else DEFAULT_RESTFREQ
    pattern = args.pattern if args.pattern else DEFAULT_PATTERN

    input, convention, direction, specsys, smooth = args.input, args.convention, args.direction, args.specsys, args.smooth

    # if args.batch or os.path.isdir(args.input):
    if args.batch:
        # Batch processing
        batch_convert_directory(input,
                              rest_freq=rest_freq,
                              convention=convention,
                              specsys=specsys,
                              pattern=pattern,
                              smooth=smooth,
                              direction=direction)
    else:
        # Single file processing
        convert_cube(input,
                    rest_freq=rest_freq,
                    convention=convention,
                    specsys=specsys,
                    smooth=smooth,
                    direction=direction)

if __name__ == "__main__":
    main()

