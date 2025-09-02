#!/usr/bin/env python3
"""
Convert SoFiA mask to CASA-compatible format.
Reads SoFiA output, extracts beam info from mom0, creates binary mask for tclean.
"""

import numpy as np
import os
from astropy.io import fits
from astropy.wcs import WCS
import argparse
import glob
from scipy.ndimage import label


def find_sofia_files(sofia_dir):
    """
    Find SoFiA output files in the given directory.
    
    Parameters:
    -----------
    sofia_dir : str
        Path to SoFiA output directory
        
    Returns:
    --------
    dict : Dictionary with paths to mask and mom0 files
    """
    
    # Common SoFiA output patterns
    mask_patterns = ['*_mask.fits']
    mom0_patterns = ['*_mom0.fits']
    
    files_found = {'mask': None, 'mom0': None}
    
    # Find mask file
    for pattern in mask_patterns:
        mask_files = glob.glob(os.path.join(sofia_dir, pattern))
        if mask_files:
            files_found['mask'] = mask_files[0]  # Take first match
            break
    
    # Find mom0 file
    for pattern in mom0_patterns:
        mom0_files = glob.glob(os.path.join(sofia_dir, pattern))
        if mom0_files:
            files_found['mom0'] = mom0_files[0]  # Take first match
            break
    
    return files_found

def extract_beam_from_mom0(mom0_file):
    """
    Extract beam parameters from SoFiA mom0 FITS file.
    
    Parameters:
    -----------
    mom0_file : str
        Path to mom0 FITS file
        
    Returns:
    --------
    dict : Beam parameters (bmaj, bmin, bpa in degrees)
    """
    
    print(f"Reading beam parameters from: {mom0_file}")
    
    with fits.open(mom0_file) as hdul:
        header = hdul[0].header
        
        # Try different possible beam parameter names
        beam_keys = {
            'bmaj': ['BMAJ', 'BMMAJ', 'BEAM_MAJ'],
            'bmin': ['BMIN', 'BMMIN', 'BEAM_MIN'], 
            'bpa': ['BPA', 'BMPA', 'BEAM_PA', 'BANGLE']
        }
        
        beam_params = {}
        
        for param, possible_keys in beam_keys.items():
            value = None
            for key in possible_keys:
                if key in header:
                    value = header[key]
                    print(f"Found {param}: {key} = {value}")
                    break
            
            if value is None:
                if param == 'bpa':
                    value = 0.0  # Default PA
                    print(f"Warning: No beam PA found, using default: {value}")
                else:
                    raise ValueError(f"Could not find beam parameter {param} in {mom0_file}")
            
            beam_params[param] = value
    
    # Convert to degrees if needed (some files store in arcsec)
    for param in ['bmaj', 'bmin']:
        if beam_params[param] > 1.0:  # Likely in arcsec
            print(f"Converting {param} from arcsec to degrees: {beam_params[param]} -> {beam_params[param]/3600.0}")
            beam_params[param] = beam_params[param] / 3600.0
    
    return beam_params

def convert_sofia_mask_to_casa(mask_file, mom0_file, output_file, central_mask=True,
                              binary_threshold=0, source_ids=None):
    """
    Convert SoFiA mask to CASA-compatible binary mask.
    
    Parameters:
    -----------
    mask_file : str
        Path to SoFiA mask FITS file
    mom0_file : str
        Path to mom0 file (for beam parameters)
    output_file : str
        Output CASA-ready mask file
    binary_threshold : float
        Threshold for creating binary mask (0 = any non-zero value)
    source_ids : list, optional
        List of specific source IDs to include in mask
    """
    
    print(f"Converting SoFiA mask: {mask_file}")
    print(f"Output file: {output_file}")
    
    # Extract beam parameters from mom0
    beam_params = extract_beam_from_mom0(mom0_file)
    
    # Read the SoFiA mask
    with fits.open(mask_file) as hdul:
        mask_data = hdul[0].data
        mask_header = hdul[0].header.copy()
        mask_data = mask_data.astype(np.int32)  # Ensure integer type for labels
        
        print(f"Original mask shape: {mask_data.shape}")
        print(f"Original mask data type: {mask_data.dtype}")
        print(f"Mask value range: {np.min(mask_data)} to {np.max(mask_data)}")
        
        # Get unique source IDs
        unique_ids = np.unique(mask_data)
        unique_ids = unique_ids[unique_ids > 0]  # Remove background (0)
        print(f"Found {len(unique_ids)} sources with IDs: {unique_ids}")

        if central_mask:
            
            Nz, Ny, Nx = mask_data.shape
            centre = (Nz//2, Ny//2, Nx//2)
            # radius = 5
            # mask_data[0:(Nz//2)-radius,:,:] = 0
            # mask_data[(Nz//2)+radius:,:,:] = 0
            structure = np.ones((3,3,3), dtype=int)   # 26-connected neighbourhood
            labels, n_islands = label(mask_data, structure=structure)

            cz, cy, cx = centre
            target_label = labels[cz, cy, cx]

            if target_label == 0:
                cz += 1
                target_label = labels[cz, cy, cx]
                if target_label == 0:
                    cz-=2
                    target_label = labels[cz, cy, cx]
                    if target_label == 0:
                        raise ValueError("The centre voxel is not masked – check centre coordinates!")
            
            central_mask = (labels == target_label).astype(int)  # 1 = central object

            mask_data = central_mask
        
        # Create binary mask
        if source_ids is not None:
            # Include only specified source IDs
            print(f"Including only source IDs: {source_ids}")
            binary_mask = np.zeros_like(mask_data, dtype=np.float32)
            for src_id in source_ids:
                binary_mask[mask_data == src_id] = 1.0
        else:
            # Include all sources above threshold
            if binary_threshold > 0:
                binary_mask = (mask_data > binary_threshold).astype(np.float32)
            else:
                binary_mask = (mask_data > 0).astype(np.float32)
        
        # Count masked pixels
        n_masked = np.sum(binary_mask > 0)
        n_total = binary_mask.size
        fraction_masked = n_masked / n_total * 100
        
        print(f"Binary mask statistics:")
        print(f"  Masked pixels: {n_masked:,}")
        print(f"  Total pixels: {n_total:,}")
        print(f"  Fraction masked: {fraction_masked:.2f}%")
    
    # Update header with beam parameters
    print("Adding beam parameters to header:")
    mask_header['BMAJ'] = (beam_params['bmaj'], 'Beam major axis (deg)')
    mask_header['BMIN'] = (beam_params['bmin'], 'Beam minor axis (deg)')
    mask_header['BPA'] = (beam_params['bpa'], 'Beam position angle (deg)')
    
    # Add CASA-specific keywords
    mask_header['BTYPE'] = 'Intensity'
    mask_header['BUNIT'] = ' '  # Empty for masks
    
    # Add creation info
    mask_header['HISTORY'] = 'Converted from SoFiA mask for CASA tclean'
    mask_header['HISTORY'] = f'Original mask: {os.path.basename(mask_file)}'
    mask_header['HISTORY'] = f'Beam from: {os.path.basename(mom0_file)}'
    
    for param, value in beam_params.items():
        print(f"  {param.upper()}: {value:.6f} deg")
    
    # Write the binary mask
    hdu = fits.PrimaryHDU(data=binary_mask, header=mask_header)
    hdu.writeto(output_file, overwrite=True)
    
    print(f"Successfully created CASA mask: {output_file}")
    
    return {
        'output_file': output_file,
        'beam_params': beam_params,
        'n_sources': len(unique_ids),
        'n_masked_pixels': n_masked,
        'fraction_masked': fraction_masked
    }

def analyze_sofia_mask(mask_file):
    """
    Analyze SoFiA mask to understand source distribution.
    
    Parameters:
    -----------
    mask_file : str
        Path to SoFiA mask file
        
    Returns:
    --------
    dict : Analysis results
    """
    
    print(f"\n=== Analyzing SoFiA mask: {mask_file} ===")
    
    with fits.open(mask_file) as hdul:
        mask_data = hdul[0].data
        header = hdul[0].header
        
        # Basic info
        print(f"Mask dimensions: {mask_data.shape}")
        print(f"Data type: {mask_data.dtype}")
        
        # Source statistics
        unique_ids = np.unique(mask_data)
        background_pixels = np.sum(mask_data == 0)
        source_pixels = np.sum(mask_data > 0)
        
        print(f"\nSource statistics:")
        print(f"  Background pixels (0): {background_pixels:,}")
        print(f"  Source pixels (>0): {source_pixels:,}")
        print(f"  Total pixels: {mask_data.size:,}")
        print(f"  Source fraction: {source_pixels/mask_data.size*100:.2f}%")
        
        # Individual sources
        source_ids = unique_ids[unique_ids > 0]
        print(f"\nFound {len(source_ids)} individual sources:")
        
        source_info = []
        for src_id in source_ids[:20]:  # Show first 20
            src_pixels = np.sum(mask_data == src_id)
            source_info.append((src_id, src_pixels))
            
        source_info.sort(key=lambda x: x[1], reverse=True)  # Sort by size
        
        print("  Source ID | Pixels")
        print("  ----------|-------")
        for src_id, pixels in source_info:
            print(f"  {src_id:8d} | {pixels:6d}")
        
        if len(source_ids) > 20:
            print(f"  ... and {len(source_ids)-20} more sources")
        
        # Coordinate info
        if 'CRVAL1' in header and 'CRVAL2' in header:
            print(f"\nCoordinate reference:")
            print(f"  RA center: {header['CRVAL1']:.6f} deg")
            print(f"  Dec center: {header['CRVAL2']:.6f} deg")
            if 'CDELT1' in header:
                print(f"  Pixel scale: {abs(header['CDELT1'])*3600:.2f} arcsec/pixel")
        
        return {
            'n_sources': len(source_ids),
            'source_ids': source_ids,
            'source_pixels': source_pixels,
            'background_pixels': background_pixels,
            'source_info': source_info
        }

def main():
    """Main function with command line interface."""
    
    parser = argparse.ArgumentParser(
        description='Convert SoFiA mask to CASA-compatible format'
    )
    parser.add_argument('sofia_dir', help='SoFiA output directory')
    parser.add_argument('-o', '--output', help='Output mask file', 
                       default='casa_mask.fits')
    parser.add_argument('-t', '--threshold', type=float, default=0,
                       help='Binary threshold (0 = any non-zero value)')
    parser.add_argument('--source-ids', type=int, nargs='+', 
                       help='Specific source IDs to include')
    parser.add_argument('-a', '--analyze', action='store_true',
                       help='Analyze mask before conversion')
    
    args = parser.parse_args()

    print(f"Searching for SoFiA files in: {args.sofia_dir}")
    sofia_files = find_sofia_files(args.sofia_dir)
    
    if sofia_files['mask'] is None:
        raise FileNotFoundError(f"No mask file found in {args.sofia_dir}")
    if sofia_files['mom0'] is None:
        raise FileNotFoundError(f"No mom0 file found in {args.sofia_dir}")
    
    mask_file = sofia_files['mask']
    mom0_file = sofia_files['mom0']

    print(f"Using mask file: {mask_file}")
    print(f"Using mom0 file: {mom0_file}")
    
    # Analyze mask if requested
    if args.analyze:
        analyze_sofia_mask(mask_file)
    
    # Convert mask
    result = convert_sofia_mask_to_casa(
        mask_file=mask_file,
        mom0_file=mom0_file,
        output_file=args.output,
        binary_threshold=args.threshold,
        source_ids=args.source_ids
    )
    
    print(f"\n=== Conversion Summary ===")
    print(f"Input mask: {mask_file}")
    print(f"Output mask: {result['output_file']}")
    print(f"Sources found: {result['n_sources']}")
    print(f"Pixels masked: {result['n_masked_pixels']:,} ({result['fraction_masked']:.2f}%)")
    print(f"Beam parameters:")
    for param, value in result['beam_params'].items():
        print(f"  {param.upper()}: {value:.6f} deg ({value*3600:.2f} arcsec)")

    return result


if __name__ == '__main__':
    main()
