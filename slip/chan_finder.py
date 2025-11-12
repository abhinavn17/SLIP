#!/usr/bin/env python3

import argparse
import numpy as np
import sys
from casatools import table

# Speed of light (km/s)
C_KMS = 299792.458

# Default HI rest frequency (Hz)
DEFAULT_RESTFREQ = 1420.40575178e6


def compute_velocity(freq, rest_freq, definition):
    if definition == 'optical' or definition == 'opt':
        return C_KMS * (rest_freq - freq) / freq
    elif definition == 'radio' or definition == 'rad':
        return C_KMS * (rest_freq - freq) / rest_freq
    elif definition == 'relativistic' or definition == 'rel':
        ratio = rest_freq / freq
        return C_KMS * (ratio**2 - 1) / (ratio**2 + 1)
    else:
        raise ValueError("Definition must be 'optical', 'radio', or 'relativistic'. instead it is " + str(definition))


def parse_cutout_input(cut_str, velocc, chan_width):
    cut_str = str(cut_str).strip().lower()
    try:
        if cut_str.endswith(('km/s', 'kms', 'm/s', 'ms')):
            if 'k' in cut_str:
                if cut_str.endswith('km/s'):
                    cut_value = float(cut_str[:-4])
                elif cut_str.endswith('kms'):
                    cut_value = float(cut_str[:-3])
            elif 'm' in cut_str and 'k' not in cut_str:
                if cut_str.endswith('m/s'):
                    cut_value = float(cut_str[:-3]) / 1000.0
                elif cut_str.endswith('ms'):
                    cut_value = float(cut_str[:-2]) / 1000.0
            else:
                raise ValueError("Invalid velocity unit.")
            chan_offset = int(np.ceil(cut_value / np.abs(np.diff(velocc).mean())))
        elif cut_str.endswith(('hz', 'khz', 'mhz', 'ghz')):
            if cut_str.endswith('ghz'):
                cut_value = float(cut_str[:-3]) * 1e9
            elif cut_str.endswith('mhz'):
                cut_value = float(cut_str[:-3]) * 1e6
            elif cut_str.endswith('khz'):
                cut_value = float(cut_str[:-3]) * 1e3
            elif cut_str.endswith('hz') and ('k' not in cut_str and 'm' not in cut_str and 'g' not in cut_str):
                cut_value = float(cut_str[:-2])
            else:
                raise ValueError("Invalid frequency unit.")
            chan_offset = int(np.ceil(cut_value / chan_width))
        else:
            chan_offset = int(float(cut_str))
            print("Defaulting to channel number interpretation.")
    except ValueError:
        chan_offset = None
        print("Error parsing cutout input. Ensure correct format and units." \
        "\nCorrect examples: 100km/s, 0.1MHz, 50Hz, 200m/s, 2GHz, 500kHz or 200.")

    return chan_offset


def find_channel_range(msfile, velocity, definition='optical', cutout=None, rest_freq=None):
    if rest_freq is None:
        rest_freq = DEFAULT_RESTFREQ

    # Open Measurement Set and read spectral window info
    data_table = table(msfile + '/SPECTRAL_WINDOW', readonly=True)
    nfreq = data_table.getcol('CHAN_FREQ')
    nchan = int(data_table.getcol('NUM_CHAN')[0])
    data_table.close()

    freqbeg, freqend = nfreq[0][0], nfreq[-1][0]
    freqq = np.linspace(freqbeg, freqend, nchan)
    chan_width = np.abs(np.diff(freqq).mean())

    velocc = compute_velocity(freqq, rest_freq, definition)
    chan_no = np.argmin(np.abs(velocc - velocity))
    found_vel = velocc[chan_no]
    found_freq = freqq[chan_no]

    print(f"Channel number corresponding to {velocity} km/s "
          f"({definition} definition): {chan_no}")
    print(f"Found velocity at this channel: {found_vel:.4f} km/s")
    print(f"Found frequency at this channel: {found_freq/1e6:.6f} MHz")

    if cutout:
        try:
            chan_offset = parse_cutout_input(cutout, velocc, chan_width)
            cutout_range = (max(0, chan_no - chan_offset),
                            min(nchan - 1, chan_no + chan_offset - 1))
            print(f"Cutout channel range: {cutout_range[0]}–{cutout_range[1]}")
            return (0, nchan - 1, cutout_range[0], cutout_range[1])
        except ValueError as e:
            print(f"Error: {e}")
            return (0, nchan - 1, None, None)
    else:
        print("No cutout selected.")
        return (0, nchan - 1, None, None)

def main():
    parser = argparse.ArgumentParser(
        description="Find the channel corresponding to a given velocity in a CASA Measurement Set."
    )

    parser.add_argument("msfile", help="Path to the input Measurement Set (*.ms).")
    # parser.add_argument("-v","--velocity", type=float, help="Target velocity in km/s.")
    parser.add_argument("-v","--velocity", nargs = 2, type = str, default = ['optical', '0.0'], help='Velocity definition and central line velocity in km/s (e.g., optical 1000.0)')

    # parser.add_argument(
    #     "-d", "--definition", choices=["optical", "opt", "radio", "rad", "relativistic", "rel"],
    #     default="optical", help="Velocity convention opt[ical]/rad[io]/rel[ativistic] (default: optical)."
    # )
    parser.add_argument(
        "-c", "--cutout", type=str, default=None,
        help="Cutout bandwidth around the found channel (e.g., 100km/s or 0.1MHz). \nTotal bandwidth of the output will be double this value (e.g. 200km/s or 0.2MHz)"
    )

    # Frequency shortcuts
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--hi', action='store_const', dest='rest_freq',
                       const=1420.40575178e6, help='Use HI 21cm rest frequency (1420.40575178 MHz).')
    group.add_argument('--oh1665', action='store_const', dest='rest_freq',
                       const=1665.4018e6, help='Use OH 1665 MHz rest frequency.')
    group.add_argument('--oh1667', action='store_const', dest='rest_freq',
                       const=1667.3590e6, help='Use OH 1667 MHz rest frequency.')
    group.add_argument('--rest-freq', type=float, help='Custom rest frequency in Hz.')

    args = parser.parse_args()
    vel = args.velocity
    vdef = vel[0]
    velocity = float(vel[1].strip())

    rest_freq = args.rest_freq if args.rest_freq else DEFAULT_RESTFREQ

    find_channel_range(args.msfile, velocity, vdef, args.cutout, rest_freq)

if __name__ == "__main__":
    main()
