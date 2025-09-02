try:
    from casampi import MPIEnvironment, MPICommandClient
except ImportError:
    casampi = None

import sys
from casatasks import tclean
import casaviewer

if __name__ == '__main__':

    ms_uvsub, out_tclean, field, datacolumn, spw_cube, outframe, veltype, restfreq, tcell, tuvrange, tuvtaper, imsize, weighting, niter, cycleniter, nsigma, robust, width = sys.argv[1:]

    timsize = [int(imsize), int(imsize)]
    niter = int(niter)
    cycleniter = int(cycleniter)
    nsigma = float(nsigma)
    robust = float(robust)
    width = float(width)

    twidth = str(width) + 'km/s'

    tclean(vis=ms_uvsub, imagename=out_tclean, field=field, spw=spw_cube, datacolumn=datacolumn, specmode='cube', outframe=outframe, veltype=veltype, restfreq=restfreq, cell=tcell, uvrange=tuvrange, uvtaper=tuvtaper, deconvolver='multiscale', scales = [0, 12, 24, 49, 98 ,195], imsize=timsize, weighting=weighting, niter=niter, cycleniter=cycleniter, interactive=False, parallel=True, usemask='user', mask = 'Crab_mask_2.crtf', nsigma=nsigma, robust=robust, restoringbeam='common', savemodel='none')