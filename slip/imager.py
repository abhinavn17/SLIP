try:
    from casampi import MPIEnvironment, MPICommandClient, MPIInterface
except ImportError:
    casampi = None

import sys
from casatasks import tclean, casalog

# client = MPICommandClient.MPICommandClient()

# client.set_log_mode('redirect')

if __name__ == '__main__':

    ms_uvsub, out_tclean, field, datacolumn, spw_cube, outframe, veltype, restfreq, tcell, tuvrange, tuvtaper, imsize, weighting, niter, cycleniter, nsigma, robust, width, deconvolver, usemask, mask, nmajor = sys.argv[1:]

    timsize = [int(imsize), int(imsize)]
    niter = int(niter)
    cycleniter = int(cycleniter)
    nsigma = float(nsigma)
    robust = float(robust)
    nmajor = int(nmajor)

    tclean(vis=ms_uvsub, imagename=out_tclean, field=field, spw=spw_cube, datacolumn=datacolumn, specmode='cube', outframe=outframe,  veltype=veltype,restfreq=restfreq, cell=tcell, uvrange=tuvrange, uvtaper=tuvtaper, deconvolver=deconvolver, imsize=timsize, weighting=weighting, niter=niter, cycleniter=cycleniter, interactive=False, parallel=True, usemask=usemask, mask= mask, robust=robust, restoringbeam='common', width=width, savemodel='none', nmajor=nmajor, nsigma=nsigma)
