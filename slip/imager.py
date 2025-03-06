try:
    from casampi import MPIEnvironment, MPICommandClient
except ImportError:
    casampi = None

import sys
from casatasks import tclean

if __name__ == '__main__':

    ms_uvsub, out_tclean, field, spw_cube, outframe, veltype, restfreq, tcell, tuvrange, tuvtaper, imsize, weighting, niter, cycleniter, nsigma, robust, width = sys.argv[1:]

    timsize = [int(imsize), int(imsize)]
    niter = int(niter)
    cycleniter = int(cycleniter)
    nsigma = float(nsigma)
    robust = float(robust)
    width = float(width)

    twidth = str(width) + 'km/s'

    if weighting != 'natural':
        if width == -999.:
            tclean(vis=ms_uvsub, imagename=out_tclean, field=field, spw=spw_cube, datacolumn='data', specmode='cube', outframe=outframe, veltype=veltype, restfreq=restfreq, cell=tcell, uvrange=tuvrange, uvtaper=tuvtaper, deconvolver='clark', imsize=timsize, weighting=weighting, niter=niter, cycleniter=cycleniter, interactive=False, parallel=True, usemask='auto-multithresh', nsigma=nsigma, robust=robust, restoringbeam='common', savemodel='none')
        else:
            tclean(vis=ms_uvsub, imagename=out_tclean, field=field, spw=spw_cube, datacolumn='data', specmode='cube', outframe=outframe,  veltype=veltype,restfreq=restfreq, cell=tcell, uvrange=tuvrange, uvtaper=tuvtaper, deconvolver='clark', imsize=timsize, weighting=weighting, niter=niter, cycleniter=cycleniter, interactive=False, parallel=True, usemask='auto-multithresh', nsigma=nsigma, robust=robust, restoringbeam='common', width=twidth, savemodel='none')
    else:
        if width == -999.:
            tclean(vis=ms_uvsub, imagename=out_tclean, field=field, spw=spw_cube, datacolumn='data', specmode='cube', outframe=outframe,  veltype=veltype,restfreq=restfreq, cell=tcell, uvrange=tuvrange, uvtaper=tuvtaper, deconvolver='clark', imsize=timsize, weighting=weighting, niter=niter, cycleniter=cycleniter, interactive=False, parallel=True, usemask='auto-multithresh', nsigma=nsigma, restoringbeam='common', savemodel='none')
        else:
            tclean(vis=ms_uvsub, imagename=out_tclean, field=field, spw=spw_cube, datacolumn='data', specmode='cube', outframe=outframe, veltype=veltype, restfreq=restfreq, cell=tcell, uvrange=tuvrange, uvtaper=tuvtaper, deconvolver='clark', imsize=timsize, weighting=weighting, niter=niter, cycleniter=cycleniter, interactive=False, parallel=True, usemask='auto-multithresh', nsigma=nsigma, restoringbeam='common', width=twidth, savemodel='none')