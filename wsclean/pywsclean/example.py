from pywsclean import WSClean

wsc=WSClean()
wsc.width=2048
wsc.height=2048
wsc.scale='5asec'
wsc.datacolumn='DATA'

wsc.set_uniform_weighting()
# Alternatively: wsc.set_natural_weighting()

# This predicts from wsclean-model.fits into the MODEL_DATA column
# (note that the column is always called MODEL_DATA, the datacolumn parameter
#  only sets the column used in the 'image' command)
wsc.predict(['myset.ms'], 'wsclean')

# This makes wsclean-dirty.fits (et al) from the DATA column
wsc.image(['myset.ms'], 'wsclean')
