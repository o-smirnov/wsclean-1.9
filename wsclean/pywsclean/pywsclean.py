"""
A Python interface to WSClean
This wrapper can be used to call the (C++) WSClean imager.
"""

class WSClean(object):
	"""The Python interface to WSClean
	"""
	
	datacolumn=''
	"""The column used for imaging; empty means CORRECTED_DATA if it exists, otherwise
	use DATA."""
	
	width=1024
	"""Image width"""
	
	height=1024
	"""Image height"""
	
	scale='1asec'
	"""Pixel scale of image. Units can e.g. be deg, amin, asec, masec. There should not
	be a space between the number and its unit."""
	
	niter=0
	"""Number of clean or moresane iterations"""
	
	gain=-1
	"""Gain per minor iteration. -1 means use WSClean's default."""
	
	mgain=-1
	"""Gain per major iteration. -1 means use WSClean's default."""
	
	__weightpar=''
	
	def __init__(self):
		return;
	
	def image(self, msnames, nameprefix):
		"""Run WSClean to make an image on the specified list of measurement sets"""
		plist=self.__get_parameterlist(nameprefix)
		import os
		msnamelist=' '.join(msnames)
		cmd='wsclean '+str(plist)+' '+msnamelist
		print cmd
		os.system(cmd)
		return;
	
	def predict(self, msnames, nameprefix):
		"""Run WSClean to predict"""
		plist=self.__get_parameterlist(nameprefix)
		import os
		msnamelist=' '.join(msnames)
		cmd='wsclean -predict '+str(plist)+' '+msnamelist
		print cmd
		os.system(cmd)
		return;
	
	def set_uniform_weighting(self):
		self.__weightpar='-weight uniform'
		
	def set_natural_weighting(self):
		self.__weightpar='-weight natural'
		
	def set_briggs_weighting(self, robustness):
		self.__weightpar='-weight briggs '+str(robustness)
	
	def __get_parameterlist(self, prefixname):
		plist = '-size '+str(self.width)+' '+str(self.height)+' -scale '+str(self.scale);
		if self.datacolumn!='':
			plist += ' -datacolumn '+self.datacolumn;
		if self.__weightpar!='':
			plist += ' '+self.__weightpar;
		if self.niter!=0:
			plist += ' -niter '+str(self.niter);
		if self.gain!=-1:
			plist += ' -gain '+str(self.gain);
		if self.mgain!=-1:
			plist += ' -mgain '+str(self.mgain);
		if prefixname!='':
			plist += ' -name '+prefixname;
		return plist;
