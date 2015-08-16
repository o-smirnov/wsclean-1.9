#ifndef MULTISCALE_DECONVOLUTION_H
#define MULTISCALE_DECONVOLUTION_H

#include <memory>
#include <string>

#include "../uvector.h"

#include "deconvolutionalgorithm.h"
#include "dynamicset.h"

#include "../multiscale/multiscalealgorithm.h"

#include "../wsclean/imagingtable.h"

class MultiScaleDeconvolution : public UntypedDeconvolutionAlgorithm
{
public:
	MultiScaleDeconvolution(class ImageBufferAllocator& allocator, double beamSize, double pixelScaleX, double pixelScaleY) :
		_allocator(allocator),
		_beamSizeInPixels(beamSize / std::max(pixelScaleX, pixelScaleY))
	{ }
	
	virtual void ExecuteMajorIteration(DynamicSet& dataImage, DynamicSet& modelImage, const ao::uvector<const double*>& psfImages, size_t width, size_t height, bool& reachedMajorThreshold)
	{
		MultiScaleAlgorithm algorithm(_allocator, width, height, _beamSizeInPixels, _threshold, _subtractionGain, _stopGain, _cleanBorderRatio, _allowNegativeComponents, _multiscaleScaleBias, _threadCount);
		
		if(_cleanMask != 0)
			algorithm.SetCleanMask(_cleanMask);
		
		algorithm.PerformMajorIteration(_iterationNumber, MaxNIter(), modelImage, dataImage, psfImages, reachedMajorThreshold);
	}
	
private:
	class ImageBufferAllocator& _allocator;
	double _beamSizeInPixels;
};

#endif
