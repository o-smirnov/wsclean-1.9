#ifndef MORESANE_H
#define MORESANE_H

#include <string>

#include "cleanalgorithm.h"
#include "imageset.h"

class MoreSane : public TypedCleanAlgorithm<clean_algorithms::SingleImageSet>
{
	public:
		MoreSane(const std::string& moreSaneLocation, const std::string& moresaneArguments) : _moresaneLocation(moreSaneLocation), _moresaneArguments(moresaneArguments)
		{ }
		
    virtual void ExecuteMajorIteration(ImageSet& dataImage, ImageSet& modelImage, std::vector<double*> psfImages, size_t width, size_t height, bool& reachedMajorThreshold)
		{
      _allocator = dataImage.Allocator();
			ExecuteMajorIteration(dataImage.GetImage(0), modelImage.GetImage(0), psfImages[0], width, height, reachedMajorThreshold);
		}
		
		void ExecuteMajorIteration(double* dataImage, double* modelImage, const double* psfImage, size_t width, size_t height, bool& reachedMajorThreshold);
	private:
		const std::string _moresaneLocation, _moresaneArguments;
		ImageBufferAllocator<double>* _allocator;
};

#endif
