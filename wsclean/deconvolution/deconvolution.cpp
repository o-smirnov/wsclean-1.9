#include "deconvolution.h"

#include "deconvolutionalgorithm.h"
#include "joinedclean.h"
#include "simpleclean.h"
#include "moresane.h"
#include "multiscaleclean.h"

#include "../casamaskreader.h"
#include "../wsclean/imagingtable.h"

Deconvolution::Deconvolution() :
	_threshold(0.0), _gain(0.1), _mGain(1.0),
	_nIter(0),
	_allowNegative(true), 
	_stopOnNegative(false),
	_multiscale(false),
	_multiscaleThresholdBias(0.7), _multiscaleScaleBias(0.6),
	_cleanBorderRatio(0.05),
	_fitsMask(), _casaMask(),
	_useMoreSane(false),
	_moreSaneLocation(), _moreSaneArgs()
{
}

Deconvolution::~Deconvolution()
{
	FreeDeconvolutionAlgorithms();
}

void Deconvolution::Perform(const class ImagingTable& groupTable, bool& reachedMajorThreshold, size_t majorIterationNr)
{
	std::cout << std::flush << " == Cleaning (" << majorIterationNr << ") ==\n";
	
	if(_summedCount != 1)
	{
		if(_squaredCount == 4)
			performJoinedPolFreqClean<4>(reachedMajorThreshold, majorIterationNr);
		else if(_squaredCount == 2)
			performJoinedPolFreqClean<2>(reachedMajorThreshold, majorIterationNr);
		else // if(_squaredCount == 1)
			performJoinedPolFreqClean<1>(reachedMajorThreshold, majorIterationNr);
	}
	else if(_squaredCount != 1) {
		size_t currentChannelIndex = groupTable.Front().outputChannelIndex;
		if(_squaredCount == 4)
			performJoinedPolClean<4>(currentChannelIndex, reachedMajorThreshold, majorIterationNr);
		else // if(_squaredCount == 2)
			performJoinedPolClean<2>(currentChannelIndex, reachedMajorThreshold, majorIterationNr);
	}
	else {
		performSimpleClean(groupTable[0].outputChannelIndex, reachedMajorThreshold, majorIterationNr, groupTable[0].polarization);
	}
}

void Deconvolution::performSimpleClean(size_t currentChannelIndex, bool& reachedMajorThreshold, size_t majorIterationNr, PolarizationEnum polarization)
{
	deconvolution::SingleImageSet
		residualImage(_imgWidth*_imgHeight, *_imageAllocator),
		modelImage(_imgWidth*_imgHeight, *_imageAllocator);
	ImageBufferAllocator<double>::Ptr psfImage;
	_imageAllocator->Allocate(_imgWidth*_imgHeight, psfImage);
		
	_residualImages->Load(residualImage.Data(), polarization, currentChannelIndex, false);
	_modelImages->Load(modelImage.Data(), polarization, currentChannelIndex, false);
	_psfImages->Load(psfImage.data(), _psfPolarization, currentChannelIndex, false);
	
	std::vector<double*> psfs(1, psfImage.data());
	TypedDeconvolutionAlgorithm<deconvolution::SingleImageSet>& tAlgorithm =
		static_cast<TypedDeconvolutionAlgorithm<deconvolution::SingleImageSet>&>(*_cleanAlgorithm);
	tAlgorithm.ExecuteMajorIteration(residualImage, modelImage, psfs, _imgWidth, _imgHeight, reachedMajorThreshold);
	
	_modelImages->Store(modelImage.Data(), polarization, currentChannelIndex, false);
	_residualImages->Store(residualImage.Data(), polarization, currentChannelIndex, false);
}

template<size_t PolCount>
void Deconvolution::performJoinedPolClean(size_t currentChannelIndex, bool& reachedMajorThreshold, size_t majorIterationNr)
{
	typename JoinedClean<deconvolution::PolarizedImageSet<PolCount>>::ImageSet
		modelSet(_imgWidth*_imgHeight, *_imageAllocator),
		residualSet(_imgWidth*_imgHeight, *_imageAllocator);
	
	ImageBufferAllocator<double>::Ptr psfImage;
	_imageAllocator->Allocate(_imgWidth*_imgHeight, psfImage);
	_psfImages->Load(psfImage.data(), _psfPolarization, currentChannelIndex, false);
	
	modelSet.Load(*_modelImages, _polarizations, currentChannelIndex);
	residualSet.Load(*_residualImages, _polarizations, currentChannelIndex);
	
	std::vector<double*> psfImages(1, psfImage.data());
	static_cast<TypedDeconvolutionAlgorithm<deconvolution::PolarizedImageSet<PolCount>>&>(*_cleanAlgorithm).ExecuteMajorIteration(residualSet, modelSet, psfImages, _imgWidth, _imgHeight, reachedMajorThreshold);
	
	psfImage.reset();
	modelSet.Store(*_modelImages, _polarizations, currentChannelIndex);
	residualSet.Store(*_residualImages, _polarizations, currentChannelIndex);
}

template<size_t PolCount>
void Deconvolution::performJoinedPolFreqClean(bool& reachedMajorThreshold, size_t majorIterationNr)
{
	typename JoinedClean<deconvolution::MultiImageSet<deconvolution::PolarizedImageSet<PolCount>>>::ImageSet
		modelSet(_imgWidth*_imgHeight, _summedCount, *_imageAllocator),
		residualSet(_imgWidth*_imgHeight, _summedCount, *_imageAllocator);
	
	std::vector<ImageBufferAllocator<double>::Ptr> psfImagePtrs(_summedCount);
	std::vector<double*> psfImages(_summedCount);
	for(size_t ch=0; ch!=_summedCount; ++ch)
	{
		_imageAllocator->Allocate(_imgWidth*_imgHeight, psfImagePtrs[ch]);
		_psfImages->Load(psfImagePtrs[ch].data(), _psfPolarization, ch, false);
		psfImages[ch] = psfImagePtrs[ch].data();
		
		modelSet.Load(*_modelImages, _polarizations, ch);
		residualSet.Load(*_residualImages, _polarizations, ch);
	}
	static_cast<TypedDeconvolutionAlgorithm<deconvolution::MultiImageSet<deconvolution::PolarizedImageSet<PolCount>>>&>(*_cleanAlgorithm).ExecuteMajorIteration(residualSet, modelSet, psfImages, _imgWidth, _imgHeight, reachedMajorThreshold);
	
	for(size_t ch=0; ch!=_summedCount; ++ch)
	{
		psfImagePtrs[ch].reset();
		modelSet.Store(*_modelImages, _polarizations, ch);
		residualSet.Store(*_residualImages, _polarizations, ch);
	}
}

void Deconvolution::FreeDeconvolutionAlgorithms()
{
	_cleanAlgorithm.reset();
}

void Deconvolution::InitializeDeconvolutionAlgorithm(const ImagingTable& groupTable, PolarizationEnum psfPolarization, class ImageBufferAllocator<double>* imageAllocator, size_t imgWidth, size_t imgHeight, double pixelScaleX, double pixelScaleY, size_t outputChannels, double beamSize, size_t threadCount)
{
	_imageAllocator = imageAllocator;
	_imgWidth = imgWidth;
	_imgHeight = imgHeight;
	_psfPolarization = psfPolarization;
	FreeDeconvolutionAlgorithms();
	
	_summedCount = groupTable.SquaredGroupCount();
	if(_summedCount == 0)
		throw std::runtime_error("Nothing to clean");
	ImagingTable firstSquaredGroup = groupTable.GetSquaredGroup(0);
	_squaredCount = firstSquaredGroup.EntryCount();
	_polarizations.clear();
	for(size_t p=0; p!=_squaredCount; ++p)
	{
		if(_polarizations.count(firstSquaredGroup[p].polarization) != 0)
			throw std::runtime_error("Two equal polarizations were given to deconvolution algorithm within a single olarized group");
		else
			_polarizations.insert(firstSquaredGroup[p].polarization);
	}
	
	if(_useMoreSane)
	{
		_cleanAlgorithm.reset(new MoreSane(_moreSaneLocation, _moreSaneArgs));
	}
	else if(_squaredCount != 1)
	{
		if(_squaredCount != 2 && _squaredCount != 4)
			throw std::runtime_error("Joined polarization cleaning was requested, but can't find a compatible set of 2 or 4 pols to clean");
		bool hasXY = _polarizations.count(Polarization::XY)!=0;
		bool hasYX = _polarizations.count(Polarization::YX)!=0;
		if((hasXY && !hasYX) || (hasYX && !hasXY))
			throw std::runtime_error("Cannot jointly clean polarization XY or YX without cleaning both.");
			
		if(_summedCount != 1)
		{
			if(_multiscale)
			{
				if(_squaredCount == 4)
				{
					_cleanAlgorithm.reset(
					new MultiScaleClean
					<deconvolution::MultiImageSet
					<deconvolution::PolarizedImageSet<4>>>(beamSize, pixelScaleX, pixelScaleY));
				}
				else {
					_cleanAlgorithm.reset(
					new MultiScaleClean
					<deconvolution::MultiImageSet<deconvolution::PolarizedImageSet<2>>>(beamSize, pixelScaleX, pixelScaleY));
				}
			}
			else {
				if(_squaredCount == 4)
					_cleanAlgorithm.reset(new JoinedClean<deconvolution::MultiImageSet<deconvolution::PolarizedImageSet<4>>>());
				else
					_cleanAlgorithm.reset(new JoinedClean<deconvolution::MultiImageSet<deconvolution::PolarizedImageSet<2>>>());
			}
		}
		else {
			if(_multiscale)
			{
				if(_squaredCount == 4)
					_cleanAlgorithm.reset(new MultiScaleClean<deconvolution::PolarizedImageSet<4>>(beamSize, pixelScaleX, pixelScaleY));
				else
					_cleanAlgorithm.reset(new MultiScaleClean<deconvolution::PolarizedImageSet<2>>(beamSize, pixelScaleX, pixelScaleY));
			}
			else
			{
				if(_squaredCount == 4)
					_cleanAlgorithm.reset(new JoinedClean<deconvolution::PolarizedImageSet<4>>());
				else
					_cleanAlgorithm.reset(new JoinedClean<deconvolution::PolarizedImageSet<2>>());
			}
		}
	}
	else { // squaredCount == 1
		if(_summedCount != 1)
		{
			if(_multiscale)
				_cleanAlgorithm.reset(new MultiScaleClean<deconvolution::MultiImageSet<deconvolution::SingleImageSet>>(beamSize, pixelScaleX, pixelScaleY));
			else
				_cleanAlgorithm.reset(new JoinedClean<deconvolution::MultiImageSet<deconvolution::SingleImageSet>>());
		}
		else {
			if(_multiscale)
				_cleanAlgorithm.reset(new MultiScaleClean<deconvolution::SingleImageSet>(beamSize, pixelScaleX, pixelScaleY));
			else
				_cleanAlgorithm.reset(new SimpleClean());
		}
	}
	
	_cleanAlgorithm->SetMaxNIter(_nIter);
	_cleanAlgorithm->SetThreshold(_threshold);
	_cleanAlgorithm->SetSubtractionGain(_gain);
	_cleanAlgorithm->SetStopGain(_mGain);
	_cleanAlgorithm->SetCleanBorderRatio(_cleanBorderRatio);
	_cleanAlgorithm->SetAllowNegativeComponents(_allowNegative);
	_cleanAlgorithm->SetStopOnNegativeComponents(_stopOnNegative);
	_cleanAlgorithm->SetThreadCount(threadCount);
	_cleanAlgorithm->SetMultiscaleScaleBias(_multiscaleScaleBias);
	_cleanAlgorithm->SetMultiscaleThresholdBias(_multiscaleThresholdBias);
	
	if(!_fitsMask.empty())
	{
		if(_cleanMask.empty())
		{
			std::cout << "Reading mask '" << _fitsMask << "'...\n";
			FitsReader maskReader(_fitsMask);
			if(maskReader.ImageWidth() != _imgWidth || maskReader.ImageHeight() != _imgHeight)
				throw std::runtime_error("Specified Fits file mask did not have same dimensions as output image!");
			ao::uvector<float> maskData(_imgWidth*_imgHeight);
			maskReader.Read(maskData.data());
			_cleanMask.assign(_imgWidth*_imgHeight, false);
			for(size_t i=0; i!=_imgWidth*_imgHeight; ++i)
				_cleanMask[i] = maskData[i]!=0.0;
		}
		_cleanAlgorithm->SetCleanMask(_cleanMask.data());
	}
	else if(!_casaMask.empty())
	{
		if(_cleanMask.empty())
		{
			std::cout << "Reading CASA mask '" << _casaMask << "'...\n";
			_cleanMask.assign(_imgWidth*_imgHeight, false);
			CasaMaskReader maskReader(_casaMask);
			if(maskReader.Width() != _imgWidth || maskReader.Height() != _imgHeight)
				throw std::runtime_error("Specified CASA mask did not have same dimensions as output image!");
			maskReader.Read(_cleanMask.data());
		}
		_cleanAlgorithm->SetCleanMask(_cleanMask.data());
	}
}
