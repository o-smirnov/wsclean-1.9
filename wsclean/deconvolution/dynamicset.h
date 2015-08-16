#ifndef DYNAMIC_SET_H
#define DYNAMIC_SET_H

#include "../uvector.h"
#include "../wsclean/imagingtable.h"
#include "../wsclean/imagebufferallocator.h"

#include <vector>
#include <map>

class DynamicSet
{
public:
	DynamicSet(const ImagingTable* table, ImageBufferAllocator& allocator) :
		_images(table->EntryCount()),
		_imageSize(0),
		_imagingTable(*table),
		_imageIndexToPSFIndex(table->EntryCount()),
		_allocator(allocator)
	{
		initializeIndices();
	}
	
	DynamicSet(const ImagingTable* table, ImageBufferAllocator& allocator, size_t width, size_t height) :
		_images(table->EntryCount()),
		_imageSize(width*height),
		_imagingTable(*table),
		_imageIndexToPSFIndex(table->EntryCount()),
		_allocator(allocator)
	{
		initializeIndices();
		AllocateImages();
	}
	
	void AllocateImages()
	{
		for(std::vector<ImageBufferAllocator::Ptr>::iterator img=_images.begin();
				img!=_images.end(); ++img)
		{
			_allocator.Allocate(_imageSize, *img);
		}
	}
	
	void AllocateImages(size_t width, size_t height)
	{
		_imageSize = width*height;
		for(std::vector<ImageBufferAllocator::Ptr>::iterator img=_images.begin();
				img!=_images.end(); ++img)
		{
			_allocator.Allocate(_imageSize, *img);
		}
	}
	
	bool IsAllocated() const
	{
		return _imageSize!=0;
	}
	
	ImageBufferAllocator& Allocator() const
	{
		return _allocator;
	}
	
	void GetIntegrated(double* dest, double* scratch) const
	{
		double nImagesAdded = 0.0;
		for(size_t sqIndex = 0; sqIndex!=_imagingTable.SquaredGroupCount(); ++sqIndex)
		{
			ImagingTable subTable = _imagingTable.GetSquaredGroup(sqIndex);
			nImagesAdded += sqrt(subTable.EntryCount());
			if(subTable.EntryCount() == 1)
			{
				const ImagingTableEntry& entry = subTable[0];
				size_t imageIndex = _tableIndexToImageIndex.find(entry.index)->second;
				assign(scratch, _images[imageIndex]);
			}
			else {
				for(size_t eIndex = 0; eIndex!=subTable.EntryCount(); ++eIndex)
				{
					const ImagingTableEntry& entry = subTable[eIndex];
					size_t imageIndex = _tableIndexToImageIndex.find(entry.index)->second;
					if(eIndex == 0)
					{
						assign(scratch, _images[0]);
						square(scratch);
					}
					else {
						addSquared(scratch, _images[imageIndex].data());
					}
				}
				squareRoot(scratch);
			}
			
			if(sqIndex == 0)
				assign(dest, scratch);
			else
				add(dest, scratch);
		}
		if(nImagesAdded > 0.0)
			multiply(dest, 1.0/nImagesAdded);
		else
			assign(dest, 0.0);
	}
	
	void GetIntegratedPSF(double* dest, const ao::uvector<const double*>& psfs)
	{
		memcpy(dest, psfs[0], sizeof(double) * _imageSize);
		for(size_t sqIndex = 1; sqIndex!=_imagingTable.SquaredGroupCount(); ++sqIndex)
		{
			add(dest, psfs[sqIndex]);
		}
		multiply(dest, 1.0/double(_imagingTable.SquaredGroupCount()));
	}
	
	size_t PSFCount() const { return _imagingTable.SquaredGroupCount(); }
	
	DynamicSet& operator=(double val)
	{
		for(size_t i=0; i!=size(); ++i)
			assign(_images[i].data(), val);
		return *this;
	}
	
	double* operator[](size_t index)
	{
		return _images[index].data();
	}
	
	const double* operator[](size_t index) const
	{
		return _images[index].data();
	}
	
	size_t size() const { return _images.size(); }
	
	size_t PSFIndex(size_t imageIndex) const { return _imageIndexToPSFIndex[imageIndex]; }
	
	const ImagingTable& Table() const { return _imagingTable; }
	
	DynamicSet* CreateTrimmed(size_t x1, size_t y1, size_t x2, size_t y2, size_t oldWidth) const
	{
		std::unique_ptr<DynamicSet> p(new DynamicSet(&_imagingTable, _allocator, x2-x1, y2-y1));
		for(size_t i=0; i!=_images.size(); ++i)
		{
			copySmallerPart(_images[i].data(), p->_images[i].data(), x1, y1, x2, y2, oldWidth);
		}
		return p.release();
	}
	
	DynamicSet& operator*=(double factor)
	{
		for(size_t i=0; i!=size(); ++i)
			multiply(_images[i].data(), factor);
		return *this;
	}
	
	DynamicSet& operator+=(const DynamicSet& other)
	{
		for(size_t i=0; i!=size(); ++i)
			add(_images[i].data(), other._images[i].data());
		return *this;
	}
	
	void FactorAdd(DynamicSet& rhs, double factor)
	{
		for(size_t i=0; i!=size(); ++i)
			addFactor(_images[i].data(), rhs._images[i].data(), factor);
	}
private:
	void assign(double* lhs, const double* rhs) const
	{
		memcpy(lhs, rhs, sizeof(double) * _imageSize);
	}
	
	void assign(double* lhs, const ImageBufferAllocator::Ptr& rhs) const
	{
		memcpy(lhs, rhs.data(), sizeof(double) * _imageSize);
	}
	
	void assign(double* image, double value) const
	{
		for(size_t i=0; i!=_imageSize; ++i)
			image[i] = value;
	}
	
	void add(double* lhs, const double* rhs) const
	{
		for(size_t i=0; i!=_imageSize; ++i)
			lhs[i] += rhs[i];
	}
	
	void square(double* image) const
	{
		for(size_t i=0; i!=_imageSize; ++i)
			image[i] *= image[i];
	}
	
	void squareRoot(double* image) const
	{
		for(size_t i=0; i!=_imageSize; ++i)
			image[i] = sqrt(image[i]);
	}
	
	void addSquared(double* lhs, const double* rhs) const
	{
		for(size_t i=0; i!=_imageSize; ++i)
			lhs[i] += rhs[i]*rhs[i];
	}
	
	void addFactor(double* lhs, const double* rhs, double factor) const
	{
		for(size_t i=0; i!=_imageSize; ++i)
			lhs[i] += rhs[i] * factor;
	}
	
	void multiply(double* image, double fact) const
	{
		for(size_t i=0; i!=_imageSize; ++i)
			image[i] *= fact;
	}
	
	void initializeIndices()
	{
		for(size_t i=0; i!=_imagingTable.EntryCount(); ++i)
		{
			_tableIndexToImageIndex.insert(
				std::make_pair(_imagingTable[i].index, i));
		}
		for(size_t sqIndex = 0; sqIndex!=_imagingTable.SquaredGroupCount(); ++sqIndex)
		{
			ImagingTable subTable = _imagingTable.GetSquaredGroup(sqIndex);
			for(size_t eIndex = 0; eIndex!=subTable.EntryCount(); ++eIndex)
			{
				const ImagingTableEntry& entry = subTable[eIndex];
				size_t imageIndex = _tableIndexToImageIndex.find(entry.index)->second;
				_imageIndexToPSFIndex[imageIndex] = sqIndex;
			}
		}
	}
	
	void copySmallerPart(const double* input, double* output, size_t x1, size_t y1, size_t x2, size_t y2, size_t oldWidth) const
	{
		size_t newWidth = x2 - x1;
		for(size_t y=y1; y!=y2; ++y)
		{
			const double* oldPtr = &input[y*oldWidth];
			double* newPtr = &output[(y-y1)*newWidth];
			for(size_t x=x1; x!=x2; ++x)
			{
				newPtr[x - x1] = oldPtr[x];
			}
		}
	}
	
	std::vector<ImageBufferAllocator::Ptr> _images;
	size_t _imageSize;
	const ImagingTable& _imagingTable;
	std::map<size_t, size_t> _tableIndexToImageIndex;
	ao::uvector<size_t> _imageIndexToPSFIndex;
	ImageBufferAllocator& _allocator;
};

#endif
