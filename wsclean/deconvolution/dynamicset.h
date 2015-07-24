#ifndef DYNAMIC_SET_H
#define DYNAMIC_SET_H

#include "../uvector.h"
#include "../wsclean/imagingtable.h"

#include <vector>
#include <map>

class DynamicSet
{
public:
	DynamicSet(const ImagingTable* table) :
		_images(table->EntryCount()),
		_imagingTable(*table),
		_imageIndexToPSFIndex(table->EntryCount())
	{
		initializeIndices();
	}
	
	DynamicSet(const ImagingTable* table, size_t width, size_t height) :
		_images(table->EntryCount()),
		_imagingTable(*table),
		_imageIndexToPSFIndex(table->EntryCount())
	{
		initializeIndices();
		AllocateImages(width, height);
	}
	
	void AllocateImages(size_t width, size_t height)
	{
		for(std::vector<ao::uvector<double>>::iterator img=_images.begin();
				img!=_images.end(); ++img)
		{
			img->resize(width*height);
		}
	}
	
	void GetIntegrated(ao::uvector<double>& dest, ao::uvector<double>& scratch) const
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
				scratch = _images[imageIndex];
			}
			else {
				for(size_t eIndex = 0; eIndex!=subTable.EntryCount(); ++eIndex)
				{
					const ImagingTableEntry& entry = subTable[eIndex];
					size_t imageIndex = _tableIndexToImageIndex.find(entry.index)->second;
					if(eIndex == 0)
					{
						scratch = _images[0];
						square(scratch);
					}
					else {
						addSquared(scratch, _images[imageIndex]);
					}
				}
				squareRoot(scratch);
			}
			
			if(sqIndex == 0)
				dest = scratch;
			else
				add(dest, scratch);
		}
		if(nImagesAdded > 0.0)
			multiply(dest, 1.0/nImagesAdded);
		else
			dest.assign(dest.size(), 0.0);
	}
	
	DynamicSet& operator=(double val)
	{
		for(size_t i=0; i!=size(); ++i)
			_images[i].assign(_images[i].size(), val);
		return *this;
	}
	
	ao::uvector<double>& operator[](size_t index)
	{
		return _images[index];
	}
	
	const ao::uvector<double>& operator[](size_t index) const
	{
		return _images[index];
	}
	
	size_t size() const { return _images.size(); }
	
	size_t PSFIndex(size_t imageIndex) const { return _imageIndexToPSFIndex[imageIndex]; }
	
	const ImagingTable& Table() const { return _imagingTable; }
	
	DynamicSet* CreateTrimmed(size_t x1, size_t y1, size_t x2, size_t y2, size_t oldWidth) const
	{
		std::unique_ptr<DynamicSet> p(new DynamicSet(&_imagingTable));
		for(size_t i=0; i!=_images.size(); ++i)
		{
			copySmallerPart(_images[i], p->_images[i], x1, y1, x2, y2, oldWidth);
		}
		return p.release();
	}
	
	DynamicSet& operator*=(double factor)
	{
		for(size_t i=0; i!=size(); ++i)
			multiply(_images[i], factor);
		return *this;
	}
	
	DynamicSet& operator+=(const DynamicSet& other)
	{
		for(size_t i=0; i!=size(); ++i)
			add(_images[i], other._images[i]);
		return *this;
	}
	
	void FactorAdd(DynamicSet& rhs, double factor)
	{
		for(size_t i=0; i!=size(); ++i)
			addFactor(_images[i], rhs._images[i], factor);
	}
private:
	static void square(ao::uvector<double>& image)
	{
		for(size_t i=0; i!=image.size(); ++i)
			image[i] *= image[i];
	}
	
	static void squareRoot(ao::uvector<double>& image)
	{
		for(size_t i=0; i!=image.size(); ++i)
			image[i] = sqrt(image[i]);
	}
	
	static void add(ao::uvector<double>& lhs, const ao::uvector<double>& rhs)
	{
		for(size_t i=0; i!=lhs.size(); ++i)
			lhs[i] += rhs[i];
	}
	
	static void addSquared(ao::uvector<double>& lhs, const ao::uvector<double>& rhs)
	{
		for(size_t i=0; i!=lhs.size(); ++i)
			lhs[i] += rhs[i]*rhs[i];
	}
	
	static void addFactor(ao::uvector<double>& lhs, const ao::uvector<double>& rhs, double factor)
	{
		for(size_t i=0; i!=lhs.size(); ++i)
			lhs[i] += rhs[i] * factor;
	}
	
	static void multiply(ao::uvector<double>& image, double fact)
	{
		for(size_t i=0; i!=image.size(); ++i)
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
	
	void copySmallerPart(const ao::uvector<double>& input, ao::uvector<double>& output, size_t x1, size_t y1, size_t x2, size_t y2, size_t oldWidth) const
	{
		size_t newWidth = x2 - x1;
		output.resize(newWidth * (y2-y1));
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
	
	std::vector<ao::uvector<double>> _images;
	const ImagingTable& _imagingTable;
	std::map<size_t, size_t> _tableIndexToImageIndex;
	ao::uvector<size_t> _imageIndexToPSFIndex;
};

#endif
