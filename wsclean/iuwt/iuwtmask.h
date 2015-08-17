#ifndef IUWT_MASK_H
#define IUWT_MASK_H

#include <vector>
#include <string>
#include <sstream>
#include <limits>

#include "../uvector.h"

class IUWTMask
{
public:
	IUWTMask(int scaleCount, size_t width, size_t height) :
		_masks(scaleCount), _width(width), _height(height)
	{
		for(int i=0; i!=scaleCount; ++i)
			_masks[i].assign(width*height, false);
	}
	
	IUWTMask* CreateTrimmed(size_t x1, size_t y1, size_t x2, size_t y2) const
	{
		std::unique_ptr<IUWTMask> p(new IUWTMask(ScaleCount(), x2-x1, y2-y1));
		
		for(size_t i=0; i!=_masks.size(); ++i)
			copySmallerPart(_masks[i], p->_masks[i], x1, y1, x2, y2);
		
		return p.release();
	}
	
	ao::uvector<bool>& operator[](int scale)
	{
		return _masks[scale];
	}
	const ao::uvector<bool>& operator[](int scale) const
	{
		return _masks[scale];
	}
	std::string Summary() const
	{
		std::ostringstream str;
		str << "IUWTMask with " << _masks.size() << " scale masks.\n";
		for(size_t i=0; i!=_masks.size(); ++i)
		{
			size_t count = 0;
			for(size_t j=0; j!=_masks[i].size(); ++j)
			{
				if(_masks[i][j]) ++count;
			}
			str << "Scale " << i << ": " << count << '\n';
		}
		return str.str();
	}
	std::string Summary(const class IUWTDecomposition& iuwt) const;
	
	size_t HighestActiveScale() const
	{
		for(int m=_masks.size()-1; m!=-1; --m)
		{
			for(size_t i=0; i!=_masks[m].size(); ++i)
			{
				if(_masks[m][i]) return m;
			}
		}
		return 0; // avoid compiler warning
	}
	
	void TransferScaleUpwards(size_t fromScale)
	{
		if(fromScale>0)
		{
			size_t toScale = fromScale-1;
			for(size_t i=0; i!=_masks[fromScale].size(); ++i)
			{
				if(_masks[fromScale][i])
					_masks[toScale][i] = true;
			}
		}
	}
	size_t ScaleCount() const { return _masks.size(); }
private:
	std::vector<ao::uvector<bool>> _masks;
	size_t _width, _height;

		void copySmallerPart(const ao::uvector<bool>& input, ao::uvector<bool>& output, size_t x1, size_t y1, size_t x2, size_t y2) const
	{
		size_t newWidth = x2 - x1;
		output.resize(newWidth * (y2-y1));
		for(size_t y=y1; y!=y2; ++y)
		{
			const bool* oldPtr = &input[y*_width];
			bool* newPtr = &output[(y-y1)*newWidth];
			for(size_t x=x1; x!=x2; ++x)
				newPtr[x - x1] = oldPtr[x];
		}
	}
};

#endif
