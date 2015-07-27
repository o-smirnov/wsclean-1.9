#include "casamaskreader.h"

#include <casacore/tables/Tables/Table.h>
#include <casacore/tables/Tables/ArrayColumn.h>

#include "uvector.h"

CasaMaskReader::CasaMaskReader(const std::string& path) : _path(path)
{
	casacore::Table table(path);
	casacore::ROArrayColumn<float> mapCol(table, "map");
	casacore::IPosition shape = mapCol.shape(0);
	_width = shape(0);
	_height = shape(1);
	_nPolarizations = shape(2);
	_nChannels = shape(3);
}

void CasaMaskReader::Read(bool* mask)
{
	casacore::Table table(_path);
	casacore::ROArrayColumn<float> mapCol(table, "map");
	casacore::Array<float> data(mapCol.get(0));
	for(size_t i=0; i!=_width*_height; ++i)
		mask[i] = false;
	casacore::Array<float>::contiter iter = data.cbegin();
	bool* maskPtr = mask;
	for(size_t j=0; j!=_nChannels*_nPolarizations; ++j)
	{
		for(size_t y=0; y!=_height; ++y)
		{
			for(size_t x=0; x!=_width; ++x)
			{
				*maskPtr = *maskPtr || (*iter!=0.0);
				++iter;
				++maskPtr;
			}
		}
	}
}
