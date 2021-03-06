#ifndef WSCLEAN_IMAGING_TABLE_H
#define WSCLEAN_IMAGING_TABLE_H

#include <string>
#include <vector>

#include "../uvector.h"
#include "../polarizationenum.h"

#include "../msproviders/partitionedms.h"

struct ImagingTableEntry
{
public:
	struct MSBandInfo
	{
		size_t bandIndex;
		size_t partIndex;
		size_t channelRangeStart, channelRangeEnd;
	};
	
	struct MSInfo
	{
		std::vector<MSBandInfo> bands;
	};
	
	ImagingTableEntry();
	
	size_t index;
	
	/**
	 * Note that mses might have overlapping frequencies.
	 */
	double lowestFrequency, highestFrequency;
	double minBandFrequency, maxBandFrequency;
	
	PolarizationEnum polarization;
	
	size_t outputChannelIndex;
	
	size_t outputTimestepIndex;
	
	std::vector<MSInfo> msData;
	
	/**
	 * The group of entries with equal squaredDeconvolutionIndex
	 * should be 'joinedly' deconvolved by adding their squared powers
	 * together. Normally, all the polarizations from a single
	 * (output)channel / timestep form such a group.
	 */
	size_t squaredDeconvolutionIndex;
	
	/**
	 * Entries with equal joinedGroupIndex are joinedly deconvolved.
	 * Such a group of entries can be further split up in 'squared'
	 * deconvolution groups.
	 */
	size_t joinedGroupIndex;
	
	/**
	 * A normal inversion results in '1' image. However, an XY
	 * imaging run results in 2 (real and imaginary), while an
	 * YX imaging run results in 0, as it is added to XY.
	 */
	size_t imageCount;
	
	std::string tmpFilePrefix;
	
	std::string ToString();
};

class ImagingTable
{
public:
	size_t IndependentGroupCount() const
	{
		return _independentGroupLookup.size();
	}
	
	ImagingTable GetIndependentGroup(size_t index) const;
	
	size_t SquaredGroupCount() const
	{
		return _squaredGroupLookup.size();
	}
	
	ImagingTable GetSquaredGroup(size_t index) const;
	
	size_t EntryCount() const
	{
		return _entries.size();
	}
	
	ImagingTableEntry& operator[](size_t index)
	{
		return _entries[index];
	}
	const ImagingTableEntry& operator[](size_t index) const
	{
		return _entries[index];
	}
	
	size_t ImageCount() const
	{
		return _imageLookup.size();
	}
	
	void GetImageInfo(size_t imageIndex, bool& isImaginary);
	
	void Clear() { _entries.clear(); }
	
	ImagingTableEntry& AddEntry()
	{
		_entries.push_back(ImagingTableEntry());
		return _entries.back();
	}
	
	void Update()
	{
		updateIndependentGroupLookup();
		updateSquaredGroupLookup();
		updateImageLookup();
	}
	
	void Print();
	
	ImagingTableEntry& Front() { return _entries.front(); }
	const ImagingTableEntry& Front() const { return _entries.front(); }
	
private:
	void printIndependentGroup(bool isFinal);
	void updateIndependentGroupLookup();
	void updateSquaredGroupLookup();
	void updateImageLookup();
	
	std::vector<ImagingTableEntry> _entries;
	
	std::vector<std::vector<ImagingTableEntry*>> _independentGroupLookup;
	std::vector<std::vector<ImagingTableEntry*>> _squaredGroupLookup;
	std::vector<std::pair<ImagingTableEntry*,bool>> _imageLookup;
};

#endif
