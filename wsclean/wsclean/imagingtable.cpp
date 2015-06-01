#include "imagingtable.h"
#include <map>

ImagingTableEntry::ImagingTableEntry() :
	lowestFrequency(0.0),
	highestFrequency(0.0),
	polarization(Polarization::StokesI),
	outputChannelIndex(0),
	outputTimestepIndex(0),
	msData(),
	squaredDeconvolutionIndex(0),
	joinedGroupIndex(0),
	imageCount(0),
	tmpFilePrefix()
{
}


ImagingTable ImagingTable::GetIndependentGroup(size_t index) const
{
	ImagingTable table;
	const std::vector<ImagingTableEntry*>& entries = _independentGroupLookup[index];
	for(std::vector<ImagingTableEntry*>::const_iterator e=entries.begin();
			e!=entries.end(); ++e)
	{
		table._entries.push_back(**e);
	}
	table.Update();
	return table;
}

ImagingTable ImagingTable::GetSquaredGroup(size_t index) const
{
	ImagingTable table;
	const std::vector<ImagingTableEntry*>& entries = _squaredGroupLookup[index];
	for(std::vector<ImagingTableEntry*>::const_iterator e=entries.begin();
			e!=entries.end(); ++e)
	{
		table._entries.push_back(**e);
	}
	table.Update();
	return table;
}

void ImagingTable::Print()
{
	std::cout << "=== IMAGING TABLE ===\n"
		"       # Pol Ch JG ²G Freq(MHz)\n";
	for(size_t i=0; i!=IndependentGroupCount(); ++i)
	{
		bool isLastGroup = ((i+1) == IndependentGroupCount());
		std::cout << "| Independent group:\n";
		GetIndependentGroup(i).printIndependentGroup(isLastGroup);
		if(!isLastGroup)
			std::cout << "|\n";
	}
}

void ImagingTable::printIndependentGroup(bool isFinal)
{
	for(size_t i=0; i!=SquaredGroupCount(); ++i)
	{
		bool isSecondFinal = (i+1)==SquaredGroupCount();
		ImagingTable table = GetSquaredGroup(i);
		for(size_t j=0; j!=table._entries.size(); ++j)
		{
			if(i == 0 && j==0)
				std::cout << "+-";
			else if(isFinal)
				std::cout << "  ";
			else
				std::cout << "| ";
			
			if(j == 0)
				std::cout << "+-";
			else if(isSecondFinal)
				std::cout << "  ";
			else
				std::cout << "| ";
			std::cout << "J-";
			std::cout << table._entries[j].ToString() << '\n';
		}
		if(isFinal && isSecondFinal)
			std::cout << '\n';
		else if(isFinal)
			std::cout << "  |\n";
		else if(isSecondFinal)
			std::cout << "|\n";
		else
			std::cout << "| |\n";
	}
}

string ImagingTableEntry::ToString()
{
	std::ostringstream str;
	if(index < 10) str << ' ';
	str << index << ' ';
	std::string polStr = Polarization::TypeToShortString(polarization);
	if(polStr.size() < 2) str << ' ';
	str << polStr << "  ";
	if(outputChannelIndex < 10) str << ' ';
	str << outputChannelIndex
		<< "  " << joinedGroupIndex
		<< "  " << squaredDeconvolutionIndex
		<< "  " << round(lowestFrequency*1e-6)
		<< "-" << round(highestFrequency*1e-6);
	
	return str.str();
}

void ImagingTable::updateIndependentGroupLookup()
{
	_independentGroupLookup.clear();
	std::map<size_t, size_t> indexToLookupSet;
	for(std::vector<ImagingTableEntry>::iterator e=_entries.begin();
			e!=_entries.end(); ++e)
	{
		size_t groupIndex = e->joinedGroupIndex;
		if(indexToLookupSet.count(groupIndex) == 0)
		{
			indexToLookupSet.insert(std::make_pair(groupIndex, _independentGroupLookup.size()));
			_independentGroupLookup.push_back(std::vector<ImagingTableEntry*>(1, &*e));
		}
		else {
			_independentGroupLookup[indexToLookupSet[groupIndex]].push_back(&*e);
		}
	}
}

void ImagingTable::updateSquaredGroupLookup()
{
	_squaredGroupLookup.clear();
	std::map<size_t, size_t> indexToLookupSet;
	for(std::vector<ImagingTableEntry>::iterator e=_entries.begin();
			e!=_entries.end(); ++e)
	{
		size_t groupIndex = e->squaredDeconvolutionIndex;
		if(indexToLookupSet.count(groupIndex) == 0)
		{
			indexToLookupSet.insert(std::make_pair(groupIndex, _squaredGroupLookup.size()));
			_squaredGroupLookup.push_back(std::vector<ImagingTableEntry*>(1, &*e));
		}
		else {
			_squaredGroupLookup[indexToLookupSet[groupIndex]].push_back(&*e);
		}
	}
}

void ImagingTable::updateImageLookup()
{
	_imageLookup.clear();
	for(std::vector<ImagingTableEntry>::iterator e=_entries.begin();
			e!=_entries.end(); ++e)
	{
		bool isImaginary = false; 
		for(size_t i=0; i!=e->imageCount; ++i)
		{
			_imageLookup.push_back(std::make_pair(&*e, isImaginary));
			isImaginary = true;
		}
	}
}
