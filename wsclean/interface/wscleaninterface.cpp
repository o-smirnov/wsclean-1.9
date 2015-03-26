#define WSCLEAN_NO_MAIN
#include "../wscleanmain.cpp"
#undef WSCLEAN_NO_MAIN

#include "wscleaninterface.h"
#include "../banddata.h"

#include <string>

#include <ms/MeasurementSets/MeasurementSet.h>

struct WSCleanUserData
{
	std::string msPath;
	int width;
	int height;
	double pixelScaleX;
	double pixelScaleY;
	std::string extraParameters;
	
	std::string dataColumn;
	bool hasAImage, hasAtImage;
};

template<typename T>
std::string str(T i)
{
	std::ostringstream s;
	s << i;
	return s.str();
}

void wsclean_main(const std::vector<std::string>& parms)
{
	char** argv = new char*[parms.size()];
	for(size_t i=0; i!=parms.size(); ++i)
	{
		std::cout << parms[i] << ' ';
		argv[i] = new char[parms[i].size()+1];
		memcpy(argv[i], parms[i].c_str(), parms[i].size()+1);
	}

	wsclean_main(parms.size(), argv);
	
	for(size_t i=0; i!=parms.size(); ++i)
		delete[] argv[i];
	delete[] argv;
}

void wsclean_initialize(
	void** userData,
	const struct purify_domain_info* domain_info,
	struct purify_domain_data_format* data_info
)
{
	WSCleanUserData* wscUserData = new WSCleanUserData();
	wscUserData->msPath = domain_info->msPath;
	wscUserData->width = domain_info->imageWidth;
	wscUserData->height = domain_info->imageHeight;
	wscUserData->pixelScaleX = domain_info->pixelScaleX;
	wscUserData->pixelScaleY = domain_info->pixelScaleY;
	wscUserData->extraParameters = domain_info->extraParameters;
	wscUserData->hasAImage = false;
	wscUserData->hasAtImage = false;
	(*userData) = static_cast<void*>(wscUserData);
	
	// Number of vis is nchannels x selected nrows; calculate both.
	// (Assuming Stokes I polarization for now)
	casa::MeasurementSet ms(wscUserData->msPath);
	casa::ROScalarColumn<int> a1Col(ms, casa::MeasurementSet::columnName(casa::MSMainEnums::ANTENNA1));
	casa::ROScalarColumn<int> a2Col(ms, casa::MeasurementSet::columnName(casa::MSMainEnums::ANTENNA2));
	BandData bandData(ms.spectralWindow());
	size_t nChannel = bandData.ChannelCount();
	size_t selectedRows = 0;
	for(size_t row=0; row!=ms.nrow(); ++row)
	{
		if(a1Col(row) != a2Col(row))
			++selectedRows;
	}
	
	data_info->data_size = selectedRows * nChannel;
	
	bool hasCorrected = ms.tableDesc().isColumn("CORRECTED_DATA");
	if(hasCorrected) {
		std::cout << "First measurement set has corrected data: tasks will be applied on the corrected data column.\n";
		wscUserData->dataColumn = "CORRECTED_DATA";
	} else {
		std::cout << "No corrected data in first measurement set: tasks will be applied on the data column.\n";
		wscUserData->dataColumn = "DATA";
	}
}

void wsclean_deinitialize(void* userData)
{
	WSCleanUserData* wscUserData = static_cast<WSCleanUserData*>(userData);
	delete wscUserData;
}

void wsclean_read(void* userData, DCOMPLEX* data, double* weights)
{
	WSCleanUserData* wscUserData = static_cast<WSCleanUserData*>(userData);
	casa::MeasurementSet ms(wscUserData->msPath);
	BandData bandData(ms.spectralWindow());
	size_t nChannels = bandData.ChannelCount();
	
	casa::ROScalarColumn<int> a1Col(ms, casa::MeasurementSet::columnName(casa::MSMainEnums::ANTENNA1));
	casa::ROScalarColumn<int> a2Col(ms, casa::MeasurementSet::columnName(casa::MSMainEnums::ANTENNA2));
	
	casa::ROArrayColumn<casa::Complex> dataCol(ms, wscUserData->dataColumn);
	casa::ROArrayColumn<float> weightCol(ms, casa::MeasurementSet::columnName(casa::MSMainEnums::WEIGHT_SPECTRUM));
	casa::ROArrayColumn<bool> flagCol(ms, casa::MeasurementSet::columnName(casa::MSMainEnums::FLAG));
	
	DCOMPLEX* dataPtr = data;
	double* weightPtr = weights;
	casa::IPosition shape = dataCol.shape(0);
	size_t polarizationCount = shape[0];
	casa::Array<casa::Complex> dataArr(shape);
	casa::Array<bool> flagArr(shape);
	casa::Array<float> weightArr(shape);
	for(size_t row=0; row!=ms.nrow(); ++row)
	{
		if(a1Col(row) != a2Col(row))
		{
			dataCol.get(row, dataArr);
			flagCol.get(row, flagArr);
			weightCol.get(row, weightArr);
			
			casa::Array<casa::Complex>::const_contiter di = dataArr.cbegin();
			casa::Array<bool>::const_contiter fi = flagArr.cbegin();
			casa::Array<float>::const_contiter wi = weightArr.cbegin();
			
			for(size_t ch=0; ch!=nChannels; ++ch)
			{
				// TODO this only works for XX/YY and LL/RR pol, but not if
				// MS contains IQUV
				std::complex<double> val = 0.5*(std::complex<double>(*di) + std::complex<double>(*(di+polarizationCount-1)));
				double weight = 0.5*(double(*wi) + (*wi+polarizationCount-1));
				bool flag = *fi || *(fi+polarizationCount-1);
				
				dataPtr[ch] = val;
				weightPtr[ch] = flag ? 0.0 : weight;
				
				di += polarizationCount;
				fi += polarizationCount;
				wi += polarizationCount;
			}
			dataPtr += nChannels;
			weightPtr += nChannels;
		}
	}
}

void wsclean_write(void* userData, const double* image)
{
	WSCleanUserData* wscUserData = static_cast<WSCleanUserData*>(userData);
	FitsWriter writer;
	writer.SetImageDimensions(wscUserData->width, wscUserData->height, wscUserData->pixelScaleX, wscUserData->pixelScaleY);
	if(wscUserData->hasAtImage)
	{
		FitsReader reader("tmp-operator-At-image.fits");
		writer = FitsWriter(reader);
	}
	writer.Write("purify-wsclean-model.fits", image);
}

void getCommandLine(std::vector<std::string>& commandline, const WSCleanUserData& userData)
{
	commandline.push_back("wsclean");
	commandline.push_back("-size");
	commandline.push_back(str(userData.width));
	commandline.push_back(str(userData.height));
	commandline.push_back("-scale");
	commandline.push_back(Angle::ToNiceString(userData.pixelScaleX));
	if(!userData.extraParameters.empty())
	{
		size_t pos = 0;
		size_t nextPos = userData.extraParameters.find(' ', 0);
		while(nextPos!=std::string::npos)
		{
			commandline.push_back(userData.extraParameters.substr(pos, nextPos-pos));
			pos = nextPos+1;
			nextPos = userData.extraParameters.find(' ', pos);
		}
		commandline.push_back(userData.extraParameters.substr(pos));
	}
	if(userData.pixelScaleX != userData.pixelScaleY)
		throw std::runtime_error("pixelscaleX should be equal to pixelscaleY for WSClean");
}

// Go from image to visibilities
// dataIn :  double[] of size width*height
// dataOut : complex double[] of size nvis: nchannels x nbaselines x ntimesteps
void wsclean_operator_A(
	void* dataIn, void* dataOut,
	void* userData)
{
	WSCleanUserData* wscUserData = static_cast<WSCleanUserData*>(userData);
	
	// Write dataIn to a fits file
	FitsWriter writer;
	writer.SetImageDimensions(wscUserData->width, wscUserData->height, wscUserData->pixelScaleX, wscUserData->pixelScaleY);
	writer.Write("tmp-operator-A-model.fits", static_cast<double*>(dataIn));
	wscUserData->hasAImage = true;
	
	// Run WSClean -predict (creates/fills new column MODEL_DATA)
	std::vector<std::string> commandline;
	getCommandLine(commandline, *wscUserData);
	commandline.push_back("-name");
	commandline.push_back("tmp-operator-A");
	commandline.push_back("-predict");
	commandline.push_back(wscUserData->msPath);
	wsclean_main(commandline);
	
	// Read MODEL_DATA into dataOut
	casa::MeasurementSet ms(wscUserData->msPath);
	BandData bandData(ms.spectralWindow());
	size_t nChannels = bandData.ChannelCount();
	
	casa::ROScalarColumn<int> a1Col(ms, casa::MeasurementSet::columnName(casa::MSMainEnums::ANTENNA1));
	casa::ROScalarColumn<int> a2Col(ms, casa::MeasurementSet::columnName(casa::MSMainEnums::ANTENNA2));
	
	casa::ROArrayColumn<casa::Complex> dataCol(ms,  casa::MeasurementSet::columnName(casa::MSMainEnums::MODEL_DATA));
	
	DCOMPLEX* dataPtr = (DCOMPLEX*) dataOut;
	casa::IPosition shape = dataCol.shape(0);
	size_t polarizationCount = shape[0];
	casa::Array<casa::Complex> dataArr(shape);
	for(size_t row=0; row!=ms.nrow(); ++row)
	{
		if(a1Col(row) != a2Col(row))
		{
			dataCol.get(row, dataArr);
			casa::Array<casa::Complex>::contiter di = dataArr.cbegin();
			
			for(size_t ch=0; ch!=nChannels; ++ch)
			{
				dataPtr[ch] = 0.5 * (std::complex<double>(*di) + std::complex<double>(*(di+polarizationCount-1)));
				di += polarizationCount;
			}
			
			dataPtr += nChannels;
		}
	}
}

// Go from visibilities to image
void wsclean_operator_At(
	void* dataIn, void* dataOut,
	void* userData)
{
	// Write dataIn to the MODEL_DATA column
	WSCleanUserData* wscUserData = static_cast<WSCleanUserData*>(userData);
	casa::MeasurementSet ms(wscUserData->msPath, casa::Table::Update);
	BandData bandData(ms.spectralWindow());
	size_t nChannels = bandData.ChannelCount();
	
	casa::ROScalarColumn<int> a1Col(ms, casa::MeasurementSet::columnName(casa::MSMainEnums::ANTENNA1));
	casa::ROScalarColumn<int> a2Col(ms, casa::MeasurementSet::columnName(casa::MSMainEnums::ANTENNA2));
	
	casa::ArrayColumn<casa::Complex> dataCol(ms,  casa::MeasurementSet::columnName(casa::MSMainEnums::MODEL_DATA));
	
	DCOMPLEX* dataPtr = (DCOMPLEX*) dataIn;
	casa::IPosition shape = dataCol.shape(0);
	size_t polarizationCount = shape[0];
	casa::Array<casa::Complex> dataArr(shape);
	for(size_t row=0; row!=ms.nrow(); ++row)
	{
		if(a1Col(row) != a2Col(row))
		{
			dataCol.get(row, dataArr);
			casa::Array<casa::Complex>::contiter di = dataArr.cbegin();
			
			for(size_t ch=0; ch!=nChannels; ++ch)
			{
				*di = std::complex<float>(dataPtr[ch]);
				*(di+polarizationCount-1) = std::complex<float>(dataPtr[ch]);
				di += polarizationCount;
			}
			
			dataCol.put(row, dataArr);
			dataPtr += nChannels;
		}
	}
	
	// Run WSClean to create dirty image
	std::vector<std::string> commandline;
	getCommandLine(commandline, *wscUserData);
	commandline.push_back("-name");
	commandline.push_back("tmp-operator-At");
	commandline.push_back(wscUserData->msPath);
	wsclean_main(commandline);
	wscUserData->hasAtImage = true;
	
	// Read dirty image and store in dataOut
	FitsReader reader("tmp-operator-At-image.fits");
	reader.Read(static_cast<double*>(dataOut));
}
