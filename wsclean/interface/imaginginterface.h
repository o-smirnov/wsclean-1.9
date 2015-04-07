#ifndef IMAGING_INTERFACE_H
#define IMAGING_INTERFACE_H

typedef struct
{
	const char* msPath;
	int imageWidth;
	int imageHeight;
	double pixelScaleX;
	double pixelScaleY;
	const char* extraParameters;
} imaging_parameters;

typedef struct
{
	long dataSize;
	// could also hold the fact that we are dealing with complex doubles.
} imaging_data;

#endif
