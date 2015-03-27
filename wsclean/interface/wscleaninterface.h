#ifndef WSCLEAN_PURIFY_INTERFACE_H
#define WSCLEAN_PURIFY_INTERFACE_H

#include <complex.h>

typedef struct
{
	const char* msPath;
	int imageWidth;
	int imageHeight;
	double pixelScaleX;
	double pixelScaleY;
	const char* extraParameters;
} purify_domain_info;

typedef struct
{
	long data_size;
	// could also hold the fact that we are dealing with complex doubles.
} purify_domain_data_format;

#ifdef __cplusplus
#include <complex>
#define DCOMPLEX std::complex<double>
extern "C" {
#else
#define DCOMPLEX double complex
#endif

/**
 * Initialize WSClean for use as measurement operator in compressed sensing
 * application. This should be called before any other function.
 * 
 * @p userdata should be a pointer to a void pointer which will be set to
 * a structure that WSClean internally uses. The way to pass this pointer
 * is like this:
 * void* userdata;
 * wsclean_initialize(&userdata, ...)
 * ...
 * wsclean_deinitialize(userdata);
 * @p domain_info domain specific information, containing the measurement set.
 * @p data_info will be filled with info describing the data.
 */
void wsclean_initialize(
	void** userData,
	const purify_domain_info* domain_info,
	purify_domain_data_format* data_info
);

/**
 * Initializes the data array
 * @p data An already allocated array of data which will be set to the selected data
 * in the measurement set.
 * @p weights An already allocated array which will be set to the weights.
 */
void wsclean_read(void* userData, DCOMPLEX* data, double* weights);

/**
 * Write the final image out. (maybe nice for later)
 * @p image The image data.
 */
void wsclean_write(void* userData, const double* image);

/**
 * Clean up.
 */
void wsclean_deinitialize(void* userData);

void wsclean_operator_A(
	void* dataIn, void* dataOut,
	void* userData);

void wsclean_operator_At(
	void* dataIn, void* dataOut,
	void* userData);

#ifdef __cplusplus
}
#endif

#endif
