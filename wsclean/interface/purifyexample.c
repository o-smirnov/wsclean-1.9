#include "wscleaninterface.h"

#include "math.h"

#include <stdlib.h>
#include <stdio.h>

int main(int argc, char* argv[])
{
	if(argc < 2)
	{
		printf("Specify measurement set: purifyexample <ms>\n");
	}
	else {
		purify_domain_info dinfo;
		purify_domain_data_format format;

		void* userdata;
		dinfo.msPath=argv[1];
		dinfo.imageWidth = 512;
		dinfo.imageHeight = 512;
		dinfo.pixelScaleX = 2.0 * M_PI/(180.0*60.0); // 2 amin
		dinfo.pixelScaleY = 2.0 * M_PI/(180.0*60.0);
		dinfo.extraParameters="-weight natural";
		
		wsclean_initialize(&userdata, &dinfo, &format);
		
		complex double* mydata = (complex double*) malloc(format.data_size * sizeof(complex double));
		complex double* emptydata = (complex double*) malloc(format.data_size * sizeof(complex double));
		double* myweights = (double*) malloc(format.data_size * sizeof(double));
		double* myimage = (double*) malloc(dinfo.imageWidth*dinfo.imageHeight * sizeof(double));
		
		wsclean_operator_At(emptydata, myimage, userdata);
		wsclean_operator_A(myimage, mydata, userdata);
		
		wsclean_read(userdata, mydata, myweights);
		wsclean_operator_At(mydata, myimage, userdata);
		wsclean_operator_A(myimage, mydata, userdata);
		
		size_t i;
		for(i=0; i!=format.data_size; ++i)
			mydata[i] = 1.0;
		
		wsclean_operator_At(mydata, myimage, userdata);
		wsclean_operator_A(myimage, mydata, userdata);
		
		// ...do the purification magic...
		
		wsclean_write(userdata, myimage);
		
		free(myimage);
		free(mydata);
		
		wsclean_deinitialize(userdata);
	}
}
