#include "../binneduvgrid.h"
#include "../imagebufferallocator.h"

int main(int argc, char* argv[])
{
	ImageBufferAllocator allocator;
	size_t width = 1000, height = 1000;
	BinnedUVGrid grid(width, height, 1.0 / width, 1.0 / height, 1, &allocator);
	grid.PrepareWLayers(11, 1e9, -1, 1);
	
	for(size_t y=0; y!=height; ++y)
	{
		for(size_t x=0; x!=width/2; ++x)
		{
			grid.GridSample(std::complex<float>(1.0, 0.0), 1.0, int(x)-int(width/2), int(y)-int(height/2), 0.0);
		}
	}
	
	/*grid.GridSample(std::complex<float>(1.0, 0.0), 1.0, 0.0, 0.0, 0.0);
	grid.GridSample(std::complex<float>(7.0, 1.0), 2.0, -15.0, 15.0, 0.0);
	grid.GridSample(std::complex<float>(1.0, 4.0), 1.0, -15.0, 15.0, 0.0);
	
	grid.GridSample(std::complex<float>(10.0, 0.0), 1.0, 30.4, 30.0, 0.0);
	
	grid.GridSample(std::complex<float>(10.0, 0.0), 1.0, 46.0, 46.0, 0.0);
	grid.GridSample(std::complex<float>(10.0, 0.0), 1.0, 46.0, 46.2, 0.0);*/
	
	grid.FinishGridding();
	
	std::cout << "Non-zero bins: " << grid.NonZeroBinCount() << '\n';
	
	double u, v, w, weight;
	std::complex<double> val;
	size_t count = 0;
	while(grid.GetNextNonZeroBin(u, v, w, val, weight))// && count < 100)
	{
		std::cout << "uvw = (" << u << ',' << v << ',' << w << "), val=" << val << ", weight=" << weight << '\n';
		++count;
	}
}
