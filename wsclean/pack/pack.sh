if [[ "$1" == "" ]] ; then
		echo Syntax: ./pack.sh \<version\>
else
		VERSION="$1"
		WORKDIR=`pwd`
		rm -rf /tmp/wsclean /tmp/wsclean-${VERSION}
		mkdir /tmp/wsclean
		mkdir /tmp/wsclean/aocommon
		mkdir /tmp/wsclean/beam
		mkdir /tmp/wsclean/cleanalgorithms
		mkdir /tmp/wsclean/msprovider
		mkdir /tmp/wsclean/parser
		cd ..
		cp -v CMakeLists.txt areaset.* banddata.* beamevaluator.* buffered_lane.* cachedimageset.* fftconvolver.* fftresampler.* fitsreader.* fitswriter.* imagecoordinates.* imagebufferallocator.* imageweights.* inversionalgorithm.* lane.* layeredimager.* multibanddata.* nlplfitter.* matrix2x2.* model.* modelrenderer.* modelsource.* msselection.* polarizationenum.* progressbar.* radeccoord.* sourcesdf.* sourcesdfwithsamples.* spectralenergydistribution.* stopwatch.* tilebeam.* uvector.* uvwdistribution.* weightmode.* wsclean.* wscleanmain.cpp wsinversion.* /tmp/wsclean/
		cp -v aocommon/lane.h aocommon/lane_03.h aocommon/lane_11.h aocommon/uvector.h /tmp/wsclean/aocommon/
		cp -v beam/*.{h,cpp} /tmp/wsclean/beam
		cp -v cleanalgorithms/*.{h,cpp} /tmp/wsclean/cleanalgorithms/
		cp -v msprovider/*.{h,cpp} /tmp/wsclean/msprovider
		cp -v parser/*.h /tmp/wsclean/parser
		cd /tmp
		mv wsclean wsclean-${VERSION}
		tar -cjvf ${WORKDIR}/wsclean-${VERSION}.tar.bz2 wsclean-${VERSION}/
		rm -rf /tmp/wsclean-${VERSION}
		tar -xjvf ${WORKDIR}/wsclean-${VERSION}.tar.bz2
		cd /tmp/wsclean-${VERSION}
		mkdir build
		cd build
		cmake ../
		make -j 12
fi
