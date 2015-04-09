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
		mkdir /tmp/wsclean/interface
		mkdir /tmp/wsclean/lofar
		mkdir /tmp/wsclean/model
		mkdir /tmp/wsclean/msproviders
		mkdir /tmp/wsclean/parser
		mkdir /tmp/wsclean/pywsclean
		cd ..
		cp -v CMakeLists.txt angle.h areaset.* banddata.* beamevaluator.* buffered_lane.* cachedimageset.* casamaskreader.* dftpredictionalgorithm.* fftconvolver.* fftresampler.* fitsiochecker.* fitsreader.* fitswriter.* gaussianfitter.h imagecoordinates.* imagebufferallocator.* imageweights.* inversionalgorithm.* lane.* layeredimager.* multibanddata.* nlplfitter.* matrix2x2.* modelrenderer.* msselection.* polarizationenum.* progressbar.* radeccoord.* sourcesdf.* sourcesdfwithsamples.* spectralenergydistribution.* stopwatch.* tilebeam.* uvector.* uvwdistribution.* weightmode.* wsclean.* wscleanmain.cpp wsinversion.* wscversion.h /tmp/wsclean/
		cp -v aocommon/lane.h aocommon/lane_03.h aocommon/lane_11.h aocommon/uvector.h aocommon/uvector_03.h /tmp/wsclean/aocommon/
		cp -v beam/*.{h,cpp} /tmp/wsclean/beam
		cp -v cleanalgorithms/*.{h,cpp} /tmp/wsclean/cleanalgorithms/
		cp -v interface/*.{c,h,cpp} /tmp/wsclean/interface
		cp -v lofar/*.{h,cpp} /tmp/wsclean/lofar/
		cp -v model/*.{h,cpp} /tmp/wsclean/model/
		cp -v msproviders/*.{h,cpp} /tmp/wsclean/msproviders
		cp -v parser/*.h /tmp/wsclean/parser
		cp -v pywsclean/*.py /tmp/wsclean/pywsclean
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
		cat ../wscversion.h
fi
