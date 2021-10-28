

tracking : PF.o ParticleFilter.o
	/opt/rh/devtoolset-3/root/bin/g++ -g -o PF PF.o ParticleFilter.o -L /opt/MATLAB/R2016a/bin/glnxa64 -L /usr/lib64/ -Wl,-rpath,/usr/lib64/,-rpath,/opt/MATLAB/R2016a/bin/glnxa64 -l opencv_calib3d -l opencv_contrib -l opencv_core -l opencv_features2d -l opencv_flann -l opencv_highgui -l opencv_imgproc -l opencv_legacy -l opencv_ml -l opencv_objdetect -l opencv_photo -l opencv_stitching -l opencv_superres -l opencv_ts -l opencv_video -l opencv_videostab -l mat -l mx

PF.o : PF.cpp ParticleFilter.h
	/opt/rh/devtoolset-3/root/bin/g++ -g -o PF.o -c PF.cpp -I /opt/MATLAB/R2016a/extern/include/

ParticleFilter.o : ParticleFilter.cpp ParticleFilter.h
	/opt/rh/devtoolset-3/root/bin/g++ -g -o ParticleFilter.o -c ParticleFilter.cpp -I /opt/MATLAB/R2016a/extern/include/

clean :
	rm PF PF.o ParticleFilter.o
