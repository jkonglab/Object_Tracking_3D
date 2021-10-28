// PF.cpp
//

#include "ParticleFilter.h"


int main(int argc, char *argv[])
{
	MATFile *matFile;
	mxArray *mxarray;
	double *data;

	string dstPath = "../data/";
	string srcFile = dstPath + "LabeledData.mat";
	matFile = matOpen(srcFile.data(), "r");
	cout << "file opened" << endl;

	mxarray = matGetVariable(matFile, "frameNum");
	data = (double*)mxGetData(mxarray);
	int frameNum = *data;
	mxFree(data);

	mxarray = matGetVariable(matFile, "width");
	data = (double*)mxGetData(mxarray);
	int xRange = *data;
	mxFree(data);

	mxarray = matGetVariable(matFile, "height");
	data = (double*)mxGetData(mxarray);
	int yRange = *data;
	mxFree(data);

	mxarray = matGetVariable(matFile, "depth");
	data = (double*)mxGetData(mxarray);
	int zRange = *data;
	mxFree(data);

	int range[3] = { zRange, yRange, xRange };
	int total = xRange * yRange * zRange;

	cout << "parameters obtained" << endl;

	int *mData, *md;
	int n = xRange*yRange;
	const int *nd;
	vector<Mat> src(frameNum);
	for (int i = 0; i < frameNum; i++)
	{
		char variableName[10];
		sprintf(variableName, "frame%02d", i + 1);
		mxarray = matGetVariable(matFile, variableName);
		mData = (int*)mxGetData(mxarray);
		md = mData;

		src[i].create(3, range, CV_16UC1);
		ushort *pd = (ushort*)src[i].data;
		for (int j = 0; j < total; j++)
			*pd++ = *md++;

		if (i == 0)
			nd = mxGetDimensions_700(mxarray);

		mxFree(mData);
	}

	matClose(matFile);
	cout << "matrice obtained" << endl;

	int N = 500;
	double q1 = 5;
	double q2 = 1;
	ParticleFilter pf;
	pf.setThresh(0.0);
	pf.setSize(xRange, yRange, zRange);
	pf.setTrackParam(N, q1, q2);
	pf.initiate(src[0]);

	for (int j = 1; j < frameNum; j++)
	{
		cout << "frame " << j << endl;

		pf.predict();
		pf.update(src[j]);

		string dstFile = dstPath + "trackingResult.mat";
		matFile = matOpen(dstFile.data(), "w7.3");

		mxarray = mxCreateNumericArray_700(3, nd, mxINT32_CLASS, mxREAL);
		int *out = (int*) mxGetData(mxarray);
		ushort *pd = (ushort*)src[0].data;
		for (int i = 0; i < total; i++)
			*out++ = *pd++;
		matPutVariable(matFile, "frame01", mxarray);

		for (int i = 1; i < frameNum; i++)
		{
			char variableName[10];
			sprintf(variableName, "frame%02d", i + 1);
			out = (int*)mxGetData(mxarray);
			pd = (ushort*)src[i].data;
			for (int j = 0; j < total; j++)
				*out++ = *pd++;
			matPutVariable(matFile, variableName, mxarray);
		}

		matClose(matFile);
	}

	return 0;
}
