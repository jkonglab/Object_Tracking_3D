#ifndef _PARTICLEFILTER_H_
#define _PARTICLEFILTER_H_

#include <stdio.h>
#include <fstream>
#include <opencv/highgui.h>
#include <opencv/cv.h>
#include <mat.h>
#include <math.h>

using namespace std;
using namespace cv;


struct ParticleState
{
	int id;
	int x, y, z;		//position
	double vx, vy, vz;	//velocity
	double sx, sy, sz;	//sigma
	double theta;
};

class ParticleFilter
{
public:
	ParticleFilter()  { rng = RNG(getTickCount()); thresh = 0; }
	void setSize(int m_width, int m_height, int m_depth)
	{ width = m_width; height = m_height; depth = m_depth; total = width*height*depth; }
	void setTrackParam(int m_Ns, double m_q1, double m_q2)
	{ Ns = m_Ns; q1 = m_q1; q2 = m_q2; }
	void setThresh(double m_thresh) { thresh = m_thresh; }

	void initiate(const Mat &img);
	void predict();
	void update(Mat &img);

	void print3DImage(char* path, const Mat &src, Size size, int depth);
	int getObjectNum()		{ return M;  }
	vector<ParticleState> getState()	{ return Xt; }

private:
	vector<ParticleState> Xt;		//current state Xt
	vector<ParticleState> particle; //predicted particles
	int t;		//current time

	//image size
	int width;
	int height;
	int depth;
	int total;

	double thresh;	//threshold in cell detection
	int maxLabel;

	int M;		//number of objects
	int Ns;		//number of particles per object
	int N;		//number of particles in total;
	int gamma;	//ratio of two sampling method

	//parameters
	int T;			//sample interval
	double q1;		//noise level for movement
	double q2;		//noise level for shape

	RNG rng;	//generate random number
	double randGuassian(double mu, double sigma);
	void state2idx(const ParticleState state, vector<int> &idx);
	void coord2idx(const vector<vector<Point3i> > &coord, vector<vector<int> > &idx);
	double calcWeight(const vector<vector<int> > &idx, ParticleState state);
	double likelihood(const vector<int> &idx1, const vector<int> &idx2);

	int detectCell(const Mat &src, vector<vector<Point3i> > &cells);
	void initState(vector<Point3i> cell, ParticleState &state);
};

#endif
