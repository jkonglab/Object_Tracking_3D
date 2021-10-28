#include "ParticleFilter.h"



void ParticleFilter::initiate(const Mat &img)
{
	T = 1;
	int type = img.type();

	//detect cells
	vector<vector<Point3i> > cells;
	maxLabel = detectCell(img, cells);

	//initiate state
	M = cells.size();
	cout << "there are " << maxLabel << " cells detected." << endl;
	vector<ParticleState> X0(M);
	for (int i = 0; i < M; i++)
	{
		initState(cells[i], X0[i]);
		X0[i].id = img.at<ushort>(cells[i][0].z, cells[i][0].y, cells[i][0].x);
	}

	Xt = X0;
	t = 0;
	N = M * Ns;
}


void ParticleFilter::predict()
{
	int i, j, idx;
	vector<ParticleState> p(N);
	for (i = 0; i < M; i++)
	{
		int t_x = Xt[i].x;
		int t_y = Xt[i].y;
		int t_z = Xt[i].z;
		//int t_I = Xt[i].I;
		double t_vx = Xt[i].vx;
		double t_vy = Xt[i].vy;
		double t_vz = Xt[i].vz;
		double t_sx = Xt[i].sx;
		double t_sy = Xt[i].sy;
		double t_sz = Xt[i].sz;
		double t_theta = Xt[i].theta;

		/***************************************** Gaussian sampling **********************************************/
		j = 0;
		while (j < Ns)
		{
			idx = i*Ns + j;
			p[idx].x = randGuassian(t_x + t_vx*T, q1*T);
			p[idx].y = randGuassian(t_y + t_vy*T, q1*T);
			p[idx].z = randGuassian(t_z + t_vz*T, 1);
			p[idx].sx = randGuassian(t_sx, q2*T);
			p[idx].sy = randGuassian(t_sy, q2*T);
			p[idx].sz = t_sz;
			p[idx].theta = randGuassian(t_theta, q2*T);
			j++;
		}
	}
	particle = p;
}


// Update procedure reference: Vermaak 2007
void ParticleFilter::update(Mat &img)
{
	int i, j, k;
	int idx;
	double factor;

	vector<double> wc(M, 0);
	vector<double> wp(N);
	vector<ParticleState> X(M);

	vector<int> disappear;
	double tx, ty, tz, tvx, tvy, tvz, tsx, tsy, tsz, ttheta, tI;

	vector<vector<Point3i> > newCells;
	int newMax = detectCell(img, newCells);
	int newNum = newCells.size();

	vector<vector<int> > newCellIdxs;
	coord2idx(newCells, newCellIdxs);

	for (i = 0; i < M; i++)
	{
		bool flag = 0;
		for (j = 0; j < Ns; j++)
		{
			double L = calcWeight(newCellIdxs, particle[i*Ns + j]);
			if (L > 0.5)
				flag = 1;
			wp[i*Ns + j] = L;
			wc[i] += L;
		}

		tx = ty = tz = tvx = tvy = tvz = tsx = tsy = tsz = ttheta = tI = 0;
		for (j = 0; j < Ns; j++)
		{
			idx = i*Ns + j;
			factor = wp[idx] / wc[i];
			tx += particle[idx].x * factor;
			ty += particle[idx].y * factor;
			tz += particle[idx].z * factor;
			tsx += particle[idx].sx * factor;
			tsy += particle[idx].sy * factor;
			tsz += particle[idx].sz * factor;
			ttheta += particle[idx].theta * factor;
		}

		X[i].x = round(tx);
		X[i].y = round(ty);
		X[i].z = round(tz);
		X[i].vx = X[i].x - Xt[i].x;
		X[i].vy = X[i].y - Xt[i].y;
		X[i].vz = X[i].z - Xt[i].z;
		X[i].sx = tsx;
		X[i].sy = tsy;
		X[i].sz = tsz;
		X[i].theta = ttheta;
		X[i].id = Xt[i].id;
	}

	cout << "updating finished. Mapping begins." << endl;

	// draw probability map
	vector<int> value2idxMap(newMax + 1, 0);
	for (i = 0; i < newNum; i++)
		value2idxMap[img.at<ushort>(newCells[i][0].z, newCells[i][0].y, newCells[i][0].x)] = i + 1;

	vector<vector<int> > updatedIdx(M);
	for (i = 0; i < M; i++)
		state2idx(X[i], updatedIdx[i]);

	double prob;
	vector<bool> matched(M, false);
	vector<int> map(newNum, -1);
	for (i = 0; i < newNum; i++)
	{
		ParticleState tmp;
		initState(newCells[i], tmp);

		double maxProb = 0;
		int maxIdx = -1;
		for (j = 0; j < M; j++)
		{
			prob = likelihood(newCellIdxs[i], updatedIdx[j]);
			if (prob > maxProb) {
				maxProb = prob;
				maxIdx = j;
			}
		}

		if (maxIdx >= 0) {
			if (!matched[maxIdx] && maxProb > thresh) {
				map[i] = maxIdx;
				matched[maxIdx] = true;
			}
		}
	}

	// re-label
	for (i = 0; i < newNum; i++)
	{
		if (map[i] != -1)
			map[i] = X[map[i]].id;
		else
			map[i] = ++maxLabel;
	}

	map.insert(map.begin(), 0);	// don't forget 0

	ushort *pd = (ushort*)img.data;
	for (i = 0; i < total; i++)
	{
		*pd = map[value2idxMap[*pd]];
		pd++;
	}

	// generate new state
	detectCell(img, newCells);
	newNum = newCells.size();
	vector<ParticleState> newX(newNum);
	for (i = 0; i < newNum; i++)
	{
		initState(newCells[i], newX[i]);
		newX[i].id = img.at<ushort>(newCells[i][0].z, newCells[i][0].y, newCells[i][0].x);
	}

	Xt = newX;
	M = newNum;
	N = M * Ns;
	t++;
}


void ParticleFilter::print3DImage(char* path, const Mat &src, Size size, int depth)
{
	uchar* p = src.data;
	for (int i = 0; i < depth; i++)
	{
		int num = size.area();
		Mat tmp(size, CV_8UC1);
		uchar* pt = tmp.data;
		for (int j = 0; j < num; j++)
			*pt++ = *p++;

		char name[30];
		sprintf(name, "%s\\%02d.jpg", path, i);
		imwrite(name, tmp);
	}
}


double ParticleFilter::randGuassian(double mu, double sigma)
{
	return mu + rng.gaussian(sigma);
}


void ParticleFilter::state2idx(const ParticleState state, vector<int> &idx)
{
	int i, j, k;
	double x, y, z;
	idx.clear();

	double p_sx = state.sx;
	double p_sy = state.sy;
	int p_sz = round(state.sz);
	double C = cos(state.theta*0.0174);
	double S = sin(state.theta*0.0174);
	if (p_sx <= 0 || p_sy <= 0 || p_sz <= 0)
		return;

	int w = p_sx*abs(C) + p_sy*abs(S);
	int h = p_sx*abs(S) + p_sy*abs(C);
	double SX = p_sx * p_sx;
	double SY = p_sy * p_sy;
	double SZ = p_sz * p_sz;

	int left = 0 > state.x - 2 * w ? 0 : state.x - 2 * w;
	int right = width - 1 < state.x + 2 * w ? width - 1 : state.x + 2 * w;
	int front = 0 > state.y - 2 * h ? 0 : state.y - 2 * h;
	int back = height - 1 < state.y + 2 * h ? height - 1 : state.y + 2 * h;
	int top = 0 > state.z - 2 * p_sz ? 0 : state.z - 2 * p_sz;
	int bottom = depth - 1 < state.z + 2 * p_sz ? depth - 1 : state.z + 2 * p_sz;

	int roiWidth = right - left + 1;
	int roiHeight = back - front + 1;
	int roiDepth = bottom - top + 1;

	if (roiWidth <= 0 || roiHeight <= 0 || roiDepth <= 0)
		return;
	for (i = top; i <= bottom; i++)
	{
		for (j = front; j <= back; j++)
		{
			for (k = left; k <= right; k++)
			{
				x = (k - state.x) * C - (j - state.y) * S;
				y = (k - state.x) * S + (j - state.y) * C;
				z = i - state.z;
				if(x*x / SX + y*y / SY + z*z / SZ <= 1)
					idx.push_back(i * width * height + j * width + k);
			}
		}
	}
}


void ParticleFilter::coord2idx(const vector<vector<Point3i> > &coord, vector<vector<int> > &idx)
{
	int size = coord.size();
	idx.clear();
	idx.resize(size);

	Point3i pt;
	for (int i = 0; i < size; i++)
	{
		int subsize = coord[i].size();
		for (int j = 0; j < subsize; j++)
		{
			pt = coord[i][j];
			idx[i].push_back(pt.x + pt.y * width + pt.z * width * height);
		}
	}
}


double ParticleFilter::calcWeight(const vector<vector<int> > &idx, ParticleState state)
{
	vector<int> stateIdx;
	state2idx(state, stateIdx);

	if (stateIdx.size() == 0)
		return 0;

	double weight = 0;
	double tmp;
	for (int i = 0 ; i < idx.size(); i++)
	{
		tmp = likelihood(stateIdx, idx[i]);
		if (tmp > weight)
			weight = tmp;
	}

	return weight;
}


double ParticleFilter::likelihood(const vector<int> &idx1, const vector<int> &idx2)
{
	int size1 = idx1.size();
	int size2 = idx2.size();

	if (size1 == 0 || size2 == 0)
		return 0;

	if (idx1[0] > idx2[size2-1] || idx2[0] > idx1[size1-1])
		return 0;

	int a, b, unionNum;
	a = b = unionNum = 0;
	while (a < size1 && b < size2)
	{
		if (idx1[a] > idx2[b])
		{
			b++;
		}
		else if (idx1[a] < idx2[b])
		{
			a++;
		}
		else
		{
			unionNum++;
			a++;
			b++;
		}
	}

	return (double) unionNum / (size1 + size2 - unionNum);
}


int ParticleFilter::detectCell(const Mat& src, vector<vector<Point3i> > &cells)
{
	int max = 0;
	ushort* pb = (ushort*)src.data;
	for (int i = 0; i < total; i++)
	{
		if (*pb > max)
			max = *pb;
		pb++;
	}
	cout << "max = " << max << endl;

	cells.clear();
	cells.resize(max);

	pb = (ushort*)src.data;
	for (int i = 0; i < depth; i++)
	{
		for (int j = 0; j < height; j++)
		{
			for (int k = 0; k < width; k++)
			{
				if (*pb != 0)
					cells[(*pb) - 1].push_back(Point3i(k, j, i));
				pb++;
			}
		}
	}

	// remove little or empty cell
	for (int i = max - 1; i >= 0; i--)
	{
		if (cells[i].size() < 50)
			cells.erase(cells.begin() + i);
	}
	return max;
}



void ParticleFilter::initState(vector<Point3i> cell, ParticleState &state)
{
	int n = cell.size();

	//find range in z axis
	int zmin = depth;
	int zmax = 0;
	for (int i = 0; i < n; i++)
	{
		int z = cell[i].z;
		if (z < zmin)
			zmin = z;
		if (z > zmax)
			zmax = z;
	}

	//find the largest layer
	int layerNum = zmax - zmin + 1;
	vector<vector<Point3i> > layers(layerNum);
	for (int i = 0; i < n; i++)
		layers[cell[i].z - zmin].push_back(cell[i]);

	int idx = 0;
	int maxSize = 0;
	for (int i = 0; i < layerNum; i++)
	{
		if (layers[i].size() > maxSize)
		{
			maxSize = layers[i].size();
			idx = i;
		}
	}

	//draw the layer on a plane and calculate its center as state position
	int centerX = 0;
	int centerY = 0;
	Mat plane = Mat::zeros(width, height, CV_8UC1);
	for (int i = 0; i < layers[idx].size(); i++)
	{
		centerX += layers[idx][i].x;
		centerY += layers[idx][i].y;
		plane.at<uchar>(layers[idx][i].y, layers[idx][i].x) = 255;
	}
	state.x = centerX / layers[idx].size();
	state.y = centerY / layers[idx].size();
	state.z = idx + zmin;

	//find contour and fit ellipsoid
	vector<vector<Point> > contour;
	findContours(plane, contour, CV_RETR_EXTERNAL, CV_CHAIN_APPROX_NONE);

	idx = 0;
	maxSize = 0;
	for (int i = 0; i < contour.size(); i++)
	{
		if (contour[i].size() > maxSize)
		{
			maxSize = contour[i].size();
			idx = i;
		}
	}

	RotatedRect rect = fitEllipse(contour[idx]);
	state.theta = -rect.angle;
	state.sx = rect.size.width / 2.0;
	state.sy = rect.size.height / 2.0;
	state.sz = (zmax - idx - zmin > idx ? zmax - idx - zmin : idx) / 2.0;

	//set velocity to zero;
	state.vx = 0;
	state.vy = 0;
	state.vz = 0;
}
