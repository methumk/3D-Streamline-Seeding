#pragma once
#include "polyhedron.H"
#include <vector>
#include <utility>
#include <list>
#include <algorithm>
#include <iostream>
#include "glError.h"
#include "gl/glew.h"
#include "gl/freeglut.h"
#include "Project1.h"

#define EPSILON 1.0e-5

// stored x-y position, z is the scalar value
static std::vector<icVector3> crit_max;
static std::vector<icVector3> crit_min;
static std::vector<icVector3> crit_saddle;

class POLYLINE {
public:
	std::list<icVector3> vertices;
	icVector3 rgb = icVector3(1.0, 0, 0);
	double weight = 1.0;
	bool isNeighbor(const POLYLINE& line) {
		if ((vertices.front() - line.vertices.front()).length() < EPSILON ||
			(vertices.front() - line.vertices.back()).length() < EPSILON ||
			(vertices.back() - line.vertices.front()).length() < EPSILON ||
			(vertices.back() - line.vertices.back()).length() < EPSILON) {
			return true;
		}
		return false;
	}

	void merge(const POLYLINE& line) {
		if ((vertices.front() - line.vertices.front()).length() < EPSILON) {
			POLYLINE l = line;
			l.vertices.pop_front();
			for (auto i = l.vertices.begin(); i != l.vertices.end(); i++) {
				vertices.push_front(*i);
			}
		}
		else if ((vertices.front() - line.vertices.back()).length() < EPSILON) {
			POLYLINE reversel = line;
			reversel.vertices.pop_back();
			reversel.vertices.reverse();
			for (auto i = reversel.vertices.begin(); i != reversel.vertices.end(); i++) {
				vertices.push_front(*i);
			}
		}
		else if ((vertices.back() - line.vertices.front()).length() < EPSILON) {
			POLYLINE l = line;
			l.vertices.pop_front();
			for (auto i = l.vertices.begin(); i != l.vertices.end(); i++) {
				vertices.push_back(*i);
			}
		}
		else if ((vertices.back() - line.vertices.back()).length() < EPSILON) {
			POLYLINE reversel = line;
			reversel.vertices.pop_back();
			reversel.vertices.reverse();
			for (auto i = reversel.vertices.begin(); i != reversel.vertices.end(); i++) {
				vertices.push_back(*i);
			}
		}
	}

	void changeHeight(const double& height) {
		for (auto& v : vertices) {
			v.z = height;
		}
	}

	void clear() { vertices.clear(); }
};

void display_poliline(std::vector<POLYLINE>& polylines) {
	glDisable(GL_LIGHTING);
	glEnable(GL_LINE_SMOOTH);
	glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	for (auto& polyline : polylines) {
		glLineWidth(polyline.weight);

		glColor3f(polyline.rgb.entry[0], polyline.rgb.entry[1], polyline.rgb.entry[2]);
		glBegin(GL_LINE_STRIP);

		for (auto it = polyline.vertices.begin(); it != polyline.vertices.end(); it++) {
			glVertex3d(it->entry[0], it->entry[1], it->entry[2]);
		}

		glEnd();
	}

	glDisable(GL_BLEND);
	glLineWidth(1);
}

Vertex linearInterpolateByScalar(const Vertex& v0, const Vertex& v1, const double& thres) {
	double f0_f1 = v0.scalar - v1.scalar;
	Vertex r(0.0, 0.0, 0.0);
	if (std::abs(f0_f1) < 1.0e-5) {
		// interpolation is middle
		r.x = (v0.x + v1.x) / 2;
		r.y = (v0.y + v1.y) / 2;
		r.z = (v0.z + v1.z) / 2;
	}
	else {
		double t = std::abs((v0.scalar - thres) / ((v0.scalar - thres) - (v1.scalar - thres)));
		r.x = v0.x + t * (v1.x - v0.x);
		r.y = v0.y + t * (v1.y - v0.y);
		r.z = v0.z + t * (v1.z - v0.z);
	}
	return r;
}

void lookUpTable(std::vector<Vertex>& r, const Vertex& v0, const Vertex& v1, const Vertex& v2, const Vertex& v3, const double& thres) {
	r.reserve(2);
	int id = 0;

	if (v0.scalar <= thres + EPSILON) {
		id += 1;
	}
	if (v1.scalar <= thres + EPSILON) {
		id += 2;
	}
	if (v2.scalar <= thres + EPSILON) {
		id += 4;
	}
	if (v3.scalar <= thres + EPSILON) {
		id += 8;
	}

	double center = 0;
	switch (id) {
	case 0:
		break;
	case 1:
		r.push_back(linearInterpolateByScalar(v0, v1, thres));
		r.push_back(linearInterpolateByScalar(v0, v3, thres));
		break;
	case 2:
		r.push_back(linearInterpolateByScalar(v0, v1, thres));
		r.push_back(linearInterpolateByScalar(v1, v2, thres));
		break;
	case 3:
		r.push_back(linearInterpolateByScalar(v1, v2, thres));
		r.push_back(linearInterpolateByScalar(v0, v3, thres));
		break;
	case 4:
		r.push_back(linearInterpolateByScalar(v1, v2, thres));
		r.push_back(linearInterpolateByScalar(v2, v3, thres));
		break;
	case 5:
		center = v0.scalar + v1.scalar + v2.scalar + v3.scalar;
		center /= 4;
		if (center <= thres) {
			r.push_back(linearInterpolateByScalar(v0, v1, thres));
			r.push_back(linearInterpolateByScalar(v1, v2, thres));
			r.push_back(linearInterpolateByScalar(v2, v3, thres));
			r.push_back(linearInterpolateByScalar(v0, v3, thres));
		}
		else {
			r.push_back(linearInterpolateByScalar(v0, v1, thres));
			r.push_back(linearInterpolateByScalar(v0, v3, thres));
			r.push_back(linearInterpolateByScalar(v1, v2, thres));
			r.push_back(linearInterpolateByScalar(v2, v3, thres));
		}
		break;
	case 6:
		r.push_back(linearInterpolateByScalar(v0, v1, thres));
		r.push_back(linearInterpolateByScalar(v2, v3, thres));
		break;
	case 7:
		r.push_back(linearInterpolateByScalar(v2, v3, thres));
		r.push_back(linearInterpolateByScalar(v0, v3, thres));
		break;
	case 8:
		r.push_back(linearInterpolateByScalar(v2, v3, thres));
		r.push_back(linearInterpolateByScalar(v0, v3, thres));
		break;
	case 9:
		r.push_back(linearInterpolateByScalar(v0, v1, thres));
		r.push_back(linearInterpolateByScalar(v2, v3, thres));
		break;
	case 10:
		// ambiguous so get center
		center = v0.scalar + v1.scalar + v2.scalar + v3.scalar;
		center /= 4;
		if (center <= thres) {
			r.push_back(linearInterpolateByScalar(v0, v1, thres));
			r.push_back(linearInterpolateByScalar(v0, v3, thres));
			r.push_back(linearInterpolateByScalar(v1, v2, thres));
			r.push_back(linearInterpolateByScalar(v2, v3, thres));
		}
		else {
			r.push_back(linearInterpolateByScalar(v0, v1, thres));
			r.push_back(linearInterpolateByScalar(v1, v2, thres));
			r.push_back(linearInterpolateByScalar(v2, v3, thres));
			r.push_back(linearInterpolateByScalar(v0, v3, thres));
		}
		break;
	case 11:
		r.push_back(linearInterpolateByScalar(v1, v2, thres));
		r.push_back(linearInterpolateByScalar(v2, v3, thres));
		break;
	case 12:
		r.push_back(linearInterpolateByScalar(v1, v2, thres));
		r.push_back(linearInterpolateByScalar(v0, v3, thres));
		break;
	case 13:
		r.push_back(linearInterpolateByScalar(v0, v1, thres));
		r.push_back(linearInterpolateByScalar(v1, v2, thres));
		break;
	case 14:
		r.push_back(linearInterpolateByScalar(v0, v1, thres));
		r.push_back(linearInterpolateByScalar(v0, v3, thres));
		break;
	case 15:
		break;
	}
}

void marchingSquare(std::list<POLYLINE>& edges, const Polyhedron& poly, const double& thresh) {
	for (int i = 0; i < poly.nquads; ++i) {
		std::vector<Vertex> r;
		lookUpTable(r,
			*poly.qlist[i]->verts[0],
			*poly.qlist[i]->verts[1],
			*poly.qlist[i]->verts[2],
			*poly.qlist[i]->verts[3],
			thresh);

		if (r.size() > 0) {
			for (int j = 0; j < r.size() / 2; j++) {
				POLYLINE line;
				auto v0 = icVector3(
					r[j * 2].x,
					r[j * 2].y,
					r[j * 2].z);
				auto v1 = icVector3(
					r[j * 2 + 1].x,
					r[j * 2 + 1].y,
					r[j * 2 + 1].z);
				line.vertices.push_back(v0);
				line.vertices.push_back(v1);
				edges.push_back(line);
			}
		}
	}
}

void makePolylineFromEdges(std::vector<POLYLINE>& polylines, const std::list<POLYLINE>& edges) {
	polylines.reserve(edges.size());
	std::list<POLYLINE> edges_temp(edges);

	while (edges_temp.size() > 0) {
		polylines.push_back(edges_temp.front());
		edges_temp.erase(edges_temp.begin());
		int init_size = 0;

		while (init_size != edges_temp.size()) {
			init_size = edges_temp.size();
			for (auto i = edges_temp.begin(); i != edges_temp.end(); ) {
				if (polylines.back().isNeighbor(*i)) {
					polylines.back().merge(*i);
					i = edges_temp.erase(i);
				}
				else
					i++;
			}
		}
	}

}

// Draws a dot at the specified location
// x, y, z are the coordinates of the dot
// radius: radius of the dot
// R: red channel for the dot color [0,1]
// B: blue channel for the dot color [0,1]
// G: green channel of the dot color [0,1]
//void drawDot(double x, double y, double z, double radius = 0.15, float R = 0.0, float G = 0.0, float B = 0.0)
//{
//	glDisable(GL_POLYGON_OFFSET_FILL);
//	glEnable(GL_DEPTH_TEST);
//	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
//	glShadeModel(GL_SMOOTH);
//	glEnable(GL_LIGHTING);
//	glEnable(GL_LIGHT0);
//	glEnable(GL_LIGHT1);
//
//	CHECK_GL_ERROR();
//
//	GLfloat mat_diffuse[4];
//	mat_diffuse[0] = R;
//	mat_diffuse[1] = G;
//	mat_diffuse[2] = B;
//	mat_diffuse[3] = 1.0;
//
//	CHECK_GL_ERROR();
//
//	glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
//
//	CHECK_GL_ERROR();
//
//	GLUquadric* quadric = gluNewQuadric();
//	glPushMatrix();
//	glTranslated(x, y, z);
//	glColor3f(R, G, B);
//	gluSphere(quadric, radius, 16, 16);
//	glPopMatrix();
//	gluDeleteQuadric(quadric);
//}


void findCritMinMax(Vertex** const vlist, const int& idx, const double& minScalar, const double& maxScalar) {
	using std::cout;
	using std::endl;
	// length is both the row and height max lenght 1s based (cuz its square)
	int POLYLINE_LENGTH = 21;

	int numsurrounding = 0;
	// {N, W, NW, SW, S, E, NE, SE}
	// 0 = same, 1 = max, -1 = min
	// int above[8] = {0, 0, 0, 0, 0, 0, 0, 0};

	int surroundTotal = 0;
	// check N
	if (idx - POLYLINE_LENGTH >= 0) {
		numsurrounding++;
		if (vlist[idx - POLYLINE_LENGTH]->scalar > vlist[idx]->scalar) {
			// above[0] = 1;
			surroundTotal += 1;
		}
		else if (vlist[idx - POLYLINE_LENGTH]->scalar < vlist[idx]->scalar) {
			// above[0] = -1;
			surroundTotal -= 1;
		}
	}

	// check W
	if ((idx + 1) % POLYLINE_LENGTH != 0) {
		numsurrounding++;
		if (vlist[(idx + 1)]->scalar > vlist[idx]->scalar) {
			// above[1] = 1;
			surroundTotal += 1;
		}
		else if (vlist[(idx + 1)]->scalar < vlist[idx]->scalar) {
			// above[1] = -1;
			surroundTotal -= 1;
		}

		// check NW 
		if (idx - POLYLINE_LENGTH >= 0) {
			numsurrounding++;
			if (vlist[(idx + 1 - POLYLINE_LENGTH)]->scalar > vlist[idx]->scalar) {
				// above[2] = 1;
				surroundTotal += 1;
			}
			else if (vlist[(idx + 1 - POLYLINE_LENGTH)]->scalar < vlist[idx]->scalar) {
				// above[2] = -1;
				surroundTotal -= 1;
			}

		}

		// check SW 
		if (idx + POLYLINE_LENGTH < POLYLINE_LENGTH * POLYLINE_LENGTH) {
			numsurrounding++;
			if (vlist[(idx + 1 + POLYLINE_LENGTH)]->scalar > vlist[idx]->scalar) {
				// above[3] = 1;
				surroundTotal += 1;
			}
			else if (vlist[(idx + 1 + POLYLINE_LENGTH)]->scalar < vlist[idx]->scalar) {
				// above[3] = -1;
				surroundTotal -= 1;
			}
		}
	}

	// check S
	if (idx + POLYLINE_LENGTH < POLYLINE_LENGTH * POLYLINE_LENGTH) {
		numsurrounding++;
		if (vlist[idx + POLYLINE_LENGTH]->scalar > vlist[idx]->scalar) {
			// above[4] = 1;
			surroundTotal += 1;
		}
		else if (vlist[idx + POLYLINE_LENGTH]->scalar < vlist[idx]->scalar) {
			// above[4] = -1;
			surroundTotal -= 1;
		}
	}

	// check E
	if (idx % POLYLINE_LENGTH != 0) {
		numsurrounding++;
		if (vlist[idx - 1]->scalar > vlist[idx]->scalar) {
			// above[5] = 1;
			surroundTotal += 1;
		}
		else if (vlist[idx - 1]->scalar < vlist[idx]->scalar) {
			// above[5] = -1;
			surroundTotal -= 1;
		}

		// check NE
		if (idx - POLYLINE_LENGTH >= 0) {
			numsurrounding++;
			if (vlist[idx - 1 - POLYLINE_LENGTH]->scalar > vlist[idx]->scalar) {
				// above[6] = 1;
				surroundTotal += 1;
			}
			else if (vlist[idx - 1 - POLYLINE_LENGTH]->scalar < vlist[idx]->scalar) {
				// above[6] = -1;
				surroundTotal -= 1;
			}
		}

		// check SE
		if (idx + POLYLINE_LENGTH < POLYLINE_LENGTH * POLYLINE_LENGTH) {
			numsurrounding++;
			if (vlist[idx - 1 + POLYLINE_LENGTH]->scalar > vlist[idx]->scalar) {
				// above[7] = 1;
				surroundTotal += 1;
			}
			else if (vlist[idx - 1 + POLYLINE_LENGTH]->scalar < vlist[idx]->scalar) {
				// above[7] = -1;
				surroundTotal -= 1;
			}
		}
	}

	/* cout << "NS: " << numsurrounding;
	if (numsurrounding == 8){
		cout << "   Reg point   ";
	}else if (numsurrounding == 5){
		cout << "   Edge point   ";
	}else if (numsurrounding == 3){
		cout << "   VERTEX   ";
	}else{
		cout << "   WTF   ";
	} */

	cout << "TOTAL: " << surroundTotal << "   " << "Scalar: " << vlist[idx]->scalar << "   ";
	// vlist[idx]->scalar
	if (surroundTotal == 8) {
		cout << "MIN " << vlist[idx]->x << " " << vlist[idx]->y << " " << vlist[idx]->z << "\n";
		crit_min.push_back(icVector3(vlist[idx]->x, vlist[idx]->y, 2 * ((vlist[idx]->scalar - minScalar) / (maxScalar - minScalar))));
	}
	else if (surroundTotal == -8) {
		cout << "MAX " << vlist[idx]->x << " " << vlist[idx]->y << " " << vlist[idx]->z << "\n";
		crit_max.push_back(icVector3(vlist[idx]->x, vlist[idx]->y, 2 * ((vlist[idx]->scalar - minScalar) / (maxScalar - minScalar))));
	}
}

//void findSaddlePoints(const Quad* const quads, const double& minScalar, const double& maxScalar) {
//	using std::cout;
//	using std::endl;
//
//	typedef double db;
//	db x1 = quads->verts[2]->x;
//	db y1 = quads->verts[2]->y;
//	db x2 = quads->verts[0]->x;
//	db y2 = quads->verts[0]->y;
//
//	db f11 = quads->verts[2]->scalar, f21 = quads->verts[3]->scalar, f12 = quads->verts[1]->scalar, f22 = quads->verts[0]->scalar;
//	db crit_denom = f11 - f21 - f12 + f22;
//
//	// Saddle crit points
//	db x0 = (x2 * f11 - x1 * f21 - x2 * f12 + x1 * f22) / (crit_denom);
//	db y0 = (y2 * f11 - y2 * f21 - y1 * f12 + y1 * f22) / (crit_denom);
//
//
//	// Saddle point if within quads
//	if ((x2 > x0 && x1 < x0) && (y2 > y0 && y1 < y0)) {
//
//		cout << "Saddle point: " << x0 << " " << y0 << endl;
//		db interpolate_scalar = (f11 + f21 + f12 + f22) / 4;
//		crit_saddle.push_back(icVector3(x0, y0, 2 * ((interpolate_scalar - minScalar) / (maxScalar - minScalar))));
//		cout << "0: " << x2 << " " << y2 << endl;
//		cout << "1: " << x1 << " " << y2 << endl;
//		cout << "2: " << x1 << " " << y1 << endl;
//		cout << "3: " << x2 << " " << y1 << endl;
//	}
//}
//
//void critPoints(const Polyhedron& poly) {
//	using std::cout;
//	using std::endl;
//#define wa '\n'
//	typedef double db;
//
//	cout << "Num verts: " << poly.nverts << wa;
//	db max = 0, min = 0;
//	findMaxMin(max, min);
//
//	for (int i = 0; i < poly.nverts; ++i) {
//		cout << "i: " << i << "\t";
//		findCritMinMax(poly.vlist, i, min, max);
//		cout << "\n";
//	}
//
//	for (int i = 0; i < poly.nquads; ++i) {
//		findSaddlePoints(poly.qlist[i], min, max);
//	}
//}

/*
Displays min, max, saddle critical points
	- Max points: 		red
	- Min points: 		green
	- Saddle points: 	blue
*/
//void display_crit_points() {
//	for (int k = 0; k < crit_max.size(); k++)
//	{
//		drawDot(crit_max[k].x, crit_max[k].y, crit_max[k].z, 0.1499999999999999944, 1, 0, 0);
//	}
//	for (int k = 0; k < crit_min.size(); k++)
//	{
//		drawDot(crit_min[k].x, crit_min[k].y, crit_min[k].z, 0.1499999999999999944, 0, 1, 0);
//	}
//	for (int k = 0; k < crit_saddle.size(); k++)
//	{
//		drawDot(crit_saddle[k].x, crit_saddle[k].y, crit_saddle[k].z, 0.1499999999999999944, 0, 0, 1);
//	}
//}

