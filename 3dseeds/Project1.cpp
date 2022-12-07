#include "Project1.h"
#include "polyhedron.h"
#include <iostream>
#include <algorithm>
#include "GL/freeglut.h"


extern Polyhedron* poly;

typedef double db;

void HSVtoRGB(icVector3& rgb, const icVector3& hsv){
	db H = hsv.x, S = hsv.y, V = hsv.z;
	db C = S*V;
	db X = C*(1-abs(fmod(H/60.0, 2) - 1));
	db m = V-C;
	db r, g, b;

	if (H >= 0 && H < 60){
		r = C;
		g = X;
		b = 0;
	}else if (H >= 60 && H < 120){
		r = X;
		g = C;
		b = 0;
	}else if (H >= 120 && H < 180){
		r = 0;
		g = C;
		b = X;
	}else if (H >= 180 && H < 240){
		r = 0;
		g = X;
		b = C;
	}else if (H >= 240 && H < 300){
		r = X;
		g = 0;
		b = C;
	}else{
		r = C;
		g = 0;
		b = X;
	}

	rgb.x = (r+m);
	rgb.y = (g+m);
	rgb.z = (b+m);
}

void RGBtoHSV(icVector3& hsv, const icVector3& rgb){
	db R = rgb.x;
	db G = rgb.y;
	db B = rgb.z;
	db& H = hsv.x;
	db& S = hsv.y;
	db& V = hsv.z;
	
	db cmax = std::max({R, G, B});
	db cmin = std::min({R, G, B});
	db D = cmax - cmin;

	if (D == 0){
		hsv.x = 0;
	}
	else if (cmax == R){
		hsv.x = fmod(60 * ((G - B)/D) + 360, 360);
	}else if (cmax == G){
		hsv.x = fmod(60 * ((B - R)/D) + 120, 360);
	}else if (cmax == B){
		hsv.x = fmod(60 * ((R - G)/D) + 240, 360);
	}
	if (cmax == 0){
		hsv.y = 0;
	}
	else if (cmax != 0){
		hsv.y = (D/cmax);
	}
	hsv.z = cmax;
}

void findMaxMin(db& max, db& min) {
	min = INFINITY;
	max = -min;
	for (int i = 0; i < poly->nverts; ++i) {
		Vertex* vrts = poly->vlist[i];
		if (vrts->scalar < min) {
			min = vrts->scalar;
		}
		if (vrts->scalar > max) {
			max = vrts->scalar;
		}
	}
}

void q1_A_GreyScaleMap() {
	db min;
	db max;

	findMaxMin(max, min);
	std::cout << "Grayscale Map\n";
	std::cout << "Max: " << max << " Min: " << min << std::endl;

	for (int i = 0; i < poly->nverts; ++i) {
		Vertex* v = poly->vlist[i];
		db s_v = v->scalar;
		db gray = (s_v - min) / (max - min);
		v->R = v ->G = v->B = gray;
	}

	glutPostRedisplay();
}

void q1_B_BiColorMap() {
	db min;
	db max;

	findMaxMin(max, min);
	std::cout << "BiColor Map\n";
	std::cout << "Max: " << max << " Min: " << min << std::endl;

	for (int i = 0; i < poly->nverts; ++i) {
		Vertex* v = poly->vlist[i];
		db s_v = v->scalar;
		icVector3 c1(1.0, 0.0, 0.0);
		icVector3 c2(0.0, 0.0, 1.0);

		db l = (s_v - min) / (max - min);
		db r = (max - s_v) / (max - min);
		icVector3 c = c1 * l + c2 * r;
		v->R = c.x;
		v->G = c.y;
		v->B = c.z;
	}
	glutPostRedisplay();
}

void q1_C_HeatMap(){
	db min;
	db max;

	findMaxMin(max, min);
	std::cout << "Height field visualized with Uniform Color\n";
	std::cout << "Max: " << max << " Min: " << min << std::endl;

	for (int i = 0; i < poly->nverts; ++i) {
		Vertex* v = poly->vlist[i];
		db s_v = v->scalar;

		icVector3 c1(1.0, 0.0, 0.0);
		icVector3 c2(0.0, 0.0, 1.0);
		icVector3 HSVc1, HSVc2;
		RGBtoHSV(HSVc1, c1);
		RGBtoHSV(HSVc2, c2);

		db l = (s_v - min) / (max - min);
		db r = (max - s_v) / (max - min);

		icVector3 HSVc = HSVc1 * l + HSVc2 * r;
		icVector3 RGBc;
		HSVtoRGB(RGBc, HSVc);

		v->R = RGBc.x;
		v->G = RGBc.y;
		v->B = RGBc.z;
	}
	glutPostRedisplay();
}

void q2_HeightMap(){
	db min;
	db max;

	findMaxMin(max, min);
	db diff = max - min;
	std::cout << "Augmented height Field\n";
	std::cout << "Max: " << max << " Min: " << min << std::endl;

	for (int i = 0; i < poly->nverts; ++i) {
		Vertex* v = poly->vlist[i];
		db s_v = v->scalar;
		db l = ((s_v - min)/diff)*2;
		v->z = l;
	}
	glutPostRedisplay();
}


void HeighTest(){
	// Not part of the assignment
	std::cout << "Height Test\n";

	for (int i = 0; i < poly->nverts; ++i) {
		Vertex* v = poly->vlist[i];
		v->z = v->scalar;
	}
	glutPostRedisplay();
}
