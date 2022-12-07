#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <fstream>
#include <vector>

#include <netcdf.h>
#include <netcdf_aux.h>
#include <netcdf_mem.h>
#include <netcdf_dispatch.h>
#include <netcdf_meta.h>

#include "glError.h"
#include "gl/glew.h"
#include "gl/freeglut.h"
#include "ply.h"
#include "icVector.H"
#include "icMatrix.H"
#include "polyhedron.h"
#include "trackball.h"
#include "tmatrix.h"
#include "ppm.h"
#include "VectorFieldTopology.h"
#include "drawUtil.h"

Polyhedron* poly;
std::vector<PolyLine> lines;
std::vector<icVector3> points;

std::vector<POLYLINE> polylines;
std::list<Singularity> singularities;

/*scene related variables*/
const float zoomspeed = 0.9;
int win_width = 1024;
int win_height = 1024;
float aspectRatio = win_width / win_height;
const int view_mode = 0;		// 0 = othogonal, 1=perspective
const double radius_factor = 0.9;

/*
Use keys 1 to 0 to switch among different display modes.
Each display mode can be designed to show one type 
visualization result.

Predefined ones: 
display mode 1: solid rendering
display mode 2: show wireframes
display mode 3: render each quad with colors of vertices
*/
int display_mode = 1;

/*User Interaction related variabes*/
float s_old, t_old;
float rotmat[4][4];
double zoom = 1.0;
double translation[2] = { 0, 0 };
int mouse_mode = -2;	// -1 = no action, 1 = tranlate y, 2 = rotate

// IBFV related variables (Van Wijk 2002)
//https://www.win.tue.nl/~vanwijk/ibfv/
#define NPN		256 //orig was 64 -> 256
#define SCALE	4.0
#define ALPHA	8
float tmax = win_width / (SCALE * NPN);
float dmax = SCALE / win_width;
unsigned char* pixels;
int alpha = (255 * .2);
bool original_image = false;
bool flow_image = false;
bool singularity_IBVF = false;

// true allows user to select stream by holding cntrl and pressing left click
bool select_stream = false;

#define STEPSIZE 0.005
#define NUM_IMAGES 8
int image_idx = 0;
std::string images[NUM_IMAGES] = 
	{
		"../data/images/tidal.ppm",
		"../data/images/ram.ppm",
		"../data/images/hypno.ppm",
		"../data/images/relentless.ppm",
		"../data/images/art.ppm",
		"../data/images/bear.ppm",
		"../data/images/obscure.ppm",
		"../data/images/Lenna.ppm"
	};

#define NUM_VEC_FILES 8
int vec_idx = 0;
char* vectors[NUM_VEC_FILES] = 
	{
		"../data/vector_data/v1.ply",
		"../data/vector_data/v3.ply",
		"../data/vector_data/v4.ply",
		"../data/vector_data/v5.ply",
		"../data/vector_data/v6.ply",
		"../data/vector_data/v8.ply",
		"../data/vector_data/v9.ply",
		"../data/vector_data/v10.ply"
	};

// Reading 3d data, 384 x 384 x 130
static const int NDIMS = 3;
static const int NX = 384;
static const int NY = 384;
static const int NZ = 130;
//using namespace netcdf;

/******************************************************************************
Forward declaration of functions
******************************************************************************/

void init(void);
void initIBFV();

/*glut attaching functions*/
void keyboard(unsigned char key, int x, int y);
void motion(int x, int y);
void displayIBFV();
void display(void);
void mouse(int button, int state, int x, int y);
void mousewheel(int wheel, int direction, int x, int y);
void reshape(int width, int height);
void makePatternsImg(const std::string& fname);
void makePatternsImgNoise(const std::string& fname, float w);
void makePatternsImgEdges(const std::string& fname);

/*functions for element picking*/
void display_vertices(GLenum mode, Polyhedron* poly);
void display_quads(GLenum mode, Polyhedron* poly);
void display_selected_vertex(Polyhedron* poly);
void display_selected_quad(Polyhedron* poly);

/*display vis results*/
void display_polyhedron(Polyhedron* poly);

/******************************************************************************
Main program.
******************************************************************************/

int main(int argc, char* argv[])
{
	/*load mesh from ply file*/
	FILE* this_file = fopen("../data/vector_data/v1.ply", "r");
	poly = new Polyhedron(this_file);
	fclose(this_file);

	/* read nc file */
	//int dataIn[NX][NY][NZ];
	// NcFile dataFile("../data/ctbl3d.nc", NcFile::ReadOnly);
	//NcFile x = ncopen("./data/ctbl3d.nc", NC_WRITE);
	/*initialize the mesh*/
	 poly->initialize(); // initialize the mesh
	 poly->write_info();


	/*init glut and create window*/
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
	glutInitWindowPosition(20, 20);
	glutInitWindowSize(win_width, win_height);
	glutCreateWindow("Scientific Visualization");


	/*initialize openGL*/
	init();
	
	/*the render function and callback registration*/
	glutKeyboardFunc(keyboard);
	glutReshapeFunc(reshape);
	glutDisplayFunc(display);
	glutMotionFunc(motion);
	glutMouseFunc(mouse);
	glutMouseWheelFunc(mousewheel);
	
	/*event processing loop*/
	glutMainLoop();
	
	/*clear memory before exit*/
	poly->finalize();	// finalize everything
	free(pixels);
	return 0;
}

/******************************************************************************
Set projection mode
******************************************************************************/

void set_view(GLenum mode)
{
	GLfloat light_ambient0[] = { 0.3, 0.3, 0.3, 1.0 };
	GLfloat light_diffuse0[] = { 0.7, 0.7, 0.7, 1.0 };
	GLfloat light_specular0[] = { 0.0, 0.0, 0.0, 1.0 };

	GLfloat light_ambient1[] = { 0.0, 0.0, 0.0, 1.0 };
	GLfloat light_diffuse1[] = { 0.5, 0.5, 0.5, 1.0 };
	GLfloat light_specular1[] = { 0.0, 0.0, 0.0, 1.0 };

	glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient0);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse0);
	glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular0);

	glLightfv(GL_LIGHT1, GL_AMBIENT, light_ambient1);
	glLightfv(GL_LIGHT1, GL_DIFFUSE, light_diffuse1);
	glLightfv(GL_LIGHT1, GL_SPECULAR, light_specular1);


	glMatrixMode(GL_PROJECTION);
	if (mode == GL_RENDER)
		glLoadIdentity();

	if (aspectRatio >= 1.0) {
		if (view_mode == 0)
			glOrtho(-radius_factor * zoom * aspectRatio, radius_factor * zoom * aspectRatio, -radius_factor * zoom, radius_factor * zoom, -1000, 1000);
		else
			glFrustum(-radius_factor * zoom * aspectRatio, radius_factor * zoom * aspectRatio, -radius_factor * zoom, radius_factor * zoom, 0.1, 1000);
	}
	else {
		if (view_mode == 0)
			glOrtho(-radius_factor * zoom, radius_factor * zoom, -radius_factor * zoom / aspectRatio, radius_factor * zoom / aspectRatio, -1000, 1000);
		else
			glFrustum(-radius_factor * zoom, radius_factor * zoom, -radius_factor * zoom / aspectRatio, radius_factor * zoom / aspectRatio, 0.1, 1000);
	}

	GLfloat light_position[3];
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	light_position[0] = 5.5;
	light_position[1] = 0.0;
	light_position[2] = 0.0;
	glLightfv(GL_LIGHT0, GL_POSITION, light_position);
	light_position[0] = -0.1;
	light_position[1] = 0.0;
	light_position[2] = 0.0;
	glLightfv(GL_LIGHT2, GL_POSITION, light_position);
}

/******************************************************************************
Update the scene
******************************************************************************/

void set_scene(GLenum mode, Polyhedron* poly)
{
	glTranslatef(translation[0], translation[1], -3.0);

	/*multiply rotmat to current mat*/
	{
		int i, j, index = 0;

		GLfloat mat[16];

		for (i = 0; i < 4; i++)
			for (j = 0; j < 4; j++)
				mat[index++] = rotmat[i][j];

		glMultMatrixf(mat);
	}

	glScalef(0.9 / poly->radius, 0.9 / poly->radius, 0.9 / poly->radius);
	glTranslatef(-poly->center.entry[0], -poly->center.entry[1], -poly->center.entry[2]);
}

/******************************************************************************
Init scene
******************************************************************************/

void init(void) {

	mat_ident(rotmat);

	/* select clearing color */
	glClearColor(0.0, 0.0, 0.0, 0.0);  // background
	glShadeModel(GL_FLAT);
	glPolygonMode(GL_FRONT, GL_FILL);

	glDisable(GL_DITHER);
	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LESS);
	
	//set pixel storage modes
	glPixelStorei(GL_PACK_ALIGNMENT, 1);
	
	glEnable(GL_NORMALIZE);
	if (poly->orientation == 0)
		glFrontFace(GL_CW);
	else
		glFrontFace(GL_CCW);
}

/******************************************************************************
Initialize IBFV patterns
******************************************************************************/

void initIBFV()
{
	pixels = (unsigned char*)malloc(sizeof(unsigned char) * win_width * win_height * 3);
	memset(pixels, 255, sizeof(unsigned char) * win_width * win_height * 3);

	tmax = win_width / (SCALE * NPN);
	dmax = SCALE / win_width;

	int lut[256];
	int phase[NPN][NPN];
	GLubyte pat[NPN][NPN][4];
	int i, j, k;

	for (i = 0; i < 256; i++) lut[i] = i < 127 ? 0 : 255;
	for (i = 0; i < NPN; i++)
		for (j = 0; j < NPN; j++) phase[i][j] = rand() % 256;

	for (i = 0; i < NPN; i++)
	{
		for (j = 0; j < NPN; j++)
		{
			pat[i][j][0] =
				pat[i][j][1] =
				pat[i][j][2] = lut[(phase[i][j]) % 255];
			pat[i][j][3] = ALPHA;
		}
	}
	
	glNewList(1, GL_COMPILE);
	glTexImage2D(GL_TEXTURE_2D, 0, 4, NPN, NPN, 0, GL_RGBA, GL_UNSIGNED_BYTE, pat);
	glEndList();
}

/******************************************************************************
Pick objects from the scene
******************************************************************************/

int processHits(GLint hits, GLuint buffer[])
{
	unsigned int i, j;
	GLuint names, * ptr;
	double smallest_depth = 1.0e+20, current_depth;
	int seed_id = -1;
	unsigned char need_to_update;

	ptr = (GLuint*)buffer;
	for (i = 0; i < hits; i++) {  /* for each hit  */
		need_to_update = 0;
		names = *ptr;
		ptr++;

		current_depth = (double)*ptr / 0x7fffffff;
		if (current_depth < smallest_depth) {
			smallest_depth = current_depth;
			need_to_update = 1;
		}
		ptr++;
		current_depth = (double)*ptr / 0x7fffffff;
		if (current_depth < smallest_depth) {
			smallest_depth = current_depth;
			need_to_update = 1;
		}
		ptr++;
		for (j = 0; j < names; j++) {  /* for each name */
			if (need_to_update == 1)
				seed_id = *ptr - 1;
			ptr++;
		}
	}
	return seed_id;
}

/******************************************************************************
Diaplay all quads for selection
******************************************************************************/

void display_quads(GLenum mode, Polyhedron* this_poly)
{
	unsigned int i, j;
	GLfloat mat_diffuse[4];

	glEnable(GL_POLYGON_OFFSET_FILL);
	glPolygonOffset(1., 1.);
	glEnable(GL_DEPTH_TEST);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glShadeModel(GL_SMOOTH);
	glDisable(GL_LIGHTING);

	for (i = 0; i < this_poly->nquads; i++) {
		if (mode == GL_SELECT)
			glLoadName(i + 1);

		Quad* temp_q = this_poly->qlist[i];
		
		glBegin(GL_POLYGON);
		for (j = 0; j < 4; j++) {
			Vertex* temp_v = temp_q->verts[j];
			glVertex3d(temp_v->x, temp_v->y, temp_v->z);
		}
		glEnd();
	}
}

/******************************************************************************
Diaplay all vertices for selection
******************************************************************************/

void display_vertices(GLenum mode, Polyhedron* this_poly)
{
	for (int i = 0; i < this_poly->nverts; i++) {
		if (mode == GL_SELECT)
			glLoadName(i + 1);

		CHECK_GL_ERROR();

		Vertex* temp_v = this_poly->vlist[i];
		drawDot(temp_v->x, temp_v->y, temp_v->z, 0.15);
	}
	CHECK_GL_ERROR();
}

/******************************************************************************
Diaplay selected quad
******************************************************************************/

void display_selected_quad(Polyhedron* this_poly)
{
	if (this_poly->selected_quad == -1)
	{
		return;
	}

	unsigned int i, j;

	glDisable(GL_POLYGON_OFFSET_FILL);
	glEnable(GL_DEPTH_TEST);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glShadeModel(GL_SMOOTH);
	glDisable(GL_LIGHTING);

	Quad* temp_q = this_poly->qlist[this_poly->selected_quad];

	glBegin(GL_POLYGON);
	for (j = 0; j < 4; j++) {
		Vertex* temp_v = temp_q->verts[j];
		glColor3f(1.0, 0.0, 0.0);
		glVertex3d(temp_v->x, temp_v->y, temp_v->z);
	}
	glEnd();
}

/******************************************************************************
Diaplay selected vertex
******************************************************************************/

void display_selected_vertex(Polyhedron* this_poly)
{
	if (this_poly->selected_vertex == -1)
	{
		return;
	}

	Vertex* temp_v = this_poly->vlist[this_poly->selected_vertex];
	drawDot(temp_v->x, temp_v->y, temp_v->z, 0.15, 1.0, 0.0,0.0);

	CHECK_GL_ERROR();
}

/******************************************************************************
Process a keyboard action.  In particular, exit the program when an
"escape" is pressed in the window.
******************************************************************************/

void keyboard(unsigned char key, int x, int y) {
	int i;

	// clear out lines and points
	lines.clear();
	points.clear();


	switch (key) {
	case 27:	// set excape key to exit program
		poly->finalize();  // finalize_everything
		exit(0);
		break;

	case '1':	// solid color display with lighting
		display_mode = 1;
		glutPostRedisplay();
		break;

	case '2':	// wireframe display
		display_mode = 2;
		glutPostRedisplay();
		break;

	case '3':	// checkerboard display
	{
		display_mode = 3;

		double L = (poly->radius * 2) / 30;
		for (int i = 0; i < poly->nquads; i++) {
			Quad* temp_q = poly->qlist[i];
			for (int j = 0; j < 4; j++) {

				Vertex* temp_v = temp_q->verts[j];

				temp_v->R = int(temp_v->x / L) % 2 == 0 ? 1 : 0;
				temp_v->G = int(temp_v->y / L) % 2 == 0 ? 1 : 0;
				temp_v->B = 0.0;
			}
		}
		glutPostRedisplay();
	}
	break;

	case '4':	// Drawing points and lines created by the dots_and_lines_example() function
		display_mode = 4;
		dots_and_lines_example(&points, &lines);
		glutPostRedisplay();
		break;

	case '5':	// IBFV vector field display
		display_mode = 5;
		std::cout << "IBFV -\tORIGINAL: " << original_image << " FLOW: " << flow_image << std::endl;
		// initIBFV();
		reshape(win_width, win_height);
		glutPostRedisplay();
		break;
	case '6':{
		display_mode = 5;
		std::cout << "Img -\tORIGINAL: " << original_image << " FLOW: " << flow_image << std::endl;
		makePatternsImg(images[image_idx]);
		glutPostRedisplay();
		break;
	}
	case '7':{
		display_mode = 5;
		std::cout << "EDGES -\tORIGINAL: " << original_image << " FLOW: " << flow_image << std::endl;
		makePatternsImgEdges(images[image_idx]);
		glutPostRedisplay();
		break;
	}
	case '8':{
		// select_stream = true;
		// display_mode = 1;
		std::cout << "Streamlines" << std::endl;
		singularities.clear();
		POLYLINE line;
		for (int i=-5; i < 5; ++i){
			line.clear();
			streamline(line, icVector3(i*1.0, 0, 1), STEPSIZE);
			line.rgb = icVector3(1, 0, 0);
			polylines.push_back(line);

			line.clear();
			streamline(line, icVector3(0, i*1, 1), STEPSIZE);
			line.rgb = icVector3(1, 0, 0);
			polylines.push_back(line);

			line.clear();
			streamline(line, icVector3(i*1, i*1, 0), STEPSIZE);
			line.rgb = icVector3(1, 0, 0);
			polylines.push_back(line);

			line.clear();
			streamline(line, icVector3(i*1, i*-1, 0), STEPSIZE);
			line.rgb = icVector3(1, 0, 0);
			polylines.push_back(line);
		}
		glutPostRedisplay();
		break;
	}
	case '9':{
		std::cout << "Singularity" << std::endl;
		// singularity_IBVF = false;
		singularities.clear();
		extractSingularity();
		classifySingularity();
		glutPostRedisplay();
		break;
	}
	case '-':{
		std::cout << "Streamline based on seeding template" << std::endl;
		// singularity_IBVF = false;
		singularities.clear();
		extractSingularity();
		classifySingularity();
		
		float OFFSET = 1.3;

		for (const auto& singularity : singularities){
			// tuple<position, color>
			std::vector<std::tuple<icVector3, icVector3>> seeds;
			seedingTemplates(singularity, seeds, OFFSET);
			POLYLINE line;
			for (auto& seed : seeds){
				line.clear();
				streamline(line, std::get<0>(seed), STEPSIZE);
				line.rgb = std::get<1>(seed);
				polylines.push_back(line);
			}
		}

		glutPostRedisplay();
		break;
	}
	case ' ':{
		image_idx++;
		if (image_idx >= NUM_IMAGES)
			image_idx = 0;
		std::cout << "Switching to image: " << images[image_idx] << std::endl;
		
		singularities.clear();
		polylines.clear();
		glutPostRedisplay();
		break;
	}
	case '`':{
		// Switch vector fields
		vec_idx++;
		if (vec_idx >= NUM_VEC_FILES)
			vec_idx = 0;
		
		std::cout << "Switching vector file: " << vectors[vec_idx] << std::endl;
		FILE* this_file = fopen(vectors[vec_idx], "r");
		delete poly;
		poly = new Polyhedron(this_file);
		fclose(this_file);
		poly->initialize();
		poly->write_info();

		singularities.clear();
		polylines.clear();
		glutPostRedisplay();
		break;
	}

	case 'o': {
		display_mode = 5;
		original_image = !original_image;
		std::cout << "O -\tORIGINAL: " << original_image << " FLOW: " << flow_image << std::endl;
		glutPostRedisplay();
		break;
	}
	case 'p': {
		display_mode = 5;
		flow_image = !flow_image;
		std::cout << "P -\tORIGINAL: " << original_image << " FLOW: " << flow_image << std::endl;
		glutPostRedisplay();
		break;
	}
	case 'q': {
		display_mode = 1;
		select_stream = !select_stream;
		std::cout << "Select stream: " << (select_stream? "on\n" : "off\n");
		break;
	}
	case 'w': {
		// Entre streamline seeds manually
		double x= 0, y = 0;
		std::cout << "Enter x coordinate of streamline seed: ";
		std::cin >> x;
		std::cout << "\nEnter y coordinate of streamline seed: ";
		std::cin >> y;

		std::cout << "CREATING STREAMLINE at " << x << ", " << y << std::endl;
		display_mode = 1;

		POLYLINE line;
		streamline(line, icVector3(x, y, 0), STEPSIZE);
		line.rgb = icVector3(0, 0, 1);
		polylines.push_back(line);
	
		glutPostRedisplay();

		break;
	}
	case 'z': {
		std::cout << "clear\n";
		polylines.clear();
		singularities.clear();
		glutPostRedisplay();
		glutPostRedisplay();
		break;
	}
	case 'r':	// reset rotation and transformation
		mat_ident(rotmat);
		translation[0] = 0;
		translation[1] = 0;
		zoom = 1.0;
		glutPostRedisplay();
		break;
	}
	
}

/******************************************************************************
Callback function for dragging mouse
******************************************************************************/

void motion(int x, int y) {
	float r[4];
	float s, t;

	s = (2.0 * x - win_width) / win_width;
	t = (2.0 * (win_height - y) - win_height) / win_height;

	if ((s == s_old) && (t == t_old))
		return;

	switch (mouse_mode) {
	case 2:

		Quaternion rvec;

		mat_to_quat(rotmat, rvec);
		trackball(r, s_old, t_old, s, t);
		add_quats(r, rvec, rvec);
		quat_to_mat(rvec, rotmat);

		s_old = s;
		t_old = t;

		display();
		break;

	case 1:

		translation[0] += (s - s_old);
		translation[1] += (t - t_old);

		s_old = s;
		t_old = t;

		display();
		break;
	}
}

/******************************************************************************
Callback function for mouse clicks
******************************************************************************/

void mouse(int button, int state, int x, int y) {

	int key = glutGetModifiers();

	if (button == GLUT_LEFT_BUTTON || button == GLUT_RIGHT_BUTTON) {
		
		if (state == GLUT_DOWN) {
			float xsize = (float)win_width;
			float ysize = (float)win_height;

			float s = (2.0 * x - win_width) / win_width;
			float t = (2.0 * (win_height - y) - win_height) / win_height;

			s_old = s;
			t_old = t;

			/*translate*/
			if (button == GLUT_LEFT_BUTTON)
			{
				mouse_mode = 1;
			}

			/*rotate*/
			if (button == GLUT_RIGHT_BUTTON)
			{
				mouse_mode = 2;
			}
		}
		else if (state == GLUT_UP) {

			if (button == GLUT_LEFT_BUTTON && key == GLUT_ACTIVE_SHIFT) {  // build up the selection feedback mode

				/*select face*/

				GLuint selectBuf[512];
				GLint hits;
				GLint viewport[4];

				glGetIntegerv(GL_VIEWPORT, viewport);

				glSelectBuffer(win_width, selectBuf);
				(void)glRenderMode(GL_SELECT);

				glInitNames();
				glPushName(0);

				glMatrixMode(GL_PROJECTION);
				glPushMatrix();
				glLoadIdentity();

				/*create 5x5 pixel picking region near cursor location */
				gluPickMatrix((GLdouble)x, (GLdouble)(viewport[3] - y), 1.0, 1.0, viewport);

				set_view(GL_SELECT);
				set_scene(GL_SELECT, poly);
				display_quads(GL_SELECT, poly);

				glMatrixMode(GL_PROJECTION);
				glPopMatrix();
				glFlush();

				glMatrixMode(GL_MODELVIEW);

				hits = glRenderMode(GL_RENDER);
				poly->selected_quad = processHits(hits, selectBuf);
				printf("Selected quad id = %d\n", poly->selected_quad);

				if (select_stream && poly->selected_quad != -1){
					Quad* q = poly->qlist[poly->selected_quad];
					Vertex* v = q->verts[0];
					std::cout << "Vx: " << v->x << " Vy: " << v->y << " Vz: " << v->z << std::endl;
					POLYLINE line;
					streamline(line, icVector3(v->x, v->y, v->z), STEPSIZE);
					line.rgb = icVector3(0, 0, 1);
					polylines.push_back(line);
				}
				glutPostRedisplay();

				CHECK_GL_ERROR();

			}
			else if (button == GLUT_LEFT_BUTTON && key == GLUT_ACTIVE_CTRL)
			{
				/*select vertex*/

				GLuint selectBuf[512];
				GLint hits;
				GLint viewport[4];

				glGetIntegerv(GL_VIEWPORT, viewport);

				glSelectBuffer(win_width, selectBuf);
				(void)glRenderMode(GL_SELECT);

				glInitNames();
				glPushName(0);

				glMatrixMode(GL_PROJECTION);
				glPushMatrix();
				glLoadIdentity();

				/*  create 5x5 pixel picking region near cursor location */
				gluPickMatrix((GLdouble)x, (GLdouble)(viewport[3] - y), 1.0, 1.0, viewport);

				set_view(GL_SELECT);
				set_scene(GL_SELECT, poly);
				display_vertices(GL_SELECT, poly);

				glMatrixMode(GL_PROJECTION);
				glPopMatrix();
				glFlush();

				glMatrixMode(GL_MODELVIEW);

				hits = glRenderMode(GL_RENDER);
				poly->selected_vertex = processHits(hits, selectBuf);
				printf("Selected vert id = %d\n", poly->selected_vertex);

				if (select_stream && poly->selected_vertex != -1){
					Vertex* v = poly->vlist[poly->selected_vertex];
					std::cout << "Vx: " << v->x << " Vy: " << v->y << " Vz: " << v->z << std::endl;
					POLYLINE line;
					streamline(line, icVector3(v->x, v->y, v->z), STEPSIZE);
					line.rgb = icVector3(0, 0, 1);
					polylines.push_back(line);
				}

				glutPostRedisplay();

			}

			mouse_mode = -1;
		}
	}
}

/******************************************************************************
Callback function for mouse wheel scroll
******************************************************************************/

void mousewheel(int wheel, int direction, int x, int y) {
	if (direction == 1) {
		zoom *= zoomspeed;
		glutPostRedisplay();
	}
	else if (direction == -1) {
		zoom /= zoomspeed;
		glutPostRedisplay();
	}
}

/******************************************************************************
Callback function for window reshaping
******************************************************************************/

void reshape(int width, int height)
{
	win_width = width;
	win_height = height;

	aspectRatio = (float)width / (float)height;

	glViewport(0, 0, width, height);

	set_view(GL_RENDER);

	// reset IBFV pixels buffer
	free(pixels);
	initIBFV();
}

/******************************************************************************
Make patterns Image (used for Project 3)
******************************************************************************/
void makePatternsImg(const std::string& fname){
    ppm img(fname);
    GLubyte pat[NPN][NPN][4];
    int i, j;
    for (i = 0; i < NPN; ++i){
        for (j = 0; j < NPN; ++j){
            pat[i][j][0] = img.r[(NPN - i - 1) * NPN + j]; 
            pat[i][j][1] = img.g[(NPN - i - 1) * NPN + j];
            pat[i][j][2] = img.b[(NPN - i - 1) * NPN + j];
			// NOTE: CHANGE
            // pat[i][j][0] = img.r[1024*(NPN - i - 1) * 2 + j]; 
            // pat[i][j][1] = img.g[1024*(NPN - i - 1) * 2 + j];
            // pat[i][j][2] = img.b[1024*(NPN - i - 1) * 2 + j];

            pat[i][j][3] = alpha;
        }
    }

    glNewList(1, GL_COMPILE);
    glTexImage2D(GL_TEXTURE_2D, 0, 4, NPN, NPN, 0, GL_RGBA, GL_UNSIGNED_BYTE, pat);
    glEndList();
}

/******************************************************************************
Make patterns Image noise (used for Project 3)
******************************************************************************/
void makePatternsImgNoise(const std::string& fname, float w){
    ppm img(fname);
    int lut[256];
    int phase[NPN][NPN];
    int i, j, t;
    GLubyte pat[NPN][NPN][4];

    for (i=0; i < 256; ++i) lut[i] = i < 127 ? 0 : 255;
    for (i=0; i < NPN; ++i)
        for (j=0; j < NPN; ++j)
            phase[i][j] = rand() %256;
    for (i=0; i < NPN; ++i)
        for (j=0; j < NPN; ++j){
            pat[i][j][0] = img.r[(NPN - i - 1) * NPN + j] * w + lut[(phase[i][j]) % 255] * (1-w); 
            pat[i][j][1] = img.g[(NPN - i - 1) * NPN + j] * w + lut[(phase[i][j]) % 255] * (1-w);
            pat[i][j][2] = img.b[(NPN - i - 1) * NPN + j] * w + lut[(phase[i][j]) % 255] * (1-w);
            pat[i][j][3] = alpha;
        }
    glNewList(1, GL_COMPILE);
    glTexImage2D(GL_TEXTURE_2D, 0, 4, NPN, NPN, 0, GL_RGBA, GL_UNSIGNED_BYTE, pat);
    glEndList();
}

/******************************************************************************
Make patterns from edges - sobel field edge(used for Project 3)
******************************************************************************/
void makePatternsImgEdges(const std::string& fname){
    float kernelx[3][3] = {
        {-1, 0, 1},
        {-2, 0, 2},
        {-1, 0, 1}
    };
    float kernely[3][3] = {
        {-1, -2, 1},
        {0, 0, 0},
        {1, 2, 1}
    };

    ppm img(fname);
    GLubyte pat[NPN][NPN][4];
    GLubyte patO[NPN][NPN][4];

    int i, j;
    for (i=0; i < NPN; ++i){
        for (j=0; j < NPN; ++j){
            float c = 0.299 * img.r[(NPN - i - 1) * NPN + j] / 255 + 
                0.587 * img.g[(NPN - i - 1) * NPN + j] / 255 + 
                0.114 * img.b[(NPN - i - 1) * NPN + j] / 255;
            patO[i][j][0] = c * 255;
            patO[i][j][1] = c * 255;
            patO[i][j][2] = c * 255;
            patO[i][j][3] = alpha;

        }
    }

    // filter
    for (i=1; i < NPN - 1; ++i){
        for (j=1; j < NPN - 1; ++j){
            float mag0x = 0;
            float mag1x = 0;
            float mag2x = 0;
            float mag0y = 0;
            float mag1y = 0;
            float mag2y = 0;

            for (int a=0; a < 3; ++a){
                for (int b = 0; b < 3; ++b){
                    mag0x += patO[i - 1 + a][j - 1 + b][0] * kernelx[a][b];
                    mag1x += patO[i - 1 + a][j - 1 + b][1] * kernelx[a][b];
                    mag2x += patO[i - 1 + a][j - 1 + b][2] * kernelx[a][b];
                    
                    mag0y += patO[i - 1 + a][j - 1 + b][0] * kernelx[a][b];
                    mag1y += patO[i - 1 + a][j - 1 + b][1] * kernelx[a][b];
                    mag2y += patO[i - 1 + a][j - 1 + b][2] * kernelx[a][b];
                }
            }

            float v0 = std::sqrt(mag0x*mag0x + mag0y * mag0y);
            float v1 = std::sqrt(mag1x*mag1x + mag1y * mag1y);
            float v2 = std::sqrt(mag2x*mag2x + mag2y * mag2y);
            if (v0 > 255) v0 = 255;
            if (v0 < 0) v0 = 0;
            if (v1 > 255) v1 = 255;
            if (v1 < 0) v1 = 0;
            if (v2 > 255) v2 = 255;
            if (v2 < 0) v2 = 0;
            pat[i][j][0] = v0;
            pat[i][j][1] = v1;
            pat[i][j][2] = v2;
            pat[i][j][3] = alpha;
        }
    }
    glNewList(1, GL_COMPILE);
    glTexImage2D(GL_TEXTURE_2D, 0, 4, NPN, NPN, 0, GL_RGBA, GL_UNSIGNED_BYTE, pat);
    glEndList();
}


/******************************************************************************
Display IBFV vector field visualization (used for Project 3)
******************************************************************************/

void displayIBFV()
{
	glDisable(GL_LIGHTING);
	glDisable(GL_LIGHT0);
	glDisable(GL_LIGHT1);
	glDisable(GL_BLEND);

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);

	glEnable(GL_TEXTURE_2D);
	glShadeModel(GL_FLAT);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    glClearColor(0.5, 0.5, 0.5, 1.0); //background for rendering color coding and lighting
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, win_width, win_height, 0, GL_RGB, GL_UNSIGNED_BYTE, pixels);

    double modelview_matrix[16], projection_matrix[16];
    int viewport[4];
    glGetDoublev(GL_MODELVIEW_MATRIX, modelview_matrix);
    glGetDoublev(GL_PROJECTION_MATRIX, projection_matrix);
    glGetIntegerv(GL_VIEWPORT, viewport);

    for (int i=0; i < poly->nquads; ++i){
        Quad* qtemp = poly->qlist[i];

        glBegin(GL_QUADS);
        for (int j=0; j < 4; ++j){
            Vertex* vtemp = qtemp->verts[j];

            double tx, ty, dummy;
            gluProject((GLdouble)vtemp->x, (GLdouble)vtemp->y, (GLdouble)vtemp->z, 
                        modelview_matrix, projection_matrix, viewport, &tx, &ty, &dummy);

            tx = tx / win_width;
            ty = ty / win_height;

            icVector2 dp;
            // key press
            if (flow_image){
                dp = icVector2(((double)rand() / (double)(RAND_MAX)), ((double)rand() / (double)(RAND_MAX)));
                normalize(dp);
                dp *= dmax;
            }
            else{
                dp = icVector2(vtemp->vx, vtemp->vy);
                normalize(dp);
                dp *= dmax;
				dp.x *= -1;
				dp.y *= -1;
            }

            float px = tx + dp.x;
            float py = ty + dp.y;

            if (original_image){
                glTexCoord2f(tx, ty);
            }else{
                glTexCoord2f(px, py);
            }
            glVertex3d(vtemp->x, vtemp->y, vtemp->z );
        }
        glEnd();
    }

    glEnable(GL_BLEND);

    // blend in noise pattern
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();

	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();

	glTranslatef(-1.0, -1.0, 0.0);
	glScalef(2.0, 2.0, 1.0);

	glCallList(1);

	glBegin(GL_QUAD_STRIP);

	glTexCoord2f(0.0, 0.0);  glVertex2f(0.0, 0.0);
	glTexCoord2f(0.0, tmax); glVertex2f(0.0, 1.0);
	glTexCoord2f(tmax, 0.0);  glVertex2f(1.0, 0.0);
	glTexCoord2f(tmax, tmax); glVertex2f(1.0, 1.0);
	glEnd();
	glDisable(GL_BLEND);

	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();

	glMatrixMode(GL_PROJECTION);
	glPopMatrix();

	glReadPixels(0, 0, win_width, win_height, GL_RGB, GL_UNSIGNED_BYTE, pixels);

    // don't advect texture coords
    // draw mesh using pixels
    glClearColor(1, 1, 1, 1);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, win_width, win_height, 0, GL_RGB, GL_UNSIGNED_BYTE, pixels);
    for (int i=0; i < poly->nquads; ++i){
        Quad* qtemp = poly->qlist[i];
        glBegin(GL_QUADS);
        for (int j=0; j < 4; ++j){ 
            Vertex* vtemp = qtemp->verts[j];
            double tx, ty, dummy;
            gluProject((GLdouble)vtemp->x, (GLdouble)vtemp->y, (GLdouble)vtemp->z, 
                            modelview_matrix, projection_matrix, viewport, &tx, &ty, &dummy);
            tx = tx / win_width;
            ty = ty / win_height;
            glTexCoord2f(tx, ty);
            glVertex3d(vtemp->x, vtemp->y, vtemp->z);
        }
        glEnd();
    }   

    glDisable(GL_TEXTURE_2D);
    glDisable(GL_BLEND);
    glShadeModel(GL_SMOOTH);
}

/******************************************************************************
Callback function for scene display
******************************************************************************/

void display(void)
{
	glClearColor(1.0, 1.0, 1.0, 1.0);  // background for rendering color coding and lighting

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	set_view(GL_RENDER);
	set_scene(GL_RENDER, poly);

	/*display the mesh*/
	display_polyhedron(poly);

	/*display selected elements*/
	display_selected_vertex(poly);
	display_selected_quad(poly);

	display_poliline(polylines);
	// display_singularities();
	drawSingularities();

	glFlush();
	glutSwapBuffers();
	glFinish();

	CHECK_GL_ERROR();
}

/******************************************************************************
Diaplay the polygon with visualization results
******************************************************************************/

void display_polyhedron(Polyhedron* poly)
{
	unsigned int i, j;

	glEnable(GL_POLYGON_OFFSET_FILL);
	glPolygonOffset(1., 1.);

	glEnable(GL_DEPTH_TEST);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glShadeModel(GL_SMOOTH);
	switch (display_mode)
	{
	case 1:	// solid color display with lighting
	{
		select_stream = true;
		glEnable(GL_LIGHTING);
		glEnable(GL_LIGHT0);
		glEnable(GL_LIGHT1);

		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		GLfloat mat_diffuse[4] = { 0.24, 0.4, 0.47, 0.0 };
		GLfloat mat_specular[] = { 1.0, 1.0, 1.0, 1.0 };
		glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
		glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
		glMaterialf(GL_FRONT, GL_SHININESS, 50.0);

		for (int i = 0; i < poly->nquads; i++) {
			Quad* temp_q = poly->qlist[i];
			glBegin(GL_POLYGON);
			for (int j = 0; j < 4; j++) {
				Vertex* temp_v = temp_q->verts[j];
				glNormal3d(temp_v->normal.entry[0], temp_v->normal.entry[1], temp_v->normal.entry[2]);
				glVertex3d(temp_v->x, temp_v->y, temp_v->z);
			}
			glEnd();
		}
	}
	break;

	case 2:	// wireframe display
	{
		glDisable(GL_LIGHTING);
		glEnable(GL_LINE_SMOOTH);
		glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		glLineWidth(1.0);

		for (int i = 0; i < poly->nquads; i++) {
			Quad* temp_q = poly->qlist[i];

			glBegin(GL_POLYGON);
			for (int j = 0; j < 4; j++) {
				Vertex* temp_v = temp_q->verts[j];
				glNormal3d(temp_q->normal.entry[0], temp_q->normal.entry[1], temp_q->normal.entry[2]);
				glColor3f(0.0, 0.0, 0.0);
				glVertex3d(temp_v->x, temp_v->y, temp_v->z);
			}
			glEnd();
		}

		glDisable(GL_BLEND);
	}
	break;

	case 3:	// checkerboard pattern display
	{
		glDisable(GL_LIGHTING);
		for (int i = 0; i < poly->nquads; i++) {
			Quad* temp_q = poly->qlist[i];
			glBegin(GL_POLYGON);
			for (int j = 0; j < 4; j++) {
				Vertex* temp_v = temp_q->verts[j];
				glColor3f(temp_v->R, temp_v->G, temp_v->B);
				glVertex3d(temp_v->x, temp_v->y, temp_v->z);
			}
			glEnd();
		}
	}
	break;

	case 4: // points and lines drawing example
	{
		glEnable(GL_LIGHTING);
		glEnable(GL_LIGHT0);
		glEnable(GL_LIGHT1);

		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		GLfloat mat_diffuse[4] = { 0.24, 0.4, 0.47, 0.0 };
		GLfloat mat_specular[] = { 1.0, 1.0, 1.0, 1.0 };
		glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
		glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
		glMaterialf(GL_FRONT, GL_SHININESS, 50.0);

		for (int i = 0; i < poly->nquads; i++) {
			Quad* temp_q = poly->qlist[i];
			glBegin(GL_POLYGON);
			for (int j = 0; j < 4; j++) {
				Vertex* temp_v = temp_q->verts[j];
				glNormal3d(temp_v->normal.entry[0], temp_v->normal.entry[1], temp_v->normal.entry[2]);
				glVertex3d(temp_v->x, temp_v->y, temp_v->z);
			}
			glEnd();
		}

		// draw lines
		for (int k = 0; k < lines.size(); k++)
		{
			drawPolyLine(lines[k], 1.0, 1.0, 0.0, 0.0);
		}

		// draw points
		for (int k = 0; k < points.size(); k++)
		{
			icVector3 point = points[k];
			drawDot(point.x, point.y, point.z);
		}
		break;
	}
	break;

	case 5:	// IBFV vector field display
	{
		select_stream = false;
		displayIBFV();
		glutPostRedisplay();
	}
	break;

	case 6: // add your own display mode
	{
		
	}
	break;

	default:
	{
		// don't draw anything
	}

	}
}
