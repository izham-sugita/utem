/*********************************************************************
 * glut_utl.c 
 *
 * Copyright (C) 2001 AzLib Developers Group
 *
 * Written by
 *  <Kenzen TAKEUCHI> <Mitsunori HIROSE> <Ryu SASAOKA>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: glut_utl.c,v 1.5 2003/07/22 08:02:11 sasaoka Exp $ */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <GL/gl.h>
#include <GL/glut.h>
#include "rc.h"
#include "fem_struct.h"
#include "glut_utl.h"

static void initialize(void);
static void initialize_light(void);
static void display(void);
static void drawObject(FEM_ELEMENT_ARRAY elem, FEM_NODE_ARRAY node);
static void reshape(int w, int h);
static void keyboard(unsigned char key, int x, int y);
static void sp_keyboard(int key, int x, int y);
static void mouse(int button, int state, int x, int y);
static void motion(int x, int y);
static void idle(void);
static void popupMenu_1(int);
static void popupMenu_2(int);
static void popupMenu_3(int);
static void drow_fem_model(FEM_ELEMENT_ARRAY elem, FEM_NODE_ARRAY node); 
static RC initialize_object(FEM_NODE_ARRAY node);
RC show_fem_model(FEM_ELEMENT_ARRAY elem, FEM_NODE_ARRAY node);

static GLdouble spin_x = 0.0; 
static GLdouble spin_y = 0.0; 
static GLint dragFlag;
static GLint x_mem1 = 0;
static GLint y_mem1 = 0;
static GLint x_mem2 = 0;
static GLint y_mem2 = 0;
static GLdouble view_size = 0.0;
static TRANSLATION3D max_object_size = {0.0, 0.0, 0.0};
static TRANSLATION3D min_object_size = {0.0, 0.0, 0.0};

/* 光源座標 */
static GLfloat light0_position[] = {0.0, 0.0, 0.0, 1.0};
/* 環境光 */
static GLfloat light0_ambient[] = {0.6, 0.6, 0.6, 1.0};
/* 参照光 */
static GLfloat light0_diffuse[] = {0.7, 0.7, 0.7, 1.0};
/* 鏡面光 */
static GLfloat light0_specular[] = {1.0, 1.0, 1.0, 1.0};
/* 環境光の反射の設定 */
static GLfloat material_ambient[8][4] = {{R(255), G(255), B(255), 1.0},
										{R(255), G(0), B(0), 1.0},
										{R(0), G(255), B(0), 1.0},
										{R(0), G(0), B(255), 1.0},
										{R(0), G(255), B(255), 1.0},
										{R(255), G(0), B(255), 1.0},
										{R(255), G(255), B(0), 1.0}};
/* 拡散光の反射の設定 */
static GLfloat material_diffuse[8][4] = {{R(0), G(0), B(0), 1.0},
										{R(0), G(0), B(0), 1.0},
										{R(0), G(0), B(0), 1.0},
										{R(0), G(0), B(0), 1.0},
										{R(0), G(0), B(0), 1.0},
										{R(0), G(0), B(0), 1.0},
										{R(0), G(0), B(0), 1.0}};
										
/* 鏡面光の反射の設定 */
static GLfloat material_specular[8][4] = {{R(0), G(0), B(0), 1.0},
										{R(0), G(0), B(0), 1.0},
										{R(0), G(0), B(0), 1.0},
										{R(0), G(0), B(0), 1.0},
										{R(0), G(0), B(0), 1.0},
										{R(0), G(0), B(0), 1.0},
										{R(0), G(0), B(0), 1.0}};
/* ハイライトの大きさの設定 */
/* 0.0-128.0 の値を取り、大きな値ほどハイライトが小さく表示される */
static GLfloat material_shininess[] = {H(64), H(64), H(64), H(64), 
									H(64), H(64), H(64), H(64)};



static GLint material = MATERIAL_CYAN; 
static FEM_NODE_ARRAY node;
static FEM_ELEMENT_ARRAY elem;
static GLint drowtype = GL_LINE_LOOP; /* デフォルト描画タイプ：POLYGON */

RC show_fem_model(FEM_ELEMENT_ARRAY elem_org, FEM_NODE_ARRAY node_org)
{
	int submenu_1, submenu_2;

	node = node_org;
	elem = elem_org;

	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
	glutInitWindowSize(WINDOW_WIDTH, WINDOW_HEIGHT);
	glutInitWindowPosition(WINDOW_POSITION_X, WINDOW_POSITION_Y);
	glutCreateWindow("MesaGL FEM MODEL VIEWER");

	initialize_light();
	initialize();
	glutKeyboardFunc(keyboard);
	glutSpecialFunc(sp_keyboard);
	glutMouseFunc(mouse);
	glutMotionFunc(motion);
	glutIdleFunc(idle);
	glutDisplayFunc(display);
	glutReshapeFunc(reshape); 

	submenu_1 = glutCreateMenu(popupMenu_1);
	glutAddMenuEntry("Polygon", POLYGON);
	glutAddMenuEntry("Line", LINE);
	glutAddMenuEntry("Points", POINT);

	submenu_2 = glutCreateMenu(popupMenu_2);
	glutAddMenuEntry("White", MATERIAL_WHITE);
	glutAddMenuEntry("Red", MATERIAL_RED);
	glutAddMenuEntry("Green", MATERIAL_GREEN);
	glutAddMenuEntry("Blue", MATERIAL_BLUE);
	glutAddMenuEntry("Cyan", MATERIAL_CYAN);
	glutAddMenuEntry("Magenta", MATERIAL_MAGENTA);
	glutAddMenuEntry("Yellow", MATERIAL_YELLOW);
	
	glutCreateMenu(popupMenu_3);
	glutAddSubMenu( "Drow type", submenu_1);
	glutAddSubMenu( "Color", submenu_2);
	glutAddMenuEntry("ReSet", RESET);
	glutAddMenuEntry("Quit", QUIT);
	glutAttachMenu(GLUT_RIGHT_BUTTON);

	glutMainLoop(); 

	return(NORMAL_RC);
}

	
static void drow_fem_model(FEM_ELEMENT_ARRAY elem, FEM_NODE_ARRAY node)
{
	int ii1,ii2;
	int index;
	
	if(dragFlag == STATE_DOWN){
		glBegin(GL_LINE_LOOP);
			glVertex3d(max_object_size.x, max_object_size.y, min_object_size.z);
			glVertex3d(max_object_size.x, max_object_size.y, max_object_size.z);
			glVertex3d(max_object_size.x, min_object_size.y, max_object_size.z);
			glVertex3d(max_object_size.x, min_object_size.y, min_object_size.z);
		glEnd();
		glBegin(GL_LINE_LOOP);
			glVertex3d(min_object_size.x, max_object_size.y, min_object_size.z);
			glVertex3d(min_object_size.x, max_object_size.y, max_object_size.z);
			glVertex3d(min_object_size.x, min_object_size.y, max_object_size.z);
			glVertex3d(min_object_size.x, min_object_size.y, min_object_size.z);
		glEnd();
		glBegin(GL_LINE_LOOP);
			glVertex3d(min_object_size.x, max_object_size.y, min_object_size.z);
			glVertex3d(min_object_size.x, max_object_size.y, max_object_size.z);
			glVertex3d(max_object_size.x, max_object_size.y, max_object_size.z);
			glVertex3d(max_object_size.x, max_object_size.y, min_object_size.z);
		glEnd();
		glBegin(GL_LINE_LOOP);
			glVertex3d(min_object_size.x, min_object_size.y, min_object_size.z);
			glVertex3d(min_object_size.x, min_object_size.y, max_object_size.z);
			glVertex3d(max_object_size.x, min_object_size.y, max_object_size.z);
			glVertex3d(max_object_size.x, min_object_size.y, min_object_size.z);
		glEnd();
	}else{ 
		for(ii1=0;ii1<elem.size;ii1++){
			if(elem.array[ii1].label < 0)continue;
			glBegin(drowtype);
			switch(elem.array[ii1].type){
			case ELEM_TRI1:
				for(ii2=0;ii2<3;ii2++){
					index = search_fem_node_label(node,
											elem.array[ii1].node[ii2]);
					if(index < 0){
						fprintf(stderr,
						        "search_fem_node_label[%d|%d] < ", ii1, ii2);
						exit(1);
					}
					glVertex3d( node.array[index].p.x,
								node.array[index].p.y,
								node.array[index].p.z);
				}
				break;
			case ELEM_QUAD1:
				for(ii2=0;ii2<4;ii2++){
					index = search_fem_node_label(node,
											elem.array[ii1].node[ii2]);
					if(index < 0){
						fprintf(stderr,
						        "search_fem_node_label[%d|%d] < ", ii1, ii2);
						exit(1);
					}
					glVertex3d( node.array[index].p.x,
								node.array[index].p.y,
								node.array[index].p.z);
				}
				break;
			default:
				break;
			}
			glEnd();
		}
	}
}
	
/* 対象となるオブジェクト中心を原点移動、視野範囲を設定 */
/* オブジェクトサイズを設定 */
static RC initialize_object(FEM_NODE_ARRAY node)
{
	int ii1;
	TRANSLATION3D center;
	double length;
	
	center.x = 0.0;
	center.y = 0.0;
	center.z = 0.0;

	for(ii1=0;ii1<node.size;ii1++){
		center.x += node.array[ii1].p.x;
		center.y += node.array[ii1].p.y;
		center.z += node.array[ii1].p.z;
	}

	center.x /= node.size;
	center.y /= node.size;
	center.z /= node.size;

	for(ii1=0;ii1<node.size;ii1++){
		node.array[ii1].p.x -= center.x;
		node.array[ii1].p.y -= center.y;
		node.array[ii1].p.z -= center.z;
	}

	for(ii1=0;ii1<node.size;ii1++){
		length = sqrt(node.array[ii1].p.x * node.array[ii1].p.x+
						node.array[ii1].p.y * node.array[ii1].p.y+
						node.array[ii1].p.z * node.array[ii1].p.z);
		if(view_size < length) view_size  = length;

		if(max_object_size.x < node.array[ii1].p.x)
			max_object_size.x = node.array[ii1].p.x;
		if(max_object_size.y < node.array[ii1].p.y)
			max_object_size.y = node.array[ii1].p.y;
		if(max_object_size.z < node.array[ii1].p.z)
			max_object_size.z = node.array[ii1].p.z;
		if(min_object_size.x > node.array[ii1].p.x)
			min_object_size.x = node.array[ii1].p.x;
		if(min_object_size.y > node.array[ii1].p.y)
			min_object_size.y = node.array[ii1].p.y;
		if(min_object_size.z > node.array[ii1].p.z)
			min_object_size.z = node.array[ii1].p.z;
	}

	return(NORMAL_RC);
	
}

static void initialize(void)
{

	dragFlag = STATE_UP;
	glClearColor(R(0), G(0), B(0), 0.0);
	glShadeModel(GL_SMOOTH); /* 一つの頂点の色でベタ塗り */
	glMatrixMode(GL_PROJECTION); /* 射影行列操作モード */
	glLoadIdentity(); /* 初期化 */
	print_rc(stderr, initialize_object(node));
	fprintf(stderr, "view_size = %f\n", 2*view_size);
	/*gluLookAt(0, 0, 20*view_size, 0, 0, 0, 0, 1, 0);  */
	glOrtho(-view_size, view_size,
			-view_size, view_size,
			-1.2*view_size, 1.2*view_size);
/*	gluPerspective(30.0, 1.0, 0, 1500); */
/*	glOrtho(-60, 60, -60, 60, -60, 150);*/
/*	gluPerspective(45.0, 1.0, 2.0, 20.0);*/
}

static void initialize_light(void)
{
	glLightfv(GL_LIGHT0, GL_AMBIENT, light0_ambient);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, light0_diffuse);
	glLightfv(GL_LIGHT0, GL_SPECULAR, light0_specular);
	glLightModelf(GL_LIGHT_MODEL_LOCAL_VIEWER, GL_TRUE);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glEnable(GL_NORMALIZE);
	glEnable(GL_DEPTH_TEST);
}

static void reshape(int w, int h)
{
	glViewport(0, 0, w, h);
}

static void display()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); 
	glMaterialfv(GL_FRONT, GL_AMBIENT, material_ambient[material]);
	glMaterialfv(GL_FRONT, GL_DIFFUSE, material_diffuse[material]);
	glMaterialfv(GL_FRONT, GL_SPECULAR, material_specular[material]);
	glMaterialf(GL_FRONT, GL_SHININESS, material_shininess[material]);
	drawObject(elem, node);
	glFlush();
	glutSwapBuffers();
}

static void drawObject(FEM_ELEMENT_ARRAY elem, FEM_NODE_ARRAY node)
{
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glPushMatrix();
	glLightfv(GL_LIGHT0, GL_POSITION, light0_position);
	glPopMatrix();
	glPushMatrix();
	glRotatef(spin_x, 0.0, 1.0, 0.0);
	glRotatef(spin_y, 1.0, 0.0, 0.0);
	drow_fem_model(elem, node);
	glPopMatrix();
}

static void popupMenu_1(int select)
{

	if(select == POINT) drowtype = GL_POINTS;
	if(select == LINE) drowtype = GL_LINE_LOOP;
	if(select == POLYGON) drowtype = GL_POLYGON;
	glutSwapBuffers();
	
}

static void popupMenu_2(int select)
{
		material = select;
		glutSwapBuffers();
}

static void popupMenu_3(int select)
{		
	if(select == QUIT){
		fprintf(stderr, "\n");
		exit(0);
	}
	if(select == RESET){
		spin_x =0.0;
		spin_y =0.0;
		glutPostRedisplay();
	}
		glutSwapBuffers();
}

static void keyboard(unsigned char key, int x,  int y)
{
	switch(key){
	case KEYBOARD_ESC:
	case KEYBOARD_Q:
		fprintf(stderr, "\n");
		exit(0);
	case KEYBOARD_C:
		if(material == MATERIAL_YELLOW){
			material = MATERIAL_WHITE;
			break;
		}else{
			material++;
			break;
		}
	case KEYBOARD_D:
		if(drowtype == GL_POLYGON){
			drowtype = GL_LINE_LOOP;
			break;
		}else if(drowtype == GL_LINE_LOOP){
			drowtype = GL_POINTS;
			break;
		}else{
			drowtype = GL_POLYGON;
			break;
		}
	case KEYBOARD_H:
		spin_x += 5;
		break;
	case KEYBOARD_J:
		spin_y -= 5;
		break;
	case KEYBOARD_K:
		spin_y += 5;
		break;
	case KEYBOARD_L:
		spin_x -= 5;
		break;
	case KEYBOARD__H:
		spin_x += 30;
		break;
	case KEYBOARD__J:
		spin_y -= 30;
		break;
	case KEYBOARD__K:
		spin_y += 30;
		break;
	case KEYBOARD__L:
		spin_x -= 30;
		break;
	default:
		break;
	}

	if(spin_x > 360.0) spin_x -= 360.0;
	else if(spin_x < -360.0) spin_x += 360.0;
	if(spin_y > 360.0) spin_y -= 360.0;
	else if(spin_y < -360.0) spin_y += 360.0;
	fprintf(stderr, " Xangle =%+6.1f°Yangle =%+6.1f°\r", spin_x, spin_y);
	glutPostRedisplay();
}

static void sp_keyboard(int key, int x,  int y)
{
	switch(key){
	case GLUT_KEY_LEFT:
		spin_x += 5;
		break;
	case GLUT_KEY_DOWN:
		spin_y -= 5;
		break;
	case GLUT_KEY_UP:
		spin_y += 5;
		break;
	case GLUT_KEY_RIGHT:
		spin_x -= 5;
		break;
	default:
		break;
	}

	if(spin_x > 360.0) spin_x -= 360.0;
	else if(spin_x < -360.0) spin_x += 360.0;
	if(spin_y > 360.0) spin_y -= 360.0;
	else if(spin_y < -360.0) spin_y += 360.0;
	fprintf(stderr, " Xangle =%+6.1f°Yangle =%+6.1f°\r", spin_x, spin_y);
	glutPostRedisplay();
}

static void mouse(int button, int state, int x, int y)
{
	if(button == GLUT_LEFT_BUTTON){
		if(state == GLUT_DOWN){
			x_mem2 = x_mem1;
			y_mem2 = y_mem1;
			x_mem1 = x;
			y_mem1 = y;
			dragFlag = STATE_DOWN;
		}else{
			dragFlag = STATE_UP;
		}
	}
}

static void motion(int x, int y)
{
	if(dragFlag == STATE_DOWN){
		if(x_mem1 != x)
			spin_x +=300.0 * 2.0 * (x_mem1-x) / WINDOW_WIDTH;
		if(y_mem1 != y)
			spin_y +=300.0 * 2.0 * (y_mem1-y) / WINDOW_HEIGHT;
		if(spin_x > 360.0) spin_x -= 360.0;
		else if(spin_x < -360.0) spin_x += 360.0;
		if(spin_y > 360.0) spin_y -= 360.0;
		else if(spin_y < -360.0) spin_y += 360.0;
		glutPostRedisplay();
		fprintf(stderr, " Xangle =%+6.1f°Yangle =%+6.1f°\r", spin_x, spin_y);
		x_mem2 = x_mem1; 
		y_mem2 = y_mem1; 
		x_mem1 = x;
		y_mem1 = y;
		glutPostRedisplay();
	}
}


static void idle(void)
{
	if(dragFlag == STATE_UP){
		spin_x += 300.0 * 2.0 * (x_mem2 - x_mem1) / WINDOW_WIDTH;
		spin_y += 300.0 * 2.0 * (y_mem2 - y_mem1) / WINDOW_HEIGHT;
		if(spin_x > 360.0) spin_x -= 360.0;
		else if(spin_x < -360.0) spin_x += 360.0;
		if(spin_y > 360.0) spin_y -= 360.0;
		else if(spin_y < -360.0) spin_y += 360.0;
		fprintf(stderr, " Xangle =%+6.1f°Yangle =%+6.1f°\r", spin_x, spin_y);
		glutPostRedisplay();
	}
}

