/************************************************************************
	 File:        Maze.h

	 Author:
				  Stephen Chenney, schenney@cs.wisc.edu
	 Modifier
				  Yu-Chi Lai, yu-chi@cs.wisc.edu

	 Comment:
						(c) 2001-2002 Stephen Chenney, University of Wisconsin at Madison

						Class header file for Maze class. Manages the maze.


	 Platform:    Visio Studio.Net 2003 (converted to 2005)

*************************************************************************/

#ifndef _MAZE_H_
#define _MAZE_H_

#include <FL/math.h> // Use FLTK's math header because it defines M_PI
#include "Cell.h"
#include "LineSeg.h"
#include <iostream>
#include <cmath>


//************************************************************************
//
// * A class for exceptions. Used by the constructor to pass file I/O errors
//   back.
//
//************************************************************************
class MazeException {
private:
	char* message;

public:
	MazeException(const char* m = "");
	~MazeException() { delete message; };

	// Return the error message string associated with the exception.
	const char* Message(void) { return message; };
};
class Frustum {
public:
	const int X = 0;
	const int Y = 1;
	Frustum() {};
	float userX, userY;
	LineSeg leftEdge, rightEdge;


	Frustum restricted(LineSeg e) {
		Frustum f;
		f.userX = this->userX;
		f.userY = this->userY;
		float angStr, angEnd;
		float leftX, leftY, rightX, rightY;
		angStr = atan2(e.start[Y] - userY, e.start[X] - userX);
		angEnd = atan2(e.end[Y] - userY, e.end[X] - userX);
		if (angStr >= angEnd) {
			leftX = e.start[X];
			leftY = e.start[Y];
			rightX = e.end[X];
			rightY = e.end[Y];
		}
		else {
			leftX = e.end[X];
			leftY = e.end[Y];
			rightX = e.start[X];
			rightY = e.start[Y];
		}
		f.leftEdge = LineSeg(userX, userY, leftX,leftY);
		f.rightEdge = LineSeg(userX, userY, rightX,rightY);
		return f;
	}
	bool pointInTriangle(float* p, float* t1, float* t2, float* t3) {
		float r1, r2, r3;
		bool pos, neg;
		r1 = (p[X] - t2[X]) * (t1[Y] - t2[Y]) - (t1[X] - t2[X]) * (p[Y] - t2[Y]);
		r2 = (p[X] - t3[X]) * (t2[Y] - t3[Y]) - (t2[X] - t3[X]) * (p[Y] - t3[Y]);
		r3 = (p[X] - t1[X]) * (t3[Y] - t1[Y]) - (t3[X] - t1[X]) * (p[Y] - t1[Y]);

		pos = (r1 > 0) || (r2 > 0) || (r3 > 0);
		neg = (r1 < 0) || (r2 < 0) || (r3 < 0);

		return !(pos && neg);
	}
	float adjustDir(float dir) {
		while (dir >= 360)
			dir -= 360;
		while (dir < 0)
			dir += 360;
		return dir;
	}
	LineSeg clip(Edge& e, float dir, float fov) {
		dir = adjustDir(dir);
		LineSeg opt;
		opt.color[0] = e.color[0];
		opt.color[1] = e.color[1];
		opt.color[2] = e.color[2];
		float userPos[] = { userX,userY };
		float leftIntersection[] = { (leftEdge.end[X] - userX) * 1000 + userX,(leftEdge.end[Y] - userY) * 1000 + userY };
		float rightIntersection[] = { (rightEdge.end[X] - userX) * 1000 + userX,(rightEdge.end[Y] - userY) * 1000 + userY };
		bool startInTriangle = pointInTriangle(e.endpoints[Edge::START]->posn, userPos, leftIntersection, rightIntersection);
		bool endInTriangle = pointInTriangle(e.endpoints[Edge::END]->posn, userPos, leftIntersection, rightIntersection);

		float deltaX = this->leftEdge.end[X] - this->leftEdge.start[X];
		float deltaY = this->leftEdge.end[Y] - this->leftEdge.start[Y];
		float leftRate;
		float leftRec[2];
		leftRate = this->leftEdge.Cross_Param(LineSeg(&e));
		leftRec[X] = userX + deltaX * leftRate;
		leftRec[Y] = userY + deltaY * leftRate;

		deltaX = this->rightEdge.end[X] - this->rightEdge.start[X];
		deltaY = this->rightEdge.end[Y] - this->rightEdge.start[Y];
		float rightRate;
		float rightRec[2];
		rightRate = this->rightEdge.Cross_Param(LineSeg(&e));
		rightRec[X] = userX + deltaX * rightRate;
		rightRec[Y] = userY + deltaY * rightRate;


		if (e.endpoints[Edge::START]->posn[Y] == e.endpoints[Edge::END]->posn[Y]) {//horizonal wall
			if (userY <= e.endpoints[Edge::START]->posn[Y]) {//top horizonal wall
				if ((dir < 90) || (dir > 270)) {//視角往右
					float mX = e.endpoints[Edge::START]->posn[X] > e.endpoints[Edge::END]->posn[X] ? e.endpoints[Edge::START]->posn[X] : e.endpoints[Edge::END]->posn[X];
					float rDir = adjustDir(dir - fov / 2);
					if (rDir < 180) {//決定虛交點方位
						if (leftRec[X] < rightRec[X] && rightRec[X] > mX) {
							rightRec[X] = mX;
						}
					}
					else {
						if (leftRec[X] > rightRec[X]) {
							rightRec[X] = mX;
						}
					}
				}
				else {//視角往左
					float mM = e.endpoints[Edge::START]->posn[X] < e.endpoints[Edge::END]->posn[X] ? e.endpoints[Edge::START]->posn[X] : e.endpoints[Edge::END]->posn[X];
					float lDir = adjustDir(dir + fov / 2);
					if (lDir < 180) {//決定虛交點方位
						if (leftRec[X] < rightRec[X] && leftRec[X] < mM) {
							leftRec[X] = mM;
						}
					}
					else {
						if (leftRec[X] > rightRec[X]) {
							leftRec[X] = mM;
						}
					}
				}
			}
			else {	//bottom wall
				if ((dir < 90) || (dir > 270)) {//視角往右
					float mX = e.endpoints[Edge::START]->posn[X] > e.endpoints[Edge::END]->posn[X] ? e.endpoints[Edge::START]->posn[X] : e.endpoints[Edge::END]->posn[X];
					float lDir = adjustDir(dir + fov / 2);
					if (lDir > 180) {//決定虛交點方位
						if (leftRec[X] > rightRec[X] && leftRec[X] > mX) {
							leftRec[X] = mX;
						}
					}
					else {
						if (leftRec[X] < rightRec[X]) {
							leftRec[X] = mX;
						}
					}
				}
				else {//視角往左
					float mM = e.endpoints[Edge::START]->posn[X] < e.endpoints[Edge::END]->posn[X] ? e.endpoints[Edge::START]->posn[X] : e.endpoints[Edge::END]->posn[X];
					float rDir = adjustDir(dir - fov / 2);
					if (rDir > 180) {//決定虛交點方位
						if (leftRec[X] > rightRec[X] && rightRec[X] < mM) {
							rightRec[X] = mM;
						}
					}
					else {
						if (leftRec[X] < rightRec[X]) {
							rightRec[X] = mM;
						}
					}
				}
			}
		}


		//(dir < 90) || (dir > 270)
		if (e.endpoints[Edge::START]->posn[X] == e.endpoints[Edge::END]->posn[X]) {//vertical wall
			if (userX <= e.endpoints[Edge::START]->posn[X]) {//right vertical wall
				if (dir < 180) {//視角往上
					float mX = e.endpoints[Edge::START]->posn[Y] > e.endpoints[Edge::END]->posn[Y] ? e.endpoints[Edge::START]->posn[Y] : e.endpoints[Edge::END]->posn[Y];
					float lDir = adjustDir(dir + fov / 2);
					if ((lDir < 90) || (lDir > 270)) {//決定虛交點方位
						if (leftRec[Y] > rightRec[Y] && leftRec[Y] > mX) {
							leftRec[Y] = mX;
						}
					}
					else {
						if (leftRec[Y] < rightRec[Y]) {
							leftRec[Y] = mX;
						}
					}
				}
				else {//視角往下
					float mM = e.endpoints[Edge::START]->posn[Y] < e.endpoints[Edge::END]->posn[Y] ? e.endpoints[Edge::START]->posn[Y] : e.endpoints[Edge::END]->posn[Y];
					float rDir = adjustDir(dir - fov / 2);
					if ((rDir < 90) || (rDir > 270)) {//決定虛交點方位
						if (leftRec[Y] > rightRec[Y] && rightRec[Y] < mM) {
							rightRec[Y] = mM;
						}
					}
					else {
						if (leftRec[Y] < rightRec[Y]) {
							rightRec[Y] = mM;
						}
					}
				}
			}
			else {	//left vertical wall
				if (dir < 180) {//視角往上
					float mX = e.endpoints[Edge::START]->posn[Y] > e.endpoints[Edge::END]->posn[Y] ? e.endpoints[Edge::START]->posn[Y] : e.endpoints[Edge::END]->posn[Y];
					float rDir = adjustDir(dir - fov / 2);
					if (!((rDir < 90) || (rDir > 270))) {//決定虛交點方位
						if (leftRec[Y] < rightRec[Y] && rightRec[Y] > mX) {
							rightRec[Y] = mX;
						}
					}
					else {
						if (leftRec[Y] > rightRec[Y]) {
							rightRec[Y] = mX;
						}
					}
				}
				else {//視角往下
					float mM = e.endpoints[Edge::START]->posn[Y] < e.endpoints[Edge::END]->posn[Y] ? e.endpoints[Edge::START]->posn[Y] : e.endpoints[Edge::END]->posn[Y];
					float lDir = adjustDir(dir + fov / 2);
					if (!((lDir < 90) || (lDir > 270))) {//決定虛交點方位
						if (leftRec[Y] < rightRec[Y] && leftRec[Y] < mM) {
							leftRec[Y] = mM;
						}
					}
					else {
						if (leftRec[Y] > rightRec[Y]) {
							leftRec[Y] = mM;
						}
					}
				}
			}
		}



		if (startInTriangle && endInTriangle) {
			opt.start[X] = e.endpoints[Edge::START]->posn[X];
			opt.start[Y] = e.endpoints[Edge::START]->posn[Y];
			opt.end[X] = e.endpoints[Edge::END]->posn[X];
			opt.end[Y] = e.endpoints[Edge::END]->posn[Y];
			return opt;
		}
		if (startInTriangle ^ endInTriangle) {
			if (!startInTriangle) {
				float distanceToLeft = (e.endpoints[Edge::START]->posn[X] - leftRec[X]) * (e.endpoints[Edge::START]->posn[X] - leftRec[X]) + (e.endpoints[Edge::START]->posn[Y] - leftRec[Y]) * (e.endpoints[Edge::START]->posn[Y] - leftRec[Y]);
				float distanceToRight = (e.endpoints[Edge::START]->posn[X] - rightRec[X]) * (e.endpoints[Edge::START]->posn[X] - rightRec[X]) + (e.endpoints[Edge::START]->posn[Y] - rightRec[Y]) * (e.endpoints[Edge::START]->posn[Y] - rightRec[Y]);
				if (distanceToLeft <= distanceToRight) {
					opt.start[X] = leftRec[X];
					opt.start[Y] = leftRec[Y];
					opt.end[X] = e.endpoints[Edge::END]->posn[X];
					opt.end[Y] = e.endpoints[Edge::END]->posn[Y];
				}
				else {
					opt.start[X] = rightRec[X];
					opt.start[Y] = rightRec[Y];
					opt.end[X] = e.endpoints[Edge::END]->posn[X];
					opt.end[Y] = e.endpoints[Edge::END]->posn[Y];
				}
			}
			else if (!endInTriangle) {
				float distanceToLeft = (e.endpoints[Edge::END]->posn[X] - leftRec[X]) * (e.endpoints[Edge::END]->posn[X] - leftRec[X]) + (e.endpoints[Edge::END]->posn[Y] - leftRec[Y]) * (e.endpoints[Edge::END]->posn[Y] - leftRec[Y]);
				float distanceToRight = (e.endpoints[Edge::END]->posn[X] - rightRec[X]) * (e.endpoints[Edge::END]->posn[X] - rightRec[X]) + (e.endpoints[Edge::END]->posn[Y] - rightRec[Y]) * (e.endpoints[Edge::END]->posn[Y] - rightRec[Y]);
				if (distanceToLeft <= distanceToRight) {
					opt.start[X] = e.endpoints[Edge::START]->posn[X];
					opt.start[Y] = e.endpoints[Edge::START]->posn[Y];
					opt.end[X] = leftRec[X];
					opt.end[Y] = leftRec[Y];
				}
				else {
					opt.start[X] = e.endpoints[Edge::START]->posn[X];
					opt.start[Y] = e.endpoints[Edge::START]->posn[Y];
					opt.end[X] = rightRec[X];
					opt.end[Y] = rightRec[Y];
				}
			}
			return opt;
		}
		if ((!startInTriangle) && (!endInTriangle)) {




			bool leftLegal = false;		//與視錐線比較方位 避免取到背後的
			bool rightLegal = false;

			if ((leftIntersection[X] >= leftRec[X] && leftRec[X] >= userX) || (leftIntersection[X] <= leftRec[X] && leftRec[X] <= userX)) {
				if ((leftIntersection[Y] >= leftRec[Y] && leftRec[Y] >= userY) || (leftIntersection[Y] <= leftRec[Y] && leftRec[Y] <= userY)) {
					leftLegal = true;
				}
			}
			if ((rightIntersection[X] >= rightRec[X] && rightRec[X] >= userX) || (rightIntersection[X] <= rightRec[X] && rightRec[X] <= userX)) {
				if ((rightIntersection[Y] >= rightRec[Y] && rightRec[Y] >= userY) || (rightIntersection[Y] <= rightRec[Y] && rightRec[Y] <= userY)) {
					rightLegal = true;
				}
			}

			if (leftLegal && rightLegal) {
				//方位排除法
				char relative_pos_start[2], relative_pos_end[2];
				relative_pos_start[0] = this->leftEdge.Point_Side(e.endpoints[Edge::START]->posn[X], e.endpoints[Edge::START]->posn[Y]);
				relative_pos_start[1] = this->rightEdge.Point_Side(e.endpoints[Edge::START]->posn[X], e.endpoints[Edge::START]->posn[Y]);
				relative_pos_end[0] = this->leftEdge.Point_Side(e.endpoints[Edge::END]->posn[X], e.endpoints[Edge::END]->posn[Y]);
				relative_pos_end[1] = this->rightEdge.Point_Side(e.endpoints[Edge::END]->posn[X], e.endpoints[Edge::END]->posn[Y]);

				bool leftTrackStart = false, rightTrackStart = false;// 紀錄相對兩條視錐線的方向
				if (relative_pos_start[0] == LineSeg::RIGHT) {
					rightTrackStart = true;
				}
				else if (relative_pos_start[0] == LineSeg::LEFT) {
					leftTrackStart = true;
				}
				if (relative_pos_start[1] == LineSeg::RIGHT) {
					rightTrackStart = true;
				}
				else if (relative_pos_start[1] == LineSeg::LEFT) {
					leftTrackStart = true;
				}


				bool leftTrackEnd = false, rightTrackEnd = false;// 紀錄相對兩條視錐線的方向
				if (relative_pos_end[0] == LineSeg::RIGHT) {
					rightTrackEnd = true;
				}
				else if (relative_pos_end[0] == LineSeg::LEFT) {
					leftTrackEnd = true;
				}
				if (relative_pos_end[1] == LineSeg::RIGHT) {
					rightTrackEnd = true;
				}
				else if (relative_pos_end[1] == LineSeg::LEFT) {
					leftTrackEnd = true;
				}
				bool startOut = (leftTrackStart ^ rightTrackStart);
				bool endOut = (leftTrackEnd ^ rightTrackEnd);
				if (startOut && endOut) {
					if (((leftTrackStart && leftTrackEnd) || (rightTrackStart && rightTrackEnd))) {	//out判定失誤 只有判斷不同邊功能 所以用確定對兩者都是同方向超界來處理
						LineSeg l;
						l.notDraw = true;
						return l;
					}
				}
				//

				opt.start[X] = leftRec[X];
				opt.start[Y] = leftRec[Y];
				opt.end[X] = rightRec[X];
				opt.end[Y] = rightRec[Y];
				return opt;
			}
			else {
				LineSeg l;
				l.notDraw = true;
				return l;
			}
		}
	}
};



//************************************************************************
//
// * The maze consists of cells, separated by edges. NOTE: The maze is defined
//   assuming that z is up, with xy forming the ground plane. This is different
//   to the OpenGL viewing assumption (which has y up and xz in the ground
//   plane). You will have to take this into account when drawing the view.
//   Also, assume that the floor of the maze is at z = -1, and the ceiling is
//   at z = 1.
//
//************************************************************************
class Maze {

public:
	// The first constructor takes the number of cells in the x and y 
	// directions, and the cell size in each dimension. This constructor
	// creates a random maze, and returns it.
	Maze(const int num_x, const int num_y,
		const float size_x, const float size_y);

	// The second constructor takes a maze file name to load. It may throw
	// exceptions of the MazeException class if there is an error.
	Maze(const char* f);

	~Maze(void);

	float* lookatmatrix;
	float* perspectivematrix;

public:
	void matrixto2d(float vec[4]) {
		float tmp[4];
		for (int i = 0; i < 4; ++i) {
			float sum = 0;
			for (int j = 0; j < 4; ++j) {
				sum += lookatmatrix[i + (j * 4)] * vec[j];
			}
			tmp[i] = sum;
		}
		for (int i = 0; i < 4; ++i)
			vec[i] = tmp[i];

		for (int i = 0; i < 4; ++i) {
			float sum = 0;
			for (int j = 0; j < 4; ++j) {
				sum += perspectivematrix[i + (j * 4)] * vec[j];
			}
			tmp[i] = sum;
		}
		for (int i = 0; i < 4; ++i)
			vec[i] = tmp[i];
	}

	// Set the viewer's location 
	void	Set_View_Posn(float x, float y, float z);

	// Set the angle in which the viewer is looking.
	void	Set_View_Dir(const float);

	// Set the horizontal field of view.
	void	Set_View_FOV(const float);

	// Move the viewer's position. This method will do collision detection
	// between the viewer's location and the walls of the maze and prevent
	// the viewer from passing through walls.
	void	Move_View_Posn(const float dx, const float dy, const float dz);

	// Draws the map view of the maze. It is passed the minimum and maximum
	// corners of the window in which to draw.
	void	Draw_Map(int, int, int, int);

	// Draws the viewer's cell and its neighbors in the map view of the maze.
	// It is passed the minimum and maximum corners of the window in which
	// to draw.
	void	Draw_Neighbors(int, int, int, int);

	// Draws the frustum on the map view of the maze. It is passed the
	// minimum and maximum corners of the window in which to draw.
	void	Draw_Frustum(int, int, int, int);

	//Draw walls
	void	Draw_Wall(const float[2], const float[2], const float[3]);

	void    initFrustum(Frustum&);
	void	Draw_Cell(Cell* c, Frustum f);


	// Draws the first-person view of the maze. It is passed the focal distance.
	// THIS IS THE FUINCTION YOU SHOULD MODIFY.
	void	Draw_View(const float, const float);

	// Save the maze to a file of the given name.
	bool	Save(const char*);

	// Functions to convert between degrees and radians.
	static double   To_Radians(double deg) { return deg / 180.0 * M_PI; };
	static double   To_Degrees(double rad) { return rad * 180.0 / M_PI; };
private:
	// Functions used when creating or loading a maze.

	// Randomly generate the edge's opaque and transparency for an empty maze
	void    Build_Connectivity(const int, const int, const float, const float);
	// Grow a maze by removing candidate edges until all the cells are
	// connected. The edges are not actually removed, they are just made
	// transparent.
	void    Build_Maze(void);
	void    Set_Extents(void);
	void    Find_View_Cell(Cell*);

private:
	Cell* view_cell;// The cell that currently contains the view
									  // point. You will need to use this.
	unsigned int    frame_num;	// The frame number we are currently drawing.
										// It isn't necessary, but you might find it
										// helpful for debugging or something.

	static const float	BUFFER;	// The viewer must be at least this far inside
											// an exterior wall of the maze.
											// Not implemented

	float	min_xp;	// The minimum x location of any vertex in the maze.
	float	min_yp;	// The minimum y location of any vertex in the maze.
	float	max_xp;	// The maximum x location of any vertex in the maze.
	float	max_yp;	// The maximum y location of any vertex in the maze.

public:
	static const char	X; // Used to index into the viewer's position
	static const char	Y;
	static const char	Z;

	int		num_vertices;	// The number of vertices in the maze
	Vertex** vertices;		// An array of pointers to the vertices.

	int		num_edges;		// The number of edges in the maze.
	Edge** edges;			// An array of pointers to the edges.

	int		num_cells;     // The number of cells in the maze
	Cell** cells;       // An array of pointers to the cells.

	float		viewer_posn[3];	// The x,y location of the viewer.
	float		viewer_dir;			// The direction in which the viewer is
										// looking. Measured in degrees about the z
										// axis, in the usual way.
	float		viewer_fov;			// The horizontal field of view, in degrees.



};


#endif

