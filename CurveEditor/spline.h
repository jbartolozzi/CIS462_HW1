#ifndef SPLINE_H
#define SPLINE_H

#include <math.h>
#include <QList>
#include <QMatrix>
#include <QVector3D>
#include <QDebug>

#include "edge.h"
#include "node.h"
#include "graphwidget.h"
#include "enums.h"
#include "SvzAlgebra.h"
#include "SvzMatrix.h"

class Node;
class Edge;
class GraphWidget;

class Spline
{
public:
	// My added functions
	float factorial(float i);
	float choose(float n, float i);
	float bernstein(float n, float i, float t);
	matrix<float> createHermiteClampedMatrix1();
	matrix<float> createHermiteClampedMatrix2();
	matrix<float> createHermiteNaturalMatrix1();
	matrix<float> createHermiteNaturalMatrix2();

    Spline(GraphWidget *graphWidget);

    // Member functions
    //Prepare data for interpolation:1.Calculate the control points
    void CalcControlPoints();
    //Interpolation selection
    void Interpolate();
    //Interpolation using CatmullRom method
    void InterpBernstein();
    //Interpolation using de Casteljau method
    void InterpCasteljau();
    //Interpolation using matrix form
    void InterpMatrix();
    //Interpolation using BSpline
    void InterpBSpline();
    //Interpolation using Hermite with "clamped" end points constraints.
    void InterpHermiteClamped();
    //Interpolation using Hermite with "natural" end points constraints.
    void InterpHermiteNatural();
    // Clear all the data.
    void Clear();
    // Add a new InterpPoint into the spline.
    void AddPoint(Node* newNode);
	// Add a new InterpPoint into the spline after the specifiec index location.
    void AddPoint(Node* newNode, int& index);
	// Delete an InterpPoint from the spline.
    void DeletePoint(int& index);
	// Clear the start and end points
	void ClearStartPoint();
	void ClearEndPoint();
    // Calculate the start and end points which are two pseudo interpolation points in Bezier spline and Hermite spline (clamped).
    void CalcStartPoint();
    void CalcEndPoint();
    // Paint
    void Paint();
    void PaintPseudeLines();
	// Set Point Visiblilty.
	void SetVisible(QList<Node *> nodes, bool visible);
	void SetVisible(QList<Edge *> edges, bool visible);
	void SetVisible(bool visible);

    // Member variables
    //Total number of interpolation points
    int totalPoints;
    //Interpolation style, refer to enum {CATMULLROM, CASTELJAU, MATRIX, BSPLINE, HERMITE}.
    Style style;
    //Two end points that helps to shape the curve.
    Node* startPoint;
    Node* endPoint;
    //Vector of interpolation points.
    QList<Node *> interpPoints;
    //Vector of control points generated from all interpolation points.
    QList<Node *> controlPoints;
    //Vector of points that displays the curve.
    QList<Node *> curvePoints;
    //Vector of Lamda.
    vector<float> L;

    QList<Edge *> curveEdges;
    QList<Edge *> pseudoEdges;
    QList<Edge *> controlEdges;

    bool showControlPoints;
	bool showCurvePoints;
    bool showPseudoPoints;

private:
    GraphWidget *graph;

    void deleteNodes(QList<Node *> &nodes);
    void deleteEdges(QList<Edge *> &edges);

    // Calculate Edges for painting
    void calcCurveEdges();
    void calcControlEdges();
    void calcPseudoEdges();

	void HermiteHelper(matrix<float>&, matrix<float>&, matrix<float>&);
};

#endif // SPLINE_H
