#include "spline.h"

const int Division = 20;

Spline::Spline(GraphWidget *graphWidget)
:graph(graphWidget)
{
	deleteNodes(interpPoints);
	deleteNodes(controlPoints);
	deleteNodes(curvePoints);
    
	deleteEdges(pseudoEdges);
	deleteEdges(curveEdges);
	deleteEdges(controlEdges);
    
	startPoint = NULL;
	endPoint = NULL;
	totalPoints = 0;
	style = BezierBernstein;
	showControlPoints = false;
	showCurvePoints = false; // LOOK
	showPseudoPoints = true;
}

void Spline::SetVisible(QList<Node *> nodes, bool visible)
{
	foreach(Node * i, nodes)
	{
		i->setVisible(visible);
		i->setEnabled(visible);
	}
}
void Spline::SetVisible(QList<Edge *> edges, bool visible)
{
	foreach(Edge * i, edges)
	{
		i->setVisible(visible);
		i->setEnabled(visible);
	}
}
void Spline::SetVisible(bool visible)
{
	if (this->startPoint && this->endPoint)
	{
		this->startPoint->setVisible(visible);
		this->startPoint->setEnabled(visible);
		this->endPoint->setVisible(visible);
		this->endPoint->setEnabled(visible);
		if (!visible)
		{
			this->graph->scene()->removeItem(startPoint);
			this->graph->scene()->removeItem(endPoint);
		}
		else
		{
			this->graph->scene()->addItem(startPoint);
			this->graph->scene()->addItem(endPoint);
		}
	}
}

void Spline::Clear()
{
	deleteNodes(interpPoints);
	deleteNodes(controlPoints);
	deleteNodes(curvePoints);
    
	deleteEdges(pseudoEdges);
	deleteEdges(curveEdges);
	deleteEdges(controlEdges);
    
	if (startPoint)
	{
		this->graph->scene()->removeItem(startPoint);
		delete startPoint;
		startPoint = NULL;
	}
	if (endPoint)
	{
		this->graph->scene()->removeItem(endPoint);
		delete endPoint;
		endPoint = NULL;
	}
    
	totalPoints = 0;
}

void Spline::PaintPseudeLines()
{
	if (startPoint && endPoint)
	{
		Interpolate();
		deleteEdges(pseudoEdges);
		calcPseudoEdges();
		Paint();
	}
}

void Spline::deleteNodes(QList<Node *> &nodes)
{
	foreach(Node * i, nodes)
	{
		delete i;
	}
	nodes.clear();
}

void Spline::deleteEdges(QList<Edge *> &edges)
{
	foreach(Edge * i, edges)
	{
		delete i;
	}
	edges.clear();
}

void Spline::Paint()
{
	deleteEdges(curveEdges);
	deleteEdges(pseudoEdges);
	deleteEdges(controlEdges);
    
	// Curve Points & Edges
    // LOOK
	if(showCurvePoints)
	{
		foreach(Node * i, curvePoints)
		{
			this->graph->scene()->addItem(i);
		}
	}
	calcCurveEdges();
    
	// Control Points & Edges
	if (showControlPoints)
	{
		foreach(Node * i, controlPoints)
		{
			this->graph->scene()->addItem(i);
		}
		calcControlEdges();
	}
    
	if (showPseudoPoints)
	{
		calcPseudoEdges();
		this->SetVisible(true);
	}
	else
	{
		this->SetVisible(false);
	}
    
	this->graph->repaint(-1000, -1000, 2000, 2000);
}

void Spline::AddPoint(Node *newNode)
{
	interpPoints.append(newNode);
	graph->scene()->addItem(newNode);
	totalPoints++;
    
	CalcStartPoint();
	CalcEndPoint();
}

// LOOK (Adding intermediate control points)
void Spline::AddPoint(Node *newNode, int& index)
{
	// Add point as above but using the index value
}

// LOOK (Adding intermediate control points)
void Spline::DeletePoint(int& index)
{
	// Remove Node
    
	// Remove the start and end points if points less than 2
	if (totalPoints<=2)
	{
		if (startPoint)
		{
			ClearStartPoint();
		}
		if (endPoint)
		{
			ClearEndPoint();
		}
	}
	else
	{
		// Calculate the first and last points if those get deleted
		if (index == 0)
		{
			ClearStartPoint();
			CalcStartPoint();
		}
		if(index==totalPoints)
		{
			ClearEndPoint();
			CalcEndPoint();
		}
	}
    
	// Recalculate start and end points if those are the only points
	if(totalPoints==2)
	{
		ClearStartPoint();
		ClearEndPoint();
		CalcStartPoint();
		CalcEndPoint();
	}
}

void Spline::ClearStartPoint()
{
	this->graph->scene()->removeItem(startPoint);
	delete startPoint;
	startPoint = NULL;
}

void Spline::ClearEndPoint()
{
	this->graph->scene()->removeItem(endPoint);
	delete endPoint;
	endPoint = NULL;
}

void Spline::CalcStartPoint()
{
	if (totalPoints == 2 || (totalPoints>2 && !startPoint))
	{
		Node* p1 = interpPoints[0];
		Node* p2 = interpPoints[1];
        
		QPointF line12 = p1->pos() - p2->pos();
		double length12 = sqrt(line12.x()*line12.x() + line12.y()*line12.y());
        
		startPoint = new Node(this->graph, PseudoPoints);
		startPoint->setPos(p1->pos() + line12 / length12 * 50);
		graph->scene()->addItem(startPoint);
	}
}

void Spline::CalcEndPoint()
{
	if (totalPoints < 2)
		return;
	if (!endPoint)
	{
		endPoint = new Node(this->graph, PseudoPoints);
		graph->scene()->addItem(endPoint);
	}
    
	Node* pN = interpPoints[totalPoints - 1];
	Node* pNm1 = interpPoints[totalPoints - 2];
    
	QPointF lineNNm1 = pN->pos() - pNm1->pos();
	double lengthNNm1 = sqrt(lineNNm1.x()*lineNNm1.x() + lineNNm1.y()*lineNNm1.y());
    
	endPoint->setPos(pN->pos() + lineNNm1 / lengthNNm1 * 50);
}

void Spline::calcCurveEdges()
{
	deleteEdges(curveEdges);
	for (int i = 0; i < curvePoints.size() - 1; i++)
	{
		Edge * tmp = new Edge(curvePoints[i]->pos(), curvePoints[i+1]->pos());
		curveEdges.append(tmp);
		graph->scene()->addItem(tmp);
	}
}

void Spline::calcControlEdges()
{
	deleteEdges(controlEdges);
    
	switch(style)
	{
        case BezierBernstein:
        case BezierDeCasteljau:
        case BezierMatrix:
            if (this->controlPoints.size() == 2*(totalPoints-1))
            {
                for (int i = 0; i < this->totalPoints - 1; i ++)
                {
                    Node * n1 = interpPoints[i];
                    Node * n2 = interpPoints[i+1];
                    Node * c1 = controlPoints[2*i];
                    Node * c2 = controlPoints[2*i+1];
                    Edge * e1 = new Edge(n1->pos(), c1->pos(), 1);
                    Edge * e2 = new Edge(c1->pos(), c2->pos(), 1);
                    Edge * e3 = new Edge(c2->pos(), n2->pos(), 1);
                    controlEdges.append(e1);
                    controlEdges.append(e2);
                    controlEdges.append(e3);
                    graph->scene()->addItem(e1);
                    graph->scene()->addItem(e2);
                    graph->scene()->addItem(e3);
                }
            }
            break;
        case BSpline:
            for (int i = 0; i < this->controlPoints.size() - 1; i++)
            {
                Node * c1 = controlPoints[i];
                Node * c2 = controlPoints[i+1];
                
                Edge * e1 = new Edge(c1->pos(), c2->pos(), 1);
                controlEdges.append(e1);
                graph->scene()->addItem(e1);
            }
            break;
        case HermiteClamped:
        case HermiteNatural:
            if (this->controlPoints.size() == totalPoints)
            {
                for (int i = 0; i < this->totalPoints; i++)
                {
                    Node * c1 = interpPoints[i];
                    Node * c2 = controlPoints[i];
                    
                    c2->setPosition(c1->getPosition() + c2->getPosition());
                    
                    Edge * e1 = new Edge(c1->pos(), c2->pos(), 1);
                    controlEdges.append(e1);
                    graph->scene()->addItem(e1);
                }
            }
            break;
	}
}

void Spline::calcPseudoEdges()
{
	if (totalPoints >= 2)
	{
		deleteEdges(pseudoEdges);
		// start
		Edge * start = new Edge(startPoint->pos(), ((Node *)interpPoints.first())->pos(), 2);
		pseudoEdges.append(start);
		graph->scene()->addItem(start);
		// end
		Edge * end = new Edge(endPoint->pos(), ((Node *)interpPoints.last())->pos(), 2);
		pseudoEdges.append(end);
		graph->scene()->addItem(end);
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Start your code here.
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Spline::Interpolate()
{
	//Clear the old curve points
	deleteNodes(curvePoints);
	//Calculate control points.
	CalcControlPoints();
	//Depending on the selected style, interpolate the curve.
	switch(style){
        case BezierBernstein:		InterpBernstein();       break;
        case BezierDeCasteljau:		InterpCasteljau();      break;
        case BezierMatrix:          InterpMatrix();         break;
        case BSpline:               InterpBSpline();        break;
        case HermiteClamped:        InterpHermiteClamped(); break;
        case HermiteNatural:        InterpHermiteNatural(); break;
	}
}


//////////////////////////////////////////////////////////////////////////
// Calculate the control points
//////////////////////////////////////////////////////////////////////////
// This function utilizes the following member variables
// interpPoints	- type: QList<Node *>
//				  discription: stores all the interpolation points
// startPoint	- type: Node*
//				  discription: stores the tangent to the first interpolation point
// endPoint		- type: Node*
//				  discription: stores the tangent to the second interpolation point
// totalPoints	- type: int
//				  discription: total number of interpolation points
// This function modifies the following member variables
// controlPoints- type: QList<Node *>
//				  discription: stores all the control points for curve interpolation and display.
//                             For Bezier curve, between very pair of consecutive interpolation points,
//                             there should be two control points. These four points determins the curve interpolation.
//                             For B-Spline, there should be totalPoints + 2 control points calculated from Ac = p.
//
// Hint: To implement Hermite and B-Spline, you need to write functions to create the A matrix as in the handouts.
//       Then you solve a linear system Ac = p, where p is the interpolation points vector and c are the control points.
//       We have provided you with a data structure to store and solve the linear system.
//       Below is an example code, read the understand it.
//
//	matrix<float> A(3,3);
//  matrix<float> c(3,1);
//  matrix<float> p(3,1);
//  A(0,0) = 1.0; A(0,1) = 0.0; A(0,2) = 0.0;
//  A(1,0) = 0.0; A(1,1) = 1.0; A(1,2) = 0.0;
//  A(2,0) = 0.0; A(2,1) = 0.0; A(2,2) = 1.0;
//  p(0,0) = 1.0; p(1,0) = 2.0; p(2,0) = 3.0;
//  c = A.Solve(p);
//
//  The result in c is c(0,0) = 1.0; c(1,0) = 2.0; c(3,0) = 3.0, which satisfies Ac = p.
// Hint2: To get the position which is a vec3 in the class Node*, we need to use the method Node::getPosition().
// Example:
// vec3 startPointPos = startPoint->getPosition();
//
// Hint3: To new a Node*, add it into the List, and set the position of it, you need to do like this:
// Example:
// Node* ctrl = new Node(graph, ControlPoints);
// ctrl->setPosition(/*some vec3*/);
// controlPoints.append(ctrl);
void Spline::CalcControlPoints()
{
	// Clear previous controlPoints
	deleteNodes(controlPoints);
	// Based on style use appropriate ways to compute control points
	switch(style)
	{
        case BezierBernstein:
        case BezierDeCasteljau:
        case BezierMatrix:
            // ADD YOUR CODE HERE
			// interp point 0 - start point
			if (totalPoints == 0 || totalPoints == 1) {
				break;
			}
			for (int i = 0; i  < totalPoints; i++) {
				if (i == 0) {
					Node* ctrl = new Node(graph,ControlPoints);
					ctrl->setPosition((interpPoints[0]->getPosition() - startPoint->getPosition()) + interpPoints[0]->getPosition());
					controlPoints.append(ctrl);
				}
				else if (i == totalPoints - 1) {
					Node* ctrl = new Node(graph,ControlPoints);
					ctrl->setPosition((interpPoints[totalPoints - 1]->getPosition() - endPoint->getPosition()) + interpPoints[totalPoints - 1]->getPosition());
					controlPoints.append(ctrl);
				}
				else {
					Node* ctrl2 = new Node(graph,ControlPoints);
					ctrl2->setPosition(interpPoints[i]->getPosition() - ((interpPoints[i+1]->getPosition() - interpPoints[i-1]->getPosition()) / 4));
					controlPoints.append(ctrl2);

					Node* ctrl = new Node(graph,ControlPoints);
					ctrl->setPosition(interpPoints[i]->getPosition() + ((interpPoints[i+1]->getPosition() - interpPoints[i-1]->getPosition()) / 4));
					controlPoints.append(ctrl);
				}
			}
            break;
        case BSpline:
            // make sure there are at least 2 totalPoints.
            // ADD YOUR CODE HERE
            break;
            // In the case for Hermite, you want to implement both "clamped" and "natural" versions.
        case HermiteClamped:
			if (interpPoints.length() >= 2) {
				matrix<float> m1 = createHermiteClampedMatrix1();
				matrix<float> m2 = createHermiteClampedMatrix2();
				matrix<float> m3 = m1.Solve(m2);
				for (int i = 0; i < interpPoints.length(); i++) {
					Node* ctrl = new Node(graph,ControlPoints);
					vec3 pos(vec3(m3(i,0),m3(i,1),m3(i,2)));
					ctrl->setPosition(pos);
					controlPoints.append(ctrl);
				}
			}
            // make sure there are at least 2 totalPoints.
            // ADD YOUR POOP HERE
            break;
        case HermiteNatural:
			if (interpPoints.length() >= 2) {
				matrix<float> m1 = createHermiteNaturalMatrix1();
				matrix<float> m2 = createHermiteNaturalMatrix2();
				matrix<float> m3 = m1.Solve(m2);
				for (int i = 0; i < interpPoints.length(); i++) {
					Node* ctrl = new Node(graph,ControlPoints);
					vec3 pos(vec3(m3(i,0),m3(i,1),m3(i,2)));
					ctrl->setPosition(pos);
					controlPoints.append(ctrl);
				}
			}
            break;
	}
}

matrix<float> Spline::createHermiteClampedMatrix2() {
	int n = interpPoints.length();
	vec3 p;
	matrix<float> output(n,3);
	for(int i = 0; i < n; i++) {
		if (i == 0) {
			p = (interpPoints[0]->getPosition() - startPoint->getPosition());
			output(i,0) = -p[0];
			output(i,1) = -p[1];
			output(i,2) = -p[2];
		}
		else if (i == n-1) {
			p = (endPoint->getPosition() - interpPoints[n-1]->getPosition());
			output(i,0) = p[0];
			output(i,1) = p[1];
			output(i,2) = p[2];
		}
		else {
			p = (3*(interpPoints[i+1]->getPosition() - interpPoints[i-1]->getPosition()));
			output(i,0) = p[0];
			output(i,1) = p[1];
			output(i,2) = p[2];
		}
	}
	return output;
}

matrix<float> Spline::createHermiteClampedMatrix1() {
	int n = interpPoints.length();
	matrix<float> output(n,n);
	for (int i = 0; i < n; i++) {
		for(int j = 0; j < n; j++) {
			if (i == 0 && j == 0) {
				output(i,j) = 1;
			}
			else if (i == (n-1) && j == (n-1)) {
				output(i,j) = 1;
			}
			else if (i == j) {
				output(i,j+1) = 1;
				output(i,j) = 4;
				output(i,j-1) = 1;
			}
			else {
				output(i,j) = 0;
			}
		}
	}
	return output;
}

matrix<float> Spline::createHermiteNaturalMatrix1() {
	int n = interpPoints.length();
	matrix<float> output(n,n);
	for (int i = 0; i < n; i++) {
		for(int j = 0; j < n; j++) {
			if (i == 0 && j == 0) {
				output(i,j) = 2;
				output(i,j+1) = 1;
			}
			else if (i == (n-1) && j == (n-1)) {
				output(i,j) = 2;
				output(i,j-1) = 1;
			}
			else if (i == j) {
				output(i,j+1) = 1;
				output(i,j) = 4;
				output(i,j-1) = 1;
			}
			else {
				output(i,j) = 0;
			}
		}
	}
	return output;
}

matrix<float> Spline::createHermiteNaturalMatrix2() {
	int n = interpPoints.length();
	vec3 p;
	matrix<float> output(n,3);
	for(int i = 0; i < n; i++) {
		if (i == 0) {
			p = (3*(interpPoints[1]->getPosition() - interpPoints[0]->getPosition()));
			output(i,0) = p[0];
			output(i,1) = p[1];
			output(i,2) = p[2];
		}
		else if (i == n-1) {
			p = (3*(interpPoints[n-1]->getPosition() - interpPoints[n-2]->getPosition()));
			output(i,0) = p[0];
			output(i,1) = p[1];
			output(i,2) = p[2];
		}
		else {
			p = (3*(interpPoints[i+1]->getPosition() - interpPoints[i-1]->getPosition()));
			output(i,0) = p[0];
			output(i,1) = p[1];
			output(i,2) = p[2];
		}
	}
	return output;
}

//////////////////////////////////////////////////////////////////////////
// Cubic Berstein Bezier Spline
//////////////////////////////////////////////////////////////////////////
// This function utilizes the following member variables
// interpPoints	- type: QList<Node *>
//				  discription: stores all the interpolation points
// controlPoints- type: QList<Node *>
//				  discription: stores the control points that helps to determine the curve.
//                             Between very pair of consecutive interpolation points,there should be two control points.
//                             These four points determins the curve interpolation.
// startPoint	- type: Node*
//				  discription: stores the tangent to the first interpolation point
// endPoint		- type: Node*
//				  discription: stores the tangent to the second interpolation point
// totalPoints	- type: int
//				  discription: total number of interpolation points
// This function modifies the following member variables
// curvePoints	- type: QList<Node *>
//				  discription: stores all the points that form the curve, including all interpolation points
//
// Hint: To get the position which is a vec3 in the class Node*, we need to use the method Node::getPosition().
// Example:
// vec3 b = interpPoints[i]->getPosition();
//
// Hint2: To new a Node*, add it into the List, and set the position of it, you need to do like this:
// Example:
// Node* newNode = new Node(graph, CurvePoints);
// newNode->setPosition(/*some vec3*/);
// curvePoints.append(newNode);
void Spline::InterpBernstein()
{
	for(int i = 0; i < interpPoints.length() - 1; i ++) {
		for( int t = 0; t < 20; t++) {
			vec3 point = interpPoints[i]->getPosition() * bernstein(3,0,(float)t/20.f);
			vec3 point2 = controlPoints[i*2]->getPosition() * bernstein(3,1,(float)t/20.f);
			vec3 point3 = controlPoints[i*2 + 1]->getPosition() * bernstein(3,2,(float)t/20.f);
			vec3 point4 = interpPoints[i+1]->getPosition() * bernstein(3,3,(float)t/20.f);

			Node* newNode = new Node(graph, CurvePoints);
			newNode->setPosition(point + point2 + point3 + point4);
			curvePoints.append(newNode);
		}
	}
}

//////////////////////////////////////////////////////////////////////////
// Cubic de Casteljau Bezier Spline
//////////////////////////////////////////////////////////////////////////
// This function utilizes the following member variables
// interpPoints	- type: QList<Node *>
//				  discription: stores all the interpolation points
// controlPoints- type: QList<Node *>
//				  discription: stores the control points that helps to determine the curve.
//                             Between very pair of consecutive interpolation points,there should be two control points.
//                             These four points determins the curve interpolation.
// startPoint	- type: Node*
//				  discription: stores the tangent to the first interpolation point
// endPoint		- type: Node*
//				  discription: stores the tangent to the second interpolation point
// totalPoints	- type: int
//				  discription: total number of interpolation points
// This function modifies the following member variables
// curvePoints	- type: QList<Node *>
//				  discription: stores all the points that form the curve, including all interpolation points
//
// Hint: To get the position which is a vec3 in the class Node*, we need to use the method Node::getPosition().
// Example:
// vec3 b = interpPoints[i]->getPosition();
//
// Hint2: To new a Node*, add it into the List, and set the position of it, you need to do like this:
// Example:
// Node* newNode = new Node(graph, CurvePoints);
// newNode->setPosition(/*some vec3*/);
// curvePoints.append(newNode);
void Spline::InterpCasteljau()
{
	// ADD YOU CODE HERE
	for(int i = 0; i < interpPoints.length() - 1; i ++) {
		for( int j = 0; j < 20; j++) {
			float t = (float)j / 20.f;
			vec3 a = interpPoints[i]->getPosition() + ((controlPoints[2*i]->getPosition() - interpPoints[i]->getPosition()) * t);
			vec3 b = controlPoints[2*i]->getPosition() + ((controlPoints[2*i+1]->getPosition() - controlPoints[2*i]->getPosition())*t);
			vec3 c = controlPoints[2*i+1]->getPosition() + ((interpPoints[i+1]->getPosition() - controlPoints[2*i+1]->getPosition())*t);

			vec3 d = a + ((b-a)*t);
			vec3 e = b + ((c-b)*t);

			vec3 f = d + ((e-d)*t);

			Node* newNode = new Node(graph, CurvePoints);
			newNode->setPosition(f);
			curvePoints.append(newNode);
		}
	}
}


//////////////////////////////////////////////////////////////////////////
// Interpolation Matrix Bezier Spline
//////////////////////////////////////////////////////////////////////////
// This function utilizes the following member variables
// interpPoints	- type: QList<Node *>
//				  discription: stores all the interpolation points
// controlPoints- type: QList<Node *>
//				  discription: stores the control points that helps to determine the curve.
//                             Between very pair of consecutive interpolation points,there should be two control points.
//                             These four points determins the curve interpolation.
// startPoint	- type: Node*
//				  discription: stores the tangent to the first interpolation point
// endPoint		- type: Node*
//				  discription: stores the tangent to the second interpolation point
// totalPoints	- type: int
//				  discription: total number of interpolation points
// This function modifies the following member variables
// curvePoints	- type: QList<Node *>
//				  discription: stores all the points that form the curve, including all interpolation points
//
// Hint: To get the position which is a vec3 in the class Node*, we need to use the method Node::getPosition().
// Example:
// vec3 b = interpPoints[i]->getPosition();
//
// Hint2: To new a Node*, add it into the List, and set the position of it, you need to do like this:
// Example:
// Node* newNode = new Node(graph, CurvePoints);
// newNode->setPosition(/*some vec3*/);
// curvePoints.append(newNode);
// Hint3: If you have any problem in vec3 vec4 mat4 operations, see the SvzAlgebra.h and svzalgebra.cpp to find some specific method out.
// For Example:
// vec3 operator * (const mat4& a, const vec3& v); has been defined in SvzAlgebra.h
// That means you can do operations like this:
// vec3 A;
// mat4 M;
// vec3 B = M*A;
// (It automatically transfer the A and B vectors to homogeneous vectors)
void Spline::InterpMatrix()
{
	// ADD YOUR CODE HERE.
	vec4 v1 = vec4(-1,3,-3,1);
	vec4 v2 = vec4(3,-6,3,0);
	vec4 v3 = vec4(-3,3,0,0);
	vec4 v4 = vec4(1,0,0,0);
	mat4 bezMatrix(v1,v2,v3,v4);

	for(int i = 0; i < interpPoints.length() - 1; i ++) {
		for( int j = 0; j < 20; j++) {
			float t = (float)j / 20.f;

			vec4 vU(t*t*t,t*t,t,1);

			vec4 p1(interpPoints[i]->getPosition(),1);
			vec4 p2(controlPoints[2*i]->getPosition(),1);
			vec4 p3(controlPoints[2*i+1]->getPosition(),1);
			vec4 p4(interpPoints[i+1]->getPosition(),1);

			mat4 vG(p1,p2,p3,p4);

			vec4 final = vU*bezMatrix*vG;
			vec3 final3(final[0],final[1],final[2]);

			Node* newNode = new Node(graph, CurvePoints);
			newNode->setPosition(final3);
			curvePoints.append(newNode);
		}
	}

}

//////////////////////////////////////////////////////////////////////////
// BSpline curve
//////////////////////////////////////////////////////////////////////////
// This function utilizes the following member variables
// interpPoints	- type: QList<Node *>
//				  discription: stores all the interpolation points
// controlPoints- type: QList<Node *>
//				  discription: stores the control points that helps to determine the curve.
//                             Between very pair of consecutive interpolation points,there should be two control points.
//                             These four points determins the curve interpolation.
// startPoint	- type: Node*
//				  discription: stores the tangent to the first interpolation point
// endPoint		- type: Node*
//				  discription: stores the tangent to the second interpolation point
// totalPoints	- type: int
//				  discription: total number of interpolation points
// This function modifies the following member variables
// curvePoints	- type: QList<Node *>
//				  discription: stores all the points that form the curve, including all interpolation points
void Spline::InterpBSpline()
{
	// ADD YOUR CODE HERE
}

//////////////////////////////////////////////////////////////////////////
// Hermite Curve ('Clamped' condition)
//////////////////////////////////////////////////////////////////////////
// This function utilizes the following member variables
// interpPoints	- type: QList<Node *>
//				  discription: stores all the interpolation points
// controlPoints- type: QList<Node *>
//				  discription: stores the control points that helps to determine the curve.
//                             Between very pair of consecutive interpolation points,there should be two control points.
//                             These four points determins the curve interpolation.
// startPoint	- type: Node*
//				  discription: stores the tangent to the first interpolation point
// endPoint		- type: Node*
//				  discription: stores the tangent to the second interpolation point
// totalPoints	- type: int
//				  discription: total number of interpolation points
// This function modifies the following member variables
// curvePoints	- type: QList<Node *>
//				  discription: stores all the points that form the curve, including all interpolation points
void Spline::InterpHermiteClamped()
{
	vec4 v1 = vec4(2,-2,1,1);
	vec4 v2 = vec4(-3,3,-2,-1);
	vec4 v3 = vec4(0,0,1,0);
	vec4 v4 = vec4(1,0,0,0);
	mat4 bezMatrix(v1,v2,v3,v4);

	for(int i = 0; i < interpPoints.length() - 1; i ++) {
		for( int j = 0; j < 20; j++) {
			float t = (float)j / 20.f;
			vec4 vU(t*t*t,t*t,t,1);

			vec4 p1(interpPoints[i]->getPosition(),1);
			vec4 p2(interpPoints[i+1]->getPosition(),1);
			vec4 p3(controlPoints[i]->getPosition(),1);
			vec4 p4(controlPoints[i+1]->getPosition(),1);

			mat4 vG(p1,p2,p3,p4);

			vec4 final = vU*bezMatrix*vG;
			vec3 final3(final[0],final[1],final[2]);

			Node* newNode = new Node(graph, CurvePoints);
			newNode->setPosition(final3);
			curvePoints.append(newNode);

		}
	}
}

//////////////////////////////////////////////////////////////////////////
// Hermite Curve ('Natural' condition)
//////////////////////////////////////////////////////////////////////////
// This function utilizes the following member variables
// interpPoints	- type: QList<Node *>
//				  discription: stores all the interpolation points
// controlPoints- type: QList<Node *>
//				  discription: stores the control points that helps to determine the curve.
//                             Between very pair of consecutive interpolation points,there should be two control points.
//                             These four points determins the curve interpolation.
// startPoint	- type: Node*
//				  discription: stores the tangent to the first interpolation point
// endPoint		- type: Node*
//				  discription: stores the tangent to the second interpolation point
// totalPoints	- type: int
//				  discription: total number of interpolation points
// This function modifies the following member variables
// curvePoints	- type: QList<Node *>
//				  discription: stores all the points that form the curve, including all interpolation points 
// Hint: It is exactly the same like the method above.
void Spline::InterpHermiteNatural()
{
	InterpHermiteClamped();
}

float Spline::bernstein(float n, float i, float t) {
	float a = choose(n,i);
	float b = pow(t,i);
	float c = pow(1-t,n-i);
	return  a*b*c;
}

float Spline::choose(float n, float i) {
	int a = factorial(n);
	int b = factorial(i) * factorial(n-i);
	return a/b;
}

float Spline::factorial(float input) {
	int out = input;
	if (input == 0) {
		return 1;
	}
	else {
		for (int i = input-1; i > 0; i--) {
			out *= i;
		}
		return out;
	}
}