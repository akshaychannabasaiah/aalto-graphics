#include "curve.h"
#include "extra.h"
#ifdef WIN32
#include <windows.h>
#endif
#include <GL/gl.h>
using namespace std;
using namespace FW;

namespace {

	// Approximately equal to.  We don't want to use == because of
	// precision issues with floating point.
	inline bool approx(const Vec3f& lhs, const Vec3f& rhs) {
		const float eps = 1e-8f;
		return (lhs - rhs).lenSqr() < eps;
	}


	// This is the core routine of the curve evaluation code. Unlike
	// evalBezier, this is only designed to work on 4 control points.
	// Furthermore, it requires you to specify an initial binormal
	// Binit, which is iteratively propagated throughout the curve as
	// the curvepoints are generated. Any other function that creates
	// cubic splines can use this function by a corresponding change
	// of basis.
	Curve coreBezier(const Vec3f& p0,
		const Vec3f& p1,
		const Vec3f& p2,
		const Vec3f& p3,
		const Vec3f& Binit,
		unsigned steps) {

		Curve R(steps + 1);

		// YOUR CODE HERE (R1): build the basis matrix and loop the given number of steps,
		// computing points on the spline

		Mat4f B;
		// 
		B.setCol(0, Vec4f(1.0, 0.0, 0, 0));
		B.setCol(1, Vec4f(-3.0, 3.0, 0, 0));
		B.setCol(2, Vec4f(3.0, -6.0, 3.0, 0));
		B.setCol(3, Vec4f(-1.0, 3.0, -3.0, 1));



		Mat4f G;
		G.setCol(0, Vec4f(p0, 0)); G.setCol(1, Vec4f(p1, 0));
		G.setCol(2, Vec4f(p2, 0)); G.setCol(3, Vec4f(p3, 0));

		Mat4f GB = G * B;
		for (unsigned i = 0; i <= steps; ++i) {
			float t = (1.0 / steps) * (float)i;

			Vec4f powBasis(1.0, t, pow(t, 2), pow(t, 3));
			Vec4f newP = GB * powBasis;
			R[i].V = newP.getXYZ();
		}

		return R;
	}

} // namespace

Curve coreBezier(const Vec3f& p0,
	const Vec3f& p1,
	const Vec3f& p2,
	const Vec3f& p3,
	const Vec3f& Binit,
	const float begin, const float end, const float errorbound, const float minstep) {

	// YOUR CODE HERE(EXTRA): Adaptive tessellation

	return Curve();
}

// the P argument holds the control points and steps gives the amount of uniform tessellation.
// the rest of the arguments are for the adaptive tessellation extra.
Curve evalBezier(const vector<Vec3f>& P, unsigned steps, bool adaptive, float errorbound, float minstep) {
	// Check
	if (P.size() < 4 || P.size() % 3 != 1) {
		cerr << "evalBezier must be called with 3n+1 control points." << endl;
		_CrtDbgBreak();
		exit(0);
	}

	// YOUR CODE HERE (R1):
	// You should implement this function so that it returns a Curve
	// (e.g., a vector<CurvePoint>).  The variable "steps" tells you
	// the number of points to generate on each piece of the spline.
	// At least, that's how the sample solution is implemented and how
	// the SWP files are written.  But you are free to interpret this
	// variable however you want, so long as you can control the
	// "resolution" of the discretized spline curve with it.

	int n = (P.size() - 1) / 3;
	Curve newCurve(n*steps + 1);
	Curve R(steps + 1);
	for (unsigned i = 0; i < n; i++) {
		R = coreBezier(P[3 * i], P[3 * i + 1], P[3 * i + 2], P[3 * i + 3], Vec3f(0, 0, 0), steps);
		for (unsigned s = 0; s < steps; s++) {
			newCurve[i*steps + s] = R[s];
		}
	}
	newCurve[n*steps] = R[steps];

	// EXTRA CREDIT NOTE:
	// Also compute the other Vec3fs for each CurvePoint: T, N, B.
	// A matrix [N, B, T] should be unit and orthogonal.
	// Also note that you may assume that all Bezier curves that you
	// receive have G1 continuity. The T, N and B vectors will not
	// have to be defined at points where this does not hold.

	cerr << "\t>>> evalBezier has been called with the following input:" << endl;

	cerr << "\t>>> Control points (type vector<Vec3f>): " << endl;
	for (unsigned i = 0; i < P.size(); ++i) {
		cerr << "\t>>> "; printTranspose(P[i]); cerr << endl;
	}


	return newCurve;
}

// the P argument holds the control points and steps gives the amount of uniform tessellation.
// the rest of the arguments are for the adaptive tessellation extra.
Curve evalBspline(const vector<Vec3f>& P, unsigned steps, bool adaptive, float errorbound, float minstep) {
	// Check
	if (P.size() < 4) {
		cerr << "evalBspline must be called with 4 or more control points." << endl;
		exit(0);
	}

	// YOUR CODE HERE (R2):
	// We suggest you implement this function via a change of basis from
	// B-spline to Bezier.  That way, you can just call your evalBezier function.

	cerr << "\t>>> evalBSpline has been called with the following input:" << endl;

	cerr << "\t>>> Control points (type vector< Vec3f >): " << endl;
	for (unsigned i = 0; i < P.size(); ++i) {
		cerr << "\t>>> "; printTranspose(P[i]); cerr << endl;
	}

	unsigned n = P.size() - 3;


	Mat4f Bbz, Bbs;
	Bbz.setCol(0, Vec4f(1.0, 0.0, 0, 0));
	Bbz.setCol(1, Vec4f(-3.0, 3.0, 0, 0));
	Bbz.setCol(2, Vec4f(3.0, -6.0, 3.0, 0));
	Bbz.setCol(3, Vec4f(-1.0, 3.0, -3.0, 1));



	Bbs.setCol(0, Vec4f(1.0, 4.0, 1.0, 0.0) / 6.0);
	Bbs.setCol(1, Vec4f(-3.0, 0.0, 3.0, 0.0) / 6.0);
	Bbs.setCol(2, Vec4f(3.0, -6.0, 3.0, 0.0) / 6.0);
	Bbs.setCol(3, Vec4f(-1.0, 3.0, -3.0, 1.0) / 6.0);

	Mat4f Bs2Bz = Bbs * Bbz.inverted();
	Mat4f G;

	Curve newCuvre(n*steps + 1);
	Curve R(steps + 1);
	for (unsigned i = 0; i < n; i++) {
		G.setCol(0, Vec4f(P[i], 0));
		G.setCol(1, Vec4f(P[i + 1], 0));
		G.setCol(2, Vec4f(P[i + 2], 0));
		G.setCol(3, Vec4f(P[i + 3], 0));

		Mat4f G1 = G * Bs2Bz;
		R = coreBezier(Vec4f(G1.getCol(0)).getXYZ(), Vec4f(G1.getCol(1)).getXYZ(),
			Vec4f(G1.getCol(2)).getXYZ(), Vec4f(G1.getCol(3)).getXYZ(), Vec3f(0, 0, 0), steps);

		for (unsigned s = 0; s < steps; s++) {
			newCuvre[i*steps + s] = R[s];
		}
	}
	newCuvre[n*steps] = R[steps];

	return newCuvre;
}

Curve evalCircle(float radius, unsigned steps) {
	// This is a sample function on how to properly initialize a Curve
	// (which is a vector<CurvePoint>).

	// Preallocate a curve with steps+1 CurvePoints
	Curve R(steps + 1);

	// Fill it in counterclockwise
	for (unsigned i = 0; i <= steps; ++i) {
		// step from 0 to 2pi
		float t = 2.0f * (float)M_PI * float(i) / steps;

		// Initialize position
		// We're pivoting counterclockwise around the y-axis
		R[i].V = radius * Vec3f(FW::cos(t), FW::sin(t), 0);

		// Tangent vector is first derivative
		R[i].T = Vec3f(-FW::sin(t), FW::cos(t), 0);

		// Normal vector is second derivative
		R[i].N = Vec3f(-FW::cos(t), -FW::sin(t), 0);

		// Finally, binormal is facing up.
		R[i].B = Vec3f(0, 0, 1);
	}

	return R;
}

void drawCurve(const Curve& curve, float framesize) {
	// Save current state of OpenGL
	glPushAttrib(GL_ALL_ATTRIB_BITS);

	// Setup for line drawing
	glDisable(GL_LIGHTING);
	glColor4f(1, 1, 1, 1);
	glLineWidth(1);

	if (framesize >= 0) {
		// Draw curve
		glBegin(GL_LINE_STRIP);
		for (unsigned i = 0; i < curve.size(); ++i) {
			glVertex(curve[i].V);
		}
		glEnd();
	}

	glLineWidth(1);

	// Draw coordinate frames if framesize nonzero
	if (framesize != 0.0f) {
		framesize = FW::abs(framesize);
		Mat4f M;

		for (unsigned i = 0; i < curve.size(); ++i) {
			M.setCol(0, Vec4f(curve[i].N, 0));
			M.setCol(1, Vec4f(curve[i].B, 0));
			M.setCol(2, Vec4f(curve[i].T, 0));
			M.setCol(3, Vec4f(curve[i].V, 1));

			glPushMatrix();
			glMultMatrixf(M.getPtr());
			glScaled(framesize, framesize, framesize);
			glBegin(GL_LINES);
			glColor3f(1, 0, 0); glVertex3d(0, 0, 0); glVertex3d(1, 0, 0);
			glColor3f(0, 1, 0); glVertex3d(0, 0, 0); glVertex3d(0, 1, 0);
			glColor3f(0, 0, 1); glVertex3d(0, 0, 0); glVertex3d(0, 0, 1);
			glEnd();
			glPopMatrix();
		}
	}

	// Pop state
	glPopAttrib();
}

