package make

import (
	. "github.com/alexozer/verb/internal"

	"github.com/alexozer/verb"
	"github.com/ungerik/go3d/float64/mat4"
	"github.com/ungerik/go3d/float64/vec3"
)

// Generate a surface by translating a profile curve along a rail curve
//
// **params**
// + profile NurbsCurveData
// + rail NurbsCurveData
//
// **returns**
// + NurbsSurfaceData object
func SweptSurface(profile, rail *verb.NurbsCurve) *verb.NurbsSurface {
	railKnots := rail.Knots()

	startu := railKnots[0]
	endu := railKnots[len(railKnots)-1]
	pt0 := rail.Point(startu)

	numSamples := 2 * len(rail.ControlPoints())
	span := (endu - startu) / (float64(numSamples) - 1)

	crvs := make([]*verb.NurbsCurve, numSamples)

	for i := range crvs {
		pt := rail.Point(startu + float64(i)*span)
		pt.Sub(&pt0)

		mat := mat4.Ident
		mat.SetTranslation(&pt)
		crvs[i] = profile.Transform(&mat)
	}

	return loftedSurface(crvs, nil)
}

// Construct a NurbsSurface by lofting between a collection of curves
//
// **params**
// + A collection of curves
//
// **returns**
// + A new NurbsSurface
func loftedSurface(curves []*verb.NurbsCurve, _degreeV *int) *verb.NurbsSurface {
	curves = verb.UnifyCurveKnotVectors(curves)

	// degree
	degreeU := curves[0].Degree()

	var degreeV int
	if _degreeV == nil {
		degreeV = 3
	} else {
		degreeV = *_degreeV
	}

	if degreeV > len(curves)-1 {
		degreeV = len(curves) - 1
	}

	// knots
	knotsU := curves[0].Knots()
	var knotsV KnotVec

	crvCtrlPts := make([][]vec3.T, len(curves))
	for i := range crvCtrlPts {
		crvCtrlPts[i] = curves[i].ControlPoints()
	}

	controlPoints := make([][]vec3.T, len(crvCtrlPts[0]))
	weights := make([][]float64, len(controlPoints))
	for i := range crvCtrlPts[0] {
		// extract the ith control pt of each curve
		points := make([]vec3.T, len(crvCtrlPts))
		for j, ctrlPts := range crvCtrlPts {
			points[j] = ctrlPts[i]
		}

		// construct an interpolating curve using this list
		_, controlPoints[i], weights[i], knotsV =
			interpCurve(points, degreeV, nil, nil)
		// computing knotsV is redundant
	}

	return verb.NewNurbsSurfaceUnchecked(degreeU, degreeV, controlPoints, weights, knotsU, knotsV)
}

func NurbsCurveByPoints(points []vec3.T, _degree *int) *verb.NurbsCurve {
	var degree int
	if _degree == nil {
		degree = 3
	} else {
		degree = *_degree
	}

	return verb.NewNurbsCurveUnchecked(interpCurve(points, degree, nil, nil))
}

func interpCurve(points []vec3.T, degree int, _startTangent, _endTangent *vec3.T) (deg int, controlPoints []vec3.T, weights []float64, knots KnotVec) {
	// 0) build knot vector for curve by normalized chord length
	// 1) construct effective basis function in square matrix (W)
	// 2) construct set of coordinattes to interpolate vector (p)
	// 3) set of control points (c)

	// Wc = p

	// 4) solve for c in all 3 dimensions

	if len(points) < degree+1 {
		panic("Must supply at least degree + 1 points")
	}

	us := make([]float64, len(points))
	for i := 1; i < len(points); i++ {
		chord := vec3.Distance(&points[1], &points[i-1])
		us[i] = us[i-1] + chord
	}

	// normalize
	max := us[len(us)-1]
	for i := range us {
		us[i] /= max
	}

	// we need two more control points, two more knots

	hasTangents := _startTangent != nil && _endTangent != nil
	var start, end int
	if hasTangents {
		end = len(us) - degree + 1
	} else {
		start = 1
		end = len(us) - degree
	}

	knots = make(KnotVec, 2*(degree+1)+(end-start))
	middleKnots := knots[degree+1 : len(knots)-(degree+1)]

	for i := range middleKnots {
		var weightSums float64
		for j := 0; j < degree; j++ {
			weightSums += us[i+j]
		}

		middleKnots[i] = (1.0 / float64(degree) * weightSums)
	}

	for i := (degree + 1) + len(middleKnots); i < len(knots); i++ {
		knots[i] = 1
	}

	// build matrix of basis function coeffs (TODO: use sparse rep)
	oldA := make(Matrix, len(us)+2)
	A := oldA[1 : len(oldA)-1]

	var n, ld int
	if hasTangents {
		n = len(points) + 1
		ld = len(points) - (degree - 1)
	} else {
		n = len(points) - 1
		ld = len(points) - (degree + 1)
	}

	for i, u := range us {
		span := knots.SpanGivenN(n, degree, u)
		basisFuncs := BasisFunctionsGivenKnotSpanIndex(span, u, degree, knots)

		row := make([]float64, ld+len(basisFuncs))

		ls := span - degree
		copy(row[ls:], basisFuncs)
		A[i] = row
	}

	if hasTangents {
		tanRow0 := make([]float64, len(A[0]))
		tanRow1 := make([]float64, len(A[0]))

		tanRow0[0], tanRow0[1] = -1, 1
		tanRow1[len(tanRow1)-2], tanRow1[len(tanRow1)-1] = -1, 1

		A = oldA
		A[0], A[1] = A[1], tanRow0
		A[len(A)-1] = tanRow1
	}

	// for each dimension, solve
	xs := make(Matrix, 3)

	mult1 := (1 - knots[len(knots)-degree-2]) / float64(degree)
	mult0 := knots[degree+1] / float64(degree)

	for i := range xs {
		var b []float64

		if !hasTangents {
			b = make([]float64, len(points))
			for j := range b {
				b[j] = points[j][i]
			}
		} else {
			b = make([]float64, len(points)+2)

			// insert the tangents at the second and second to last index
			b[0] = points[0][i]
			b[1] = mult0 * (*_startTangent)[i]
			for j := 1; j < len(points)-1; j++ {
				b[j+1] = points[j][i]
			}
			b[len(b)-2] = mult1 * (*_endTangent)[i]
			b[len(b)-1] = points[len(points)-1][i]
		}

		xs[i] = A.Solve(b)
	}

	controlPoints = make([]vec3.T, len(xs[0]))
	for i := range xs[0] {
		controlPoints[i][0] = xs[0][i]
		controlPoints[i][1] = xs[1][i]
		controlPoints[i][2] = xs[2][i]
	}

	weights = make([]float64, len(controlPoints))
	for i := range weights {
		weights[i] = 1
	}

	deg = degree

	return
}
