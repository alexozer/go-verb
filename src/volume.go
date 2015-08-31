package verb

import (
	. "github.com/alexozer/verb/internal"
	"github.com/ungerik/go3d/float64/vec3"
)

type UVW [3]float64

type volume struct {
	DegreeU int // integer degree in u direction
	DegreeV int // integer degree in v direction
	DegreeW int // integer degree in w direction

	KnotsU KnotVec // array of nondecreasing knot values in u direction
	KnotsV KnotVec // array of nondecreasing knot values in v direction
	KnotsW KnotVec // array of nondecreasing knot values in w direction

	// 3d array of control points, where rows are the u dir, and columns run along the positive v direction,
	// and where each control point is an array of length (dim)
	ControlPoints [][][]vec3.T
}

//
// Compute a point in a non-uniform, non-rational B spline volume
//
// **params**
// + VolumeData
// + u parameter at which to evaluate the volume point
// + v parameter at which to evaluate the volume point
// + w parameter at which to evaluate the volume point
//
// **returns**
// + a point represented by an array of length (dim)
func (this *volume) Point(uvw UVW) vec3.T {
	n := len(this.KnotsU) - this.DegreeU - 2
	m := len(this.KnotsV) - this.DegreeV - 2
	l := len(this.KnotsW) - this.DegreeW - 2

	return this.PointGivenNML(n, m, l, uvw)
}

//
// Compute a point in a non-uniform, non-rational B spline volume
//
// **params**
// + VolumeData
// + u parameter at which to evaluate the volume point
// + v parameter at which to evaluate the volume point
// + w parameter at which to evaluate the volume point
//
// **returns**
// + a point represented by an array of length (dim)
func (this *volume) PointGivenNML(n, m, l int, uvw UVW) vec3.T {
	if !areValidRelations(this.DegreeU, len(this.ControlPoints), len(this.KnotsU)) ||
		!areValidRelations(this.DegreeV, len(this.ControlPoints[0]), len(this.KnotsV)) ||
		!areValidRelations(this.DegreeW, len(this.ControlPoints[0][0]), len(this.KnotsW)) {
		panic("Invalid relations between control points and knot vector")
	}

	controlPoints := this.ControlPoints
	degreeU, degreeV, degreeW := this.DegreeU, this.DegreeV, this.DegreeW
	knotsU, knotsV, knotsW := this.KnotsU, this.KnotsV, this.KnotsW

	knotSpanIndexU := knotsU.SpanGivenN(n, degreeU, uvw[0])
	knotSpanIndexV := knotsV.SpanGivenN(m, degreeV, uvw[1])
	knotSpanIndexW := knotsW.SpanGivenN(l, degreeW, uvw[2])

	uBasisVals := BasisFunctionsGivenKnotSpanIndex(knotSpanIndexU, uvw[0], degreeU, knotsU)
	vBasisVals := BasisFunctionsGivenKnotSpanIndex(knotSpanIndexV, uvw[0], degreeV, knotsV)
	wBasisVals := BasisFunctionsGivenKnotSpanIndex(knotSpanIndexV, uvw[0], degreeW, knotsW)

	uind := knotSpanIndexU - degreeU
	var position, temp, temp2 vec3.T

	for i := 0; i <= degreeW; i++ {
		temp2 = vec3.Zero
		wind := knotSpanIndexW - degreeW + i

		for j := 0; j <= degreeV; j++ {
			temp = vec3.Zero
			vind := knotSpanIndexV - degreeV + j

			for k := 0; k <= degreeU; k++ {
				scaled := controlPoints[uind+k][vind][wind].Scaled(uBasisVals[k])
				temp.Add(&scaled)
			}

			// add weighted contribution of u isoline
			scaled := temp.Scaled(vBasisVals[j])
			temp2.Add(&scaled)
		}

		// add weighted contribution from uv isosurfaces
		scaled := temp2.Scaled(wBasisVals[i])
		position.Add(&scaled)
	}

	return position
}

// Confirm the relations between degree (p), number of control points(n+1), and the number of knots (m+1)
// via The NURBS Book (section 3.2, Second Edition)
//
// **params**
// + integer degree
// + integer number of control points
// + integer length of the knot Array (including duplicate knots)
//
// **returns**
// + whether the values are correct
func areValidRelations(degree, numControlPoints, knotsLength int) bool {
	return numControlPoints+degree+1 == knotsLength
}
