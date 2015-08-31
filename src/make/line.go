package make

import (
	"github.com/alexozer/verb"
	"github.com/ungerik/go3d/float64/vec3"
)

func Line(first, last *vec3.T) *verb.NurbsCurve {
	return Polyline([]vec3.T{*first, *last})
}

// Generate the control points, weights, and knots of a polyline curve
//
// **params**
// + array of points in curve
//
// **returns**
// + a NurbsCurveData object representing a NURBS curve
func Polyline(pts []vec3.T) *verb.NurbsCurve {
	knots := make([]float64, len(pts)+1)

	var lsum float64
	for i := 0; i < len(pts)-1; i++ {
		lsum += vec3.Distance(&pts[i], &pts[i+1])
		knots[i+2] = lsum
	}
	knots[len(knots)-1] = lsum

	// normalize the knot array
	for i := range knots {
		knots[i] /= lsum
	}

	weights := make([]float64, len(pts))
	for i := range weights {
		weights[i] = 1
	}

	return verb.NewNurbsCurveUnchecked(1, pts, weights, knots)
}
