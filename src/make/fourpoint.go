package make

import (
	"github.com/alexozer/verb"
	"github.com/ungerik/go3d/float64/vec3"
)

// Generate the control points, weights, and knots of a surface defined by 4 points
//
// **params**
// + first point in counter-clockwise form
// + second point in counter-clockwise form
// + third point in counter-clockwise form
// + forth point in counter-clockwise form
//
// **returns**
// + NurbsSurfaceData object
func FourPointSurface(p1, p2, p3, p4 *vec3.T, _degree *int) *verb.NurbsSurface {
	var degree int
	if _degree == nil {
		degree = 3
	} else {
		degree = *_degree
	}
	degreeFloat := float64(degree)

	pts := make([][]vec3.T, degree+1)
	for i := range pts {
		iFloat := float64(i)

		row := make([]vec3.T, degree+1)
		for j := range row {
			l := 1 - iFloat/degreeFloat
			p1p2 := vec3.Interpolate(p1, p2, l)
			p4p3 := vec3.Interpolate(p4, p3, l)

			row[j] = vec3.Interpolate(&p1p2, &p4p3, 1-float64(j)/degreeFloat)
		}

		pts[i] = row
	}

	// build uniform weights
	weightRow := make([]float64, degree+1)
	for i := range weightRow {
		weightRow[i] = 1
	}
	weights := make([][]float64, degree+1)
	for i := range weights {
		weights[i] = weightRow
	}

	knots := make([]float64, 2*(degree+1))
	for i := degree + 1; i < len(knots); i++ {
		knots[i] = 1
	}

	return verb.NewNurbsSurfaceUnchecked(degree, degree, pts, weights, knots, knots)
}
