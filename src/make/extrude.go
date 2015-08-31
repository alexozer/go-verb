package make

import (
	"github.com/alexozer/verb"
	"github.com/ungerik/go3d/float64/vec3"
)

// Generate the control points, weights, and knots of an extruded surface
//
// **params**
// + axis of the extrusion
// + length of the extrusion
// + a NurbsCurveData object representing a NURBS surface
//
// **returns**
// + an object with the following properties: controlPoints, weights, knots, degree
func ExtrudedSurface(axis *vec3.T, length float64, profile *verb.NurbsCurve) *verb.NurbsSurface {
	profControlPoints := profile.ControlPoints()
	profWeights := profile.Weights()

	controlPoints, weights := make([][]vec3.T, 3), make([][]float64, 3)
	for i := range controlPoints {
		controlPoints[i] = make([]vec3.T, len(profControlPoints))
		weights[i] = make([]float64, len(profWeights))
	}

	translation := axis.Scaled(length)
	halfTranslation := translation.Scaled(0.5)

	// original control points
	for j := range profControlPoints {
		controlPoints[2][j] = profControlPoints[j]
		controlPoints[1][j] = vec3.Add(&halfTranslation, &profControlPoints[j])
		controlPoints[0][j] = vec3.Add(&translation, &profControlPoints[j])

		weights[0][j] = profWeights[j]
		weights[1][j] = profWeights[j]
		weights[2][j] = profWeights[j]
	}

	return verb.NewNurbsSurfaceUnchecked(
		2, profile.Degree(),
		controlPoints, weights,
		[]float64{0, 0, 0, 1, 1, 1}, profile.Knots(),
	)
}

// Generate the control points, weights, and knots of a cylinder
//
// **params**
// + normalized axis of cylinder
// + xaxis in plane of cylinder
// + position of base of cylinder
// + height from base to top
// + radius of the cylinder
//
// **returns**
// + an object with the following properties: controlPoints, weights, knotsU, knotsV, degreeU, degreeV
func CylindricalSurface(axis, xaxis *vec3.T, base *vec3.T, height, radius float64) *verb.NurbsSurface {
	yaxis := vec3.Cross(axis, xaxis)
	circ := Circle(base, xaxis, &yaxis, radius)

	return ExtrudedSurface(axis, height, circ)

}
