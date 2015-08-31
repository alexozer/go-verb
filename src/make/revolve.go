package make

import (
	"math"

	"github.com/alexozer/verb"
	"github.com/alexozer/verb/internal"
	"github.com/alexozer/verb/intersect"
	"github.com/ungerik/go3d/float64/vec3"
)

// Generate the control points, weights, and knots of a revolved surface
// (Corresponds to Algorithm A7.1 from Piegl & Tiller)
//
// **params**
// + center of the rotation axis
// + axis of the rotation axis
// + angle to revolve around axis
// + degree of the generatrix
// + control points of the generatrix
// + weights of the generatrix
//
// **returns**
// + an object with the following properties: controlPoints, weights, knots, degree
func RevolvedSurface(profile *verb.NurbsCurve, center *vec3.T, axis *vec3.T, theta float64) *verb.NurbsSurface {
	prof_controlPoints := profile.ControlPoints()
	prof_weights := profile.Weights()

	var narcs int
	var knotsU []float64

	switch {
	case theta <= math.Pi/2:
		{ // less than 90
			narcs = 1
			knotsU = make([]float64, 6+2*(narcs-1))
		}
	case theta <= math.Pi:
		{ // between 90 and 180
			narcs = 2
			knotsU = make([]float64, 6+2*(narcs-1))
			knotsU[3], knotsU[4] = 0.5, 0.5
		}
	case theta <= 3*math.Pi/2:
		{ // between 180 and 270
			narcs = 3
			knotsU = make([]float64, 6+2*(narcs-1))
			knotsU[3], knotsU[4] = 1/3, 1/3
			knotsU[5], knotsU[6] = 2/3, 2/3
		}
	default:
		{ // between 270 and 360
			narcs = 4
			knotsU = make([]float64, 6+2*(narcs-1))
			knotsU[3], knotsU[4] = 1/4, 1/4
			knotsU[5], knotsU[6] = 1/2, 1/2
			knotsU[7], knotsU[8] = 3/4, 3/4
		}
	}

	dtheta := theta / float64(narcs) // divide the interval into several points
	j := 3 + 2*(narcs-1)

	// initialize the start and end knots
	// keep in mind that we only return the knot vector for thes
	for i := 0; i < 3; i++ {
		knotsU[j+i] = 1
	}

	// do some initialization
	wm := math.Cos(dtheta / 2)
	sines, cosines := make([]float64, narcs+1), make([]float64, narcs+1)

	controlPoints := make([][]vec3.T, 2*narcs+1)
	for i := range controlPoints {
		controlPoints[i] = make([]vec3.T, len(prof_controlPoints))
	}

	weights := make([][]float64, 2*narcs+1)
	for i := range weights {
		weights[i] = make([]float64, len(prof_controlPoints))
	}

	// initialize the sines and cosines
	var angle float64
	for i := 1; i <= narcs; i++ {
		angle += dtheta
		cosines[i] = math.Cos(angle)
		sines[i] = math.Sin(angle)
	}

	// for each pt in the generatrix
	// i.e. for each row of the 2d knot vectors
	for j := range prof_controlPoints {
		// get the closest point of the generatrix point on the axis
		O := rayClosestPoint(prof_controlPoints[j], center, axis)
		// X is the vector from the axis to generatrix control pt
		X := vec3.Sub(&prof_controlPoints[j], &O)
		// radius at that height
		r := X.Length()
		// Y is perpendicular to X and axis, and complete the coordinate system
		Y := vec3.Cross(axis, &X)

		if r > internal.Epsilon {
			X.Scale(1 / r)
			Y.Scale(1 / r)
		}

		// the first row of controlPoints and weights is just the generatrix
		controlPoints[0][j] = prof_controlPoints[j]
		P0 := prof_controlPoints[j]
		weights[0][j] = prof_weights[j]

		// store T0 as the Y vector
		var T0 = Y
		var index int

		// proceed around the circle
		for i := 1; i <= narcs; i++ {
			// O + r * cos(theta) * X + r * sin(theta) * Y
			// rotated generatrix pt
			var P2 vec3.T
			if r == 0 {
				P2 = O
			} else {
				xCompon := X.Scaled(r * cosines[i])
				yCompon := Y.Scaled(r * sines[i])
				offset := xCompon.Add(&yCompon)
				P2 = vec3.Add(&O, offset)
			}

			controlPoints[index+2][j] = P2
			weights[index+2][j] = prof_weights[j]

			// construct the vector tangent to the rotation
			temp0 := Y.Scaled(cosines[i])
			temp1 := X.Scaled(sines[i])
			T2 := temp0.Sub(&temp1)

			// construct the next control pt
			if r == 0 {
				controlPoints[index+1][j] = O
			} else {
				T0Norm := T0.Normalized()
				T2Norm := T2.Normalized()
				inters := intersect.Rays(&P0, &T0Norm, &P2, &T2Norm)

				T0Scaled := T0.Scaled(inters.U0)
				P1 := T0Scaled.Add(&P0)

				controlPoints[index+1][j] = *P1
			}

			weights[index+1][j] = wm * prof_weights[j]

			index += 2

			if i < narcs {
				P0 = P2
				T0 = *T2
			}
		}
	}

	return verb.NewNurbsSurfaceUnchecked(2, profile.Degree(), controlPoints, weights, knotsU, profile.Knots())
}

//
// Generate the control points, weights, and knots of a sphere
//
// **params**
// + the center of the sphere
// + normalized axis of sphere
// + vector perpendicular to axis of sphere, starting the rotation of the sphere
// + radius of the sphere
//
// **returns**
// + an object with the following properties: controlPoints, weights, knotsU, knotsV, degreeU, degreeV
//
func SphericalSurface(center *vec3.T, axis, xaxis *vec3.T, radius float64) *verb.NurbsSurface {
	invAxis := axis.Inverted()
	arc := Arc(center, &invAxis, xaxis, radius, 0, math.Pi)

	return RevolvedSurface(arc, center, axis, 2*math.Pi)
}

//
// Generate the control points, weights, and knots of a cone
//
// **params**
// + normalized axis of cone
// + position of base of cone
// + height from base to tip
// + radius at the base of the cone
//
// **returns**
// + an object with the following properties: controlPoints, weights, knots, degree
//
func ConicalSurface(axis, xaxis *vec3.T, base *vec3.T, height, radius float64) *verb.NurbsSurface {
	angle := 2 * math.Pi
	profDegree := 1
	heightCompon := axis.Scaled(height)
	radiusCompon := xaxis.Scaled(radius)
	profCtrlPts := []vec3.T{vec3.Add(base, &heightCompon), vec3.Add(base, &radiusCompon)}
	profKnots := []float64{0, 0, 1, 1}
	profWeights := []float64{1, 1}
	prof := verb.NewNurbsCurveUnchecked(profDegree, profCtrlPts, profWeights, profKnots)

	return RevolvedSurface(prof, base, axis, angle)
}
