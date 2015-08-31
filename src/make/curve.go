package make

import (
	"math"

	"github.com/alexozer/verb"
	"github.com/alexozer/verb/intersect"
	"github.com/ungerik/go3d/float64/vec3"
)

// Generate the control points, weights, and knots of an arbitrary arc
// (Corresponds to Algorithm A7.1 from Piegl & Tiller)
//
// **params**
// + the center of the arc
// + the xaxis of the arc
// + orthogonal yaxis of the arc
// + radius of the arc
// + start angle of the arc, between 0 and 2pi
// + end angle of the arc, between 0 and 2pi, greater than the start angle
//
// **returns**
// + a NurbsCurveData object representing a NURBS curve
func Arc(center *vec3.T, xaxis, yaxis *vec3.T, radius float64, startAngle, endAngle float64) *verb.NurbsCurve {
	xaxisScaled, yaxisScaled := xaxis.Scaled(radius), yaxis.Scaled(radius)
	return EllipseArc(center, &xaxisScaled, &yaxisScaled, startAngle, endAngle)
}

// Create a circle
//
// **params**
// + Length 3 array representing the center of the circle
// + Length 3 array representing the xaxis
// + Length 3 array representing the perpendicular yaxis
// + Radius of the circle
func Circle(center *vec3.T, xaxis, yaxis *vec3.T, radius float64) *verb.NurbsCurve {
	return Arc(center, xaxis, yaxis, radius, 0, 2*math.Pi)
}

func Ellipse(center *vec3.T, xaxis, yaxis *vec3.T) *verb.NurbsCurve {
	return EllipseArc(center, xaxis, yaxis, 0, 2*math.Pi)
}

// Generate the control points, weights, and knots of an elliptical arc
//
// **params**
// + the center
// + the scaled x axis
// + the scaled y axis
// + start angle of the ellipse arc, between 0 and 2pi, where 0 points at the xaxis
// + end angle of the arc, between 0 and 2pi, greater than the start angle
//
// **returns**
// + a NurbsCurveData object representing a NURBS curve
func EllipseArc(center *vec3.T, xaxis, yaxis *vec3.T, startAngle, endAngle float64) *verb.NurbsCurve {
	xradius, yradius := xaxis.Length(), yaxis.Length()

	xaxisNorm, yaxisNorm := xaxis.Normalized(), yaxis.Normalized()

	// if the end angle is less than the start angle, do a circle
	if endAngle < startAngle {
		endAngle = 2.0*math.Pi + startAngle
	}

	theta := endAngle - startAngle

	// how many arcs?
	var numArcs int
	if theta <= math.Pi/2 {
		numArcs = 1
	} else {
		if theta <= math.Pi {
			numArcs = 2
		} else if theta <= 3*math.Pi/2 {
			numArcs = 3
		} else {
			numArcs = 4
		}
	}

	dtheta := theta / float64(numArcs)
	w1 := math.Cos(dtheta / 2)

	xCompon := xaxisNorm.Scaled(xradius * math.Cos(startAngle))
	yCompon := yaxisNorm.Scaled(yradius * math.Sin(startAngle))
	P0 := vec3.Add(&xCompon, &yCompon)

	temp0 := yaxisNorm.Scaled(math.Cos(startAngle))
	temp1 := xaxisNorm.Scaled(math.Sin(startAngle))
	T0 := vec3.Sub(&temp0, &temp1)

	controlPoints := make([]vec3.T, 2*numArcs+1)
	knots := make([]float64, 2*numArcs+3)
	index := 0
	angle := startAngle
	weights := make([]float64, numArcs*2)

	controlPoints[0] = P0
	weights[0] = 1.0

	for i := 1; i <= numArcs; i++ {
		angle += dtheta
		xCompon = xaxisNorm.Scaled(xradius * math.Cos(angle))
		yCompon = yaxisNorm.Scaled(yradius * math.Sin(angle))
		offset := vec3.Add(&xCompon, &yCompon)
		P2 := vec3.Add(center, &offset)

		weights[index+2] = 1
		controlPoints[index+2] = P2

		temp0 := yaxisNorm.Scaled(math.Cos(angle))
		temp1 := xaxisNorm.Scaled(math.Sin(angle))
		T2 := vec3.Sub(&temp0, &temp1)

		T0Norm := T0.Normalized()
		T2Norm := T2.Normalized()
		inters := intersect.Rays(&P0, &T0Norm, &P2, &T2Norm)

		T0Scaled := T0.Scaled(inters.U0)
		P1 := vec3.Add(&P0, &T0Scaled)

		weights[index+1] = w1
		controlPoints[index+1] = P1

		index += 2

		if i < numArcs {
			P0 = P2
			T0 = T2
		}
	}

	j := 2*numArcs + 1

	for i := 0; i < 3; i++ {
		knots[i] = 0.0
		knots[i+j] = 1.0
	}

	switch numArcs {
	case 2:
		knots[3] = 0.5
		knots[4] = 0.5
	case 3:
		knots[3] = 1 / 3
		knots[4] = 1 / 3

		knots[5] = 2 / 3
		knots[6] = 2 / 3
	case 4:
		knots[3] = 0.25
		knots[4] = 0.25

		knots[5] = 0.5
		knots[6] = 0.5

		knots[7] = 0.75
		knots[8] = 0.75
	}

	return verb.NewNurbsCurveUnchecked(2, controlPoints, weights, knots)
}

// generate the control points, weights, and knots for a bezier curve of any degree
//
// **params**
// + first point in counter-clockwise form
// + second point in counter-clockwise form
// + third point in counter-clockwise form
// + forth point in counter-clockwise form
//
// **returns**
// + nurbssurfacedata object
func BezierCurve(controlPoints []vec3.T) *verb.NurbsCurve {
	degree := len(controlPoints) - 1

	// build uniform weights
	weights := make([]float64, len(controlPoints))
	for i := range weights {
		weights[i] = 1
	}

	knots := make([]float64, 2*degree+2)
	for i := range knots[len(knots)/2:] {
		knots[i] = 1
	}

	return verb.NewNurbsCurveUnchecked(degree, controlPoints, weights, knots)
}
