package internal

import "github.com/ungerik/go3d/float64/vec3"

type HomoPoint struct {
	Vec3 vec3.T
	W    float64
}

func (this *HomoPoint) Add(pt *HomoPoint) *HomoPoint {
	this.Vec3.Add(&pt.Vec3)
	this.W += pt.W

	return this
}

func (this *HomoPoint) Scale(scale float64) *HomoPoint {
	this.Vec3.Scale(scale)
	this.W *= scale

	return this
}

func Homogenized(pt vec3.T, w float64) HomoPoint {
	return HomoPoint{pt.Scaled(w), w}
}

// Transform a 1d array of points into their homogeneous equivalents
//
// **params**
// + 1d array of control points, (actually a 2d array of size (m x dim) )
// + array of control point weights, the same size as the array of control points (m x 1)
//
// **returns**
// + 1d array of control points where each point is (wi*pi, wi) where wi
// i the ith control point weight and pi is the ith control point,
// hence the dimension of the point is dim + 1
func Homogenize1d(pts []vec3.T, weights []float64) []HomoPoint {
	homoPts := make([]HomoPoint, 0, len(pts))
	for i, pt := range pts {
		homoPts = append(homoPts, Homogenized(pt, weights[i]))
	}

	return homoPts
}

// **params**
// + 2d array of control points, (actually a 3d array of size m x n x dim)
// + array of control point weights, the same size as the control points array (m x n x 1)
//
// **returns**
// + 1d array of control points where each point is (wi*pi, wi) where wi
// i the ith control point weight and pi is the ith control point, the size is
// (m x n x dim+1)
func Homogenize2d(pts [][]vec3.T, weights [][]float64) [][]HomoPoint {
	homoPts := make([][]HomoPoint, len(pts))
	for i := range homoPts {
		homoPts[i] = Homogenize1d(pts[i], weights[i])
	}

	return homoPts
}

// Dehomogenize a point
//
// **params**
// + a point represented by an array (wi*pi, wi) with length (dim+1)
//
// **returns**
// + a point represented by an array pi with length (dim)
func (this *HomoPoint) Dehomogenized() vec3.T {
	return this.Vec3.Scaled(1 / this.W)
}

// Dehomogenize an array of points
//
// **params**
// + array of points represented by an array (wi*pi, wi) with length (dim+1)
//
// **returns**
// + an array of points, each of length dim
func Dehomogenize1d(homoPoints []HomoPoint) []vec3.T {
	result := make([]vec3.T, 0, len(homoPoints))
	for _, homoPt := range homoPoints {
		result = append(result, homoPt.Dehomogenized())
	}

	return result
}

// Dehomogenize a 2d array of pts
//
// **params**
// + array of arrays of points represented by an array (wi*pi, wi) with length (dim+1)
//
// **returns**
// + array of arrays of points, each of length dim
func Dehomogenize2d(homoPoints [][]HomoPoint) [][]vec3.T {
	result := make([][]vec3.T, len(homoPoints))
	for i := range result {
		result[i] = Dehomogenize1d(homoPoints[i])
	}

	return result
}

// Obtain the weight from a collection of points in homogeneous space, assuming all
// are the same dimension
//
// **params**
// + array of points represented by an array (wi*pi, wi) with length (dim+1)
//
// **returns**
// + a point represented by an array pi with length (dim)
func Weight1d(homoPoints []HomoPoint) (weights []float64) {
	weights = make([]float64, len(homoPoints))
	for i := range weights {
		weights[i] = homoPoints[i].W
	}

	return
}

// Obtain the weight from a collection of points in homogeneous space, assuming all
// are the same dimension
//
// **params**
// + array of arrays of of points represented by an array (wi*pi, wi) with length (dim+1)
//
// **returns**
// +  array of arrays of points, each represented by an array pi with length (dim)
func Weight2d(homoPoints [][]HomoPoint) (weights [][]float64) {
	weights = make([][]float64, len(homoPoints))
	for i := range weights {
		weights[i] = Weight1d(homoPoints[i])
	}

	return
}

func HomoInterpolated(hpt0, hpt1 *HomoPoint, t float64) HomoPoint {
	return HomoPoint{
		vec3.Interpolate(&hpt0.Vec3, &hpt1.Vec3, t),
		(1-t)*hpt0.W + t*hpt1.W,
	}
}
