package intersect

import (
	"math"

	"github.com/alexozer/verb"
	"github.com/ungerik/go3d/float64/vec3"
)

//
// Get min coordinate on an axis
//
// **params**
// + array of length 3 arrays of numbers representing the points
// + length 3 array of point indices for the triangle
// + index of the axis to test - 0 for x, 1 for y, 2 for z
//
// **returns**
// + the minimum coordinate
//
func minCoordOnAxis(points []vec3.T, tri *verb.Tri, axis int) float64 {
	min := math.Inf(1)

	for _, iPt := range tri {
		if coord := points[iPt][axis]; coord < min {
			min = coord
		}
	}

	return min
}

//
// Get triangle normal
//
// **params**
// + array of length 3 arrays of numbers representing the points
// + length 3 array of point indices for the triangle
//
// **returns**
// + a normal vector represented by an array of length 3
//
func TriangleNormal(points []vec3.T, tri *verb.Tri) vec3.T {
	v0 := points[tri[0]]
	v1 := points[tri[1]]
	v2 := points[tri[2]]

	v1.Sub(&v0)
	v2.Sub(&v0)
	n := vec3.Cross(&v1, &v2)

	return *n.Normalize()
}

//
// Get triangle centroid
//
// **params**
// + array of length 3 arrays of numbers representing the points
// + length 3 array of point indices for the triangle
//
// **returns**
// + a point represented by an array of length 3
//
func TriangleCentroid(points []vec3.T, tri *verb.Tri) vec3.T {
	centroid := vec3.Zero

	for _, iPt := range tri {
		for j, compon := range points[iPt] {
			centroid[j] += compon
		}
	}

	for i := range centroid {
		centroid[i] /= 3
	}

	return centroid
}
