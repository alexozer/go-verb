package internal

import "github.com/ungerik/go3d/float64/vec3"

type Ray struct {
	Origin, Dir vec3.T
}

// Find the closest point on a ray
//
// **params**
// + point to project
// + origin for ray
// + direction of ray 1, assumed normalized
//
// **returns**
// + pt
func (this Ray) ClosestPoint(pt vec3.T) vec3.T {
	o2pt := vec3.Sub(&pt, &this.Origin)
	do2ptr := vec3.Dot(&o2pt, &this.Dir)
	dirScaled := this.Dir.Scaled(do2ptr)
	proj := vec3.Add(&this.Origin, &dirScaled)

	return proj
}

// Find the distance of a point to a ray
//
// **params**
// + point to project
// + origin for ray
// + direction of ray 1, assumed normalized
//
// **returns**
// + the distance
func (this Ray) DistToPoint(pt vec3.T) float64 {
	d := this.ClosestPoint(pt)

	return vec3.Distance(&d, &pt)
}
