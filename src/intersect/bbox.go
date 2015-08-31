package intersect

import (
	"math"

	"github.com/ungerik/go3d/float64/vec3"
)

const BoundingBoxTolerance = 1e-4

// The zero value for BoundingBox is ready to use
type BoundingBox struct {
	Min, Max    vec3.T
	initialized bool
}

// Adds a point to the bounding box, expanding the bounding box if the point is outside of it.
// If the bounding box is not initialized, this method has that side effect.
//
// **params**
// + A length-n array of numbers
//
// **returns**
// + This BoundingBox for chaining
func (this *BoundingBox) Add(point *vec3.T) *BoundingBox {
	if !this.initialized {
		copy(this.Min[:], point[:])
		copy(this.Max[:], point[:])
		this.initialized = true

		return this
	}

	for i, val := range point[:] {
		if val > this.Max[i] {
			this.Max[i] = val
		}
		if val < this.Min[i] {
			this.Min[i] = val
		}
	}

	return this

}

// Asynchronously add an array of points to the bounding box
//
// **params**
// + An array of length-n array of numbers
//
// **returns**
// + this BoundingBox for chaining
func (this *BoundingBox) AddRange(points []vec3.T) *BoundingBox {
	if points != nil {
		for _, pt := range points {
			this.Add(&pt)
		}
	}

	return this
}

// Determines if point is contained in the bounding box
//
// **params**
// + the point
// + the tolerance
//
// **returns**
// + true if the two intervals overlap, otherwise false
func (this *BoundingBox) Contains(point *vec3.T, tol float64) bool {
	if !this.initialized {
		return false
	}

	return this.Intersects(new(BoundingBox).Add(point), tol)
}

// Determines if two intervals on the real number line intersect
//
// **params**
// + Beginning of first interval
// + End of first interval
// + Beginning of second interval
// + End of second interval
//
// **returns**
// + true if the two intervals overlap, otherwise false
func (this *BoundingBox) intervalsOverlap(a1, a2, b1, b2 float64, tol float64) bool {
	if tol < -0.5 {
		tol = BoundingBoxTolerance
	}

	x1, x2 := math.Min(a1, a2)-tol, math.Max(a1, a2)+tol
	y1, y2 := math.Min(b1, b2)-tol, math.Max(b1, b2)+tol

	return (x1 >= y1 && x1 <= y2) || (x2 >= y1 && x2 <= y2) || (y1 >= x1 && y1 <= x2) || (y2 >= x1 && y2 <= x2)
}

// Determines if this bounding box intersects with another
//
// **params**
// + BoundingBox to check for intersection with this one
//
// **returns**
// +  true if the two bounding boxes intersect, otherwise false

func (this *BoundingBox) Intersects(bb *BoundingBox, tol float64) bool {
	if !this.initialized || !bb.initialized {
		return false
	}

	a1, a2, b1, b2 := this.Min, this.Max, bb.Min, bb.Max

	for i := 0; i < len(this.Min); i++ {
		if !this.intervalsOverlap(a1[i], a2[i], b1[i], b2[i], tol) {
			return false
		}
	}

	return true
}

// Clear the bounding box, leaving it in an uninitialized state.  Call add, addRange in order to
// initialize
//
// **returns**
// + this BoundingBox for chaining

func (this *BoundingBox) Clear() *BoundingBox {
	this.initialized = false
	return this
}

// Get longest axis of bounding box
//
// **returns**
// + Index of longest axis
func (this *BoundingBox) LongestAxis() int {
	id, max := 0, 0.0

	for i := range this.Min {
		l := this.AxisLength(i)
		if l > max {
			max = l
			id = i
		}
	}

	return id
}

// Get length of given axis.
//
// **params**
// + Index of axis to inspect (between 0 and 2)
//
// **returns**
// + Length of the given axis.  If axis is out of bounds, returns 0.
func (this *BoundingBox) AxisLength(i int) float64 {
	if i < 0 || i > len(this.Min)-1 {
		return 0
	}
	return math.Abs(this.Min[i] - this.Max[i])
}

// Compute the boolean intersection of this with another axis-aligned bounding box.  If the two
// bounding boxes do not intersect, returns null.
//
// **params**
// + BoundingBox to intersect with
//
// **returns**
// + The bounding box formed by the intersection or null if there is no intersection.
func (this *BoundingBox) Intersect(bb *BoundingBox, tol float64) *BoundingBox {
	if !this.initialized {
		return nil
	}
	a1, a2, b1, b2 := this.Min, this.Max, bb.Min, bb.Max

	if !this.Intersects(bb, tol) {
		return nil
	}

	maxbb, minbb := new(vec3.T), new(vec3.T)

	for i := range this.Min {
		maxbb[i] = math.Min(a2[i], b2[i])
		minbb[i] = math.Max(a1[i], b1[i])
	}

	return new(BoundingBox).Add(minbb).Add(maxbb)
}
