package verb

import (
	. "github.com/alexozer/verb/internal"
	"github.com/ungerik/go3d/float64/vec3"
)

func distToSegment(a, b, c *vec3.T) float64 {
	// check if ac is zero length
	acv := vec3.Sub(c, a)
	acl := acv.Length()

	// subtract b from a
	var bma = vec3.Sub(b, a)

	if acl < Tolerance {
		return bma.Length()
	}

	// normalize acv
	acv.Normalize()

	// project b - a to acv = p
	p := vec3.Dot(&bma, &acv)

	// multiply ac by d = acd
	acv.Scale(p).Add(a)

	return vec3.Distance(&acv, b)
}

// Determine if three points form a straight line within a given tolerance for their 2 * squared area
//
//          * p2
//         / \
//        /   \
//       /     \
//      /       \
//     * p1 ---- * p3
//
// The area metric is 2 * the squared norm of the cross product of two edges, requiring no square roots and no divisions
//
// **params**
// + p1
// + p2
// + p3
// + The tolerance
//
// **returns**
// + Whether the triangle passes the test
//
func threePointsAreCollinear(p1, p2, p3 *vec3.T, tol float64) bool {
	// find the area of the triangle without using a square root
	p2mp1 := vec3.Sub(p2, p1)
	p3mp1 := vec3.Sub(p3, p1)
	norm := vec3.Cross(&p2mp1, &p3mp1)
	area := vec3.Dot(&norm, &norm)

	return area < tol
}

// Find the closest point on a segment
//
// **params**
// + point to project
// + first point of segment
// + second point of segment
// + first param of segment
// + second param of segment
//
// **returns**
// + *Object* with u and pt properties
func segmentClosestPoint(pt, segpt0, segpt1 *vec3.T, u0, u1 float64) CurvePoint {
	dif := vec3.Sub(segpt1, segpt0)
	l := dif.Length()

	if l < Epsilon {
		return CurvePoint{u0, *segpt0}
	}

	o := segpt0
	r := dif.Normalize()
	o2pt := vec3.Sub(pt, o)
	do2ptr := vec3.Dot(&o2pt, r)

	if do2ptr < 0 {
		return CurvePoint{u0, *segpt0}
	} else if do2ptr > l {
		return CurvePoint{u1, *segpt1}
	}

	return CurvePoint{
		u0 + (u1-u0)*do2ptr/l,
		vec3.Add(o, r.Scale(do2ptr)),
	}
}
