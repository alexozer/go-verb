package intersect

import (
	"github.com/alexozer/verb"
	"github.com/ungerik/go3d/vec3"
)

type (
	CurveCurveIntersection struct {
		Point0, Point1 *vec3.T
		U0, U1         float64
	}

	CurveCurveIntersectionOptions struct {
		SampleTol, Tol float64
	}

	CurveSurfaceIntersection struct {
		U            float64
		UV           verb.UV
		CurvePoint   *vec3.T
		SurfacePoint *vec3.T
	}

	PolylineMeshIntersection struct {
		Point                    *vec3.T
		U                        float64
		UV                       verb.UV
		PolylineIndex, FaceIndex int
	}

	SurfaceSurfaceIntersectionPoint struct {
		UV0, UV1 verb.UV
		Point    *vec3.T
		Dist     float64
	}

	TriSegmentIntersection struct {
		Point *vec3.T // where the intersection took place
		S     float64 // the u param where u is the axis from v0 to v1
		T     float64 // the v param where v is the axis from v0 to v2
		P     float64 // the parameter along the segment
	}
)
