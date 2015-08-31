package intersect

import "errors"

var errNotImplemented = errors.New("Not implemented")

func Meshes(...interface{}) *Mesh {
	panic(errNotImplemented)
	return nil
}

func Rays(...interface{}) *CurveCurveIntersection {
	panic(errNotImplemented)
	return nil
}
