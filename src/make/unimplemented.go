package make

import (
	"errors"

	"github.com/ungerik/go3d/float64/vec3"
)

var errNotImplemented = errors.New("Not implemented")

func rayClosestPoint(...interface{}) vec3.T {
	panic(errNotImplemented)
	return vec3.Zero
}
