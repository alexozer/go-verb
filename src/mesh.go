package verb

import "github.com/ungerik/go3d/float64/vec3"

type Tri [3]int

type Mesh struct {
	Faces   []Tri
	Points  []vec3.T
	Normals []vec3.T
	UVs     []UV
}

func newMesh() *Mesh {
	return &Mesh{
		Faces:   make([]Tri, 0),
		Points:  make([]vec3.T, 0),
		Normals: make([]vec3.T, 0),
		UVs:     make([]UV, 0),
	}
}
