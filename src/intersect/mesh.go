package intersect

import (
	"math"
	"sort"

	"github.com/alexozer/verb"
	"github.com/ungerik/go3d/float64/vec3"
)

type Mesh verb.Mesh

// Form axis-aligned bounding box from triangles of mesh
//
// **params**
// + a mesh
// + face indices of the mesh to include in the bounding box
//
// **returns**
// + a BoundingBox containing the mesh
//
func (this *Mesh) BoundingBox(faceIndices []int) BoundingBox {
	bb := BoundingBox{}

	for _, iFace := range faceIndices {
		for _, iPt := range this.Faces[iFace] {
			bb.Add(&this.Points[iPt])
		}
	}

	return bb
}

//
// Make tree of axis aligned bounding boxes
//
// **params**
// + array of length 3 arrays of numbers representing the points
// + array of length 3 arrays of number representing the triangles
// + array of numbers representing the relevant triangles to use to form aabb
//
// **returns**
// + a point represented by an array of length (dim)
//
func (this *Mesh) BoundingBoxTree(faceIndices []int) BoundingBoxNode {
	bbox := this.BoundingBox(faceIndices)

	if len(faceIndices) == 1 {
		return NewBoundingBoxLeaf(bbox, faceIndices[0])
	}

	sortedIndices := this.SortedTrianglesOnLongestAxis(bbox, faceIndices)

	halfLen := len(sortedIndices) / 2
	leftIndices := sortedIndices[:halfLen]
	rightIndices := sortedIndices[halfLen:]

	return NewBoundingBoxInnerNode(bbox, [2]BoundingBoxNode{
		this.BoundingBoxTree(leftIndices),
		this.BoundingBoxTree(rightIndices),
	})
}

type faceCoord struct {
	FaceIndex int
	Coord     float64
}

type minFaceCoordMap []faceCoord

// Implements sort.Interface
func (this minFaceCoordMap) Len() int {
	return len(this)
}

func (this minFaceCoordMap) Swap(i, j int) {
	this[i], this[j] = this[j], this[i]
}

func (this minFaceCoordMap) Less(i, j int) bool {
	return this[i].Coord < this[j].Coord
}

//
// Sort particular faces of a mesh on the longest axis
//
// **params**
// + bounding box containing the faces
// + the mesh it self
// + the indices of the mesh faces to inspect
//
// **returns**
// + a point represented by an array of length (dim)
//
func (this *Mesh) SortedTrianglesOnLongestAxis(bbox BoundingBox, faceIndices []int) (sortedFaceIndices []int) {
	longAxis := bbox.LongestAxis()

	minCoords := make(minFaceCoordMap, len(faceIndices))
	for i, faceIndex := range faceIndices {
		triMin := minCoordOnAxis(this.Points, &this.Faces[faceIndex], longAxis)
		minCoords[i] = faceCoord{faceIndex, triMin}
	}

	sort.Sort(minCoords)

	sortedFaceIndices = make([]int, len(minCoords))
	for i, faceCoord := range minCoords {
		sortedFaceIndices[i] = faceCoord.FaceIndex
	}

	return
}

type MeshIntersectionPoint struct {
	UV0, UV1               verb.UV
	Point                  vec3.T
	FaceIndex0, FaceIndex1 int

	// Opp, Adj, and Visited are fine with default zero-values
	// (tags to navigate a segment structure)
	Opp, Adj *MeshIntersectionPoint
	Visited  bool
}

// TODO
// Intersect two meshes, yielding a list of polylines
//
// **params**
// + MeshData for the first mesh
// + MeshData for the latter
//
// **returns**
// + array of array of MeshIntersectionPoints
func (this *Mesh) IntersectSlices(min, max, step float64) [][][]*MeshIntersectionPoint {
	bbtree := newLazyMeshBoundingBoxTree(this, nil)
	bb := bbtree.BoundingBox()

	x0, y0 := bb.Min[0], bb.Min[1]
	x1, y1 := bb.Max[0], bb.Max[1]

	numSlices := int(math.Ceil(math.Abs(max-min) / step))
	slices := make([]*Mesh, numSlices)

	var i int
	for z := min; z <= max; z += step {
		pts := []vec3.T{
			vec3.T{x0, y0, z},
			vec3.T{x1, y0, z},
			vec3.T{x1, y1, z},
			vec3.T{x0, y1, z},
		}

		uvs := []verb.UV{{0, 0}, {1, 0}, {0, 1}, {1, 1}}

		faces := []verb.Tri{{0, 1, 2}, {0, 2, 3}}
		plane := Mesh{faces, pts, nil, uvs}

		slices[i] = Meshes(this, plane, bbtree)
		i++
	}

	return nil
}

func (this *Mesh) TriangleUvFromPoint(faceIndex int, f *vec3.T) verb.UV {
	tri := this.Faces[faceIndex]

	p0 := this.Points[tri[0]]
	p1 := this.Points[tri[1]]
	p2 := this.Points[tri[2]]

	uv0 := this.UVs[tri[0]]
	uv1 := this.UVs[tri[1]]
	uv2 := this.UVs[tri[2]]

	f0 := vec3.Sub(&p0, f)
	f1 := vec3.Sub(&p1, f)
	f2 := vec3.Sub(&p2, f)

	// calculate the areas and factors (order of parameters doesn't matter):
	p1.Sub(&p0)
	p2.Sub(&p0)
	aVec := vec3.Cross(&p1, &p2)
	a := aVec.Length()

	a0Vec := vec3.Cross(&f1, &f2)
	a1Vec := vec3.Cross(&f2, &f0)
	a2Vec := vec3.Cross(&f0, &f1)

	a0 := a0Vec.Length() / a
	a1 := a1Vec.Length() / a
	a2 := a2Vec.Length() / a

	// find the uv corresponding to point f (uv1/uv2/uv3 are associated to p1/p2/p3):
	return verb.UV{
		a0*uv0[0] + a1*uv1[0] + a2*uv2[0],
		a0*uv0[1] + a1*uv1[1] + a2*uv2[1],
	}
}
