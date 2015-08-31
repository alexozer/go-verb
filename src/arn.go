package verb

import (
	"math"

	. "github.com/alexozer/verb/internal"
	"github.com/ungerik/go3d/float64/vec3"
)

type adaptiveRefinementOptions struct {
	NormTol            float64
	MinDepth           int
	MaxDepth           int
	Refine             bool
	MinDivsU, MinDivsV int
}

var defaultAdaptiveRefinementOptions = adaptiveRefinementOptions{
	NormTol:  8.5e-2,
	MinDepth: 0,
	MaxDepth: 10,
	Refine:   true,
	MinDivsU: 1,
	MinDivsV: 1,
}

type adaptiveRefinementNode struct {
	srf                   *NurbsSurface
	Neighbors             [4]*adaptiveRefinementNode
	children              [2]*adaptiveRefinementNode
	corners, midPoints    [4]*SurfacePoint
	splitVert, splitHoriz bool
	horizontal            bool
	uv05                  UV
	tolerance             float64
}

// If nil passed for corners, they are constructed from the surface
func newAdaptiveRefinementNode(srf *NurbsSurface, corners *[4]*SurfacePoint, neighbors *[4]*adaptiveRefinementNode) *adaptiveRefinementNode {

	//
	// Structure of the child nodes
	// in the adaptive refinement tree
	//
	//  v
	//  ^
	//  |
	//  +--> u
	//
	//                        neighbors[2]
	//
	//                (u0,v1)---(u05,v1)---(u1,v1)
	//                  |           |          |
	//                  |     3     |     2    |
	//                  |           |          |
	// neighbors[3]   (u0,v05)--(u05,v05)--(u1,v05)   neighbors[1]
	//                  |           |          |
	//                  |     0     |     1    |
	//                  |           |          |
	//                (u0,v0)---(u05,v0)---(u1,v0)
	//
	//                        neighbors[0]
	//

	this := new(adaptiveRefinementNode)

	this.srf = srf

	if neighbors == nil {
		this.Neighbors = [4]*adaptiveRefinementNode{}
	} else {
		this.Neighbors = *neighbors
	}

	// if no corners, we need to construct initial corners from the surface
	if corners == nil {
		u0 := srf.knotsU[0]
		u1 := srf.knotsU[len(srf.knotsU)-1]
		v0 := srf.knotsV[0]
		v1 := srf.knotsV[len(srf.knotsV)-1]

		this.corners = [4]*SurfacePoint{
			NewSurfacePoint(UV{u0, v0}),
			NewSurfacePoint(UV{u1, v0}),
			NewSurfacePoint(UV{u1, v1}),
			NewSurfacePoint(UV{u0, v1}),
		}
	} else {
		this.corners = *corners
	}

	return this
}

func (this *adaptiveRefinementNode) IsLeaf() bool {
	return this.children[0] == nil && this.children[1] == nil
}

func (this *adaptiveRefinementNode) Center() *SurfacePoint {
	return this.EvalSrf(this.uv05)
}

func (this *adaptiveRefinementNode) EvalCorners() {
	// eval the center
	this.uv05 = UV{(this.corners[0].UV[0] + this.corners[2].UV[0]) / 2,
		(this.corners[0].UV[1] + this.corners[2].UV[1]) / 2}

	// eval all of the corners
	for i, c := range this.corners {
		// if it's not already evaluated
		if c.Point == nil {
			// evaluate it
			this.corners[i] = this.EvalSrf(c.UV)
		}
	}
}

func (this *adaptiveRefinementNode) EvalSrf(uv UV) *SurfacePoint {
	derivs := this.srf.Derivatives(uv, 1)
	pt := derivs[0][0]
	norm := vec3.Cross(&derivs[0][1], &derivs[1][0])

	degen := true
	for _, compon := range norm {
		if math.Abs(compon) > this.tolerance {
			degen = false
			break
		}
	}

	if !degen {
		norm.Normalize()
	}

	return &SurfacePoint{uv, &pt, &norm, -1, degen}
}

func (this *adaptiveRefinementNode) EdgeCorners(edgeIndex int) []*SurfacePoint {

	// if its a leaf, there are no children to obtain uvs from
	if this.IsLeaf() {
		return []*SurfacePoint{this.corners[edgeIndex]}
	}

	if this.horizontal {
		switch edgeIndex {
		case 0:
			return this.children[0].EdgeCorners(0)
		case 1:
			return append(this.children[0].EdgeCorners(1), this.children[1].EdgeCorners(1)...)
		case 2:
			return this.children[1].EdgeCorners(2)
		case 3:
			return append(this.children[1].EdgeCorners(3), this.children[0].EdgeCorners(3)...)
		}

	}

	// vertical case
	switch edgeIndex {
	case 0:
		return append(this.children[0].EdgeCorners(0), this.children[1].EdgeCorners(0)...)
	case 1:
		return this.children[1].EdgeCorners(1)
	case 2:
		return append(this.children[1].EdgeCorners(2), this.children[0].EdgeCorners(2)...)
	case 3:
		return this.children[0].EdgeCorners(3)
	}

	return nil
}

func (this *adaptiveRefinementNode) GetAllCorners(edgeIndex int) []*SurfacePoint {
	baseArr := []*SurfacePoint{this.corners[edgeIndex]}

	if this.Neighbors[edgeIndex] == nil {
		return baseArr
	}

	// get opposite edges uvs
	corners := this.Neighbors[edgeIndex].EdgeCorners((edgeIndex + 2) % 4)

	funcIndex := edgeIndex % 2

	e := Epsilon

	// range clipping functions
	rangeFuncMap := []func(c *SurfacePoint) bool{
		func(c *SurfacePoint) bool {
			return c.UV[0] > this.corners[0].UV[0]+e && c.UV[0] < this.corners[2].UV[0]-e
		},
		func(c *SurfacePoint) bool {
			return c.UV[1] > this.corners[0].UV[1]+e && c.UV[1] < this.corners[2].UV[1]-e
		},
	}

	// clip the range of uvs to match this one
	cornerCopy := make([]*SurfacePoint, 0, len(corners))
	filterFunc := rangeFuncMap[funcIndex]
	for i := len(corners) - 1; i >= 0; i-- {
		if filterFunc(corners[i]) {
			cornerCopy = append(cornerCopy, corners[i])
		}
	}

	return append(baseArr, cornerCopy...)
}

//public function midpoint( index ){
func (this *adaptiveRefinementNode) midpoint(index int) *SurfacePoint {
	if this.midPoints[index] != nil {
		return this.midPoints[index]
	}

	switch index {
	case 0:
		this.midPoints[0] = this.EvalSrf(UV{this.uv05[0], this.corners[0].UV[1]})
	case 1:
		this.midPoints[1] = this.EvalSrf(UV{this.corners[1].UV[0], this.uv05[1]})
	case 2:
		this.midPoints[2] = this.EvalSrf(UV{this.uv05[0], this.corners[2].UV[1]})
	case 3:
		this.midPoints[3] = this.EvalSrf(UV{this.corners[0].UV[0], this.uv05[1]})
	}

	return this.midPoints[index]

}

func (this *adaptiveRefinementNode) HasBadNormals() bool {
	return this.corners[0].Degen || this.corners[1].Degen || this.corners[2].Degen || this.corners[3].Degen
}

func (this *adaptiveRefinementNode) FixNormals() {
	l := len(this.corners)
	for i := range this.corners {
		if this.corners[i].Degen {

			// get neighbors
			v1 := this.corners[(i+1)%l]
			v2 := this.corners[(i+3)%l]

			// correct the normal
			if v1.Degen {
				this.corners[i].Normal = v2.Normal
			} else {
				this.corners[i].Normal = v1.Normal
			}
		}
	}
}

func (this *adaptiveRefinementNode) ShouldDivide(options *adaptiveRefinementOptions, currentDepth int) bool {
	if currentDepth < options.MinDepth {
		return true
	}
	if currentDepth >= options.MaxDepth {
		return false
	}

	if this.HasBadNormals() {
		this.FixNormals()
		// don't divide any further when encountering a degenerate normal
		return false
	}

	normDiff01 := vec3.Sub(this.corners[0].Normal, this.corners[1].Normal)
	normDiff23 := vec3.Sub(this.corners[2].Normal, this.corners[3].Normal)
	this.splitVert = normDiff01.LengthSqr() > options.NormTol || normDiff23.LengthSqr() > options.NormTol

	normDiff12 := vec3.Sub(this.corners[1].Normal, this.corners[2].Normal)
	normDiff30 := vec3.Sub(this.corners[3].Normal, this.corners[0].Normal)
	this.splitHoriz = normDiff12.LengthSqr() > options.NormTol || normDiff30.LengthSqr() > options.NormTol

	if this.splitVert || this.splitHoriz {
		return true
	}

	center := this.Center()

	for _, corner := range this.corners {
		diffVec := vec3.Sub(center.Normal, corner.Normal)
		if diffVec.LengthSqr() > options.NormTol {
			return true
		}
	}
	return false
}

func (this *adaptiveRefinementNode) Divide(options *adaptiveRefinementOptions) {
	if options == nil {
		options = &defaultAdaptiveRefinementOptions
	}

	this.divide(options, 0, true)
}

func (this *adaptiveRefinementNode) divide(options *adaptiveRefinementOptions, currentDepth int, horiz bool) {
	this.EvalCorners()

	if !this.ShouldDivide(options, currentDepth) {
		return
	}

	currentDepth++

	// is the quad flat in one dir and curved in the other?
	if this.splitVert && !this.splitHoriz {
		horiz = false
	} else if !this.splitVert && this.splitHoriz {
		horiz = true
	}

	this.horizontal = horiz

	if this.horizontal {
		var bott = [4]*SurfacePoint{this.corners[0], this.corners[1], this.midPoints[1], this.midPoints[3]}
		var top = [4]*SurfacePoint{this.midPoints[3], this.midPoints[1], this.corners[2], this.corners[3]}

		this.children = [2]*adaptiveRefinementNode{
			newAdaptiveRefinementNode(
				this.srf, &bott,

				// assign neighbors to bottom node
				&[4]*adaptiveRefinementNode{
					this.Neighbors[0], this.Neighbors[1],
					this.children[1], this.Neighbors[3],
				},
			),

			newAdaptiveRefinementNode(
				this.srf, &top,

				// assign neighbors to top node
				&[4]*adaptiveRefinementNode{
					this.children[0], this.Neighbors[1],
					this.Neighbors[2], this.Neighbors[3],
				},
			),
		}
	} else {
		left := [4]*SurfacePoint{
			this.corners[0], this.midPoints[0],
			this.midPoints[2], this.corners[3],
		}
		right := [4]*SurfacePoint{
			this.midPoints[0], this.corners[1],
			this.corners[2], this.midPoints[2],
		}

		this.children = [2]*adaptiveRefinementNode{
			newAdaptiveRefinementNode(
				this.srf, &left,

				&[4]*adaptiveRefinementNode{
					this.Neighbors[0], this.children[1],
					this.Neighbors[2], this.Neighbors[3],
				},
			),

			newAdaptiveRefinementNode(
				this.srf, &right,

				&[4]*adaptiveRefinementNode{
					this.Neighbors[0], this.Neighbors[1],
					this.Neighbors[2], this.children[0],
				},
			),
		}
	}

	// divide all children recursively
	for _, child := range this.children {
		child.divide(options, currentDepth, !horiz)
	}

}

func (this *adaptiveRefinementNode) Triangulate(mesh *Mesh) *Mesh {
	if mesh == nil {
		mesh = newMesh()
	}

	if this.IsLeaf() {
		return this.TriangulateLeaf(mesh)
	}

	// recurse on the children
	for _, x := range this.children {
		if x == nil {
			break
		}
		x.Triangulate(mesh)
	}

	return mesh
}

func (this *adaptiveRefinementNode) TriangulateLeaf(mesh *Mesh) *Mesh {
	baseIndex := len(mesh.Points)
	uvs := make([]*SurfacePoint, 0)
	ids := make([]int, 0)
	splitid := 0

	// enumerate all uvs in counter clockwise direction
	for i := 0; i < 4; i++ {

		edgeCorners := this.GetAllCorners(i)

		// this is the vertex that is split
		if len(edgeCorners) == 2 {
			splitid = i + 1
		}

		for _, corn := range edgeCorners {
			uvs = append(uvs, corn)
		}
	}

	for _, corner := range uvs {

		// if the id is defined, we can just push it and continue
		if corner.Id != -1 {
			ids = append(ids, corner.Id)
			continue
		}

		mesh.UVs = append(mesh.UVs, corner.UV)
		mesh.Points = append(mesh.Points, *corner.Point)
		mesh.Normals = append(mesh.Normals, *corner.Normal)

		corner.Id = baseIndex
		ids = append(ids, baseIndex)

		baseIndex++
	}

	if len(uvs) == 4 {

		// if the number of points is 4, we're just doing a
		// rectangle - just build the basic triangulated square
		mesh.Faces = append(mesh.Faces, Tri{ids[0], ids[3], ids[1]})
		mesh.Faces = append(mesh.Faces, Tri{ids[3], ids[2], ids[1]})

		// all done
		return mesh
	}

	if len(uvs) == 5 {

		// use the splitcorner to triangulate
		il := len(ids)

		// there will be 3 triangles
		mesh.Faces = append(mesh.Faces, Tri{
			ids[splitid],
			ids[(splitid+1)%il],
			ids[(splitid+2)%il],
		}, Tri{
			ids[(splitid+4)%il],
			ids[(splitid+3)%il],
			ids[splitid],
		}, Tri{
			ids[splitid],
			ids[(splitid+2)%il],
			ids[(splitid+3)%il],
		})

		return mesh
	}

	// make point at center of face
	center := this.Center()

	mesh.UVs = append(mesh.UVs, center.UV)
	mesh.Points = append(mesh.Points, *center.Point)
	mesh.Normals = append(mesh.Normals, *center.Normal)

	// get index
	centerIndex := len(mesh.Points) - 1

	// build triangle fan from center
	j := len(uvs) - 1

	for i := 0; i < len(uvs); i++ {
		mesh.Faces = append(mesh.Faces, Tri{centerIndex, ids[j], ids[i]})
		j = i
	}

	return mesh
}

type adaptiveRefinementNodeTree []*adaptiveRefinementNode

func (this *adaptiveRefinementNodeTree) Triangulate() *Mesh {
	// triangulate all of the nodes of the tree
	mesh := newMesh()
	for _, x := range *this {
		x.Triangulate(mesh)
	}

	return mesh
}
