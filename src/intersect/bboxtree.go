package intersect

import (
	"math"
	"math/rand"

	"github.com/alexozer/verb"
	. "github.com/alexozer/verb/internal"
)

type BoundingBoxNode interface {
	BoundingBox() BoundingBox
}

type BoundingBoxInnerNode struct {
	bbox     BoundingBox
	children [2]BoundingBoxNode
}

func NewBoundingBoxInnerNode(bbox BoundingBox, children [2]BoundingBoxNode) *BoundingBoxInnerNode {
	return &BoundingBoxInnerNode{bbox, children}
}

func (this *BoundingBoxInnerNode) BoundingBox() BoundingBox {
	return this.bbox
}

type BoundingBoxLeaf struct {
	bbox      BoundingBox
	faceIndex int
}

func NewBoundingBoxLeaf(bbox BoundingBox, faceIndex int) *BoundingBoxLeaf {
	return &BoundingBoxLeaf{bbox, faceIndex}
}

func (this *BoundingBoxLeaf) BoundingBox() BoundingBox {
	return this.bbox
}

type BoundingBoxTree interface {
	BoundingBoxNode
	Empty() bool
	Indivisible(tolerance float64) bool
	Split() (BoundingBoxTree, BoundingBoxTree)
	Yield() interface{}
}

type lazyCurveCoundingBoxTree struct {
	curve   *verb.NurbsCurve
	knots   KnotVec
	knotTol float64
}

func newLazyCurveBoundingBoxTree(curve *verb.NurbsCurve, knotTol *float64) *lazyCurveCoundingBoxTree {
	this := &lazyCurveCoundingBoxTree{
		curve: curve,
		knots: KnotVec(curve.Knots()),
	}

	if knotTol == nil {
		this.knotTol = this.knots.Domain() / 64
	} else {
		this.knotTol = *knotTol
	}

	return this
}

func (this *lazyCurveCoundingBoxTree) Split() (BoundingBoxTree, BoundingBoxTree) {
	min := this.knots[0]
	max := this.knots[len(this.knots)-1]
	dom := max - min

	crv0, crv1 := this.curve.Split((max+min)/2.0 + dom*0.1*rand.Float64())

	return newLazyCurveBoundingBoxTree(crv0, &this.knotTol),
		newLazyCurveBoundingBoxTree(crv1, &this.knotTol)
}

func (this *lazyCurveCoundingBoxTree) BoundingBox() BoundingBox {
	var bbox BoundingBox
	pts := this.curve.ControlPoints()
	bbox.AddRange(pts)

	return bbox
}

func (this *lazyCurveCoundingBoxTree) Yield() interface{} {
	return this.curve
}

func (this *lazyCurveCoundingBoxTree) Indivisible(tolerance float64) bool {
	return this.knots.Domain() < this.knotTol
}

func (this *lazyCurveCoundingBoxTree) Empty() bool {
	return false
}

type lazyPolylineBoundingBoxTree struct {
	min, max    int
	polyline    Polyline
	boundingBox BoundingBox
}

func newLazyPolylineBoundingBoxTree(polyline Polyline, min, max int) *lazyPolylineBoundingBoxTree {
	return &lazyPolylineBoundingBoxTree{
		min:      min,
		max:      max,
		polyline: polyline,
	}
}

func (this *lazyPolylineBoundingBoxTree) Split() (BoundingBoxTree, BoundingBoxTree) {
	min := this.min
	max := this.max

	pivot := min + int(math.Ceil((float64(max)-float64(min))/2.0))

	return newLazyPolylineBoundingBoxTree(this.polyline, min, pivot),
		newLazyPolylineBoundingBoxTree(this.polyline, pivot, max)
}

func (this *lazyPolylineBoundingBoxTree) BoundingBox() BoundingBox {
	var bbox BoundingBox
	bbox.AddRange(this.polyline.Points)
	return bbox
}

func (this *lazyPolylineBoundingBoxTree) Yield() interface{} {
	return this.min
}

func (this *lazyPolylineBoundingBoxTree) Indivisible(tolerance float64) bool {
	return this.max-this.min == 1
}

func (this *lazyPolylineBoundingBoxTree) Empty() bool {
	return this.max-this.min == 0
}

type lazySurfaceBoundingBoxTree struct {
	surface            *verb.NurbsSurface
	knotsU, knotsV     KnotVec
	boundingBox        *BoundingBox
	splitV             bool
	knotTolU, knotTolV float64
}

func newLazySurfaceBoundingBoxTree(surface *verb.NurbsSurface, splitV bool, knotTolU, knotTolV *float64) *lazySurfaceBoundingBoxTree {
	this := lazySurfaceBoundingBoxTree{
		surface: surface,
		knotsU:  surface.KnotsU(),
		knotsV:  surface.KnotsV(),
		splitV:  splitV, // default: false
	}

	if knotTolU == nil {
		this.knotTolU = this.knotsU.Domain() / 16
	} else {
		this.knotTolU = *knotTolU
	}
	if knotTolV == nil {
		this.knotTolV = this.knotsV.Domain() / 16
	} else {
		this.knotTolV = *knotTolV
	}

	return &this
}

func (this *lazySurfaceBoundingBoxTree) Split() (BoundingBoxTree, BoundingBoxTree) {
	var min, max float64

	if this.splitV {
		min = this.knotsV[0]
		max = this.knotsV[len(this.knotsV)-1]
	} else {
		min = this.knotsU[0]
		max = this.knotsU[len(this.knotsU)-1]
	}

	//dom := max - min
	pivot := (min + max) / 2.0 //+ dom*0.01*rand.Float64()

	srf0, srf1 := this.surface.Split(pivot, this.splitV)

	return newLazySurfaceBoundingBoxTree(srf0, !this.splitV, &this.knotTolU, &this.knotTolV),
		newLazySurfaceBoundingBoxTree(srf1, !this.splitV, &this.knotTolU, &this.knotTolV)
}

func (this *lazySurfaceBoundingBoxTree) BoundingBox() BoundingBox {
	var bbox BoundingBox

	for _, row := range this.surface.ControlPoints() {
		bbox.AddRange(row)
	}

	return bbox
}

func (this *lazySurfaceBoundingBoxTree) Yield() interface{} {
	return this.surface
}

func (this *lazySurfaceBoundingBoxTree) Indivisible(tolerance float64) bool {
	if this.splitV {
		return this.knotsV.Domain() < this.knotTolV
	} else {
		return this.knotsU.Domain() < this.knotTolU
	}

	return this.knotsV.Domain() < this.knotTolV &&
		this.knotsU.Domain() < this.knotTolU
}

func (this *lazySurfaceBoundingBoxTree) Empty() bool {
	return false
}

type surfaceBoundingBoxTree struct {
	children       [2]BoundingBoxTree
	surface        *verb.NurbsSurface
	knotsU, knotsV KnotVec
	boundingBox    *BoundingBox
}

func newSurfaceBoundingBoxTree(surface *verb.NurbsSurface, splitV bool, knotTolU, knotTolV *float64) *surfaceBoundingBoxTree {
	this := new(surfaceBoundingBoxTree)
	this.surface = surface

	knotsU, knotsV := KnotVec(surface.KnotsU()), KnotVec(surface.KnotsV())

	var _knotTolU float64
	if knotTolU == nil {
		_knotTolU = knotsU.Domain() / 16.0
	} else {
		_knotTolU = *knotTolU
	}

	var _knotTolV float64
	if knotTolV == nil {
		_knotTolV = knotsV.Domain() / 16.0
	} else {
		_knotTolV = *knotTolV
	}

	var divisible bool
	if splitV {
		divisible = knotsV.Domain() > _knotTolV
	} else {
		divisible = knotsU.Domain() > _knotTolU
	}

	if !divisible {
		return this
	}

	var min, max float64
	if splitV {
		min = knotsV[0]
		max = knotsV[len(knotsV)-1]
	} else {
		min = knotsU[0]
		max = knotsU[len(knotsU)-1]
	}

	dom := max - min
	pivot := (min+max)/2.0 + dom*0.1*rand.Float64()

	srf0, srf1 := this.surface.Split(pivot, splitV)

	this.children = [2]BoundingBoxTree{
		newSurfaceBoundingBoxTree(srf0, !splitV, &_knotTolU, &_knotTolV),
		newSurfaceBoundingBoxTree(srf1, !splitV, &_knotTolU, &_knotTolV),
	}

	return this
}

func (this *surfaceBoundingBoxTree) Split() (BoundingBoxTree, BoundingBoxTree) {
	return this.children[0], this.children[1]
}

func (this *surfaceBoundingBoxTree) BoundingBox() BoundingBox {
	var bbox BoundingBox

	for _, row := range this.surface.ControlPoints() {
		bbox.AddRange(row)
	}

	return bbox
}

func (this *surfaceBoundingBoxTree) Yield() interface{} {
	return this.surface
}

func (this *surfaceBoundingBoxTree) Indivisible(tolerance float64) bool {
	return this.children[0] == nil && this.children[1] == nil
}

func (this *surfaceBoundingBoxTree) Empty() bool {
	return false
}

type lazyMeshBoundingBoxTree struct {
	mesh            *Mesh
	faceIndices     []int
	boundingBox     BoundingBox
	bboxInitialized bool
}

func newLazyMeshBoundingBoxTree(mesh *Mesh, faceIndices []int) *lazyMeshBoundingBoxTree {
	this := new(lazyMeshBoundingBoxTree)

	this.mesh = mesh

	if faceIndices == nil {
		faceIndices = make([]int, len(mesh.Faces))
		for i := range faceIndices {
			faceIndices[i] = i
		}
	}
	this.faceIndices = faceIndices

	return this
}

func (this *lazyMeshBoundingBoxTree) BoundingBox() BoundingBox {
	if !this.bboxInitialized {
		this.boundingBox = this.mesh.BoundingBox(this.faceIndices)
		this.bboxInitialized = true
	}

	return this.boundingBox
}

func (this *lazyMeshBoundingBoxTree) Split() (BoundingBoxTree, BoundingBoxTree) {
	as := this.mesh.SortedTrianglesOnLongestAxis(this.BoundingBox(), this.faceIndices)

	halfLen := len(as) / 2
	l, r := as[:halfLen], as[halfLen:]

	return newLazyMeshBoundingBoxTree(this.mesh, l),
		newLazyMeshBoundingBoxTree(this.mesh, r)
}

func (this *lazyMeshBoundingBoxTree) Yield() interface{} {
	return this.faceIndices[0]
}

func (this *lazyMeshBoundingBoxTree) Indivisible(tolerance float64) bool {
	return len(this.faceIndices) == 1
}

func (this *lazyMeshBoundingBoxTree) Empty() bool {
	return len(this.faceIndices) == 0
}
