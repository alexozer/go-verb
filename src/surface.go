package verb

import (
	"errors"
	"math"

	. "github.com/alexozer/verb/internal"
	"github.com/ungerik/go3d/float64/mat4"
	"github.com/ungerik/go3d/float64/vec3"
)

type UV [2]float64

type SurfacePoint struct {
	UV            UV
	Point, Normal *vec3.T
	Id            int
	Degen         bool // normally false
}

func NewSurfacePoint(uv UV) *SurfacePoint {
	return &SurfacePoint{
		UV: uv,
		Id: -1,
	}
}

type NurbsSurface struct {
	// integer degree of surface in u direction
	degreeU int

	// integer degree of surface in v direction
	degreeV int

	// 2d array of control points, the vertical direction (u) increases from top to bottom, the v direction from left to right,
	// and where each control point is an array of length (dim)
	controlPoints [][]HomoPoint

	// array of nondecreasing knot values in u direction
	knotsU KnotVec

	// array of nondecreasing knot values in v direction
	knotsV KnotVec
}

func NewNurbsSurfaceUnchecked(degreeU, degreeV int, controlPoints [][]vec3.T, weights [][]float64, knotsU, knotsV []float64) *NurbsSurface {
	return &NurbsSurface{
		degreeU, degreeV,
		Homogenize2d(controlPoints, weights),
		KnotVec(knotsU).Clone(), KnotVec(knotsV).Clone(),
	}
}

func NewNurbsSurface(degreeU, degreeV int, controlPoints [][]vec3.T, weights [][]float64, knotsU, knotsV []float64) (*NurbsSurface, error) {
	this := NewNurbsSurfaceUnchecked(degreeU, degreeV, controlPoints, weights, knotsU, knotsV)
	if err := this.check(); err != nil {
		return nil, err
	}

	return this, nil
}

func (this *NurbsSurface) DegreeU() int {
	return this.degreeU
}

func (this *NurbsSurface) DegreeV() int {
	return this.degreeV
}

func (this *NurbsSurface) ControlPoints() [][]vec3.T {
	return Dehomogenize2d(this.controlPoints)
}

func (this *NurbsSurface) Weights() [][]float64 {
	return Weight2d(this.controlPoints)
}

func (this *NurbsSurface) KnotsU() []float64 {
	return []float64(this.knotsU.Clone())
}

func (this *NurbsSurface) KnotsV() []float64 {
	return []float64(this.knotsV.Clone())
}

// uDir normally true
func (this *NurbsSurface) isClosed(uDir bool) bool {
	var cpts [][]HomoPoint

	if uDir {
		cpts = this.controlPoints
	} else {
		cpts := make([][]HomoPoint, len(this.controlPoints))
		for i := range cpts {
			cpts[i] = make([]HomoPoint, len(this.controlPoints[0]))
			copy(cpts[i], this.controlPoints[i])
		}

		cpts = transposed(cpts)
	}

	for i := range cpts[0] {
		// TODO there's probably a more efficient, equally effective way
		first, last := cpts[0][i], cpts[len(cpts)-1][i]
		dist := math.Sqrt(
			vec3.SquareDistance(&first.Vec3, &last.Vec3) +
				(first.W-last.W)*(first.W-last.W),
		)
		if dist >= Epsilon {
			return false
		}
	}

	return true
}

func (this *NurbsSurface) DomainU() (min, max float64) {
	min = this.knotsU[0]
	max = this.knotsU[len(this.knotsU)-1]
	return
}

func (this *NurbsSurface) DomainV() (min, max float64) {
	min = this.knotsV[0]
	max = this.knotsV[len(this.knotsV)-1]
	return
}

func (this *NurbsSurface) ClosestPoint(p vec3.T) vec3.T {
	uv := this.ClosestParam(p)
	return this.Point(uv)
}

func (this *NurbsSurface) ClosestParam(p vec3.T) UV {
	// for surfaces, we try to minimize the following:
	//
	// f = Su(u,v) * r = 0
	// g = Sv(u,v) * r = 0
	//
	//  where r = S(u,v) - P
	//
	// Again, this requires newton iteration, but this time our objective function is vector valued
	//
	//    J d = k
	//
	//      d =   [ u* - u, v* - v ]
	//		k = - [ f(u,v), g(u,v) ]
	//		J =
	//          |Su|^2   +  Suu * r       Su*Sv  +  Suv * r
	//		     Su*Sv   +  Svu * r      |Sv|^2  +  Svv * r
	//
	//
	// 	we have similar halting conditions:
	//
	//  point coincidence
	//
	//		|S(u,v) - p| < e1
	//
	//  cosine
	//
	//   |Su(u,v)*(S(u,v) - P)|
	//   ----------------------  < e2
	//   |Su(u,v)| |S(u,v) - P|
	//
	//   |Sv(u,v)*(S(u,v) - P)|
	//   ----------------------  < e2
	//   |Sv(u,v)| |S(u,v) - P|
	//
	//  1) first check 2 & 3
	// 	2) if at least one of these is not, compute new value, otherwise halt
	// 	3) ensure the parameter stays within range
	// 			* if not closed, don't allow outside of range a-b
	// 			* if closed (e.g. circle), allow to move back to beginning
	//  4)  if |(u* - u)C'(u)| < e1, halt
	//

	maxits := 5
	var i int
	var e [][]vec3.T
	eps1, eps2 := 0.0001, 0.0005
	var dif vec3.T
	minu, maxu := this.knotsU[0], this.knotsU[len(this.knotsU)-1]
	minv, maxv := this.knotsV[0], this.knotsV[len(this.knotsV)-1]
	closedu, closedv := this.isClosed(true), this.isClosed(false)
	var cuv UV

	// TODO divide surface instead of a full on tessellation

	// approximate closest point with tessellation
	tess := this.tessellateAdaptive(&defaultAdaptiveRefinementOptions)

	dmin := math.MaxFloat64

	for i, x := range tess.Points {
		d := vec3.SquareDistance(&p, &x)

		if d < dmin {
			dmin = d
			cuv = tess.UVs[i]
		}
	}

	f := func(uv UV) [][]vec3.T {
		return this.Derivatives(uv, 2)
	}

	n := func(uv UV, e [][]vec3.T, r vec3.T) UV {
		// f = Su(u,v) * r = 0
		// g = Sv(u,v) * r = 0

		Su, Sv := e[1][0], e[0][1]
		Suu, Svv := e[2][0], e[0][2]
		Suv, Svu := e[1][1], e[1][1]

		f := vec3.Dot(&Su, &r)
		g := vec3.Dot(&Sv, &r)

		k := [2]float64{-f, -g}

		J00 := vec3.Dot(&Su, &Su) + vec3.Dot(&Suu, &r)
		J01 := vec3.Dot(&Su, &Sv) + vec3.Dot(&Suv, &r)
		J10 := vec3.Dot(&Su, &Sv) + vec3.Dot(&Svu, &r)
		J11 := vec3.Dot(&Sv, &Sv) + vec3.Dot(&Svv, &r)

		//J := [2][2]float64{{J00, J01}, {J10, J11}}
		//J := Mat2{J00, J01, J10, J11}

		//    d =   [ u* - u, v* - v ]
		//		k = - [ f(u,v), g(u,v) ]
		//		J =
		//          |Su|^2   +  Suu * r       Su*Sv  +  Suv * r
		//		     Su*Sv   +  Svu * r      |Sv|^2  +  Svv * r
		//

		//d := J.Solve(k)
		x, y := Mat2Solve(J00, J01, J10, J11, k[0], k[1])

		//return UV{d[0] + uv[0], d[1] + uv[1]}
		return UV{x + uv[0], y + uv[1]}
	}

	for i < maxits {
		e = f(cuv)

		//  point coincidence
		//
		//		|S(u,v) - p| < e1
		c1v := vec3.Distance(&e[0][0], &p)

		//
		//  cosine
		//
		//   |Su(u,v)*(S(u,v) - P)|
		//   ----------------------  < e2
		//   |Su(u,v)| |S(u,v) - P|
		//
		//   |Sv(u,v)*(S(u,v) - P)|
		//   ----------------------  < e2
		//   |Sv(u,v)| |S(u,v) - P|
		//
		c2an := vec3.Dot(&e[1][0], &dif)
		c2ad := e[1][0].Length() * c1v

		c2bn := vec3.Dot(&e[0][1], &dif)
		c2bd := e[0][1].Length() * c1v

		c2av := c2an / c2ad
		c2bv := c2bn / c2bd

		c1 := c1v < eps1
		c2a := c2av < eps2
		c2b := c2bv < eps2

		// if all of the tolerance are met, we're done
		if c1 && c2a && c2b {
			return cuv
		}

		// otherwise, take a step
		ct := n(cuv, e, dif)

		// correct for exceeding bounds
		if ct[0] < minu {
			if closedu {
				ct = UV{maxu - (ct[0] - minu), ct[1]}
			} else {
				ct = UV{minu + Epsilon, ct[1]}
			}
		} else if ct[0] > maxu {
			if closedu {
				ct = UV{minu + (ct[0] - maxu), ct[1]}
			} else {
				ct = UV{maxu - Epsilon, ct[1]}
			}
		}

		if ct[1] < minv {
			if closedv {
				ct = UV{ct[0], maxv - (ct[1] - minv)}
			} else {
				ct = UV{ct[0], minv + Epsilon}
			}
		} else if ct[1] > maxv {
			if closedv {
				ct = UV{ct[0], minv + (ct[0] - maxv)}
			} else {
				ct = UV{ct[0], maxv - Epsilon}
			}
		}

		// if |(u* - u) C'(u)| < e1, halt
		c3v0v := e[1][0].Scaled(ct[0] - cuv[0])
		c3v0 := c3v0v.Length()
		c3v1v := e[0][1].Scaled(ct[1] - cuv[1])
		c3v1 := c3v1v.Length()

		if c3v0+c3v1 < eps1 {
			return cuv
		}

		cuv = ct
		i++

	}

	return cuv
}

// Validate a NurbsSurfaceData object
//
// **params**
// + The data object
//
// **returns**
// + The original, unmodified data
func (this *NurbsSurface) check() error {
	if this.controlPoints == nil {
		return errors.New("Control points array cannot be nil!")
	}

	if this.degreeU < 1 {
		return errors.New("degreeU must be greater than 1!")
	}
	if this.degreeV < 1 {
		return errors.New("degreeV must be greater than 1!")
	}

	if this.knotsU == nil {
		return errors.New("knotsU cannot be nil!")
	}
	if this.knotsV == nil {
		return errors.New("knotsV cannot be nil!")
	}

	if len(this.knotsU) != len(this.controlPoints)+this.degreeU+1 {
		return errors.New("len(controlPointsU) + degreeU + 1 must equal len(knotsU)!")
	}
	if len(this.knotsV) != len(this.controlPoints[0])+this.degreeV+1 {
		return errors.New("len(controlPointsV) + degreeV + 1 must equal len(knotsV)!")
	}

	if !this.knotsU.IsValid(this.degreeU) || !this.knotsV.IsValid(this.degreeV) {
		return errors.New("Invalid knot vector format! Should begin with degree + 1 repeats and end with degree + 1 repeats!")
	}

	return nil
}

//
// Tessellate a NURBS surface on equal spaced intervals in the parametric domain
//
// **params**
// + NurbsSurfaceData object
// + number of divisions in the u direction
// + number of divisions in the v direction
//
// **returns**
// + MeshData object
//
func (this *NurbsSurface) tessellateNaive(divsU, divsV int) *Mesh {
	if divsU < 1 {
		divsU = 1
	}
	if divsV < 1 {
		divsV = 1
	}

	//degreeU, degreeV := this.DegreeU, this.DegreeV
	//controlPoints := this.ControlPoints
	knotsU, knotsV := this.knotsU, this.knotsV

	uSpan := knotsU[len(knotsU)-1] - knotsU[0]
	vSpan := knotsV[len(knotsV)-1] - knotsV[0]

	spanU := uSpan / float64(divsU)
	spanV := vSpan / float64(divsV)

	numPoints := (divsU + 1) * (divsV + 1)
	points := make([]vec3.T, numPoints)
	uvs := make([]UV, numPoints)
	normals := make([]vec3.T, 0, numPoints)

	var counter int
	for i := 0; i <= divsU; i++ {
		for j := 0; j <= divsV; j++ {
			uv := UV{float64(i) * spanU, float64(j) * spanV}
			uvs[counter] = uv

			derivs := this.Derivatives(uv, 1)
			pt := derivs[0][0]
			points[counter] = pt

			normal := vec3.Cross(&derivs[1][0], &derivs[0][1])
			normals = append(normals, *normal.Normalize())

			counter++
		}
	}

	faces := make([]Tri, 0, 2*divsU*divsV)

	for i := 0; i < divsU; i++ {
		for j := 0; j < divsV; j++ {
			ai := i*(divsV+1) + j
			bi := (i+1)*(divsV+1) + j
			ci := bi + 1
			di := ai + 1
			abc := Tri{ai, bi, ci}
			acd := Tri{ai, ci, di}

			faces = append(faces, abc, acd)
		}
	}

	return &Mesh{faces, points, normals, uvs}
}

//
// Divide a NURBS surface int equal spaced intervals in the parametric domain as AdaptiveRefinementNodes
//
// **params**
// + NurbsSurfaceData object
// + SurfaceDivideOptions object
//
// **returns**
// + MeshData object
//
func (this *NurbsSurface) adaptiveDivisions(options *adaptiveRefinementOptions) []*adaptiveRefinementNode {
	if options == nil {
		options = &defaultAdaptiveRefinementOptions
	}

	minU := (len(this.controlPoints) - 1) * 2
	minV := (len(this.controlPoints[0]) - 1) * 2

	var divsU, divsV int
	if options.MinDivsU > minU {
		divsU = options.MinDivsU
	} else {
		divsU = minU
	}
	if options.MinDivsU > minV {
		divsV = options.MinDivsV
	} else {
		divsV = minV
	}

	// get necessary intervals
	umax := this.knotsU[len(this.knotsU)-1]
	umin := this.knotsU[0]
	vmax := this.knotsV[len(this.knotsV)-1]
	vmin := this.knotsV[0]

	du := (umax - umin) / float64(divsU)
	dv := (vmax - vmin) / float64(divsV)

	pts := make([][]*SurfacePoint, divsV+1)

	// 1) evaluate all of the corners
	for i := range pts {
		ptrow := make([]*SurfacePoint, divsU+1)
		for j := range ptrow {
			uv := UV{umin + du*float64(j), vmin + dv*float64(i)}

			// todo: make this faster by specifying n,m
			ds := this.Derivatives(uv, 1)

			norm := vec3.Cross(&ds[0][1], &ds[1][0])
			norm.Normalize()

			degen := true
			for _, compon := range norm {
				if math.Abs(compon) > Tolerance {
					degen = false
					break
				}
			}
			ptrow[j] = &SurfacePoint{uv, &ds[0][0], &norm, -1, degen}
		}
		pts[i] = ptrow
	}

	divs := make([]*adaptiveRefinementNode, divsU*divsV)

	// 2) make all of the nodes
	var divsI int
	for i := 0; i < divsV; i++ {
		for j := 0; j < divsU; j++ {
			corners := [4]*SurfacePoint{
				pts[divsV-i-1][j],
				pts[divsV-i-1][j+1],
				pts[divsV-i][j+1],
				pts[divsV-i][j],
			}

			divs[divsI] = newAdaptiveRefinementNode(this, &corners, nil)
			divsI++
		}
	}

	if !options.Refine {
		return divs
	}

	// 3) assign all of the neighbors and divide
	for i := 0; i < divsV; i++ {
		for j := 0; j < divsU; j++ {
			ci := i*divsU + j
			n := north(ci, i, j, divsU, divsV, divs)
			e := east(ci, i, j, divsU, divsV, divs)
			s := south(ci, i, j, divsU, divsV, divs)
			w := west(ci, i, j, divsU, divsV, divs)

			divs[ci].Neighbors = [4]*adaptiveRefinementNode{s, e, n, w}
			divs[ci].Divide(options)
		}
	}

	return divs
}

func north(index, i, j, divsU, divsV int, divs []*adaptiveRefinementNode) *adaptiveRefinementNode {
	if i == 0 {
		return nil
	}
	return divs[index-divsU]
}

func south(index, i, j, divsU, divsV int, divs []*adaptiveRefinementNode) *adaptiveRefinementNode {
	if i == divsV-1 {
		return nil
	}
	return divs[index+divsU]
}

func east(index, i, j, divsU, divsV int, divs []*adaptiveRefinementNode) *adaptiveRefinementNode {
	if j == divsU-1 {
		return nil
	}
	return divs[index+1]
}

func west(index, i, j, divsU, divsV int, divs []*adaptiveRefinementNode) *adaptiveRefinementNode {
	if j == 0 {
		return nil
	}
	return divs[index-1]
}

// Tessellate the surface
//
// **params**
// + an AdaptiveRefinementOptions object
//
// **returns**
// + A MeshData object
func (this *NurbsSurface) Tessellate() *Mesh {
	return this.tessellateAdaptive(nil)
}

func (this *NurbsSurface) tessellateAdaptive(options *adaptiveRefinementOptions) *Mesh {
	if options == nil {
		options = &defaultAdaptiveRefinementOptions
	}

	// adaptive divide
	arrTrees := adaptiveRefinementNodeTree(this.adaptiveDivisions(options))

	// triangulation
	return arrTrees.Triangulate()
}

func (this *NurbsSurface) Reverse(useV bool) *NurbsSurface {
	reversed := NurbsSurface{degreeU: this.degreeU, degreeV: this.degreeV}

	if useV {
		reversed.knotsU = append(KnotVec(nil), this.knotsU...)
		reversed.knotsV = this.knotsV.Reversed()

		reversed.controlPoints = make([][]HomoPoint, len(this.controlPoints))

		for i := range reversed.controlPoints {
			col := make([]HomoPoint, 0, len(this.controlPoints[i]))

			for j := len(col) - 1; j >= 0; j-- {
				col = append(col, this.controlPoints[i][j])
			}

			reversed.controlPoints[i] = col
		}
	} else {
		reversed.knotsU = this.knotsU.Reversed()
		reversed.knotsV = append(KnotVec(nil), this.knotsV...)

		cpts := make([][]HomoPoint, 0, len(this.controlPoints))

		for i := len(this.controlPoints) - 1; i >= 0; i-- {
			cpts = append(cpts, append([]HomoPoint(nil), this.controlPoints[i]...))
		}
	}

	return &reversed
}

func (this *NurbsSurface) Transform(mat *mat4.T) *NurbsSurface {
	pts := Dehomogenize2d(this.controlPoints)

	for i := range pts {
		for j := range pts[0] {
			pts[i][j] = mat.MulVec3(&pts[i][j])
		}
	}

	return &NurbsSurface{
		this.degreeU,
		this.degreeV,

		Homogenize2d(pts, Weight2d(this.controlPoints)),

		append(KnotVec(nil), this.knotsU...),
		append(KnotVec(nil), this.knotsV...),
	}
}

func (this *NurbsSurface) knotRefine(knotsToInsert KnotVec, useV bool) *NurbsSurface {
	// TODO: make this faster by taking advantage of repeat computations in every row
	// 			 i.e. no reason to recompute the knot vectors on every row

	var newPts [][]HomoPoint
	var knots KnotVec
	var degree int
	var ctrlPts [][]HomoPoint

	// u dir
	if !useV {
		ctrlPts = transposed(this.controlPoints)
		knots = this.knotsU
		degree = this.degreeU
		// v dir
	} else {
		ctrlPts = this.controlPoints
		knots = this.knotsV
		degree = this.degreeV
	}

	// do knot refinement on every row
	baseCurve := NurbsCurve{degree: degree, knots: knots}
	var c *NurbsCurve
	for _, cptrow := range ctrlPts {
		baseCurve.controlPoints = cptrow
		c = baseCurve.knotRefine(knotsToInsert)
		newPts = append(newPts, c.controlPoints)
	}

	newknots := c.knots

	// u dir
	if !useV {
		newPts = transposed(newPts)
		return &NurbsSurface{
			this.degreeU, this.degreeV,
			newPts,
			newknots, append(KnotVec(nil), this.knotsV...),
		}
		// v dir
	} else {
		return &NurbsSurface{
			this.degreeU, this.degreeV,
			newPts,
			append(KnotVec(nil), this.knotsU...), newknots,
		}
	}
}

func (this *NurbsSurface) Split(u float64, useV bool) (*NurbsSurface, *NurbsSurface) {
	var (
		knots         KnotVec
		degree        int
		controlPoints [][]HomoPoint
	)

	if !useV {
		controlPoints = transposed(this.controlPoints)
		knots = this.knotsU
		degree = this.degreeU
	} else {
		controlPoints = this.controlPoints
		knots = this.knotsV
		degree = this.degreeV
	}

	knotsToInsert := make([]float64, degree+1)
	for i := range knotsToInsert {
		knotsToInsert[i] = u
	}

	newpts0, newpts1 := make([][]HomoPoint, 0), make([][]HomoPoint, 0)

	s := knots.Span(degree, u)
	var res *NurbsCurve
	baseCurve := NurbsCurve{degree: degree, knots: knots}

	for _, cps := range controlPoints {
		baseCurve.controlPoints = cps
		res = baseCurve.knotRefine(knotsToInsert)

		newpts0 = append(newpts0, res.controlPoints[:s+1])
		newpts1 = append(newpts1, res.controlPoints[s+1:])
	}

	knots0, knots1 := res.knots[:s+degree+2], res.knots[s+1:]

	if !useV {
		newpts0 = transposed(newpts0)
		newpts1 = transposed(newpts1)

		return &NurbsSurface{
				degree, this.degreeV,
				newpts0,
				knots0, append(KnotVec(nil), this.knotsV...),
			}, &NurbsSurface{
				degree, this.degreeV,
				newpts1,
				knots1, append(KnotVec(nil), this.knotsV...),
			}
	}

	// v dir
	return &NurbsSurface{
			this.degreeU, degree,
			newpts0,
			append(KnotVec(nil), this.knotsU...), knots0,
		}, &NurbsSurface{
			this.degreeU, degree,
			newpts1,
			append(KnotVec(nil), this.knotsU...), knots1,
		}
}

func transposed(mat [][]HomoPoint) (result [][]HomoPoint) {
	result = make([][]HomoPoint, len(mat))
	for i := range result {
		result[i] = make([]HomoPoint, len(mat[0]))
	}

	for col := 0; col < len(mat)-1; col++ {
		for row := col + 1; row < len(mat[0]); row++ {
			result[col][row], result[row][col] = result[row][col], result[col][row]
		}
	}

	return
}

// Compute the derivatives at a point on a NURBS surface
//
// **params**
// + NurbsSurfaceData object representing the surface
// + u parameter
// + v parameter
//
// **returns**
// + a Vector represented by an array of length (dim)
func (this *NurbsSurface) Normal(uv UV) vec3.T {
	derivs := this.Derivatives(uv, 1)
	return vec3.Cross(&derivs[1][0], &derivs[0][1])
}

// Compute the derivatives at a point on a NURBS surface
//
// **params**
// + NurbsSurfaceData object representing the surface
// + number of derivatives to evaluate
// + u parameter at which to evaluate the derivatives
// + v parameter at which to evaluate the derivatives
//
// **returns**
// + a point represented by an array of length (dim)
func (this *NurbsSurface) Derivatives(uv UV, numDerivs int) [][]vec3.T {
	ders := this.nonRationalDerivatives(uv, numDerivs)
	wders := Weight2d(ders)
	skl := make([][]vec3.T, numDerivs+1)

	for k := 0; k <= numDerivs; k++ {
		for l := 0; l <= numDerivs-k; l++ {
			v := ders[k][l].Vec3

			for j := 1; j <= l; j++ {
				scaled := skl[k][l-j].Scaled(binomial(l, j) * wders[0][j])
				v.Sub(&scaled)
			}

			for i := 1; i <= k; i++ {
				scaled := skl[k-i][l].Scaled(binomial(k, i) * wders[i][0])
				v.Sub(&scaled)

				var v2 vec3.T

				for j := 1; j <= l; j++ {
					scaled := skl[k-i][l-j].Scaled(binomial(l, j) * wders[i][j])
					v2.Add(&scaled)
				}

				scaled = v2.Scaled(binomial(k, i))
				v.Sub(&scaled)
			}

			v.Scale(1 / wders[0][0])
			skl[k][l] = v
		}
	}

	return skl
}

//
// Compute a point on a NURBS surface
//
// **params**
// + integer degree of surface in u direction
// + array of nondecreasing knot values in u direction
// + integer degree of surface in v direction
// + array of nondecreasing knot values in v direction
// + 3d array of control points (tensor), top to bottom is increasing u direction, left to right is increasing v direction,
// and where each control point is an array of length (dim+1)
// + u parameter at which to evaluate the surface point
// + v parameter at which to evaluate the surface point
//
// **returns**
// + a point represented by an array of length (dim)
func (this *NurbsSurface) Point(uv UV) vec3.T {
	homoPt := this.nonRationalPoint(uv)
	return homoPt.Dehomogenized()
}

// Compute the derivatives on a non-uniform, non-rational B spline surface
//
// **params**
// + NurbsSurfaceData object representing the surface
// + number of derivatives to evaluate
// + u parameter at which to evaluate the derivatives
// + v parameter at which to evaluate the derivatives
//
// **returns**
// + a 2d jagged array representing the derivatives - u derivatives increase by row, v by column
func (this *NurbsSurface) nonRationalDerivatives(uv UV, numDerivs int) [][]HomoPoint {
	n := len(this.knotsU) - this.degreeU - 2
	m := len(this.knotsV) - this.degreeV - 2

	return this.nonRationalDerivativesGivenNM(n, m, uv, numDerivs)
}

// Compute the derivatives on a non-uniform, non-rational B spline surface
// (corresponds to algorithm 3.6 from The NURBS book, Piegl & Tiller 2nd edition)
//
// **params**
// + integer number of basis functions in u dir - 1 = knotsU.length - degreeU - 2
// + integer number of basis functions in v dir - 1 = knotsU.length - degreeU - 2
// + NurbsSurfaceData object representing the surface
// + u parameter at which to evaluate the derivatives
// + v parameter at which to evaluate the derivatives
//
// **returns**
// + a 2d jagged array representing the derivatives - u derivatives increase by row, v by column
func (this *NurbsSurface) nonRationalDerivativesGivenNM(n, m int, uv UV, numDerivs int) [][]HomoPoint {
	degreeU := this.degreeU
	degreeV := this.degreeV
	controlPoints := this.controlPoints
	knotsU := this.knotsU
	knotsV := this.knotsV

	if !areValidRelations(degreeU, len(controlPoints), len(knotsU)) ||
		!areValidRelations(degreeV, len(controlPoints[0]), len(knotsV)) {
		panic("Invalid relations between control points, knot vector, and n")
	}

	var du, dv int
	if numDerivs < degreeU {
		du = numDerivs
	} else {
		du = degreeU
	}
	if numDerivs < degreeV {
		dv = numDerivs
	} else {
		dv = degreeV
	}

	skl := make([][]HomoPoint, du+1)
	for i := range skl {
		skl[i] = make([]HomoPoint, dv+1)
	}

	knotSpanIndexU := knotsU.SpanGivenN(n, degreeU, uv[0])
	knotSpanIndexV := knotsV.SpanGivenN(m, degreeV, uv[1])
	uders := DerivativeBasisFunctionsGivenNI(knotSpanIndexU, uv[0], degreeU, n, knotsU)
	vders := DerivativeBasisFunctionsGivenNI(knotSpanIndexV, uv[1], degreeV, m, knotsV)
	temp := make([]HomoPoint, degreeV+1)
	var dd int

	for k := 0; k <= du; k++ {
		for s := range temp {
			temp[s] = HomoPoint{}

			for r := 0; r <= degreeU; r++ {
				scaled := controlPoints[knotSpanIndexU-degreeU+r][knotSpanIndexV-degreeV+s]
				scaled.Scale(uders[k][r])
				temp[s].Add(&scaled)
			}
		}

		nk := numDerivs - k
		if nk < dv {
			dd = nk
		} else {
			dd = dv
		}

		for l := 0; l <= dd; l++ {
			skl[k][l] = HomoPoint{}

			for s := 0; s <= degreeV; s++ {
				scaled := temp[s].Scale(vders[l][s])
				skl[k][l].Add(scaled)
			}
		}
	}

	return skl
}

// Compute a point on a non-uniform, non-rational B-spline surface
//
// **params**
// + NurbsSurfaceData object representing the surface
// + u parameter at which to evaluate the surface point
// + v parameter at which to evaluate the surface point
//
// **returns**
// + a point represented by an array of length (dim)
func (this *NurbsSurface) nonRationalPoint(uv UV) HomoPoint {
	n := len(this.knotsU) - this.degreeU - 2
	m := len(this.knotsV) - this.degreeV - 2

	return this.nonRationalPointGivenNM(n, m, uv)
}

// Compute a point on a non-uniform, non-rational B spline surface
// (corresponds to algorithm 3.5 from The NURBS book, Piegl & Tiller 2nd edition)
//
// **params**
// + integer number of basis functions in u dir - 1 = knotsU.length - degreeU - 2
// + integer number of basis functions in v dir - 1 = knotsV.length - degreeV - 2
// + NurbsSurfaceData object representing the surface
// + u parameter at which to evaluate the surface point
// + v parameter at which to evaluate the surface point
//
// **returns**
// + a point represented by an array of length (dim)
func (this *NurbsSurface) nonRationalPointGivenNM(n, m int, uv UV) HomoPoint {
	degreeU := this.degreeU
	degreeV := this.degreeV
	controlPoints := this.controlPoints
	knotsU := this.knotsU
	knotsV := this.knotsV

	if !areValidRelations(degreeU, len(controlPoints), len(knotsU)) ||
		!areValidRelations(degreeV, len(controlPoints[0]), len(knotsV)) {
		panic("Invalid relations between control points, knot vector, and n")
	}

	knotSpanIndexU := knotsU.SpanGivenN(n, degreeU, uv[0])
	knotSpanIndexV := knotsV.SpanGivenN(m, degreeV, uv[1])
	uBasisVals := BasisFunctionsGivenKnotSpanIndex(knotSpanIndexU, uv[0], degreeU, knotsU)
	vBasisVals := BasisFunctionsGivenKnotSpanIndex(knotSpanIndexV, uv[1], degreeV, knotsV)
	uind := knotSpanIndexU - degreeU
	vind := knotSpanIndexV
	var position HomoPoint

	for l := 0; l <= degreeV; l++ {
		temp := HomoPoint{}
		vind = knotSpanIndexV - degreeV + l

		// sample u isoline
		for k := 0; k <= degreeU; k++ {
			scaled := controlPoints[uind+k][vind]
			scaled.Scale(uBasisVals[k])
			temp.Add(&scaled)
		}

		// add point from u isoline
		temp.Scale(vBasisVals[l])
		position.Add(&temp)
	}

	return position
}

// Extract the boundary curves from a surface
//
// **returns**
// + an array containing 4 elements, first 2 curves in the V direction, then 2 curves in the U direction
func (this *NurbsSurface) Boundaries() []*NurbsCurve {
	return []*NurbsCurve{
		this.Isocurve(this.knotsU[0], false),
		this.Isocurve(this.knotsU[len(this.knotsU)-1], false),
		this.Isocurve(this.knotsV[0], true),
		this.Isocurve(this.knotsV[len(this.knotsV)-1], true),
	}
}

func (this *NurbsSurface) Isocurve(u float64, useV bool) *NurbsCurve {
	var knots KnotVec
	if useV {
		knots = this.knotsV
	} else {
		knots = this.knotsU
	}

	var degree int
	if useV {
		degree = this.degreeV
	} else {
		degree = this.degreeU
	}

	knotMults := knots.Multiplicities()

	// if the knot already exists in the array, don't make duplicates
	reqKnotIndex := -1
	for i, knotMult := range knotMults {
		if math.Abs(u-knotMult.Knot) < Epsilon {
			reqKnotIndex = i
			break
		}
	}

	numKnotsToInsert := degree + 1
	if reqKnotIndex >= 0 {
		numKnotsToInsert -= knotMults[reqKnotIndex].Mult
	}

	// insert the knots
	var newSrf *NurbsSurface
	if numKnotsToInsert > 0 {
		newKnots := make(KnotVec, numKnotsToInsert)
		for i := range newKnots {
			newKnots[i] = u
		}

		newSrf = this.knotRefine(newKnots, useV)
	} else {
		newSrf = this
	}

	// obtain the correct index of control points to extract
	span := knots.Span(degree, u)

	if math.Abs(u-knots[0]) < Epsilon {
		span = 0
	} else if math.Abs(u-knots[len(knots)-1]) < Epsilon {
		if useV {
			span = len(newSrf.controlPoints[0])
		} else {
			span = len(newSrf.controlPoints) - 1
		}
	}

	if useV {
		controlPoints := make([]HomoPoint, 0, len(newSrf.controlPoints))
		for _, row := range newSrf.controlPoints {
			controlPoints = append(controlPoints, row[span])
		}

		return &NurbsCurve{newSrf.degreeU, controlPoints, newSrf.knotsU}
	}

	return &NurbsCurve{newSrf.degreeV, newSrf.controlPoints[span], newSrf.knotsV}
}
