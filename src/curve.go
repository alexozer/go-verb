package verb

import (
	"errors"
	"math"
	"math/rand"

	. "github.com/alexozer/verb/internal"

	"github.com/ungerik/go3d/float64/mat4"
	"github.com/ungerik/go3d/float64/vec3"
)

type (
	CurvePoint struct {
		U  float64
		Pt vec3.T
	}
)

type NurbsCurve struct {
	// degree of curve
	degree int

	// slice of control points, each a homogeneous coordinate
	controlPoints []HomoPoint

	// slice of nondecreasing knot values
	knots KnotVec
}

func NewNurbsCurve(degree int, controlPoints []vec3.T, weights []float64, knots []float64) (*NurbsCurve, error) {
	this := NewNurbsCurveUnchecked(degree, controlPoints, weights, knots)
	if err := this.check(); err != nil {
		return nil, err
	}

	return this, nil
}

func NewNurbsCurveUnchecked(degree int, controlPoints []vec3.T, weights []float64, knots []float64) *NurbsCurve {
	return &NurbsCurve{degree, Homogenize1d(controlPoints, weights), KnotVec(knots).Clone()}
}

func (this *NurbsCurve) Degree() int {
	return this.degree
}

func (this *NurbsCurve) ControlPoints() []vec3.T {
	return Dehomogenize1d(this.controlPoints)
}

func (this *NurbsCurve) Weights() []float64 {
	return Weight1d(this.controlPoints)
}

func (this *NurbsCurve) Knots() []float64 {
	return []float64(this.knots.Clone())
}

// clone() is not exported because NurbsCurve is immutable to the client,
// so there's no point in making a deep copy.
// Should only be used when control points and knots can't be shared
func (this *NurbsCurve) clone() *NurbsCurve {
	return &NurbsCurve{
		degree:        this.degree,
		controlPoints: append([]HomoPoint(nil), this.controlPoints...),
		knots:         this.knots.Clone(),
	}
}

// Determine the valid domain of the curve
//
// **returns**
// + An array representing the high and end point of the domain of the curve
func (this *NurbsCurve) Domain() (min, max float64) {
	min = this.knots[0]
	max = this.knots[len(this.knots)-1]
	return
}

// Split a curve into two parts
//
// **params**
// + NurbsCurveData object representing the curve
// + location to split the curve
//
// **returns**
// + *Array* two new curves, defined by degree, knots, and control points
//
func (this *NurbsCurve) Split(u float64) (*NurbsCurve, *NurbsCurve) {
	degree, knots := this.degree, this.knots

	knotsToInsert := make(KnotVec, degree+1)
	for i := range knotsToInsert {
		knotsToInsert[i] = u
	}
	res := this.knotRefine(knotsToInsert)

	s := knots.Span(degree, u)

	knots0 := res.knots[:s+degree+2 : s+degree+2]
	knots1 := res.knots[s+1:]

	cpts0 := res.controlPoints[:s+1 : s+1]
	cpts1 := res.controlPoints[s+1:]

	return &NurbsCurve{degree, cpts0, knots0}, &NurbsCurve{degree, cpts1, knots1}
}

// Insert a collection of knots on a curve
//
// Corresponds to Algorithm A5.4 (Piegl & Tiller)
//
// **params**
// + NurbsCurveData object representing the curve
// + array of knots to insert
//
// **returns**
// +  NurbsCurveData object representing the curve
//
func (this *NurbsCurve) knotRefine(knotsToInsert KnotVec) *NurbsCurve {
	if len(knotsToInsert) == 0 {
		return this.clone()
	}

	degree := this.degree
	controlPoints := this.controlPoints
	knots := this.knots

	n := len(controlPoints) - 1
	m := n + degree + 1
	r := len(knotsToInsert) - 1
	a := knots.Span(degree, knotsToInsert[0])
	b := knots.Span(degree, knotsToInsert[r])

	// TODO can we use slices instead of maps?
	controlPointsPost := make(map[int]HomoPoint)
	knotsPost := make(map[int]float64)

	// new control pts
	for i := 0; i <= a-degree; i++ {
		controlPointsPost[i] = controlPoints[i]
	}

	for i := b - 1; i <= n; i++ {
		controlPointsPost[i+r+1] = controlPoints[i]
	}

	// new knot vector
	for i := 0; i <= a; i++ {
		knotsPost[i] = knots[i]
	}

	for i := b + degree; i <= m; i++ {
		knotsPost[i+r+1] = knots[i]
	}

	i := b + degree - 1
	k := b + degree + r
	j := r

	for j >= 0 {
		for knotsToInsert[j] <= knots[i] && i > a {
			controlPointsPost[k-degree-1] = controlPoints[i-degree-1]
			knotsPost[k] = knots[i]
			k = k - 1
			i = i - 1
		}

		controlPointsPost[k-degree-1] = controlPointsPost[k-degree]

		for l := 1; l <= degree; l++ {
			ind := k - degree + l
			alfa := knotsPost[k+l] - knotsToInsert[j]

			if math.Abs(alfa) < Epsilon {
				controlPointsPost[ind-1] = controlPointsPost[ind]
			} else {
				alfa /= knotsPost[k+l] - knots[i-degree+l]

				oldP, newP := controlPointsPost[ind], controlPointsPost[ind-1]
				controlPointsPost[ind-1] = HomoPoint{
					vec3.Interpolate(&oldP.Vec3, &newP.Vec3, alfa),
					controlPointsPost[ind-1].W,
				}
			}
		}

		knotsPost[k] = knotsToInsert[j]
		k = k - 1

		j--
	}

	var maxCptsI int
	for i := range controlPointsPost {
		if i > maxCptsI {
			maxCptsI = i
		}
	}
	finalCpts := make([]HomoPoint, maxCptsI+1)
	for i, cpt := range controlPointsPost {
		finalCpts[i] = cpt
	}

	var maxKnotsI int
	for i := range knotsPost {
		if i > maxKnotsI {
			maxKnotsI = i
		}
	}
	finalKnots := make(KnotVec, maxKnotsI+1)
	for i, knot := range knotsPost {
		finalKnots[i] = knot
	}

	return &NurbsCurve{degree, finalCpts, finalKnots}
}

func (this *NurbsCurve) ClosestPoint(p vec3.T) vec3.T {
	return this.Point(this.ClosestParam(p))
}

func (this *NurbsCurve) ClosestParam(p vec3.T) float64 {
	//  We want to solve:
	//
	//   C'(u) * ( C(u) - P ) = 0 = f(u)
	//
	//  C(u) is the curve, p is the point, * is a dot product
	//
	// We'll use newton's method:
	//
	// 	 u* = u - f / f'
	//
	// We use the product rule in order to form the derivative, f':
	//
	//	f' = C"(u) * ( C(u) - p ) + C'(u) * C'(u)
	//
	// What is the conversion criteria? (Piegl & Tiller suggest)
	//
	// |C(u) - p| < e1
	//
	// |C'(u)*(C(u) - P)|
	// ------------------  < e2
	// |C'(u)| |C(u) - P|
	//
	//  1) first check 2 & 3
	// 	2) if at least one of these is not, compute new value, otherwise halt
	// 	3) ensure the parameter stays within range
	// 			* if not closed, don't allow outside of range a-b
	// 			* if closed (e.g. circle), allow to move back to beginning
	//  4)  if |(u* - u)C'(u)| < e1, halt
	//

	min := math.MaxFloat64
	var u float64

	pts := this.regularSample(len(this.controlPoints) * this.degree)

	for i := 0; i < len(pts)-1; i++ {
		u0, u1 := pts[i].U, pts[i+1].U

		p0 := pts[i].Pt
		p1 := pts[i+1].Pt

		proj := segmentClosestPoint(&p, &p0, &p1, u0, u1)
		dv := vec3.Sub(&p, &proj.Pt)
		d := dv.Length()

		if d < min {
			min = d
			u = proj.U
		}
	}

	maxits := 5
	var i int
	var e []vec3.T
	eps1, eps2 := 0.0001, 0.0005
	var dif vec3.T
	minu, maxu := this.knots[0], this.knots[len(this.knots)-1]

	firstCtrlPt := this.controlPoints[0].Dehomogenized()
	lastCtrlPt := this.controlPoints[len(this.controlPoints)-1].Dehomogenized()
	closed := vec3.SquareDistance(&firstCtrlPt, &lastCtrlPt) < Epsilon

	cu := u

	f := func(u float64) []vec3.T {
		return this.Derivatives(u, 2)
	}

	n := func(u float64, e []vec3.T, d vec3.T) float64 {
		//   C'(u) * ( C(u) - P ) = 0 = f(u)
		f := vec3.Dot(&e[1], &d)

		//	f' = C"(u) * ( C(u) - p ) + C'(u) * C'(u)
		s0 := vec3.Dot(&e[2], &d)
		s1 := vec3.Dot(&e[1], &e[1])
		df := s0 + s1

		return u - f/df
	}

	for i < maxits {
		e = f(cu)
		dif = vec3.Sub(&e[0], &p)

		// |C(u) - p| < e1
		c1v := dif.Length()

		// C'(u) * (C(u) - P)
		// ------------------ < e2
		// |C'(u)| |C(u) - P|
		c2n := vec3.Dot(&e[1], &dif)
		c2d := e[1].Length() * c1v

		c2v := c2n / c2d

		c1 := c1v < eps1
		c2 := math.Abs(c2v) < eps2

		// if both tolerances are met
		if c1 && c2 {
			return cu
		}

		ct := n(cu, e, dif)

		// are we outside of the bounds of the curve?
		if ct < minu {
			if closed {
				ct = maxu - (ct - minu)
			} else {
				ct = minu
			}
		} else if ct > maxu {
			if closed {
				ct = minu + (ct - maxu)
			} else {
				ct = maxu
			}
		}

		// will our next step force us out of the curve?
		c3vv := e[1].Scaled(ct - cu)
		c3v := c3vv.Length()

		if c3v < eps1 {
			return cu
		}

		cu = ct
		i++

	}

	return cu
}

// Determine the parameter of the curve at the given arc length
//
// **params**
// + The arc length at which to determine the parameter
//
// **returns**
// + The length of the curve at the given parameter
func (this *NurbsCurve) ParamAtArcLength(length float64) (u float64) {
	return this.paramAtArcLength(length, nil, nil, nil)
}

// tol is optional; nil uses the default tolerance
func (this *NurbsCurve) paramAtArcLength(length float64, tol *float64, beziers []*NurbsCurve, bezierLengths []float64) float64 {
	if length < Epsilon {
		return this.knots[0]
	}

	var crvs []*NurbsCurve
	if beziers != nil {
		crvs = beziers
	} else {
		crvs = this.decomposeIntoBeziers()
	}
	var i int
	//cc := crvs[i]
	cl := -Epsilon
	var bezier_lengths []float64
	if bezierLengths != nil {
		bezier_lengths = bezierLengths
	} else {
		bezier_lengths = make([]float64, 0)
	}

	// iterate through the curves consuming the bezier's, summing their length along the way
	for cl < length && i < len(crvs) {
		if i < len(bezier_lengths) {
			bezier_lengths[i] = bezier_lengths[i]
		} else {
			bezier_lengths[i] = this.bezierArcLength(nil, nil)
		}

		cl += bezier_lengths[i]

		if length < cl+Epsilon {
			return this.bezierParamAtArcLength(length, tol, &bezier_lengths[i])
		}

		i++
	}

	return -1
}

//
// Get the curve parameter at an arc length
//
// **params**
// + NurbsCurveData object representing the curve
// + the arc length to find the parameter
// + the tolerance - increasing the tolerance can make this computation quite expensive
// + the total length of the curve, if already computed
//
// **returns**
// + the parameter
func (this *NurbsCurve) bezierParamAtArcLength(length float64, tol, totalLength *float64) float64 {
	if length < 0 {
		return this.knots[0]
	}

	// we compute the whole length. if desired length is outside of that, give up
	var totalLen float64
	if totalLength != nil {
		totalLen = *totalLength
	} else {
		totalLen = this.bezierArcLength(nil, nil)
	}

	if length > totalLen {
		return this.knots[len(this.knots)-1]
	}

	// divide & conquer
	// TODO can we use derivative?
	startP, startL := this.knots[0], 0.0
	endP, endL := this.knots[len(this.knots)-1], totalLen
	var midP, midL float64

	var _tol float64
	if tol != nil {
		_tol = *tol
	} else {
		_tol = Tolerance * 2
	}

	for endL-startL > _tol {
		midP = (startP + endP) / 2
		midL = this.bezierArcLength(&midP, nil)

		if midL > length {
			endP = midP
			endL = midL
		} else {
			startP = midP
			startL = midL
		}

	}

	return (startP + endP) / 2
}

// Determine the arc length of the curve
//
// **returns**
// + The length of the curve
func (this *NurbsCurve) Length() float64 {
	return this.arcLength(nil, nil)
}

// Determine the arc length of the curve at the given parameter
//
// **params**
// + The parameter at which to evaluate
//
// **returns**
// + The length of the curve at the given parameter
func (this *NurbsCurve) LengthAtParam(u float64) float64 {
	return this.arcLength(&u, nil)
}

//
// Approximate the length of a rational curve by gaussian quadrature - assumes a smooth curve
//
// **params**
// + NurbsCurveData object representing the curve
// + the parameter at which to approximate the length
// + the degree of gaussian quadrature to perform - a higher number yields a more exact result
//
// **returns**
// + the approximate length
//
func (this *NurbsCurve) arcLength(u *float64, gaussDegIncrease *int) float64 {
	var _u float64
	if u == nil {
		_u = this.knots[len(this.knots)-1]
	} else {
		_u = *u
	}

	crvs := this.decomposeIntoBeziers()
	var i int
	cc := crvs[0]
	var sum float64

	for i < len(crvs) && cc.knots[0]+Epsilon < _u {
		param := math.Min(cc.knots[len(cc.knots)-1], _u)
		sum += this.bezierArcLength(&param, gaussDegIncrease)
		i++
		cc = crvs[i]
	}

	return sum
}

//
// Approximate the length of a rational bezier curve by gaussian quadrature - assumes a smooth curve
//
// **params**
// + NurbsCurveData object representing the curve
// + the parameter at which to approximate the length
// + the degree of gaussian quadrature to perform - a higher number yields a more exact result
//
// **returns**
// + the approximate length
//
func (this *NurbsCurve) bezierArcLength(u *float64, gaussDegIncrease *int) float64 {
	var _gaussDegIncrease int
	if gaussDegIncrease != nil {
		_gaussDegIncrease = *gaussDegIncrease
	} else {
		_gaussDegIncrease = 16
	}

	var _u float64
	if u == nil {
		_u = this.knots[len(this.knots)-1]
	} else {
		_u = *u
	}

	z := (_u - this.knots[0]) / 2
	var sum float64
	gaussDeg := this.degree + _gaussDegIncrease
	var cu float64
	var tan []vec3.T

	for i := 0; i < gaussDeg; i++ {
		cu = z*tValues[gaussDeg][i] + z + this.knots[0]
		tan = this.Derivatives(cu, 1)

		sum += cValues[gaussDeg][i] * tan[1].Length()

	}

	return z * sum
}

// Validate a NurbsCurveData object
//
// **params**
// + The data object
//
// **returns**
// + The original, unmodified data
func (this *NurbsCurve) check() error {
	if this.controlPoints == nil {
		return errors.New("Control points array cannot be nil!")
	}

	if this.degree < 1 {
		return errors.New("Degree must be greater than 1!")
	}

	if this.knots == nil {
		return errors.New("Knots cannot be nil!")
	}

	if len(this.knots) != len(this.controlPoints)+this.degree+1 {
		return errors.New("len(controlPoints) + degree + 1 must equal len(knots)!")
	}

	if !this.knots.IsValid(this.degree) {
		return errors.New("Invalid knot vector format! Should begin with degree + 1 repeats and end with degree + 1 repeats!")
	}

	return nil
}

type CurveLengthSample struct {
	U, Len float64
}

func (this *NurbsCurve) DivideByEqualArcLength(num int) []CurveLengthSample {
	tlen := this.arcLength(nil, nil)
	inc := tlen / float64(num)

	return this.DivideByArcLength(inc)
}

func (this *NurbsCurve) DivideByArcLength(l float64) []CurveLengthSample {
	crvs := this.decomposeIntoBeziers()

	crvlens := make([]float64, len(crvs))
	for i, crv := range crvs {
		crvlens[i] = crv.bezierArcLength(nil, nil)
	}

	var totlen float64
	for _, length := range crvlens {
		totlen += length
	}

	pts := []CurveLengthSample{{this.knots[0], 0}}

	if l > totlen {
		return pts
	}

	inc := l
	lc := inc
	var runsum, runsum1 float64
	var u float64

	for i, crvlen := range crvlens {
		runsum += crvlen

		for lc < runsum+Epsilon {
			tol := Tolerance
			u = crvs[i].bezierParamAtArcLength(lc-runsum1, &tol, &crvlen)

			pts = append(pts, CurveLengthSample{u, lc})
			lc += inc

		}

		runsum1 += crvlen
	}

	return pts
}

//
// Sample a NURBS curve at equally spaced parametric intervals
//
// **params**
// + NurbsCurveData object
// + integer number of samples
// + whether to prefix the point with the parameter
//
// **returns**
// + an array of points, prepended by the point param if required
//
func (this *NurbsCurve) regularSample(numSamples int) []CurvePoint {
	return this.regularSampleRange(
		this.knots[0], this.knots[len(this.knots)-1],
		numSamples,
	)
}

//
// Sample a range of a NURBS curve at equally spaced parametric intervals
//
// **params**
// + NurbsCurveData object
// + start parameter for sampling
// + end parameter for sampling
// + integer number of samples
// + whether to prefix the point with the parameter
//
// **returns**
// + an dictionary of parameter - point pairs
//
func (this *NurbsCurve) regularSampleRange(start, end float64, numSamples int) []CurvePoint {
	if numSamples < 1 {
		numSamples = 2
	}

	samples := make([]CurvePoint, numSamples)
	span := (end - start) / float64(numSamples-1)
	var u float64

	for i := range samples {
		u = start + span*float64(i)

		samples[i] = CurvePoint{u, this.Point(u)}
	}

	return samples
}

func (this *NurbsCurve) Tessellate() []CurvePoint {
	return this.adaptiveSample(nil, nil)
}

//
// Sample a NURBS curve over its entire domain, corresponds to http://ariel.chronotext.org/dd/defigueiredo93adaptive.pdf
//
// **params**
// + NurbsCurveData object
// + tol for the adaptive scheme
// + whether to prefix the point with the parameter
//
// **returns**
// + an array of dim + 1 length where the first element is the param where it was sampled and the remaining the pt
//
func (this *NurbsCurve) adaptiveSample(tol *float64, includeU *bool) []CurvePoint {
	var _tol float64
	if tol == nil {
		_tol = 1e-6
	} else {
		_tol = *tol
	}

	var _includeU bool
	if includeU != nil {
		_includeU = *includeU
	}

	// if degree is 1, just return the dehomogenized control points
	if this.degree == 1 {
		samples := make([]CurvePoint, len(this.controlPoints))
		// the first element of each array is the parameter
		for i := range this.controlPoints {
			samples[i] = CurvePoint{this.knots[i+1], this.controlPoints[i].Dehomogenized()}
		}
		return samples
	}

	return this.adaptiveSampleRange(this.knots[0], this.knots[len(this.knots)-1], _tol, _includeU)
}

//
// Sample a NURBS curve at 3 points, facilitating adaptive sampling
//
// **params**
// + NurbsCurveData object
// + start parameter for sampling
// + end parameter for sampling
// + whether to prefix the point with the parameter
//
// **returns**
// + an array of dim + 1 length where the first element is the param where it was sampled and the remaining the pt
//
func (this *NurbsCurve) adaptiveSampleRange(start, end, tol float64, includeU bool) []CurvePoint {
	// sample curve at three pts
	p1, p3 := this.Point(start), this.Point(end)
	t := 0.5 + 0.2*rand.Float64()
	mid := start + (end-start)*t
	p2 := this.Point(mid)

	// if the two end control points are coincident, the three point test will always return 0, let's split the curve
	diff := vec3.Sub(&p1, &p3)
	diff2 := vec3.Sub(&p1, &p2)

	// the first condition checks if the curve makes up a loop, if so, we will need to continue evaluation
	if (vec3.Dot(&diff, &diff) < tol && vec3.Dot(&diff2, &diff2) > tol) || !threePointsAreCollinear(&p1, &p2, &p3, tol) {

		// get the exact middle
		exact_mid := start + (end-start)*0.5

		// recurse on the two halves
		left_pts := this.adaptiveSampleRange(start, exact_mid, tol, includeU)
		right_pts := this.adaptiveSampleRange(exact_mid, end, tol, includeU)

		// concatenate the two
		leftEnd := len(left_pts) - 1
		return append(left_pts[:leftEnd:leftEnd], right_pts...)

	} else {
		return []CurvePoint{{start, p1}, {end, p3}}
	}
}

func (this *NurbsCurve) Reverse() *NurbsCurve {
	reversed := NurbsCurve{
		degree:        this.degree,
		controlPoints: make([]HomoPoint, 0, len(this.controlPoints)),
		knots:         this.knots.Reversed(),
	}

	for i := len(this.controlPoints) - 1; i >= 0; i-- {
		reversed.controlPoints = append(reversed.controlPoints, this.controlPoints[i])
	}

	return &reversed
}

func (this *NurbsCurve) elevateDegree(finalDegree int) *NurbsCurve {
	if finalDegree <= this.degree {
		return this.clone()
	}

	// args
	n := len(this.knots) - this.degree - 2
	newDegree := this.degree
	knots := this.knots
	controlPoints := this.controlPoints
	degreeInc := finalDegree - this.degree

	// intermediate values
	rows, cols := newDegree+degreeInc+1, newDegree+1
	bezalfs := make([][]float64, rows)
	for i := range bezalfs {
		bezalfs[i] = make([]float64, cols)
	}

	bpts := make([]HomoPoint, 0)
	ebpts := make([]HomoPoint, 0)
	Nextbpts := make([]HomoPoint, 0)
	//alphas := make([]float64, 0)

	m := n + newDegree + 1
	ph := finalDegree
	ph2 := ph / 2

	Qw := make([]HomoPoint, 0)
	Uh := make(KnotVec, 0)

	bezalfs[0][0] = 1
	bezalfs[ph][newDegree] = 1

	for i := 1; i <= ph2; i++ {
		inv := 1 / binomial(ph, i)
		mpi := imin(newDegree, i)
		for j := imax(0, i-degreeInc); j <= mpi; j++ {
			bezalfs[i][j] = inv * binomial(newDegree, j) * binomial(degreeInc, i-j)
		}
	}
	for i := ph2 + 1; i < ph; i++ {
		mpi := imin(newDegree, i)
		for j := imax(0, i-degreeInc); j <= mpi; j++ {
			bezalfs[i][j] = bezalfs[ph-i][newDegree-j]
		}
	}
	//mh := ph
	kind := ph + 1
	r := -1
	a := newDegree
	b := newDegree + 1
	cind := 1
	ua := knots[0]
	Qw[0] = controlPoints[0]
	for i := 0; i <= ph; i++ {
		Uh[i] = ua
	}
	for i := 0; i <= newDegree; i++ {
		bpts[i] = controlPoints[i]
	}
	for b < m {
		i := b
		for b < m && knots[b] == knots[b+1] {
			b = b + 1
		}
		mul := b - i + 1
		//mh := mh + mul + degreeInc
		ub := knots[b]
		oldr := r
		r = newDegree - mul
		// check for integer arithmetic
		var lbz int
		if oldr > 0 {
			lbz = (oldr + 2) / 2
		} else {
			lbz = 1
		}
		var rbz int
		if r > 0 {
			rbz = ph - (r+1)/2
		} else {
			rbz = ph
		}
		if r > 0 {
			numer := ub - ua
			alfs := make(map[int]float64) // TODO can this use a slice instead?
			k := newDegree
			for k > mul {
				alfs[k-mul-1] = numer / (knots[a+k] - ua) // integer arithmetic?
				k--
			}
			for j := 1; j <= r; j++ {
				save := r - j
				s := mul + j
				k := newDegree
				for k >= s {
					bpts[k] = HomoInterpolated(&bpts[k], &bpts[k-1], alfs[k-s])
					k--
				}
				Nextbpts[save] = bpts[newDegree]
			}
		}

		for i := lbz; i <= ph; i++ {
			ebpts[i] = HomoPoint{}
			mpi := imin(newDegree, i)
			for j := imax(0, i-degreeInc); j <= mpi; j++ {
				alf := bezalfs[i][j]
				bptsj := HomoPoint{bpts[j].Vec3.Scaled(alf), alf * bpts[j].W}
				ebpts[i].Vec3.Add(&bptsj.Vec3)
				ebpts[i].W += bptsj.W
			}
		}

		if oldr > 1 {
			first := kind - 2
			last := kind
			den := ub - ua
			bet := (ub - Uh[kind-1]) / den // integer arithmetic
			for tr := 1; tr < oldr; tr++ {
				i := first
				j := last
				kj := j - kind + 1
				for j-i > tr {
					if i < cind {
						alf := (ub - Uh[i]) / (ua - Uh[i]) // integer arithmetic
						Qw[i] = HomoInterpolated(&Qw[i], &Qw[i-1], alf)
					}
					if j >= lbz {
						if j-tr <= kind-ph+oldr {
							gam := (ub - Uh[j-tr]) / den
							ebpts[kj] = HomoInterpolated(&ebpts[kj], &ebpts[kj+1], gam)
						}
					} else {
						ebpts[kj] = HomoInterpolated(&ebpts[kj], &ebpts[kj+1], bet)
					}
					i = i + 1
					j = j - 1
					kj = kj - 1
				}
				first = first - 1
				last = last + 1
			}
		}

		if a != newDegree {
			for i := 0; i < ph-oldr; i++ {
				Uh[kind] = ua
				kind = kind + 1
			}
		}

		for j := lbz; j <= rbz; j++ {
			Qw[cind] = ebpts[j]
			cind = cind + 1
		}

		if b < m {
			for j := 0; j < r; j++ {
				bpts[j] = Nextbpts[j]
			}
			for j := r; j <= newDegree; j++ {
				bpts[j] = controlPoints[b-newDegree+j]
			}
			a = b
			b = b + 1
			ua = ub
		} else {
			for i := 0; i <= ph; i++ {
				Uh[kind+i] = ub
			}
		}
	}

	return &NurbsCurve{finalDegree, Qw, Uh}
}

func UnifyCurveKnotVectors(curves []*NurbsCurve) (unified []*NurbsCurve) {
	var maxDegree int
	for _, curve := range curves {
		if curve.degree > maxDegree {
			maxDegree = curve.degree
		}
	}

	// elevate all curves to the same degree
	unified = make([]*NurbsCurve, len(curves))
	for i, curve := range curves {
		if curve.degree < maxDegree {
			unified[i] = curve.elevateDegree(maxDegree)
		} else {
			unified[i] = &NurbsCurve{
				curve.degree,
				curve.controlPoints,
				curve.knots.Clone(),
			}
		}
	}

	var maxSpan float64
	for _, curve := range unified {
		min, max := curve.knots[0], curve.knots[len(curve.knots)-1]

		// shift all knot vectors to start at 0.0
		for iKnot, knot := range curve.knots {
			curve.knots[iKnot] = knot - min
		}

		// find the max knot span
		maxSpan = math.Max(maxSpan, max-min)
	}

	// scale all of the knot vectors to match
	for _, curve := range unified {
		scale := maxSpan / (curve.knots[len(curve.knots)-1] - curve.knots[0])
		for iKnot := range curve.knots {
			curve.knots[iKnot] *= scale
		}
	}

	// merge all of the knot vectors
	mergedKnotSet := Set(unified[0].knots)
	for _, curve := range unified[1:] {
		mergedKnotSet = mergedKnotSet.SortedUnion(Set(curve.knots))
	}

	// knot refinement on each curve
	for i, curve := range unified {
		rem := KnotVec(mergedKnotSet.SortedSub(Set(curve.knots)))
		unified[i] = curve.knotRefine(rem)
	}

	return
}

func (this *NurbsCurve) Transform(mat *mat4.T) *NurbsCurve {
	pts := Dehomogenize1d(this.controlPoints)

	for i := range pts {
		pts[i] = mat.MulVec3(&pts[i])
	}

	return &NurbsCurve{
		this.degree,
		Homogenize1d(pts, Weight1d(this.controlPoints)),
		this.knots,
	}
}

// Decompose a NURBS curve into a collection of bezier's.  Useful
// as each bezier fits into it's convex hull.  This is a useful starting
// point for intersection, closest point, divide & conquer algorithms
//
// **params**
// + NurbsCurveData object representing the curve
//
// **returns**
// + *Array* of NurbsCurveData objects, defined by degree, knots, and control points
func (this *NurbsCurve) decomposeIntoBeziers() []*NurbsCurve {
	degree := this.degree
	controlPoints := this.controlPoints
	knots := this.knots

	// find all of the unique knot values and their multiplicity
	// for each, increase their multiplicity to degree + 1

	knotmults := knots.Multiplicities()
	reqMult := degree + 1

	// insert the knots
	baseCurve := NurbsCurve{degree: degree}
	for _, knotmult := range knotmults {
		if knotmult.Mult < reqMult {
			knotsToInsert := make(KnotVec, reqMult-knotmult.Mult)
			for i := range knotsToInsert {
				knotsToInsert[i] = knotmult.Knot
			}
			baseCurve.knots = knots
			baseCurve.controlPoints = controlPoints
			res := baseCurve.knotRefine(knotsToInsert)

			knots = res.knots
			controlPoints = res.controlPoints
		}
	}

	//numCrvs := len(knots)/reqMult - 1
	crvKnotLength := reqMult * 2

	crvs := make([]*NurbsCurve, 0, len(controlPoints)/reqMult)

	for i := 0; i < len(controlPoints); i += reqMult {
		kts := knots[i : i+crvKnotLength : i+crvKnotLength]
		pts := controlPoints[i : i+reqMult : i+reqMult]

		crvs = append(crvs, &NurbsCurve{degree, pts, kts})
	}

	return crvs
}

// Insert a knot along a rational curve.  Note that this algorithm only works
// for r + s <= degree, where s is the initial multiplicity (number of duplicates) of the knot.
//
// Corresponds to algorithm A5.1 (Piegl & Tiller)
//
// Use the curveKnotRefine for applications like curve splitting.
//
// **params**
// + integer degree
// + array of nondecreasing knot values
// + array of control points
// + parameter at which to insert the knot
// + number of times to insert the knot
//
// **returns**
// + *Object* the new curve, defined by knots and controlPoints
//
func (this *NurbsCurve) insertKnot(u float64, r int) *NurbsCurve {
	degree := this.degree
	controlPoints := this.controlPoints
	knots := this.knots

	// numPts is num control points for the initial curve
	// k is the span on which the knots are inserted
	// s is the initial multiplicity of the knot
	// r is the number of times to insert the knot
	// controlPoints is initial set of control points

	s := 0 // assume original multiplicity is 0 - TODO add check for multiplicity in knots

	numPts := len(controlPoints)
	k := knots.Span(degree, u)                         // the span in which the knot will be inserted
	numPtsPost := numPts + r                           // a new control pt for every new knot
	controlPointsTemp := make([]HomoPoint, degree-s)   // new Array( degree - s )
	knotsPost := make(KnotVec, len(knots)+r)           // new Array( knots.length + r )  // r new knots
	controlPointsPost := make([]HomoPoint, numPtsPost) // new Array( numPtsPost )

	// new knot vector

	// insert the k knots that will not be affected
	for i := 1; i <= k; i++ {
		knotsPost[i] = knots[i]
	}

	// insert the new repeat knots
	for i := 1; i <= r; i++ {
		knotsPost[k+i] = u
	}

	// insert the rest of the knots
	for i := k + 1; i < len(knots); i++ {
		knotsPost[i+r] = knots[i]
	}

	// control point generation

	// copy the original control points before the insertion span
	for i := 0; i <= k-degree; i++ {
		controlPointsPost[i] = controlPoints[i]
	}

	// copy the original controls after the insertion span
	for i := k - s; i < numPts; i++ {
		controlPointsPost[i+r] = controlPoints[i]
	}

	// collect the affected control points in this temporary array
	for i := 0; i <= degree-s; i++ {
		controlPointsTemp[i] = controlPoints[k-degree+i]
	}

	var L int
	var alpha float64

	// insert knot r times
	for j := 1; j <= r; j++ {
		L = k - degree + j

		for i := 0; i <= degree-j-s; i++ {
			alpha = (u - knots[L+i]) / (knots[i+k+1] - knots[L+i])

			controlPointsTemp[i] = HomoInterpolated(
				&controlPointsTemp[i],
				&controlPointsTemp[i+1],
				alpha,
			)
		}

		controlPointsPost[L] = controlPointsTemp[0]
		controlPointsPost[k+r-j-s] = controlPointsTemp[degree-j-s]

	}

	// not so confident about this part
	for i := L + 1; i < k-s; i++ {
		controlPointsPost[i] = controlPointsTemp[i-L]
	}

	return &NurbsCurve{degree, controlPointsPost, knotsPost}
}

// Compute the tangent at a point on a NURBS curve
//
// **params**
// + NurbsCurveData object representing the curve
// + u parameter
// + v parameter
//
// **returns**
// + a Vector represented by an array of length (dim)
func (this *NurbsCurve) Tangent(u float64) vec3.T {
	return this.Derivatives(u, 1)[1]
}

//
// Determine the derivatives of a NURBS curve at a given parameter
//
// **params**
// + NurbsCurveData object representing the curve - the control points are in homogeneous coordinates
// + parameter on the curve at which the point is to be evaluated
// + number of derivatives to evaluate
//
// **returns**
// + a point represented by an array of length (dim)
//
func (this *NurbsCurve) Derivatives(u float64, numDerivs int) []vec3.T {
	ders := this.nonRationalDerivatives(u, numDerivs)
	ck := make([]vec3.T, 0, numDerivs+1)

	for k := 0; k <= numDerivs; k++ {
		v := ders[k].Vec3

		for i := 1; i <= k; i++ {
			scaled := ck[k-i].Scaled(binomial(k, i) * ders[i].W)
			v.Sub(&scaled)
		}
		v.Scale(1 / ders[0].W)
		ck = append(ck, v)
	}

	return ck
}

// Compute a point on a NURBS curve
//
// **params**
// + integer degree of curve
// + array of nondecreasing knot values
// + 2d array of homogeneous control points, where each control point is an array of length (dim+1)
// and form (wi*pi, wi)
// + parameter on the curve at which the point is to be evaluated
//
// **returns**
// + a point represented by an array of length (dim)
func (this *NurbsCurve) Point(u float64) vec3.T {
	homoPt := this.nonRationalPoint(u)
	return homoPt.Dehomogenized()
}

// Determine the derivatives of a non-uniform, non-rational B-spline curve at a given parameter
//
// **params**
// + NurbsCurveData object representing the curve
// + parameter on the curve at which the point is to be evaluated
// + number of derivatives to evaluate
//
// **returns**
// + a point represented by an array of length (dim)
func (this *NurbsCurve) nonRationalDerivatives(u float64, numDerivs int) []HomoPoint {
	n := len(this.knots) - this.degree - 2
	return this.nonRationalDerivativesGivenNM(n, u, numDerivs)
}

// Determine the derivatives of a non-uniform, non-rational B-spline curve at a given parameter
// (corresponds to algorithm 3.1 from The NURBS book, Piegl & Tiller 2nd edition)
//
// **params**
// + integer number of basis functions - 1 = knots.length - degree - 2
// + NurbsCurveData object representing the curve
// + parameter on the curve at which the point is to be evaluated
//
// **returns**
// + a point represented by an array of length (dim)
func (this *NurbsCurve) nonRationalDerivativesGivenNM(n int, u float64, numDerivs int) []HomoPoint {
	degree := this.degree
	controlPoints := this.controlPoints
	knots := this.knots

	if !areValidRelations(degree, len(controlPoints), len(knots)) {
		panic("Invalid relations between control points, knot vector, and n")
	}

	var du int
	if numDerivs < degree {
		du = numDerivs
	} else {
		du = degree
	}

	ck := make([]HomoPoint, du+1)
	knotSpanIndex := knots.SpanGivenN(n, degree, u)
	nders := DerivativeBasisFunctionsGivenNI(knotSpanIndex, u, degree, du, knots)

	for k := 0; k <= du; k++ {
		for j := 0; j <= degree; j++ {
			scaled := controlPoints[knotSpanIndex-degree+j]
			scaled.Scale(nders[k][j])
			ck[k].Add(&scaled)
		}
	}

	return ck
}

// Compute a point on a non-uniform, non-rational b-spline curve
//
// **params**
// + NurbsCurveData object representing the curve
// + parameter on the curve at which the point is to be evaluated
//
// **returns**
// + a point represented by an array of length (dim)
func (this *NurbsCurve) nonRationalPoint(u float64) HomoPoint {
	n := len(this.knots) - this.degree - 2
	return this.nonRationalPointGivenN(n, u)
}

// Compute a point on a non-uniform, non-rational b-spline curve
// (corresponds to algorithm 3.1 from The NURBS book, Piegl & Tiller 2nd edition)
//
// **params**
// + integer number of basis functions - 1 = knots.length - degree - 2
// + NurbsCurveData object representing the curve
// + parameter on the curve at which the point is to be evaluated
//
// **returns**
// + a point represented by an array of length (dim)
func (this *NurbsCurve) nonRationalPointGivenN(n int, u float64) HomoPoint {
	degree := this.degree
	controlPoints := this.controlPoints
	knots := this.knots

	if !areValidRelations(degree, len(controlPoints), len(knots)) {
		panic("Invalid relations between control points, knot Array, and n")
	}

	knotSpanIndex := knots.SpanGivenN(n, degree, u)
	basisValues := BasisFunctionsGivenKnotSpanIndex(knotSpanIndex, u, degree, knots)
	var position HomoPoint

	for j := 0; j <= degree; j++ {
		scaled := controlPoints[knotSpanIndex-degree+j]
		scaled.Scale(basisValues[j])
		position.Add(&scaled)
	}

	return position
}
