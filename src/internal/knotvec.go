package internal

import (
	"math"
)

type KnotVec []float64

func (this KnotVec) Clone() KnotVec {
	return append(KnotVec(nil), this...)
}

func (this KnotVec) Domain() float64 {
	return this[len(this)-1] - this[0]
}

// Find the span on the knot Array without supplying n
//
// **params**
// + integer degree of function
// + float parameter
// + array of nondecreasing knot values
//
// **returns**
// + the index of the knot span
//
func (this KnotVec) Span(degree int, u float64) int {
	m := len(this) - 1
	n := m - degree - 1

	return this.SpanGivenN(n, degree, u)
}

// Find the span on the knot Array knots of the given parameter
// (corresponds to algorithm 2.1 from The NURBS book, Piegl & Tiller 2nd edition)
//
// **params**
// + integer number of basis functions - 1 = knots.length - degree - 2
// + integer degree of function
// + parameter
// + array of nondecreasing knot values
//
// **returns**
// + the index of the knot span
//
func (this KnotVec) SpanGivenN(n int, degree int, u float64) int {
	if u >= this[n+1] {
		return n
	}

	if u < this[degree] {
		return degree
	}

	low, high := degree, n+1
	mid := (low + high) / 2

	for u < this[mid] || u >= this[mid+1] {
		if u < this[mid] {
			high = mid
		} else {
			low = mid
		}

		mid = (low + high) / 2
	}

	return mid
}

//
// Determine the multiplicities of the values in a knot vector
//
// **params**
// + array of nondecreasing knot values
//
// **returns**
// + *Array* of length 2 arrays, [knotValue, knotMultiplicity]
//
func (this KnotVec) Multiplicities() []KnotMultiplicity {
	mults := []KnotMultiplicity{{this[0], 0}}

	var currI int
	for _, knot := range this {
		if math.Abs(knot-mults[currI].Knot) > Epsilon {
			mults = append(mults, KnotMultiplicity{knot, 0})
			currI++
		}

		mults[currI].Mult++
	}

	return mults
}

func (this KnotVec) IsValid(degree int) bool {
	if len(this) == 0 {
		return false
	}

	if len(this) < (degree+1)*2 {
		return false
	}

	rep := this[0]

	for _, knot := range this[:degree+1] {
		if math.Abs(knot-rep) > Epsilon {
			return false
		}
	}

	rep = this[len(this)-1]

	for _, knot := range this[len(this)-degree-1:] {
		if math.Abs(knot-rep) > Epsilon {
			return false
		}
	}

	return this.IsNonDecreasing()
}

func (this KnotVec) IsNonDecreasing() bool {
	rep := this[0]
	for _, knot := range this[1:] {
		if knot < rep-Epsilon {
			return false
		}
		rep = knot
	}
	return true
}

func (this KnotVec) Reversed() KnotVec {
	l := make(KnotVec, len(this))
	l[0] = this[0]

	length := len(this)
	for i := 1; i < length; i++ {
		l[i] = l[i-1] + (this[length-1] - this[length-i-1])
	}

	return l
}

type KnotMultiplicity struct {
	Knot float64
	Mult int
}
