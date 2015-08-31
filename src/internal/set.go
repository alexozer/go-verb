package internal

import "math"

type Set []float64

func (this Set) SortedUnion(set Set) (merged Set) {
	merged = make(Set, 0)

	var thisI, setI int
	for thisI < len(this) || setI < len(set) {
		if thisI >= len(this) {
			merged = append(merged, set[setI])
			setI++
			continue
		} else if setI >= len(set) {
			merged = append(merged, this[thisI])
			thisI++
			continue
		}

		diff := this[thisI] - set[setI]

		if math.Abs(diff) < Epsilon {
			merged = append(merged, this[thisI])
			thisI++
			setI++
			continue
		}

		if diff > 0.0 {
			// add the smaller
			merged = append(merged, set[setI])
			setI++
			continue
		}

		// thus diff < 0.0
		merged = append(merged, this[thisI])
		thisI++

	}

	return
}

// a is superset, hence it is always longer or equal
func (this Set) SortedSub(set Set) (result Set) {
	result = make(Set, 0)

	var thisI, setI int

	for thisI < len(this) {

		if setI >= len(set) {
			result = append(result, this[thisI])
			thisI++
			continue
		}

		if math.Abs(this[thisI]-set[setI]) < Epsilon {
			thisI++
			setI++
			continue
		}

		result = append(result, this[thisI])
		thisI++
	}

	return
}
