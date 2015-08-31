package internal

// Compute the non-vanishing basis functions
//
// **params**
// + float parameter
// + integer degree of function
// + array of nondecreasing knot values
//
// **returns**
// + list of non-vanishing basis functions
//
func basisFunctions(u float64, degree int, knots KnotVec) []float64 {
	knotSpanIndex := knots.Span(degree, u)
	return BasisFunctionsGivenKnotSpanIndex(knotSpanIndex, u, degree, knots)
}

// Compute the non-vanishing basis functions
// (corresponds to algorithm 2.2 from The NURBS book, Piegl & Tiller 2nd edition)
//
// **params**
// + *Number*, integer knot span index
// + *Number*, float parameter
// + *Number*, integer degree of function
// + array of nondecreasing knot values
//
// **returns**
// + list of non-vanishing basis functions
//
func BasisFunctionsGivenKnotSpanIndex(knotSpanIndex int, u float64, degree int, knots KnotVec) []float64 {
	basisFunctions := make([]float64, degree+1)
	left := make([]float64, degree+1)
	right := make([]float64, degree+1)

	basisFunctions[0] = 1

	for j := 1; j <= degree; j++ {
		left[j] = u - knots[knotSpanIndex+1-j]
		right[j] = knots[knotSpanIndex+j] - u
		var saved float64

		for r := 0; r < j; r++ {
			temp := basisFunctions[r] / (right[r+1] + left[j-r])
			basisFunctions[r] = saved + right[r+1]*temp
			saved = left[j-r] * temp
		}

		basisFunctions[j] = saved
	}

	return basisFunctions
}

// Compute the non-vanishing basis functions and their derivatives
//
// **params**
// + float parameter
// + integer degree
// + array of nondecreasing knot values
//
// **returns**
// + 2d array of basis and derivative values of size (n+1, p+1) The nth row is the nth derivative and the first row is made up of the basis function values.
func derivativeBasisFunctions(u float64, degree int, knots KnotVec) [][]float64 {
	knotSpanIndex := knots.Span(degree, u)
	m := len(knots) - 1
	n := m - degree - 1

	return DerivativeBasisFunctionsGivenNI(knotSpanIndex, u, degree, n, knots)
}

// Compute the non-vanishing basis functions and their derivatives
// (corresponds to algorithm 2.3 from The NURBS book, Piegl & Tiller 2nd edition)
//
// **params**
// + integer knot span index
// + float parameter
// + integer degree
// + integer number of basis functions - 1 = knots.length - degree - 2
// + array of nondecreasing knot values
//
// **returns**
// + 2d array of basis and derivative values of size (n+1, p+1) The nth row is the nth derivative and the first row is made up of the basis function values.
func DerivativeBasisFunctionsGivenNI(knotSpanIndex int, u float64, p, n int, knots KnotVec) [][]float64 {
	ndu := zeros2d(p+1, p+1)

	left := make([]float64, p+1)
	right := make([]float64, p+1)

	ndu[0][0] = 1

	for j := 1; j <= p; j++ {
		left[j] = u - knots[knotSpanIndex+1-j]
		right[j] = knots[knotSpanIndex+j] - u
		var saved float64

		for r := 0; r < j; r++ {
			ndu[j][r] = right[r+1] + left[j-r]
			temp := ndu[r][j-1] / ndu[j][r]

			ndu[r][j] = saved + right[r+1]*temp
			saved = left[j-r] * temp

		}
		ndu[j][j] = saved
	}

	ders := zeros2d(n+1, p+1)

	for j := 0; j <= p; j++ {
		ders[0][j] = ndu[j][p]
	}

	a := zeros2d(2, p+1)
	var j1, j2 int

	for r := 0; r <= p; r++ {
		s1, s2 := 0, 1
		a[0][0] = 1

		for k := 1; k <= n; k++ {
			var d float64
			rk := r - k
			pk := p - k

			if r >= k {
				a[s2][0] = a[s1][0] / ndu[pk+1][rk]
				d = a[s2][0] * ndu[rk][pk]
			}

			if rk >= -1 {
				j1 = 1
			} else {
				j1 = -rk
			}

			if r-1 <= pk {
				j2 = k - 1
			} else {
				j2 = p - r
			}

			for j := j1; j <= j2; j++ {
				a[s2][j] = (a[s1][j] - a[s1][j-1]) / ndu[pk+1][rk+j]
				d += a[s2][j] * ndu[rk+j][pk]
			}

			if r <= pk {
				a[s2][k] = -a[s1][k-1] / ndu[pk+1][r]
				d += a[s2][k] * ndu[r][pk]
			}

			ders[k][r] = d

			s1, s2 = s2, s1
		}
	}

	acc := p
	for k := 1; k <= n; k++ {
		for j := 0; j <= p; j++ {
			ders[k][j] *= float64(acc)
		}
		acc *= (p - k)
	}

	return ders
}

func zeros2d(n, m int) [][]float64 {
	result := make([][]float64, n)
	for i := range result {
		result[i] = make([]float64, m)
	}

	return result
}
