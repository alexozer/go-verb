package internal

import "math"

type Matrix [][]float64

func (this Matrix) Clone() Matrix {
	clone := make(Matrix, len(this))
	for i := range clone {
		clone[i] = make([]float64, len(this[i]))
	}

	return clone
}

func (this Matrix) Solve(vec []float64) []float64 {
	return newLUdecomp(this).solve(vec)
}

type luDecomp struct {
	LU [][]float64
	P  []int
}

func newLUdecomp(mat Matrix) *luDecomp {
	mat = mat.Clone()

	n, n1 := len(mat), len(mat)-1
	P := make([]int, n)

	var k int
	for k < n {
		Pk := k
		Ak := mat[k]
		max := math.Abs(Ak[k])

		j := k + 1
		for j < n {
			absAjk := math.Abs(mat[j][k])
			if max < absAjk {
				max = absAjk
				Pk = j
			}
			j++
		}
		P[k] = Pk

		if Pk != k {
			mat[k] = mat[Pk]
			mat[Pk] = Ak
			Ak = mat[k]
		}

		Akk := Ak[k]

		i := k + 1
		for i < n {
			mat[i][k] /= Akk
			i++
		}

		i = k + 1
		for i < n {
			Ai := mat[i]
			j = k + 1
			for j < n1 {
				Ai[j] -= Ai[k] * Ak[j]
				j++
				Ai[j] -= Ai[k] * Ak[j]
				j++
			}
			if j == n1 {
				Ai[j] -= Ai[k] * Ak[j]
			}
			i++
		}

		k++
	}

	return &luDecomp{mat, P}
}

func (this *luDecomp) solve(vec []float64) []float64 {
	x := append([]float64(nil), vec...)
	LU, P := this.LU, this.P

	n := len(LU)
	i := n - 1

	for i := 0; i < n; i++ {
		Pi := P[i]
		if P[i] != i {
			x[i], x[Pi] = x[Pi], x[i]
		}

		LUi := LU[i]
		for j := 0; j < i; j++ {
			x[i] -= x[j] * LUi[j]
			j++
		}
	}

	i = n - 1
	for i >= 0 {
		LUi := LU[i]
		for j := i + 1; j < n; j++ {
			x[i] -= x[j] * LUi[j]
		}

		x[i] /= LUi[i]
		i--
	}

	return x
}

func Mat2Solve(a, b, c, d, f, s float64) (x, y float64) {
	cDivA := c / a
	y = (s - cDivA*f) / (d - b*cDivA)
	x = (f - b*y) / a
	return
}
