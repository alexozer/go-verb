package verb

var binomCache map[[2]int]float64

func init() {
	binomCache = make(map[[2]int]float64)
}

func binomial(n, k int) float64 {
	if k == 0 {
		return 1
	}

	if n == 0 || k > n {
		return 0
	}

	if k > n-k {
		k = n - k // optimization
	}

	if result, ok := binomCache[[2]int{n, k}]; ok {
		return result
	}

	nO := n
	var r float64
	for d := 1; d <= k; d++ {
		if cacheR, ok := binomCache[[2]int{nO, d}]; ok {
			n--
			r = cacheR
			continue
		}

		r *= float64(n) / float64(d)
		n--

		binomCache[[2]int{nO, d}] = r

	}

	return r
}

func binomialNoCache(n, k int) float64 {
	if k == 0 {
		return 1
	}

	if n == 0 || k > n {
		return 0
	}

	if k > n-k {
		k = n - k // optimization
	}

	r := 1.0
	for d := 1; d <= k; d++ {
		r *= float64(n) / float64(d)
		n--
	}

	return r
}
