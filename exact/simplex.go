package simplex

import (
	"errors"
	"fmt"
	"math"
	"math/big"
)

// const epsilon = 0.0001

// func equals(a, b float64) bool {
// return equaleps(a, b, epsilon)
// }

// func equaleps(a, b float64, eps float64) bool {
// return (a-b) < eps && (b-a) < eps
// }

func equals(a, b *big.Rat) bool {
	return a.Cmp(b) == 0
}

func makerats(n int) []*big.Rat {
	xs := make([]*big.Rat, n)
	for i := range xs {
		xs[i] = &big.Rat{}
	}
	return xs
}

func makevec(n int) Vec {
	return Vec(makerats(n))
}

var (
	ErrUnbounded  = errors.New("simplex: problem is unbounded")
	ErrInfeasible = errors.New("simplex: problem is infeasible")

	OSR = errors.New("simplex: optimal solution reached")
)

type Var struct {
	int
	Val *big.Rat
}

type Coef struct {
	C int
	V Var
}

type Relation int

const (
	GreaterEq Relation = iota - 1
	Equal
	LessEq
)

type Constraint struct {
	Cs []Coef
	Eq Relation
	To *big.Rat
}

func Constrain(vs ...Coef) Constraint {
	return Constraint{Cs: vs, To: &big.Rat{}}
}

func (cn Constraint) LessEq(x int64) Constraint {
	cn.Eq = LessEq
	cn.To = &big.Rat{}
	cn.To.SetInt64(x)
	return cn
}

func (cn Constraint) GreaterEq(x int64) Constraint {
	cn.Eq = GreaterEq
	cn.To = &big.Rat{}
	cn.To.SetInt64(x)
	return cn
}

func (cn Constraint) Equal(x int64) Constraint {
	cn.Eq = Equal
	cn.To = &big.Rat{}
	cn.To.SetInt64(x)
	return cn
}

type Program struct {
	c Vec // coefficients
	a Mat // constraints
	b Vec // equalities

	s Mat // surplus/slack
	r Mat // artificials

	tbl   Mat // tableau
	isMin bool
}

func (prg *Program) Var(coef int64) Var {
	tmp := big.NewRat(coef, 1)
	prg.c = append(prg.c, tmp)
	return Var{len(prg.c) - 1, &big.Rat{}}
}

func (prg *Program) Z() *big.Rat {
	lr := prg.tbl[len(prg.tbl)-1]
	z := lr[len(lr)-1]
	if prg.isMin {
		return z
	} else {
		tmp := &big.Rat{}
		return tmp.Mul(z, minusone) // return as minus for maximization, see HACK in iter.
	}
}

func (prg *Program) AddConstraints(cns ...Constraint) {
	for _, cn := range cns {
		// if rhs is negative, multiply by -1
		if cn.To.Cmp(zero) == -1 {
			for i := range cn.Cs {
				cn.Cs[i].C *= -1
			}
			switch cn.Eq {
			case GreaterEq:
				cn.Eq = LessEq
			case LessEq:
				cn.Eq = GreaterEq
			case Equal:
				panic("equality must be greater than or equal to zero")
			}
		}

		// TODO handle unrestricted variables

		// vcn := make(Vec, len(prg.c))
		vcn := makevec(len(prg.c))
		for _, v := range cn.Cs {
			vcn[v.V.int].SetInt64(int64(v.C))
		}
		prg.a = append(prg.a, vcn)
		prg.b = append(prg.b, cn.To)

		// introduce slack/surplus
		if cn.Eq != Equal {
			prg.s = append(prg.s, make([]Vec, len(prg.b)-len(prg.s))...)
			prg.s[0] = append(prg.s[0], &big.Rat{})
			for i := range prg.s {
				n := len(prg.s[0]) - len(prg.s[i])
				if n > 0 {
					// prg.s[i] = append(prg.s[i], make([]*big.Rat, n)...)
					prg.s[i] = append(prg.s[i], makerats(n)...)
				}
			}
			s := prg.s[len(prg.s)-1]
			// TODO check for nil
			s[len(s)-1].SetInt64(int64(cn.Eq))
		}

		// introduce artificials for >= and ==
		// TODO it's not necessary to blindly introduce artificials, this is already
		// accounted for in the iter method where non-basics that create an identity
		// will be used as basics initially for a feasible solution to start.
		// TODO artificials should account for cn.To being strictly positive
		if cn.Eq != LessEq {
			prg.r = append(prg.r, make([]Vec, len(prg.b)-len(prg.r))...)
			prg.r[0] = append(prg.r[0], &big.Rat{})
			for i := range prg.r {
				n := len(prg.r[0]) - len(prg.r[i])
				if n > 0 {
					// prg.r[i] = append(prg.r[i], make([]*big.Rat, n)...)
					prg.r[i] = append(prg.r[i], makerats(n)...)
				}
			}
			r := prg.r[len(prg.r)-1]
			r[len(r)-1] = big.NewRat(1, 1)
		}
	}
}

func (prg *Program) init(min bool) Mat {
	if len(prg.s) > 0 {
		ns := len(prg.b) - len(prg.s)
		for i := 0; i < ns; i++ {
			// prg.s = append(prg.s, make(Vec, len(prg.s[0])))
			prg.s = append(prg.s, makevec(len(prg.s[0])))
		}
	}

	if len(prg.r) > 0 {
		nr := len(prg.b) - len(prg.r)
		for i := 0; i < nr; i++ {
			// prg.r = append(prg.r, make(Vec, len(prg.r[0])))
			prg.r = append(prg.r, makevec(len(prg.r[0])))
		}
	}

	tbl := make(Mat, len(prg.b))
	for i := range tbl {
		tbl[i] = append(tbl[i], prg.a[i]...)
		if len(prg.s) != 0 {
			tbl[i] = append(tbl[i], prg.s[i].Copy()...)
		}
		tmp := &big.Rat{}
		tbl[i] = append(tbl[i], tmp.Set(prg.b[i]))
	}
	// tbl = append(tbl, make(Vec, len(tbl[0])))
	tbl = append(tbl, makevec(len(tbl[0])))
	n := len(tbl) - 1
	if min {
		for j, x := range prg.c {
			tmp := &big.Rat{}
			// tbl[n][j] = -x
			tbl[n][j] = tmp.Mul(x, minusone)
		}
	} else {
		copy(tbl[n], prg.c)
	}

	return tbl
}

func (prg *Program) iter() error {
	// scratch space to store pivot row multiplied by non-pivot row's pivot column.
	// dt := make(Vec, len(prg.tbl[0]))
	dt := makevec(len(prg.tbl[0]))

	// scratch space to perform theta calc
	// bnd := make(Vec, len(prg.tbl)-1)
	bnd := makevec(len(prg.tbl) - 1)
	// col := make(Vec, len(prg.tbl)-1)
	col := makevec(len(prg.tbl) - 1)

	var err error
	var pi, pj int                      // pivot indices
	var px, py = &big.Rat{}, &big.Rat{} // pivot values
	for {
		// find pivot indices
		prg.tbl.Column(&bnd, len(prg.tbl[0])-1)
		pi, pj = -1, -1
		px.SetInt64(0)
		py.SetInt64(0)
		cz := prg.tbl[len(prg.tbl)-1]
		for j, y := range cz[:len(cz)-1] {
			// don't select if ident since this var has already entered basis
			if prg.tbl.IsColIdent(j) {
				continue
			}
			// TODO near-zero cases are a common point of error during iteration but ultimately
			// it'd be worth investigating the use of something like math/big for exact results.
			if equals(y, zero) {
				y = &big.Rat{}
			}
			if y.Cmp(py) == 1 || (py.Cmp(zero) == 0 && y.Cmp(zero) == 0 && j < len(prg.c)) {
				prg.tbl.Column(&col, j)
				pj = j
				py.Set(y)
			}
		}
		if py.Cmp(zero) == -1 || pj == -1 {
			err = OSR
			break
		}

		// Theta may equal zero if there is a tie for minimum (degeneracy), this is still valid but wasteful.
		if pi, px = Theta(bnd, col); px.Cmp(zero) == -1 || pi == -1 {
			err = ErrUnbounded
			break
		}

		// divide pivot row by pivot element so pivot element == 1
		div := &big.Rat{}
		div.Set(prg.tbl[pi][pj])
		if div.Cmp(zero) != 1 {
			// if here, either pivot indices don't match pivot checks or pivot checks incorrect.
			panic(fmt.Errorf("divisor <= 0: %.2f", div))
		}
		for j := range prg.tbl[pi] {
			prg.tbl[pi][j] = prg.tbl[pi][j].Quo(prg.tbl[pi][j], div)
		}

		// adjust so pivot column is all zeros (except pivot element) and becomes an identity column.
		//
		// HACK apply same transform to calculate Cj-Zj (last row) understanding that
		// real Z equals -Z from tableau.
		for i := range prg.tbl {
			if i != pi {
				prg.tbl[pi].MulScalar(&dt, prg.tbl[i][pj])
				prg.tbl[i].Sub(&prg.tbl[i], dt)
			}
		}
	}

	if err != OSR {
		return err
	}
	return nil
}

// isBasicFeasible returns true if basic vars give feasible solution; that is,
// the solution values are greater than zero.
// TODO this isn't exactly a valid check
func (prg *Program) isBasicFeasible() bool {
	if len(prg.s) == 0 {
		return false
	}
	for _, row := range prg.s {
		for _, y := range row {
			if y.Cmp(zero) == -1 {
				return false
			}
		}
	}
	return true
}

func (prg *Program) twophase(min bool) error {
	// TODO this is all pretty lame but isolated
	prg.tbl.Transpose(&prg.tbl)
	prg.r.Transpose(&prg.r)
	for _, row := range prg.r {
		// xs := make([]*big.Rat, len(row)+1)
		xs := makerats(len(row) + 1)
		copy(xs, row)
		prg.tbl.Insert(len(prg.tbl)-1, xs)
	}
	prg.r.Transpose(&prg.r)
	prg.tbl.Transpose(&prg.tbl)

	// make new objective function for artificials
	i := len(prg.tbl) - 1
	offset := len(prg.c)
	if len(prg.s) != 0 {
		offset += len(prg.s[0])
	}
	end := len(prg.tbl[0]) - 1
	for j := range prg.tbl[i] {
		if j < offset || j == end {
			prg.tbl[i][j] = &big.Rat{}
		} else if min {
			prg.tbl[i][j] = big.NewRat(-1, 1)
		} else {
			prg.tbl[i][j] = big.NewRat(1, 1)
		}
	}

	// calculate Cj-Zj
	// HACK calculate Cj-Zj as in iter; so, fake pivot on artificials to calc
	// last row which affects sign of final Z if maximization.
	//
	// calculating Cj-Zj by hand
	//
	//       0    0         -1   -1
	//    | x1 | x2 | ... | a1 | a2
	// --------------
	// a1 |  2 |  3 |
	// a2 |  5 |  2 |
	// --------------
	//    |  7 |  5 |
	//
	// C1 = x1 = 0
	// Z1 = (a1*x1_a1) + (a2*x1_a2)
	// Z1 = (-1*2) + (-1*5) = -2 + -5 = -7
	// C1-Z1 = 0 - -7 = 7
	//
	// C2 = x1 = 0
	// Z2 = (a1*x2_a1) + (a2*x2_a2)
	// Z2 = (-1*3) + (-1*2) = -3 + -2 = -5
	// C2-Z2 = 0 - -5 = 5
	//
	// dt := make(Vec, len(prg.tbl[0]))
	dt := makevec(len(prg.tbl[0]))
	for ri, rrow := range prg.r {
		for rj, rx := range rrow {
			if rx.Cmp(one) == 0 {
				pi, pj := ri, offset+rj
				prg.tbl[pi].MulScalar(&dt, prg.tbl[i][pj])
				prg.tbl[i].Sub(&prg.tbl[i], dt)
			}
		}
	}

	// find basic feasible
	if err := prg.iter(); err == ErrUnbounded {
		// TODO return value of the artificial still in basis to provide hint to user on
		// adjusting constraint equality value so that a feasible solution is possible.
		return ErrInfeasible
	}

	// drop artificials
	prg.tbl.Transpose(&prg.tbl)
	prg.tbl = append(prg.tbl[:offset], prg.tbl[len(prg.tbl)-1])
	prg.tbl.Transpose(&prg.tbl)
	// dt = make(Vec, len(prg.tbl[0]))
	dt = makevec(len(prg.tbl[0]))

	// calc new Cj-Zj with original coefs
	for xi := range prg.tbl[i] {
		prg.tbl[i][xi] = &big.Rat{}
	}
	if min {
		for xi, x := range prg.c {
			tmp := &big.Rat{}
			prg.tbl[i][xi] = tmp.Set(x).Mul(tmp, minusone)
		}
	} else {
		copy(prg.tbl[i], prg.c)
	}
	for xi, xrow := range prg.tbl {
		for xj, xy := range xrow {
			if xy.Cmp(one) == 0 { // safe to assume first 1 seen is for non-basic of original
				prg.tbl[xi].MulScalar(&dt, prg.tbl[i][xj])
				prg.tbl[i].Sub(&prg.tbl[i], dt)
				continue
			}
		}
	}

	return prg.iter()
}

func (prg *Program) Minimize() error {
	prg.isMin = true
	prg.tbl = prg.init(true)
	if !prg.isBasicFeasible() {
		return prg.twophase(true)
	}
	return prg.iter()
}

func (prg *Program) Maximize() error {
	prg.isMin = false
	prg.tbl = prg.init(false)
	if !prg.isBasicFeasible() {
		return prg.twophase(false)
	}
	return prg.iter()
}

func (prg *Program) For(vars ...*Var) {
	j := len(prg.tbl[0]) - 1
	for _, v := range vars {
		for i, x := range prg.tbl.ColumnAlloc((*v).int) {
			if x.Cmp(one) == 0 {
				(*v).Val = prg.tbl[i][j]
				break
			}
		}
	}
}

// Theta returns the minimum result of division ops of only positive divisors.
func Theta(a, b Vec) (int, *big.Rat) {
	if len(a) != len(b) {
		panic("a Vec length does not match b Vec length")
	}

	// TODO If there is more than one row for which the minimum is achieved
	// then a dropping variable choice rule can be used to make the determination.
	// .. choosing variable with smallest subscript is easy way to avoid cycling

	j := -1
	theta := big.NewRat(math.MaxInt64, 1)
	for i, x := range b {
		if x.Cmp(zero) == 1 {
			y := &big.Rat{}
			y.Quo(a[i], x)
			if theta.Cmp(y) == 1 {
				// if y := a[i] / x; theta > y {
				j = i
				theta = y
			}
		}
	}
	return j, theta
}
