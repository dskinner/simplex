package simplex

import (
	"bytes"
	"fmt"
	"math/big"
)

var (
	minusone = big.NewRat(-1, 1)
	zero     = &big.Rat{}
	one      = big.NewRat(1, 1)
)

type Vec []*big.Rat

func (a Vec) Copy() Vec {
	b := make(Vec, len(a))
	for i := range a {
		x := &big.Rat{}
		b[i] = x.Set(a[i])
	}
	return b
}

func (a Vec) Dot(b Vec) *big.Rat {
	if len(a) != len(b) {
		panic("a Vec length does not match b Vec length")
	}
	c := &big.Rat{}
	tmp := &big.Rat{}
	for i, x := range a {
		// c += x * b[i]
		c.Add(c, tmp.Mul(x, b[i]))

	}
	return c
}

func (a Vec) Sub(c *Vec, b Vec) {
	if len(a) != len(b) {
		panic("a Vec length does not match b Vec length")
	}
	for i, x := range a {
		// (*c)[i] = x - b[i]
		(*c)[i].Sub(x, b[i])
	}
}

func (a Vec) MulScalar(b *Vec, x *big.Rat) {
	for i, v := range a {
		// (*b)[i] = v * x
		(*b)[i].Mul(v, x)
	}
}

func (a Vec) Div(b Vec) (c Vec) {
	if len(a) != len(b) {
		panic("a Vec length does not match b Vec length")
	}
	for i, x := range a {
		tmp := &big.Rat{}
		c = append(c, tmp.Quo(x, b[i]))
	}
	return
}

func (a Vec) Max() (i int, x *big.Rat) {
	for j, y := range a {
		// if j == 0 || y > x {
		if j == 0 || y.Cmp(x) == 1 {
			// i, x = j, y
			i = j
			x.Set(y)
		}
	}
	return
}

func (v *Vec) Insert(i int, x *big.Rat) {
	*v = append(*v, nil)
	copy((*v)[i+1:], (*v)[i:])
	(*v)[i] = x
}

type Mat []Vec

func (m *Mat) Insert(i int, xs []*big.Rat) {
	*m = append(*m, nil)
	copy((*m)[i+1:], (*m)[i:])
	(*m)[i] = xs
}

func (m Mat) Column(v *Vec, j int) {
	n := len(*v)
	for i, row := range m {
		if i == n {
			break
		}
		// (*v)[i] = row[j]
		(*v)[i].Set(row[j])
	}
}

func (m Mat) ColumnAlloc(j int) (v Vec) {
	for _, row := range m {
		x := &big.Rat{}
		v = append(v, x.Set(row[j]))
	}
	return
}

func (m Mat) Transpose(a *Mat) {
	*a = make(Mat, len(m[0]))
	for j := range *a {
		(*a)[j] = m.ColumnAlloc(j)
	}
}

func (a Mat) MulVec(b Vec) (c Vec) {
	if len(a[0]) != len(b) {
		panic("row length of a Math does not equal b Vec length")
	}
	for _, r := range a {
		c = append(c, r.Dot(b))
	}
	return
}

func (a Mat) IsColIdent(j int) bool {
	var ok bool
	for _, row := range a {
		y := row[j]
		if y.Cmp(one) == 0 && !ok {
			ok = true
		} else if y.Cmp(zero) != 0 {
			return false
		}
	}
	return ok
}

func (m Mat) String() string {
	var buf bytes.Buffer
	for _, row := range m {
		buf.WriteString(fmt.Sprintf("% 8v\n", row))
	}
	n := buf.Len() - 1
	if n < 0 {
		n = 0
	}
	buf.Truncate(n)
	return buf.String()
}
