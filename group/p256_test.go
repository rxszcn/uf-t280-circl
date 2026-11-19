package group

import (
	"crypto/rand"
	"math/big"
	"testing"
)

// Differential test: the Montgomery field arithmetic in p256.go must agree
// with big.Int modular arithmetic on random and edge-case inputs.
func TestP256MontField(t *testing.T) {
	p, _ := new(big.Int).SetString("FFFFFFFF00000001000000000000000000000000FFFFFFFFFFFFFFFFFFFFFFFF", 16)
	toBig := func(z *p256Mont) *big.Int {
		return new(big.Int).SetBits([]big.Word{
			big.Word(z[0]), big.Word(z[1]), big.Word(z[2]), big.Word(z[3]),
		})
	}
	for i := 0; i < 1000; i++ {
		a, _ := rand.Int(rand.Reader, p)
		b, _ := rand.Int(rand.Reader, p)
		var am, bm, cm p256Mont
		p256MontEncode(&am, a)
		p256MontEncode(&bm, b)

		if got := p256MontToBig(&am); got.Cmp(a) != 0 {
			t.Fatalf("roundtrip: got %x want %x", got, a)
		}
		p256MontMul(&cm, &am, &bm)
		want := new(big.Int).Mul(a, b)
		want.Mod(want, p)
		if got := p256MontToBig(&cm); got.Cmp(want) != 0 {
			t.Fatalf("mul: got %x want %x", got, want)
		}
		p256MontAdd(&cm, &am, &bm)
		want.Add(a, b).Mod(want, p)
		if got := p256MontToBig(&cm); got.Cmp(want) != 0 {
			t.Fatalf("add: got %x want %x", got, want)
		}
		p256MontSub(&cm, &am, &bm)
		want.Sub(a, b).Mod(want, p)
		if got := p256MontToBig(&cm); got.Cmp(want) != 0 {
			t.Fatalf("sub: got %x want %x", got, want)
		}
		p256MontNeg(&cm, &am)
		want.Neg(a).Mod(want, p)
		if got := p256MontToBig(&cm); got.Cmp(want) != 0 {
			t.Fatalf("neg: got %x want %x", got, want)
		}
		p256MontExp(&cm, &am, p256C1)
		want.Exp(a, p256C1, p)
		if got := p256MontToBig(&cm); got.Cmp(want) != 0 {
			t.Fatalf("exp: got %x want %x", got, want)
		}
		// canonical form must be fully reduced
		if raw := toBig(&cm); raw.Cmp(p) >= 0 {
			t.Fatalf("not reduced: %x", raw)
		}
	}
	// edge cases: 0, 1, p-1
	for _, v := range []*big.Int{big.NewInt(0), big.NewInt(1), new(big.Int).Sub(p, big.NewInt(1))} {
		var vm, cm p256Mont
		p256MontEncode(&vm, v)
		if got := p256MontToBig(&vm); got.Cmp(v) != 0 {
			t.Fatalf("edge roundtrip %x: got %x", v, got)
		}
		p256MontMul(&cm, &vm, &vm)
		want := new(big.Int).Mul(v, v)
		want.Mod(want, p)
		if got := p256MontToBig(&cm); got.Cmp(want) != 0 {
			t.Fatalf("edge sqr %x: got %x want %x", v, got, want)
		}
	}
}
