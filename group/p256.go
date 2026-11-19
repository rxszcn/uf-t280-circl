package group

import (
	"crypto/elliptic"
	"math/big"
	"math/bits"
)

// This file implements the Simplified SWU map (for curves with p = 3 mod 4)
// specialized to NIST P-256 using Montgomery arithmetic on 64-bit limbs.
// All field operations work on stack-allocated [4]uint64 values, so the
// map performs zero heap allocations, unlike the generic big.Int path.
//
// P-256 has p' = -p^-1 mod 2^64 = 1, which keeps the Montgomery reduction
// (CIOS-style, T < p*R on input) short. Correctness is anchored by the
// hash-to-curve test vectors in testdata/ (P256_XMD-SHA-256_SSWU_{RO,NU}_),
// which exercise this code path through HashToElement and
// HashToElementNonUniform.

// p256Mont is an element of GF(p256) in Montgomery form: it holds x*R mod p
// with R = 2^256, limbs in little-endian order, always fully reduced.
type p256Mont [4]uint64

// p256Modulus is the P-256 field order
// p = 2^256 - 2^224 + 2^192 + 2^96 - 1 in little-endian 64-bit limbs.
var p256Modulus = p256Mont{
	0xffffffffffffffff,
	0x00000000ffffffff,
	0x0000000000000000,
	0xffffffff00000001,
}

var (
	// p256R2 is 2^512 mod p, used to convert into the Montgomery domain.
	p256R2 p256Mont
	// p256One is 1 in Montgomery form, i.e. R mod p.
	p256One p256Mont
	// Curve and map constants in Montgomery form.
	p256Z  p256Mont // -10
	p256C2 p256Mont // 0x78bc71a02d89ec07214623f6d0f955072c7cc05604a5a6e23ffbf67115fa5301
	p256A  p256Mont // -3
	p256B  p256Mont // curve coefficient B
	// Exponents (plain integers, only their bits are read).
	p256C1    *big.Int // (p - 3) / 4
	p256PSub2 *big.Int // p - 2, for inversion via Fermat's little theorem
)

func init() {
	params := elliptic.P256().Params()
	p := params.P
	R := new(big.Int).Lsh(big.NewInt(1), 256)
	R.Mod(R, p)

	toMont := func(x *big.Int) (out p256Mont) {
		v := new(big.Int).Mul(x, R)
		v.Mod(v, p)
		p256MontFromBig(&out, v)
		return out
	}

	p256MontFromBig(&p256R2, new(big.Int).Exp(big.NewInt(2), big.NewInt(512), p))
	p256One = toMont(big.NewInt(1))
	p256Z = toMont(big.NewInt(-10))
	c2, _ := new(big.Int).SetString("0x78bc71a02d89ec07214623f6d0f955072c7cc05604a5a6e23ffbf67115fa5301", 0)
	p256C2 = toMont(c2)
	p256A = toMont(big.NewInt(-3))
	p256B = toMont(params.B)

	p256C1 = new(big.Int).Sub(p, big.NewInt(3))
	p256C1.Rsh(p256C1, 2)
	p256PSub2 = new(big.Int).Sub(p, big.NewInt(2))
}

// p256MontFromBig sets z to the limbs of v, requiring 0 <= v < p.
func p256MontFromBig(z *p256Mont, v *big.Int) {
	words := v.Bits()
	for i := range z {
		if i < len(words) {
			z[i] = uint64(words[i])
		} else {
			z[i] = 0
		}
	}
}

// p256MontToBig returns a fresh big.Int holding the canonical value of z.
func p256MontToBig(z *p256Mont) *big.Int {
	var plain p256Mont
	p256MontMul(&plain, z, &p256Mont{1, 0, 0, 0}) // decode: z * 1 * R^-1
	return new(big.Int).SetBits([]big.Word{
		big.Word(plain[0]), big.Word(plain[1]),
		big.Word(plain[2]), big.Word(plain[3]),
	})
}

// p256MontEncode sets z to the Montgomery form of v, requiring 0 <= v < p.
func p256MontEncode(z *p256Mont, v *big.Int) {
	var plain p256Mont
	p256MontFromBig(&plain, v)
	p256MontMul(z, &plain, &p256R2) // v * R^2 * R^-1 = v * R
}

// p256MontMul sets c = a * b * R^-1 mod p. Safe for c to alias a and/or b.
//
// Note: bits.Add64/bits.Sub64 require their carry/borrow input to be 0 or 1,
// so full-word carries are always passed as regular addends, never as the
// carry argument.
func p256MontMul(c, a, b *p256Mont) {
	// t = a * b, a 512-bit product in 8 limbs (t[8] stays 0).
	var t [9]uint64
	for i := 0; i < 4; i++ {
		var carry uint64
		for j := 0; j < 4; j++ {
			hi, lo := bits.Mul64(a[i], b[j])
			sum, c1 := bits.Add64(t[i+j], lo, 0)
			sum, c2 := bits.Add64(sum, carry, 0)
			t[i+j] = sum
			// hi <= 2^64-2, and c1 = c2 = 1 with hi = 2^64-2 would require
			// lo = 1 (m = b = 2^64-1) and carry >= 2^64, which is impossible,
			// so this addition cannot overflow.
			carry = hi + c1 + c2
		}
		t[i+4] = carry
	}

	// Montgomery reduction: for each limb, add t[i] * p shifted by i limbs,
	// which cancels t[i] since p * 1 = -1 (mod 2^64) and p' = 1.
	for i := 0; i < 4; i++ {
		m := t[i]
		var carry uint64
		for j := 0; j < 4; j++ {
			hi, lo := bits.Mul64(m, p256Modulus[j])
			sum, c1 := bits.Add64(t[i+j], lo, 0)
			sum, c2 := bits.Add64(sum, carry, 0)
			t[i+j] = sum
			carry = hi + c1 + c2
		}
		for k := i + 4; carry != 0; k++ {
			sum, c1 := bits.Add64(t[k], carry, 0)
			t[k] = sum
			carry = c1
		}
	}

	// The result (t[8], t[7..4]) is < 2p; subtract p conditionally.
	d0, b1 := bits.Sub64(t[4], p256Modulus[0], 0)
	d1, b1 := bits.Sub64(t[5], p256Modulus[1], b1)
	d2, b1 := bits.Sub64(t[6], p256Modulus[2], b1)
	d3, b1 := bits.Sub64(t[7], p256Modulus[3], b1)
	// If t[8] < b1 the subtraction borrowed out, meaning the value was
	// already < p: keep the original limbs, otherwise keep the difference.
	_, bOut := bits.Sub64(t[8], b1, 0)
	mask := uint64(0) - bOut
	c[0] = (t[4] & mask) | (d0 &^ mask)
	c[1] = (t[5] & mask) | (d1 &^ mask)
	c[2] = (t[6] & mask) | (d2 &^ mask)
	c[3] = (t[7] & mask) | (d3 &^ mask)
}

// p256MontAdd sets c = a + b mod p.
func p256MontAdd(c, a, b *p256Mont) {
	s0, c1 := bits.Add64(a[0], b[0], 0)
	s1, c1 := bits.Add64(a[1], b[1], c1)
	s2, c1 := bits.Add64(a[2], b[2], c1)
	s3, c1 := bits.Add64(a[3], b[3], c1)

	d0, b1 := bits.Sub64(s0, p256Modulus[0], 0)
	d1, b1 := bits.Sub64(s1, p256Modulus[1], b1)
	d2, b1 := bits.Sub64(s2, p256Modulus[2], b1)
	d3, b1 := bits.Sub64(s3, p256Modulus[3], b1)

	// a + b < 2p: use the difference unless it borrowed out with no carry.
	useD := c1 | (b1 ^ 1)
	mask := uint64(0) - useD
	c[0] = (d0 & mask) | (s0 &^ mask)
	c[1] = (d1 & mask) | (s1 &^ mask)
	c[2] = (d2 & mask) | (s2 &^ mask)
	c[3] = (d3 & mask) | (s3 &^ mask)
}

// p256MontSub sets c = a - b mod p.
func p256MontSub(c, a, b *p256Mont) {
	d0, b1 := bits.Sub64(a[0], b[0], 0)
	d1, b1 := bits.Sub64(a[1], b[1], b1)
	d2, b1 := bits.Sub64(a[2], b[2], b1)
	d3, b1 := bits.Sub64(a[3], b[3], b1)

	// On borrow, add the modulus back.
	mask := uint64(0) - b1
	s0, c1 := bits.Add64(d0, p256Modulus[0]&mask, 0)
	s1, c1 := bits.Add64(d1, p256Modulus[1]&mask, c1)
	s2, c1 := bits.Add64(d2, p256Modulus[2]&mask, c1)
	s3, _ := bits.Add64(d3, p256Modulus[3]&mask, c1)
	c[0], c[1], c[2], c[3] = s0, s1, s2, s3
}

// p256MontNeg sets c = -a mod p.
func p256MontNeg(c, a *p256Mont) {
	var zero p256Mont
	p256MontSub(c, &zero, a)
}

// p256MontSelect sets c = b if k is true, and c = a otherwise.
func p256MontSelect(c, a, b *p256Mont, k bool) {
	if k {
		*c = *b
	} else {
		*c = *a
	}
}

// p256MontIsZero reports whether a is zero.
func p256MontIsZero(a *p256Mont) bool {
	return (a[0] | a[1] | a[2] | a[3]) == 0
}

// p256MontSgn0 returns the least significant bit of the canonical value of a.
func p256MontSgn0(a *p256Mont) uint {
	var plain p256Mont
	p256MontMul(&plain, a, &p256Mont{1, 0, 0, 0})
	return uint(plain[0] & 1)
}

// p256MontExp sets z = x^e mod p for e > 0 by square-and-multiply.
func p256MontExp(z, x *p256Mont, e *big.Int) {
	r := p256One
	b := *x
	for i := 0; i < e.BitLen(); i++ {
		if e.Bit(i) == 1 {
			p256MontMul(&r, &r, &b)
		}
		p256MontMul(&b, &b, &b)
	}
	*z = r
}

// p256MapToCurve runs the Simplified SWU map on u (0 <= u < p) and returns
// fresh big.Int values holding the affine coordinates of the mapped point.
func p256MapToCurve(u *big.Int) (xOut, yOut *big.Int) {
	var um p256Mont
	p256MontEncode(&um, u)

	var tv1, tv2, tv3, tv4 p256Mont
	var xn, xd, x1n, x2n p256Mont
	var gx1, gxd, y1, y2 p256Mont
	var x, y p256Mont

	sqr := func(c, a *p256Mont) { p256MontMul(c, a, a) }

	sqr(&tv1, &um)                              // 1.  tv1 = u^2
	p256MontMul(&tv3, &p256Z, &tv1)             // 2.  tv3 = Z * tv1
	sqr(&tv2, &tv3)                             // 3.  tv2 = tv3^2
	p256MontAdd(&xd, &tv2, &tv3)                // 4.   xd = tv2 + tv3
	p256MontAdd(&x1n, &xd, &p256One)            // 5.  x1n = xd + 1
	p256MontMul(&x1n, &x1n, &p256B)             // 6.  x1n = x1n * B
	p256MontNeg(&tv4, &p256A)                   //
	p256MontMul(&xd, &tv4, &xd)                 // 7.   xd = -A * xd
	e1 := p256MontIsZero(&xd)                   // 8.   e1 = xd == 0
	p256MontMul(&tv4, &p256Z, &p256A)           //
	p256MontSelect(&xd, &xd, &tv4, e1)          // 9.   xd = CMOV(xd, Z * A, e1)
	sqr(&tv2, &xd)                              // 10. tv2 = xd^2
	p256MontMul(&gxd, &tv2, &xd)                // 11. gxd = tv2 * xd
	p256MontMul(&tv2, &p256A, &tv2)             // 12. tv2 = A * tv2
	sqr(&gx1, &x1n)                             // 13. gx1 = x1n^2
	p256MontAdd(&gx1, &gx1, &tv2)               // 14. gx1 = gx1 + tv2
	p256MontMul(&gx1, &gx1, &x1n)               // 15. gx1 = gx1 * x1n
	p256MontMul(&tv2, &p256B, &gxd)             // 16. tv2 = B * gxd
	p256MontAdd(&gx1, &gx1, &tv2)               // 17. gx1 = gx1 + tv2
	sqr(&tv4, &gxd)                             // 18. tv4 = gxd^2
	p256MontMul(&tv2, &gx1, &gxd)               // 19. tv2 = gx1 * gxd
	p256MontMul(&tv4, &tv4, &tv2)               // 20. tv4 = tv4 * tv2
	p256MontExp(&y1, &tv4, p256C1)              // 21.  y1 = tv4^c1
	p256MontMul(&y1, &y1, &tv2)                 // 22.  y1 = y1 * tv2
	p256MontMul(&x2n, &tv3, &x1n)               // 23. x2n = tv3 * x1n
	p256MontMul(&y2, &y1, &p256C2)              // 24.  y2 = y1 * c2
	p256MontMul(&y2, &y2, &tv1)                 // 25.  y2 = y2 * tv1
	p256MontMul(&y2, &y2, &um)                  // 26.  y2 = y2 * u
	sqr(&tv2, &y1)                              // 27. tv2 = y1^2
	p256MontMul(&tv2, &tv2, &gxd)               // 28. tv2 = tv2 * gxd
	e2 := tv2 == gx1                            // 29.  e2 = tv2 == gx1
	p256MontSelect(&xn, &x2n, &x1n, e2)         // 30.  xn = CMOV(x2n, x1n, e2)
	p256MontSelect(&y, &y2, &y1, e2)            // 31.   y = CMOV(y2, y1, e2)
	e3 := p256MontSgn0(&um) == p256MontSgn0(&y) // 32.  e3 = sgn0(u) == sgn0(y)
	p256MontNeg(&tv1, &y)                       //
	p256MontSelect(&y, &tv1, &y, e3)            // 33.   y = CMOV(-y, y, e3)
	// xd is always non-zero here (step 9 forces xd = Z * A when xd == 0),
	// so Fermat's little theorem yields the same inverse as ModInverse.
	p256MontExp(&tv1, &xd, p256PSub2) // tv1 = xd^(p-2) = xd^-1 (mod p)
	p256MontMul(&x, &xn, &tv1)        // 34. return (xn, xd, y, 1)

	return p256MontToBig(&x), p256MontToBig(&y)
}
