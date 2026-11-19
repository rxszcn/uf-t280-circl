package group_test

import (
	"encoding/json"
	"fmt"
	"path/filepath"
	"sync"
	"testing"

	"github.com/cloudflare/circl/group"
	"github.com/cloudflare/circl/internal/test"
)

func TestHashToElement(t *testing.T) {
	fileNames, err := filepath.Glob("./testdata/P*.json.gz")
	if err != nil {
		t.Fatal(err)
	}

	for _, fileName := range fileNames {
		input, err := test.ReadGzip(fileName)
		if err != nil {
			t.Fatal(err)
		}

		var v vectorSuite
		err = json.Unmarshal(input, &v)
		if err != nil {
			t.Fatal(err)
		}

		t.Run(v.Ciphersuite, v.testHashing)
	}
}

func (vs *vectorSuite) testHashing(t *testing.T) {
	var G group.Group
	switch vs.Ciphersuite[0:4] {
	case "P256":
		G = group.P256
	case "P384":
		G = group.P384
	case "P521":
		G = group.P521
	default:
		t.Fatal("non supported suite")
	}

	hashFunc := G.HashToElement
	if !vs.RandomOracle {
		hashFunc = G.HashToElementNonUniform
	}

	want := G.NewElement()
	for i, v := range vs.Vectors {
		got := hashFunc([]byte(v.Msg), []byte(vs.Dst))
		err := want.UnmarshalBinary(v.P.toBytes())
		if err != nil {
			t.Fatal(err)
		}

		if !got.IsEqual(want) {
			test.ReportError(t, got, want, i)
		}
	}
}

type vectorSuite struct {
	L           string `json:"L"`
	Z           string `json:"Z"`
	Ciphersuite string `json:"ciphersuite"`
	Curve       string `json:"curve"`
	Dst         string `json:"dst"`
	Expand      string `json:"expand"`
	Field       struct {
		M string `json:"m"`
		P string `json:"p"`
	} `json:"field"`
	Hash string `json:"hash"`
	K    string `json:"k"`
	Map  struct {
		Name string `json:"name"`
	} `json:"map"`
	RandomOracle bool     `json:"randomOracle"`
	Vectors      []vector `json:"vectors"`
}

type point struct {
	X test.HexBytes `json:"x"`
	Y test.HexBytes `json:"y"`
}

func (p point) toBytes() []byte {
	return append(append([]byte{0x04}, p.X...), p.Y...)
}

type vector struct {
	P   point    `json:"P"`
	Q0  point    `json:"Q0,omitempty"`
	Q1  point    `json:"Q1,omitempty"`
	Q   point    `json:"Q,omitempty"`
	Msg string   `json:"msg"`
	U   []string `json:"u"`
}

func BenchmarkHash(b *testing.B) {
	for _, g := range allGroups {
		name := g.(fmt.Stringer).String()
		b.Run(name+"/HashToElement", func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				g.HashToElement(nil, nil)
			}
		})
		b.Run(name+"/HashToElementNonUniform", func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				g.HashToElementNonUniform(nil, nil)
			}
		})
		b.Run(name+"/HashToScalar", func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				g.HashToScalar(nil, nil)
			}
		})
	}
}

// Allocation regression cases for the P-256 hash-to-curve hot path.
//
// Baseline before the optimization (measured with go test -bench=. -benchmem
// on this tree): HashToElement on P-256 cost 246 allocs/op. The requirement is
// to get below one quarter of the baseline, i.e. under 61.5 allocs/op, with
// the batch-aggregate reading converging in the same proportion.
const maxHashToElementAllocsP256 = 61.5

const hashAllocDST = "QUUX-V01-CS02-with-P256_XMD:SHA-256_SSWU_RO_"

// Single-call allocation threshold assertion.
func TestHashToElementAllocsP256(t *testing.T) {
	msg := []byte("allocation-threshold-input")
	allocs := testing.AllocsPerRun(200, func() {
		group.P256.HashToElement(msg, []byte(hashAllocDST))
	})
	if allocs >= maxHashToElementAllocsP256 {
		t.Fatalf("HashToElement allocs/op = %v, want < %v (baseline 246)", allocs, maxHashToElementAllocsP256)
	}
}

// Batch-aggregate allocation convergence assertion: a 1024-input batch must
// stay below 1024 times the single-call threshold.
func TestHashToElementBatchAllocsP256(t *testing.T) {
	const batchSize = 1024
	inputs := make([][]byte, batchSize)
	for i := range inputs {
		inputs[i] = []byte(fmt.Sprintf("batch-input-%d", i))
	}
	dst := []byte(hashAllocDST)
	allocs := testing.AllocsPerRun(10, func() {
		for i := 0; i < batchSize; i++ {
			group.P256.HashToElement(inputs[i], dst)
		}
	})
	if limit := batchSize * maxHashToElementAllocsP256; allocs >= limit {
		t.Fatalf("batch of %d allocs = %v, want < %v (per-input %v >= %v)",
			batchSize, allocs, limit, allocs/batchSize, maxHashToElementAllocsP256)
	}
}

// Concurrent reuse assertion: many goroutines hammering HashToElement through
// the shared workspace pool must each get results identical to sequential
// computation (no crosstalk). Run with -race for the full guarantee.
func TestHashToElementConcurrentReuse(t *testing.T) {
	const goroutines = 8
	const perGoroutine = 32

	dst := []byte(hashAllocDST)
	msgs := make([][]byte, goroutines*perGoroutine)
	want := make([]group.Element, len(msgs))
	for i := range msgs {
		msgs[i] = []byte(fmt.Sprintf("concurrent-input-%d", i))
		want[i] = group.P256.HashToElement(msgs[i], dst)
	}

	got := make([]group.Element, len(msgs))
	var wg sync.WaitGroup
	wg.Add(goroutines)
	for g := 0; g < goroutines; g++ {
		go func(g int) {
			defer wg.Done()
			for j := 0; j < perGoroutine; j++ {
				idx := g*perGoroutine + j
				got[idx] = group.P256.HashToElement(msgs[idx], dst)
			}
		}(g)
	}
	wg.Wait()

	for i := range got {
		if !got[i].IsEqual(want[i]) {
			t.Fatalf("input %q: concurrent result differs from sequential", msgs[i])
		}
	}
}

// Allocation benchmarks for the two constructed inputs: a single P-256
// hash-to-curve call, and the aggregate over a 1024-input batch.
func BenchmarkHashAllocsP256(b *testing.B) {
	dst := []byte(hashAllocDST)
	b.Run("Single", func(b *testing.B) {
		msg := []byte("allocation-benchmark-input")
		b.ReportAllocs()
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			group.P256.HashToElement(msg, dst)
		}
	})
	b.Run("Batch1024", func(b *testing.B) {
		inputs := make([][]byte, 1024)
		for i := range inputs {
			inputs[i] = []byte(fmt.Sprintf("batch-input-%d", i))
		}
		b.ReportAllocs()
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			for j := range inputs {
				group.P256.HashToElement(inputs[j], dst)
			}
		}
	})
}
