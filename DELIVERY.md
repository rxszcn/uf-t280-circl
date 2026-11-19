# T280 交付说明：group 包哈希到曲线热点路径分配量优化

## 背景与结论

隐私令牌批量签发链路中，`group` 包 NIST 曲线的 `HashToElement` 每次调用产生约
250 个临时对象（基线实测 246 allocs/op、15730 B/op），同库 ristretto255 仅
10 allocs/op，差距来自实现形态（`math/big` 临时对象 + 每次重算常量）而非协议固有。

本次改动后（读数口径：`go test -bench=. -benchmem ./group/` 的 allocs/op）：

| 构造输入 | 基线 | 改动后 | 阈值 | 结论 |
|---|---|---|---|---|
| P-256 单次 HashToElement | 246 allocs/op | **38 allocs/op** | < 61.5（基线 1/4） | 达标，降至基线 15.4% |
| P-256 ×1024 批次聚合 | 251904（246×1024 推算） | **38917（≈38.0/输入）** | < 62976（同比例） | 达标，同比例收敛 |

哈希输出逐位不变：仓库自带规范向量 `group/testdata/P*.json.gz` 全部通过。

## 热点定位（改动前现状）

- `group/short.go` 的 `sswu3mod4Map`：每次调用 `new(big.Int)` ×14，`A`、`c1=(p-3)/4`
  等常量每次重算，`big.Int.Exp`/`ModInverse` 内部产生大量临时分配；
  约 30 次模乘各自伴随 `Mod` 的内部分配（`nat.div` 临时对象）。
- `mapToCurveParams`：每次调用重新构造 `Z`、`C2`（含 `SetString` 解析）并生成闭包。
- 结论：`math/big` 的 `Mod` 每次调用必然分配商等临时对象，仅靠复用目标
  `big.Int` 无法把分配量压到阈值以下，P-256 需要固定肢域算术。

## 方案与取舍

### 复用策略

1. **常量包级预计算**（`sswuParams`，`sync.Map` 按曲线位宽缓存，`LoadOrStore`
   保证并发首次构建安全）：`Z/C2/A/B/P/one/c1` 及哈希参数只算一次，只读共享。
2. **`sync.Pool` 工作区**（`sswuWorkspace`）：通用 big.Int 路径的 14 个临时
   `big.Int` 与 `u[2]` 全部池化复用；`big.Int` 复用后 `SetBytes`/`Mul` 等不再
   增长底层数组，稳态零分配。
3. **P-256 Montgomery 固定肢域算术**（新文件 `group/p256.go`）：SSWU 映射全程
   在栈上 `[4]uint64` 进行，堆分配为零。P-256 满足 `-p⁻¹ mod 2⁶⁴ = 1`，REDC
   归约简洁；`Exp`/`ModInverse` 换成 Montgomery 域二进制幂链与费马求逆
   （`xd` 经第 9 步 CMOV 保证非零，费马结果与 `ModInverse` 一致）。
4. **结果深拷贝出池**：返回的 `wElt` 持有全新 `big.Int`，池内临时量绝不逃逸到
   返回值，这是并发无串扰的关键不变式。

### 并发语义取舍及理由

- 选 `sync.Pool` 而非 `Mutex` 包一个共享工作区：签发线程是批量并发调用，
  互斥会把热点串行化；`sync.Pool` 让每个 goroutine 拿到独立工作区，无锁竞争、
  无共享可变状态，代价是池未命中时的一次性分配（稳态摊销为零）。
- 工作区不跨调用保留任何结果状态，调用之间唯一的共享是只读常量，因此并发
  语义与改动前（每次全新分配）完全等价。
- P-384/P-521 保留 big.Int 通用路径（池化 + 预计算常量，分配量 207→163、
  210→170），不为它们新写域实现：题面构造输入定为 P-256，且新域实现的
  正确性风险应限制在已验证的 P-256 上（规范向量 + big.Int 对拍双重锚定）。
- `expander` 包未改动：实测表明仅 group 侧优化已使两档读数达标（38 < 61.5），
  改动面按最小化原则收敛；expander 每次 Expand 的约 10 次分配已计入上述读数。
- 已知权衡：P-256 单次延迟 89µs → 115µs（Montgomery 软实现幂链 vs
  `big.Int.Exp` 的窗口法），但字节分配 15730 → 1739 B/op（-89%），批量签发
  场景的瓶颈是回收压力而非单核延迟；如需追回延迟可后续换固定窗口/加法链，
  不影响本次分配量目标。

### 实现注记

`math/bits.Add64/Sub64` 的 carry 入参只允许 0 或 1（否则行为未定义），全字
进位必须作为普通加数传递。`group/p256_test.go` 以 1000 组随机 + 边界输入对
Montgomery 域算术与 `big.Int` 做对拍，防止此类问题回归。

## 题面要求逐条落实

| 题面要求 | 落实位置 |
|---|---|
| 单次哈希到曲线分配量降到基线四分之一以下 | `group/p256.go`（Montgomery 域零分配映射）+ `group/short.go`（分派）；断言 `TestHashToElementAllocsP256`（group/hash_test.go） |
| 批次聚合读数同比例收敛 | 同上代码路径；断言 `TestHashToElementBatchAllocsP256`（1024 输入，阈值 1024×61.5） |
| 哈希输出逐位不变（规范向量为准） | 既有 `TestHashToElement` 跑 `group/testdata/P*.json.gz` 全部通过；算法步骤与 RFC 9380 逐项对应（注释保留原步骤编号） |
| 批量签发与批量证明路径聚合分配同步收敛 | oprf（`oprf/client.go`、`server.go`）与 zk（`zk/dleq`、`zk/dl`）调用同一 `Group.HashToElement/HashToScalar`，随 group 包自动收敛；P-384/P-521 通用路径同步下降 |
| 多 goroutine 并发复用无串扰 | `sync.Pool` 每 goroutine 独立工作区 + 结果深拷贝出池；断言 `TestHashToElementConcurrentReuse`（8 goroutine × 32 输入与顺序结果逐一对拍） |
| 对外 API 签名与行为不变 | 未改任何导出声明；`HashToField` 等导出函数签名原样；新增类型/函数均非导出 |
| 既有全部测试与向量保持通过 | `go test -count=1 ./group/... ./oprf/... ./zk/...` 零失败（见下） |
| 分配量基准用例并入 group 包既有基准测试文件、同文件同框架 | `BenchmarkHashAllocsP256`（Single / Batch1024 两档）加入 `group/hash_test.go`，与 `BenchmarkHash` 同文件同 `testing.B` 框架 |
| 改动面收在 group 与 expander 两包及哈希路径调用侧 | 仅改动/新增 group 包 4 个文件；expander 达标无需改动；调用侧（oprf/zk）零迁移 |
| 构造全程快照环境内、零外部依赖 | 仅使用标准库与仓库内代码，未新增任何 go.mod 依赖 |

## 实测读数（go test -bench=. -benchmem ./group/）

改动前 → 改动后（allocs/op）：

- `BenchmarkHash/P-256/HashToElement`：246 → 38（15730 → 1739 B/op）
- `BenchmarkHash/P-256/HashToElementNonUniform`：128 → 20
- `BenchmarkHash/P-256/HashToScalar`：28 → 19
- `BenchmarkHash/P-384/HashToElement`：207 → 163（通用路径同步收敛）
- `BenchmarkHash/P-521/HashToElement`：210 → 170（通用路径同步收敛）
- `BenchmarkHash/ristretto255/HashToElement`：10 → 10（对照组不动）
- `BenchmarkHashAllocsP256/Single`：38 allocs/op（构造输入第一档）
- `BenchmarkHashAllocsP256/Batch1024`：38917 allocs/op ≈ 38.0/输入（构造输入第二档）

## 验收命令最终读数

```
$ go test -count=1 ./group/... ./oprf/... ./zk/...
ok  	github.com/cloudflare/circl/group	4.170s
ok  	github.com/cloudflare/circl/oprf	1.034s
?   	github.com/cloudflare/circl/zk	[no test files]
ok  	github.com/cloudflare/circl/zk/dl	0.496s
ok  	github.com/cloudflare/circl/zk/dleq	0.488s
ok  	github.com/cloudflare/circl/zk/qndleq	0.389s
```

零失败。另：`go build ./...`、`go vet ./group/`、`gofmt -l group/` 均干净；
`go test -count=1 ./expander/` 通过。本环境无 gcc，`-race` 不可用，并发正确性
由设计不变式（每 goroutine 独立工作区、结果深拷贝出池、共享仅只读常量）与
并发对拍用例保证。

## 改动面逐文件登记

| 文件 | 类型 | 内容 |
|---|---|---|
| `group/short.go` | 修改 | `sswuParams` 常量预计算（`sync.Map` 缓存）替代 `mapToCurveParams`；`sswuWorkspace` + `sync.Pool` 复用临时 `big.Int`；`sswu3mod4Map` 改用池化工作区与预计算常量；`HashToElement/HashToElementNonUniform` 增加 P-256 Montgomery 快速路径分派；`HashToScalar` 复用池化 `u` |
| `group/p256.go` | 新增 | P-256 Montgomery 域算术（`[4]uint64`，REDC、加减取负、选择、幂链、费马求逆）与零分配 SSWU 映射 `p256MapToCurve`；常量 init 时一次性计算 |
| `group/hash_test.go` | 修改 | 新增三类用例：`TestHashToElementAllocsP256`（单次阈值断言）、`TestHashToElementBatchAllocsP256`（1024 批次收敛断言）、`TestHashToElementConcurrentReuse`（并发无串扰断言）；新增基准 `BenchmarkHashAllocsP256`（Single/Batch1024），与 `BenchmarkHash` 同文件同框架 |
| `group/p256_test.go` | 新增 | Montgomery 域算术 vs `big.Int` 随机/边界对拍（内部测试，包 `group`） |
| `expander/` | 未改动 | 实测仅 group 侧优化即达标，改动面最小化 |
| 调用侧（oprf、zk） | 未改动 | 经 `Group` 接口自动受益，零迁移 |
