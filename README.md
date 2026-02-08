# TB2J -> Vampire (isotropic, Jiso-only) UCF Exporter

中文说明在前，English below.

---

## 中文

### 这是什么

一个小脚本，把 TB2J 的结果（`TB2J.pickle`）转换成 Vampire 可读的 unit cell file（UCF），并且**只输出各向同性交换** `Jiso`（Vampire `isotropic` 交换类型：每条相互作用只有一个 `Jij` 标量）。

脚本文件：`vampire-TB2J.py`

### 适用场景（为什么要用 Jiso-only）

很多人用 Vampire 做蒙特卡洛算 Tc，本质上是跑各向同性 Heisenberg：只需要一个标量 `Jij`。  
TB2J 默认导出的 Vampire 交换通常是 `tensorial`（3x3 矩阵，可能还包含各向异性/反对称项）。当你的目标只是 Tc 时，输入过于复杂往往会带来：

- 模型不一致（你想跑 isotropic，但输入包含更复杂项）
- 版本/单位/系数细节更容易踩坑，表现为磁化率/磁化强度曲线异常、Tc 很不稳定

所以这个脚本提供一个更“干净”的 isotropic 输入：只写 `Jiso`。

### 依赖

需要 Python 环境能 `import TB2J`，以及：

```bash
pip install TB2J ase numpy
```

### 输入/输出

输入：

- `--inpath TB2J_results`（目录）或 `--inpath TB2J_results/TB2J.pickle`（文件）

如何得到 `TB2J.pickle`：

- 你正常运行 TB2J（例如 `wann2J.py` / `siesta2J.py` / `abacus2J.py`）后，输出目录（默认 `TB2J_results` 或你指定的 `--output_path`）里会生成 `TB2J.pickle`。

输出：

- 一个 Vampire UCF 文件（默认：`TB2J_results/Vampire/vampire_isotropic.UCF`）
- 内容包含：
  - unit cell size / lattice vectors
  - magnetic atoms（按 `index_spin >= 0`）
  - interactions（`isotropic`，每条一列 `Jij`）
- 不写 biquadratic 段，不写 DMI，不写 3x3 张量

### 快速使用

全量导出（不做过滤）：

```bash
python vampire-TB2J.py --inpath TB2J_results --out TB2J_results/Vampire/vampire_isotropic.UCF
```

按阈值过滤（只保留强相互作用）：

```bash
# 按导出后的 Jij（Joule）过滤
python vampire-TB2J.py --inpath TB2J_results --out TB2J_results/Vampire/vampire_isotropic.UCF --cutoff 1e-23 --cutoff-unit joule

# 按 |2*Jiso|（eV）过滤
python vampire-TB2J.py --inpath TB2J_results --out TB2J_results/Vampire/vampire_isotropic.UCF --cutoff 1e-4 --cutoff-unit eV
```

### 参数说明

- `--inpath`：TB2J 结果目录或 `TB2J.pickle` 路径（默认 `TB2J_results`）
- `--out`：输出 UCF 路径（默认 `TB2J_results/Vampire/vampire_isotropic.UCF`）
- `--cutoff`：可选阈值，过滤弱相互作用（默认不启用）
- `--cutoff-unit`：`joule` 或 `eV`
  - `joule`：按最终写入的 `Jij`（Joule）过滤
  - `eV`：按 `|2*Jiso|`（eV）过滤

### 单位与系数约定（重要）

Vampire 手册中 `Jij` 的单位是 **Joules**。TB2J 的 `exchange_Jdict` 在实践中通常按 **eV** 使用。  
本脚本写入时使用：

```
Jij (Joule) = 2 * Jiso (eV) / ase.units.J
```

- `/ ase.units.J`：eV -> Joule
- `* 2`：沿用 TB2J 现有 Vampire 输出的惯例（与双计数相关）

如果你有自己明确的 convention（比如不想乘 2），可以在 `vampire-TB2J.py` 的 `_prepare_isotropic_interactions()` 里改这一行。

### 其他注意事项

- 检测到**非正交晶格**会给 warning（不强制退出）。
- 默认**不做 (i,j,R)/(j,i,-R) 去重**：保持和 TB2J `exchange_Jdict` 一致，避免额外假设。

---

## English

### What this is

A small helper script that converts TB2J results (`TB2J.pickle`) into a Vampire unit cell file (UCF) using **isotropic exchange only** (`Jiso`).  
It writes Vampire interactions as `isotropic` with a single scalar `Jij` per interaction.

Script: `vampire-TB2J.py`

### When to use (why Jiso-only)

For Monte Carlo Curie temperature (Tc) runs, many workflows want an **isotropic Heisenberg** model (scalar `Jij`).  
TB2J’s default Vampire export is typically `tensorial` (3x3) and may include additional terms. If your goal is Tc with an isotropic model, the extra complexity can make results unstable or trigger convention/unit pitfalls.

This script exports a minimal isotropic UCF (Jiso only).

### Requirements

You need TB2J importable in Python, plus:

```bash
pip install TB2J ase numpy
```

### Input / Output

Input:

- `--inpath TB2J_results` (directory) or `--inpath TB2J_results/TB2J.pickle` (file)

How to get `TB2J.pickle`:

- Run TB2J as usual (e.g. `wann2J.py` / `siesta2J.py` / `abacus2J.py`). The output directory (default `TB2J_results` or your `--output_path`) will contain `TB2J.pickle`.

Output:

- A Vampire UCF file (default: `TB2J_results/Vampire/vampire_isotropic.UCF`)
- Includes unit cell, magnetic atoms (`index_spin >= 0`), and `isotropic` interactions only

### Usage

Export everything:

```bash
python vampire-TB2J.py --inpath TB2J_results --out TB2J_results/Vampire/vampire_isotropic.UCF
```

Filter by cutoff:

```bash
# Cutoff on exported Jij (Joules)
python vampire-TB2J.py --inpath TB2J_results --out TB2J_results/Vampire/vampire_isotropic.UCF --cutoff 1e-23 --cutoff-unit joule

# Cutoff on |2*Jiso| (eV)
python vampire-TB2J.py --inpath TB2J_results --out TB2J_results/Vampire/vampire_isotropic.UCF --cutoff 1e-4 --cutoff-unit eV
```

### Convention (important)

Vampire expects `Jij` in **Joules**. TB2J `exchange_Jdict` values are typically treated as **eV** in practice.  
This script writes:

```
Jij (Joule) = 2 * Jiso (eV) / ase.units.J
```

The `*2` follows TB2J’s Vampire convention (related to double counting). If you use a different convention, adjust the conversion in `_prepare_isotropic_interactions()`.
