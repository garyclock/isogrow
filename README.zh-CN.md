# IsoGrow

**本地可复现的 SwissBioisostere 替代实现，外加一条从片段替换到分子生长、构象对齐、对接与排序的完整流水线。**

English: [README.md](README.md)

SwissBioisostere 是个网页服务：提交片段，返回替换片段。IsoGrow 用 ChEMBL 里的
**MMP（matched molecular pairs）** 在本地重建这套检索，输出格式与官方一致，然后继续往下走：
生成变体 → 3D 构象 → 对齐参考构象 → 对接 → 排序。

```text
片段 / 母核
   │ A. 等排体检索（本地 mmpdb 库，或直接读官方 SwissBioisostere CSV）
   ▼
替换片段集
   │ B. CReM mutate（上下文约束的变体）
   ▼
变体集
   │ C. ETKDG+MMFF   D. MCS 对齐   E. Meeko → PDBQT
   ▼
对齐后的 3D 构象
   │ F. 对接（默认 Vina，可插 Uni-Dock / 外部引擎）
   ▼
G. 排序表 + HTML 报告
```

## 模块

| 层 | 模块 | 作用 |
|---|---|---|
| 检索 | `isogrow.bioisostere` | mmpdb 建库与查询，输出 `Frequency`、`ΔlogP`、`ΔtPSA`、`ΔMW`、`#Better/#Similar/#Worse`、success score |
| 生长 | `isogrow.pipeline` | 上面 A→G 七阶段，纯 Python（RDKit + Meeko + Vina），无需 py2 / `sbatch` / 编译型对齐工具 |
| 数据 | `scripts/data/` | ChEMBL chemreps、SQLite、活性样本、CReM 库、Vina 二进制的下载与准备 |
| 命令行 | `isogrow` | `doctor`、`fetch`、`build`、`search`、`run`、`benchmark` |

## 安装

```bash
git clone https://github.com/garyclock/isogrow
cd isogrow
python -m venv .venv && source .venv/bin/activate      # Windows: .venv\Scripts\activate
pip install -e ".[all]"
isogrow doctor
```

> **一个必踩的坑**：PyPI 上名为 `mmpdb` 的包能装上但**没有模块**。可用的发行在
> [rdkit/mmpdb](https://github.com/rdkit/mmpdb) 仓库里——模块名 `mmpdblib`、命令行 `mmpdb`。
> 所以依赖写成 `mmpdb @ git+https://github.com/rdkit/mmpdb.git`，而不是版本号。
> 另外 `meeko` 必须带上 `gemmi` 才能 import。

## 30 秒上手（不需下载任何数据）

仓库自带 15 个分子的 MMP 系列，足以建出真实的 2 MB 库：

```bash
isogrow build --smi tests/data/tiny.smi --out data/mmpdb/tiny.mmpdb
isogrow search --db data/mmpdb/tiny.mmpdb --query '[*:1]c1cc2ccccc2o1'
```

## 用真实数据检索

```bash
python scripts/data/download_chemreps.py                 # 279 MB，全 ChEMBL 结构
isogrow build --chemreps data/chembl/chembl_37_chemreps.txt.gz \
              --out data/mmpdb/chembl50k --limit 50000 --jobs 6
isogrow search --db data/mmpdb/chembl50k.mmpdb --query '[*:1]c1cc2ccccc2o1' --topn 25
```

默认 `--num-cuts 3`，与 SwissBioisostere 一致（1 刀=侧链、2 刀=linker、3 刀=骨架）。
多附件点片段（如 `[*:1]C(=O)Oc1ccc([*:2])cc1`）必须建 2/3 刀库才查得到。

**体积警告**：3 刀库膨胀很快——5 万分子约 **6.5 GB**，15 万约 **20 GB+**，建索引时内存峰值约 9 GB。

### 活性维度（#Better / #Similar / #Worse）

```bash
python scripts/data/fetch_activities.py --target CHEMBL203 --max-rows 5000
isogrow build --chemreps data/chembl/chembl_37_chemreps.txt.gz \
              --activity-table data/chembl_sample/molecules_CHEMBL203.csv \
              --out data/mmpdb/egfr
isogrow search --db data/mmpdb/egfr.mmpdb --query '[*:1]F' --activity
```

只有需要 ChEMBL 自己的 logP 模型或 assay/靶点/PubMed 元数据时，才需要下 5.4 GB 的 SQLite 全量库。

## 完整流水线

```bash
cp config/example_local.json config/my_run.json     # 改路径与查询片段
isogrow run --config config/my_run.json
```

`config/example_official_csv.json` 是另一种输入方式：直接吃官方 CSV 导出，后续阶段完全相同。
产物按阶段落盘：`1_bioisosteres.csv`、`2_3d/`、`3_aligned/`、`4_pdbqt/`、`5_dock/`、`6_ranked.csv|html`。

变体预算按**每个候选**分配（`variants_per_candidate`），不是全局截断——否则排序会塌到频率最高的
那一个替换上。

## 与官方服务的对标

查询 `[*:1]c1cc2ccccc2o1`（苯并呋喃），5 万分子库：

| 指标 | 值 | 含义 |
|---|---|---|
| 相对官方的 precision | **72.3 %** | 本地候选里 72 % 确实是官方结果 |
| recall@50 | **64 %** | 官方前 50 基本复现 |
| 覆盖全部 1663 个候选 | 6.1 % | 受库规模限制（5 万≈ChEMBL 的 1/48） |
| Δ 误差（tPSA）中位数 / 均值 | **0.00 / 0.13** | 基本精确；3 条官方自相矛盾的行把均值拉高 |
| Δ 误差（MW）中位数 / 均值 | **0.00 / 1.81** | 同样那 3 行（例如官方给出 ΔMW −78.7，而实际只变 −4.6） |
| Δ 误差（logP）中位数 / 均值 | **0.50 / 0.51** | 系统性：官方用 ACD/Labs logP，公开 ChEMBL 里没有这个值 |

阶段 B 与原 notebook 的历史产物逐分子比对：CReM 变体集 **100 % 重合**（50/50、53/53），
RDKit 对齐 RMSD 0.06 Å。细节见 [docs/benchmark.md](docs/benchmark.md)。

## 已知局限（引用数字前请先读）

* **覆盖率随库规模线性变化**：5 万子集只能复现 240 万分子资源的约 6%；生产用请建更大库。
* **ΔlogP 比官方低 0.50**（sd 0.29，r = 0.92）。官方用 ACD/Labs logP，而 ChEMBL 自 24 版起
  `alogp` 就是 RDKit 算的、ACD 值不在公开 dump 里——我们分别建了两种库验证过，**这个差距无法用公开数据消除**。
  ΔtPSA / ΔMW 在 102 个候选里 99 个精确吻合（中位误差 0.00）——**看中位数**，少数官方异常行会把均值毁掉；
  `report.json` 两者都给。加 +0.5 常数校准可把 logP 误差从 0.51 降到 0.22，但**默认不启用**（见 docs/benchmark.md）。
* **success score 是代理公式**（净提升率 × log(对数)），官方还加权靶点/靶点类别/Murcko 骨架多样性。
* **KCOMBU 换成 RDKit MCS 对齐**、py2 AutoDockTools 换成 Meeko：分子身份不变（逐条核过重原子数），
  构象几何可能与原 KCOMBU 输出略有差别。
* **对接分数依赖引擎**：默认 Vina，Uni-Dock / 外部引擎可通过 `isogrow.pipeline.dock.BACKENDS` 接入。

## 数据策略

仓库**不含**任何 ChEMBL 数据、构建产物或官方 SwissBioisostere 导出。`tests/data/` 只有约 113 KB
夹具（15 分子 SMILES 系列、一个受体、一个参考配体）。许可与署名见 [THIRD_PARTY.md](THIRD_PARTY.md)。

## 许可

代码 MIT。使用该方法时请一并引用 SwissBioisostere 原文与 ChEMBL，见 [CITATION.cff](CITATION.cff)。
