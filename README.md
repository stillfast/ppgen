# pp-gen —— 赝势生成流水线

已实现三条命令：

| 命令 | 干什么 |
|---|---|
| **`ppgen prepare`** | 把 run 目录里的 `input.json` 翻译成可直接喂给 QE `ld1.x` 的输入文件（如 `U.in`）；扫描型工作流会一次生成**一批** |
| **`ppgen run`** | 跑 `ld1.x < U.in > U.out`（等价于参考脚本 `ppgen/temporary/run.sh`）。**两套流程**：不给 `-i` 就在**当前目录直接跑**现成的 `U.in`（不生成、不校验）；给 `-i input.json` 就**先生成、再运行**（扫描型逐点生成逐个跑） |
| **`ppgen plot`** | 按 `plot.json` 画图：波函数（`ld1.wfc`）或 AE/PS 对数导数（`ld1.dlog` / `ld1ps.dlog`） |
| **`ppgen report`** | 把一次扫描汇成 `report.md`（每个轨道最佳的 `r_c` + 该参数下各通道极点/零点的极值与是否重合），并画一张"凭什么选这个 `r_c`"的图 |

两种工作流形状：

| 工作流 | `prepare` 生成 | `run` 跑 |
|---|---|---|
| **`ld1.x`** | 一份：`<元素>.in`（全电子、生成赝势……都算，用哪套参数由 `parameters.input` 决定） | 一次 |
| `rc_scan` | 一批：`runs/<轨道>/rc_<值>/U.in`（5 轨道 × 11 点 = 55 份） | 逐个点跑，每个点在自己的目录里 |

`workflow` 只说"跑一次还是扫一批"，**具体跑什么由 `parameters` 里的 YAML 决定** ——
所以 AE（`iswitch=1`）、扫完挑出来的最佳 `r_c`（`rc_best.yml`，`iswitch=3`）用的是
**同一个** `ld1.x` 工作流，换 `parameters.input` 就行：

```json
{"workflow": "ld1.x", "parameters": {"input": "ae_radial"}, "qe_path": "..."}   // 全电子
{"workflow": "ld1.x", "parameters": {"input": "rc_best"},   "qe_path": "..."}   // 最佳 r_c 生成赝势
```

## 用法

```bash
# 单份输入（全电子径向波函数）：run -i 会先把 U.in 生成一遍再跑
ppgen prepare -i /home/liguozhou/oncv/uranium-ppgen/ae_radial/input.json
ppgen run     -i /home/liguozhou/oncv/uranium-ppgen/ae_radial/input.json
ppgen plot    -i /home/liguozhou/oncv/uranium-ppgen/ae_radial/plot.json

# 只生成、不跑（想看会生成什么）
ppgen prepare -i /home/liguozhou/oncv/uranium-ppgen/ae_radial/input.json --dry-run

# 扫描型（逐轨道扫 rcut）：一次生成 55 份输入并逐个跑
ppgen run -i /home/liguozhou/oncv/uranium-ppgen/rc_scan/input.json
ppgen run -i /home/liguozhou/oncv/uranium-ppgen/rc_scan/input.json --only 5F --rc 1.0   # 只跑一个点

# 对数导数参考半径扫描（看极点在哪个能量上最干净）
ppgen run    -i /home/liguozhou/oncv/uranium-ppgen/best/test_rlderiv/input.json
ppgen report -i /home/liguozhou/oncv/uranium-ppgen/best/test_rlderiv/output.json \
                /home/liguozhou/oncv/uranium-ppgen/best/test_rlderiv/plot_dlog.json

# 局域势芯半径扫描（找最佳的 rcloc）
ppgen run -i /home/liguozhou/oncv/uranium-ppgen/rcloc_scan/input.json
ppgen run -i /home/liguozhou/oncv/uranium-ppgen/rcloc_scan/input.json --rc 2.5   # 只跑一个点

# 能量扫描（给 occ=0 的额外投影算子找 ener）
ppgen run -i /home/liguozhou/oncv/uranium-ppgen/energy_scan/input.json
ppgen run -i /home/liguozhou/oncv/uranium-ppgen/energy_scan/input.json --only 5F-2 --ener=-0.6

# 手改过某个点的 U.in，只想原样重跑：进那个目录直接用（不生成、不覆盖）
cd /home/liguozhou/oncv/uranium-ppgen/rc_scan/runs/5F/rc_1p2 && ppgen run
```

`prepare` 的输出（产物落在 `input.json` 同目录，即 run 目录）：

```
workflow   : ld1.x
参数文件   : /home/liguozhou/oncv/pp-gen/src/pp_gen/parameters/input/ae_radial.yml
qe_path    : /home/liguozhou/oncv/work/qe/q-e-qe-7.4/bin（input.json 里是占位符，已补默认值，并写回 input.json）
ld1.x      : /home/liguozhou/oncv/work/qe/q-e-qe-7.4/bin/ld1.x
已生成     : /home/liguozhou/oncv/uranium-ppgen/ae_radial/U.in
```

`run` 的输出（真跑本机 QE 7.4）：

```
run 目录   : /home/liguozhou/oncv/uranium-ppgen/ae_radial
输入       : U.in（266 B，已按参数重新生成，覆盖了原来的同名文件）
qe_path    : /home/liguozhou/oncv/work/qe/q-e-qe-7.4/bin
ld1.x      : /home/liguozhou/oncv/work/qe/q-e-qe-7.4/bin/ld1.x
运行命令   : cd /home/liguozhou/oncv/uranium-ppgen/ae_radial && /home/liguozhou/oncv/work/qe/q-e-qe-7.4/bin/ld1.x < U.in > U.out
退出码     : 0
耗时       : 0.6 s
stdout     : U.out（5578 B，覆盖）
stderr     : U.err（170 B，覆盖，已回显到终端）
产物       : ld1.wfc（551902 B，覆盖）、vx.pot（181631 B，覆盖）
```

`输入` 与 `产物` 两行刻意写得直白，因为这两件事最常被问：

* 输入文件是**直接用的**（不给 `-i`：`直接跑：没读 input.json、没生成、没校验`）
  还是**按配置重新生成的**（给 `-i`：`已按参数重新生成，覆盖了原来的同名文件`；
  干跑时写「干跑：会按参数重新生成并覆盖，盘上现在这份 N B 没动」）；
* 每个产物是**覆盖**了原来就有的文件，还是这次**新建**的 —— 光看文件名和大小
  分不出这两件事。

### `ld1ps.wfc` 里的两半（QE 的约定）

`ld1.x` 的 `ld1ps.wfc` **一份文件里装两套轨道**（源码 `write_resultsps.f90:181-204`）：
前一半是**赝轨道** φ，后一半是把对应**全电子轨道** Ψ 乘一个相位对齐后写进来 ——
所以表头里的态标签会成对出现（U 上是 `5F 6S 6P 6D 7S | 5F 6S 6P 6D 7S`）。
`ppgen plot` 认这个约定：读出来的态名是 `5F (PS)` / `5F (AE)`，画图时
**同一态共用一个颜色、实线＝AE、虚线＝PS** —— 也就是论文 Fig. 2 左panel的画法。
副作用是没有这个约定的重复标签不再报错，而是加 `#2`、`#3` 后缀。

所以想画"赝轨道 vs 全电子轨道"，`data.wfcs` **只写 `ld1ps.wfc` 一份就够**；
再叠 `ld1.wfc` 只会多出重复的 AE 曲线（那份文件里的 AE 列与 `(AE)` 半完全一致，
个别态差一个整体符号，就是 QE 做的相位对齐）。

`ppgen run -i input.json` 这种相对写法也可以：路径会被 resolve 成绝对路径再显示
（否则报告里会出现 `run 目录 : .`、`cd .` 这种看不出到底在哪儿跑的显示）。

`plot` 的输出：

```
plot.json  : /home/liguozhou/oncv/uranium-ppgen/ae_radial/plot.json
mode       : wavefunction
数据       : ld1.wfc（label 'ae_radial'，3427 个网格点，r = 1.93e-06–53.08 Å）
           态：4F、5D、5F、6S、6P、6D、7S
量         : rR(r)（由 ylabel 'rR(r)' 推断）
横轴单位   : Å（由 xlabel 'r(Å)' 推断）
轴标题     : x='r(Å)'，y='rR(r)'
图标题     : ae_radial
坐标范围   : x 0.01–10 Å（对数轴），y -2–2
核态       : 未指定（所有态共用一套配色；要按核/价配色就写 figure.core）
已生成     : /home/liguozhou/oncv/uranium-ppgen/ae_radial/ae_radial.png（134863 B）
```

### 开关

| 命令 | 开关 | 作用 |
|---|---|---|
| 三条 | `-i/--input PATH` | 必填，配置文件（`input.json` / `plot.json`）的路径 |
| `report` | `-i/--input JSON JSON` | 必填两个：`output.json`（运行记录）+ `plot.json`（分析模板） |
| 三条 | `-o/--outdir DIR` | 目录（默认与配置文件同目录） |
| `prepare` | `-n/--dry-run` | 只把将生成的内容打到屏幕上，不写文件、不回写 `input.json` |
| `run` | `--qe PATH` | 覆盖 `qe_path`（QE 的 bin 目录或 `ld1.x` 本身） |
| `run` | `-n/--dry-run` | 只打印将要执行的命令，不跑、不写任何文件 |
| `run` | `--timeout SEC` | 超时秒数（默认 3600；`0` 表示不限） |
| `plot` | `-n/--dry-run` | 只读数据并报告将要画什么，不写图片 |
| `report` | `--median`/`--average`/`--maximum`/`--minimum` | 用哪个统计量排序（默认 `--median`） |
| `report` | `--scope scanned`/`all` | 排序范围：被扫通道（默认）/ 所有通道池化 |
| `report` | `-n/--dry-run` | 只打印 `report.md` 正文，不写文件 |
| `report` | `--no-figures` | 不给每个扫描点画图/写 `parser.json`，只写 `report.md` |

生成的 `U.in` **一个注释都没有**，只有 namelist 本身（下面这份对应当前
`ae_radial.yml`，也就是 U 全电子 + `dx=0.005` 的细网格）：

```
&input
  title       = 'U all-electron',
  atom        = 'U',
  zed         = 92,
  dft         = 'PBE',
  config      = '[Rn] 5f3 6d1 7s2',
  rel         = 1,
  iswitch     = 1,
  max_out_wfc = 7,
  xmin        = -8.0,
  dx          = 0.005,
  rmax        = 100,
/
```

`U.in` 里出现的变量就是 yml 里写的那些，一个不多一个不少 —— 所以上面这三行网格
参数是**你写进 yml 的**，不是 pp-gen 自己加的。`prepare` 之后再跑 `run`，产物
大小会随网格变：`dx=0.005` 是 3427 个网格点、`ld1.wfc` 551902 B；`dx=0.008`
（ld1.x 对 Z≥19 的自动值）是 2141 点、344856 B。

## `plot.json` 怎么写

```json
{
  "mode": "wavefunction",
  "data": {
    "wfcs": ["ld1.wfc"],
    "labels": ["ae_radial"]
  },
  "figure": {
    "xlabel": "r(Å)",
    "ylabel": "rR(r)",
    "title": "ae_radial",
    "xscale": "log",
    "xlim": [0.01, 10.0],
    "ylim": [-2.0, 2.0],
    "fig_name": "ae_radial.png"
  }
}
```

* **`mode`**：画图模式，目前只有 `wavefunction`（由 `ld1.wfc` 画径向波函数）。
* **`data.wfcs`**：要画哪些 `ld1.wfc`，**相对 `plot.json` 所在目录**解析。
  多个文件会叠加在同一张图上（对比用）。**`data.labels`**：与 `wfcs` 一一对应的
  标签，不写就用文件名（`ld1`）。
* **`figure`**：`xlabel` / `ylabel` / `title` / `xlim` / `ylim` / `fig_name`
  都按字面用；`fig_name` 不写就用 `plot.json` 的名字（`plot.png`），只取文件名、
  不能带目录（输出目录用 `-o`）。可选：`xscale` / `yscale`（`linear`/`log`）、
  `figsize`（英寸）、`dpi`（默认 180）、`quantity`、`x_unit`、`core`。
* **不认识的键直接报错**（把 `fig_name` 打成 `figname` 不该被静静忽略）。

### 对数横轴（`"xscale": "log"`）

对数轴上没有 0，所以 `xlim` 的下限必须为正。写法上不做硬性要求，但会被**修好
并说明**：

| 你写的 | 实际用的 | 报告里会写 |
|---|---|---|
| `"xlim": [0.0, 10.0]` | `0.01`–`10` Å | `xlim 下限 0 在对数轴上无效（对数轴没有 0 和负数），已换成 0.01 Å` |
| 没写 `xlim` | `0.01`–`10` Å | `横轴是对数轴但没写 xlim，按参考图取 0.01–10 Å；想改就写 figure.xlim` |
| `"xlim": [0.001, 5.0]` | 原样 | — |

（**纵轴不一样**：横轴有参考图定下的约定下限 0.01 Å，纵轴没有这种约定，
所以对数纵轴上不可用的 `ylim` 是"丢掉 + 按数据自动 + 警告"，不猜下限。）

下限的"同一段物理半径"按显示单位换算：`x_unit: "bohr"` 时 0.01 Å = 0.0189 bohr。
默认范围取的是参考图 `fig1_ae_radial_wfc_u_logx_lineary.png` 那套（0.01–10 Å）——
不这么定的话，自动缩放会把网格的 1.9e-6 Å 到 53 Å 全铺在轴上，7 个量级挤在一起。

另外：`r ≤ 0` 的网格点在对数轴上画不出来，会被跳过（真跳过了会报一句）；
纵轴仍是线性（只改横轴），要纵轴也取对数就写 `"yscale": "log"`。

### 画哪个量、横轴什么单位

优先看显式的 `figure.quantity` / `figure.x_unit`；没写就从轴标题推断，
**推断结果连同来源一起报出来**：

| `ylabel` 写成 | 画的量 |
|---|---|
| `rR(r)`、`u` | `rR(r)`（带符号，能看出节点） |
| `R(r)` | `R(r) = rR/r`（原点按 0 处理，那里 R 本身发散） |
| `4πr²R²(r)`、`rho` | 径向概率密度 `(rR)²` |

| `xlabel` 写成 | 横轴单位 |
|---|---|
| `r(Å)`、`r (ang)` | Å（`ld1.wfc` 里的 bohr × 0.529177210903） |
| `r (bohr)`、`r (a.u.)` | bohr（不换算） |

认不出来就落到默认：画 `rR(r)`、横轴用 Å，并在报告里说明是默认值。

### `ld1.wfc` 的列到底是 R(r) 还是 rR(r)？

**是 `rR(r)`**（即 `u(r) = r·R(r)`）。对着 QE 7.4 源码核过：

1. 写文件处不做任何缩放：`write_results.f90` 的 `psiaux(:,counter)=psi(:,1,i)`
   → `write_files.f90` 的 `write_wfcfile`，那里就是
   `WRITE(19,'(38f20.12)') grid%r(n), (wfc(n,ns),…)`（表头 `# r <标签…>`、
   每行 12 位小数，与实际文件吻合）。
2. U 用 `rel=1`，`scf.f90` 走 `CALL lschps(…, psi(1,1,n), nstop)`；而
   `lschps.f90` 的输出形参**就叫** `u`，它的归一化是 `sn = Σ grid%dx·r(i)·u(i)²`
   —— 对数网格上 `dr = r·dx`，也就是 **∫u²dr = 1**，正是 `u = rR` 的归一化。
3. 实测这份 `ld1.wfc`：七个态的 `∫ψ²dr` 都是 `1.000000`，`∫ψ²r²dr` 则
   0.35～18.7 各不相同；并且 `∫ψ²r dr`、`∫ψ²r²dr` 恰好等于 `ld1.x` 自己打印的
   `<r>`、`<r2>`（7S：4.0641 / 18.6879）。

于是 `ppgen plot` 直接把这些列当 `rR(r)` 用；想要 `R(r)` 或概率密度就写
`figure.quantity`（见下表），不用自己换算。

### 图的画法

* **颜色按态**：`figure.core` 给了（如 `["4F", "5D"]`）就用你参考脚本那套
  "核=蓝系 / 价=暖色"；没给就所有态共用一套定性配色，**也不声称谁是核谁是价**
  —— 核/价划分是元素相关的，不该由代码替你定。
* **线型按文件**：叠加多个 `wfcs` 时第 1 个文件实线、第 2 个虚线、第 3 点划线、
  第 4 个点线；单文件全是实线。这时图例分两块：`states`（颜色）与 `files`（线型）。
* 态的绘图顺序按 `n`/`l` 排（`7S 6D 5F 6P 4F 5D 6S` → `4F 5D 5F 6S 6P 6D 7S`），
  线色按这个顺序分配。
* 线性纵轴会画一条 `y=0` 参考线；对数纵轴画的是 `|量|`，轴标题会自动写成
  `|rR(r)|`（节点表现为向下尖峰）。
* 只有整条曲线落在 `ylim` 之外、或网格比 `xlim` 上限还短时，会额外提醒一句。
* **图内文字用英文**：matplotlib 默认字体不含中日韩字形，写了中文会被画成方框
  （真写了会给警告）。`Å` 这类拉丁扩展字符没问题。

画图需要 matplotlib（`pip install -e '.[plot]'`）；`cli` 里是懒导入，所以没装
matplotlib 也能正常用 `prepare` / `run`。

### 对数导数图（`"mode": "logfunction"`）

扫 `r_c` 时真正要看的是**对数导数曲线**：在固定半径 `r_d`（`U.in` 里的 `rlderiv`）
上解径向方程得 `L_l(E) = u_l'(r_d,E)/u_l(r_d,E)`，`ld1.x` 把 AE 侧写进
`ld1.dlog`、PS 侧写进 `ld1ps.dlog`，每行 `能量  L_{l=0}  L_{l=1} …`（Ry）。

```json
{
  "mode": "logfunction",
  "data": {"wfcs": ["ld1.dlog", "ld1ps.dlog"], "labels": ["AE", "PS"]},
  "figure": {"xlabel": "Energy(Ry)", "ylabel": "Log. der. (arbs. units)",
             "yscale": "log", "xlim": [-5, 5], "fig_name": "logfunction.png"}
}
```

* **`data.wfcs`**：任意几个 `dlog` 文件（**第一个当 AE、第二个当 PS**，极点对比就
  在这两个之间做）。`labels` 是图例里给文件的名字。
* **颜色按通道 `l`、线型按文件**（第 1 个实线、第 2 个虚线），所以同一通道的
  AE/PS 一眼能对上。
* **纵轴取对数时画 `|L|`**（轴标题自动写成 `|L|`）：`L` 在极点是发散的（这份数据
  `|L|` 从 3e-4 到 1e4，7 个量级），线性轴必然要截断。
* **对数纵轴上不可用的 `ylim` 会被丢掉并警告**，而不是替你猜一个下限：
  纵轴的合理下限取决于数据，猜出来只会误导。所以 `ylim: [-1.5, 2.0]` 这种
  （下限非正）会变成"按数据自动取值"+ 一条警告；合法但偏窄的 `ylim`
  （如两个正数）照用，只提示"裁掉了多少点"。**不写 `ylim` 就自动**，
  报告里的 `坐标范围` 会标出每一轴是"给定"还是"自动按数据"。
* **不写极值 CSV**（那张表是 `rR` 的极值，跟对数导数无关）：本模式默认不写，
  写了 `figure.csv_name` 会被忽略并提醒。
* 线性纵轴的 0 处**不画**参考线 —— 对 `L` 来说 0 是渐近线，不是参考值。

### `data.orbitals`：只画这几个轨道

```json
"data": {"wfcs": ["ld1.wfc", "ld1ps.wfc"], "labels": ["AE_radial", "PS_radial"],
         "orbitals": ["6S", "7S", "6P", "5F", "6D"]}
```

不写就全画。规则：

* 匹配按**态名**、忽略大小写：`6S` 同时选中 `6S`、`6S (PS)`、`6S (AE)`
  （`ld1ps.wfc` 那种成对表头两半都算同一个态）；
* 过滤在"读文件之后、比态之前"做 —— 所以**图、极值 CSV、报告里的"态："**
  都只剩选中的那几个，跨文件那条"态不一致"的提醒也不会再被没画出来的态干扰
  （例如筛掉 `ld1.x` 里没有赝轨道的 `4F`/`5D`，那条提醒就消失了）；
* 某个名字在所有文件里都没有 → **提醒**并忽略；全都匹配不上 → 报错（不会画出一张空图）。

### `analysis` 段：极点容差与误差统计

```json
"analysis": {
  "pole_threshold": 0.02,
  "mark_poles": true,
  "error": {"median": true, "average": true, "maximum": true}
}
```

**图上要不要竖线由你显式说了算**：`mark_poles` / `mark_zeros` 不写就是**不标**
（下面 `error` 的规矩也一样 —— 写了就按写的来，缺失的键算 `false`；只有整块
不写时才有默认值）。

* **`mark_poles`**（不写＝**不标**）：把极点**标在图上** —— 每个极点一条竖线，
  **颜色跟通道走**（与曲线同色）、**线型区分来源**（AE 点/PS 划）；未配上的那一对
  （PS 极点**和它的 AE 伙伴**）画得更粗更实，并在图顶竖排写上 **`ΔE`**
  （如 `Δ=0.180 Ry`）—— 于是"没对齐的那一对"在图上就是两根粗竖线。
  图例里会多出 `AE: pole / zero` / `PS: pole / zero` 两条说明。不写（或写
  `false`）就一条竖线都没有 —— **报告里照样给，只是不画到图上**。

  > 竖线一度太淡（`lw=1.0 / alpha=0.35` 的点线，在 8 条曲线里几乎看不见 —— 实测
  > 那一列只比邻列"深"15%）。现在配上是 `lw=1.4 / alpha=0.85`、未配上的是
  > `lw=2.0 / alpha=1.0`，同样几列实测"深"了 1.5–4 倍。
* **`mark_zeros`**（不写就**跟着 `mark_poles`**）：把 **`L(E)` 的零点**（`u'(r_d,E) = 0`，
  在对数纵轴上就是 `|L|` 向下尖到 0 的那些位置）也标出来 —— 画成**贴底的短刻度**，
  所以**从高度就能区分**：贯穿整幅的是极点、只贴底的是零点；颜色仍按通道、
  线型仍区分来源（AE 点线、PS 划线）。图例里那两条说明相应写成
  `AE: pole / zero`、`PS: pole / zero`。于是"只写 `mark_poles: true`"＝极点与
  零点都标，想只要极点就再写 `"mark_zeros": false`，想只要零点刻度就写
  `"mark_poles": false, "mark_zeros": true`。

  > 判据与极点互补、阈值同为 `amp=20`：**符号翻转且一侧 |L| 很大的是极点**
  > （`u=0`）、**翻转而两侧 |L| 都不大的是零点**（`u'=0`）；零点位置按线性插值
  > 给到 `L=0`。报告里每个通道也顺带给出零点个数（`L 零点 AE/PS 2/2`）。
* **`pole_threshold`**（Ry，默认 `0.15`，同参考脚本）：PS 的某个极点与**最近的
  AE 极点**在 x 轴（能量）上差得不超过它，就算**同一个极点**；差得更远的算
  "未配"，报告里逐条列出它到最近 AE 极点的 **x 轴距离 `ΔE`**（"否则显示 pole 的
  x 轴距离"）。窗口里一个 AE 极点都没有时写明"（窗口内无 AE 极点）"。
  另外仍会标注 `ΔE > 1 Ry` 的为**疑似鬼态**（参考脚本的判据），以及"AE 有、PS
  没有的极点"—— 反过来那一侧也看得见。
* **`error`**：要显示哪些 **`|L_PS − L_AE|`** 统计（在 `xlim` 窗口内、**逐通道**）：
  **`minimum`（最小值）/ `median`（中位数）/ `average`（平均值）/ `maximum`（最大值）**。
  **写了 `error` 就按写的来**：没列出的键算 `false` —— `{"median": true}` 就只显示
  中位数，不会冒出平均和最大；四个都写 `false`（或写个空块）就一个字都不印。
  **整块不写**才是四个都要。**只有这一套统计** —— 极点位置差不再另给平均/中位
  之类的汇总（逐对差值见 `parser.json` 的 `pairs`）。
* 两者都**只在 `mode="logfunction"` 下生效**（其它模式写了会提醒并忽略）。
  统计既进命令行报告，也**写进图的图例**（图例里用英文缩写 `med/avg/max` ——
  matplotlib 默认字体画不出中文，会变方框）。

命令行报告除了轴范围，还给出**极点诊断**（这是这个模式的价值所在）：

```
analysis   : pole_threshold=0.02 Ry，图上标极点，error: minimum、median、average、maximum（同一极点的能量容差；误差是 |L_PS − L_AE|）
极点诊断   : 第一个文件当 AE、第二个当 PS
           : l=0  AE/PS 2/2，同一点 2、未配 0（其中疑似鬼态 0）；L 零点 AE/PS 2/2   |ΔL| 最小 3.666e-09，中位 0.00041067，平均 0.9187，最大 263.82
           : l=1  AE/PS 2/2，同一点 1、未配 1（其中疑似鬼态 0）   |ΔL| 最小 5.6659e-05，中位 0.0074168，平均 10.884，最大 1981
             · 未配上的极点 4.490->4.310  Δ=0.180
           : l=2  AE/PS 1/1，同一点 1、未配 0（其中疑似鬼态 0）   |ΔL| 最小 1.3406e-07，中位 0.0028996，平均 3.2466，最大 767.33
           : l=3  AE/PS 2/2，同一点 1、未配 1（其中疑似鬼态 0）   |ΔL| 最小 0.00028888，中位 0.019203，平均 16.197，最大 4846.7
             · 未配上的极点 4.070->4.390  Δ=0.320
结论       : 各通道都没有多出来的极点
```

判据与 `ppgen/script/rc_scan_transferability.py` 完全一致（已逐通道对照过）：
极点 = 相邻两点 `L` 符号翻转且一侧 `|L| > 20`；PS 的极点按与 AE 极点的距离分成
`matched`（≤0.15 Ry）、`shifted`（≤1 Ry，高能区几乎必然出现，**不是**鬼态）、
`extra`（1 Ry 内配不上 —— **疑似鬼态**）。**AE 自己也有极点**，所以判鬼态要看
"多出来的"，而不是数个数。

### `parser.json`：机器可读的分析结果

对数导数模式下，`plot` 还会在图片旁边写一份 **`parser.json`**（`figure.parser_name`
可改名；写 `""` 就不写），把这次分析整理成便于**扫完 55 个点后统一解析**的形式：

```json
{
  "plot_json": ".../plot.json",
  "generated_at": "2026-09-24T15:46:03+08:00",
  "mode": "logfunction",
  "files": [{"label": "AE", "path": ".../ld1.dlog", "n_points": 501,
             "channels": [0,1,2,3], "energy_min": -5.0, "energy_max": 5.0},
            {"label": "PS", "path": ".../ld1ps.dlog", ...}],
  "window": {"energy_min": -5.0, "energy_max": 5.0},
  "settings": {"pole_threshold": 0.02, "mark_poles": true, "mark_zeros": true,
               "error": {"minimum": true, "median": true, "average": true,
                         "maximum": true}},
  "channels": [
    {"l": 1,
     "poles": {"ae": {"n": 2, "energies": [-1.45, 4.49], "min": -1.45, "max": 4.49},
               "ps": {"n": 2, "energies": [-1.45, 4.31], "min": -1.45, "max": 4.31},
               "coincide": false,
               "pairs": [{"ae": -1.45, "ps": -1.45, "delta": 0.0},
                         {"ae": 4.49, "ps": 4.31, "delta": 0.18}]},
     "zeros": {"ae": {"n": 2, "energies": [-1.904947, 1.347777], ...},
               "ps": {...}, "coincide": false, "pairs": [...]},
     "error": {"minimum": 5.67e-05, "median": 0.0074168,
               "average": 10.884, "maximum": 1981.0}}
  ],
  "figure": ".../logfunction.png"
}
```

* **`settings`** 把 plot.json 里决定结果的参数原样记下来（窗口、`pole_threshold`、
  标不标极点/零点、要哪些误差统计），所以单看这一份就知道当时怎么算的。
* **`poles` / `zeros`** 两块结构相同：AE 与 PS 各自的极端位置（`min`/`max`，
  即"极大值/极小值"两个端点）、**是否重合**（`coincide`：个数相同且每个都能在
  `pole_threshold` 内配上；两边都空也算重合）、逐对差（`pairs` 里每对的
  `delta`，以及 `ae: null` 表示窗口里没有 AE 特征可配）。
  **不另给平均/中位之类的汇总统计**（已按要求删除）。
* **`error`** 是该通道 `|L_PS − L_AE|` 的 最小 / 中位 / 平均 / 最大。
* 干跑（`-n`）不落盘，但报告里照样给出路径；非 `logfunction` 模式写了
  `parser_name` 会被忽略并提醒。

### 极值表 CSV

`plot` 除了 PNG，还会写一份**极值表**（`rR(r)` 的每一个局部极大/极小），
默认文件名是 `<fig_name 主名>_extrema.csv` —— 你的 `ae_radial.png` 得到
`ae_radial_extrema.csv`，与你参考 run 里那份同名同格式：

```
state,n_nodes,idx,r_bohr,r_angstrom,rR
4F,0,1,2.8599281895e-01,1.5134088227e-01,1.7585306514e+00
5D,2,3,7.9708928745e-01,4.2180148597e-01,1.2443807668e+00
5D,2,2,2.7340841474e-01,1.4468150235e-01,-6.9094497348e-01
5D,2,1,8.4856933398e-02,4.4904355341e-01,4.8458155693e-01
…
```

* **列**：只要 `state, n_nodes, idx, r_bohr, r_angstrom, rR`（参考文件里的
  `class` / `kind` / `abs_rR` 都不要）。数值一律 `%.10e`，与参考文件一致。
* **排序**：先按态（`n`/`l`，与图的图例、配色同一顺序：`4F 5D 5F 6S 6P 6D 7S`），
  **态内按 `|rR|` 从大到小**。
* **`idx` 仍是该态内按 `r` 递增的序号**，不随排序重编号 —— 所以排序后 `idx`
  会跳号（上面 `5D` 是 3、2、1），这是有意的：能一眼看出"这是第几个极值"。
* **节点与极值**的算法沿用 `ppgen/script/ae_radial_wfc.py`：节点是符号变化点
  （线性插值到零点），极值过滤掉 `< 1e-8` 的数值抖动 —— `ld1.wfc` 只写 12 位
  小数，不加这道过滤，尾部会冒出一堆假极值。
* **改名 / 关掉**：`figure.csv_name`（如 `"peaks.csv"`）改名；写 `""` 就完全不
  写 CSV。多个 `wfcs` 时**每个文件一份**（`<主名>_<label>.csv`）—— 不同文件的
  同名态在同一张表里分不清。
* 态的顺序与颜色分配用的是同一套 `n`/`l` 排序，所以 CSV 的行序与图的图例一一对应。

## `ppgen run` 的两套流程

```bash
ppgen run -i .../ae_radial/input.json      # ① 先生成、再运行（按 input.json）
cd .../runs/5F/rc_0p9 && ppgen run         # ② 直接跑当前目录里的 U.in
ppgen run -i U.in                          # 同上（给的就是一份 ld1.x 输入）
```

**① 给 `-i input.json`：先生成、再运行。** 每次都按配置重新生成一遍输入
（`<元素>.in`；扫描型则是 `runs/<轨道>/rc_<值>/U.in` 每个点一份），**覆盖同名
文件**，然后逐个跑 `ld1.x`。所以手改过的 `U.in` 会被覆盖 —— 想保住手改版就用
②（或者手改 `input.json` 里对应的参数，让生成结果就是你要的）。

**② 不给 `-i`（或 `-i` 给的是 `U.in` 这种）：在当前目录直接跑。** 就是"在这个
目录里把它跑一遍"：**不读 `input.json`、不生成输入、也不做任何参数校验**（坏输入
直接由 `ld1.x` 报错）。输入文件按 `U.in` → 目录里唯一的 `*.in` → 报错让你点名 的
顺序找；`ld1.x` 取 `PPGEN_QE_PATH`（默认 `/home/liguozhou/oncv/work/qe/q-e-qe-7.4/bin`），
也可以用 `--qe PATH` 指定。`--only`/`--rc` 在直接模式下会直接报错（那些是
`input.json` 那一套）；`--dry-run`/`--timeout` 照旧可用。

两种模式跑完都会在输入旁边写 `output.json`（直接跑那份里 `"direct": true`、
`"workflow": ""`），报告里输入那行也会写明"直接跑：没读 input.json、没生成、没校验"。

## `ppgen report`：把一次扫描汇成报告

```bash
ppgen report -i output.json plot.json --median   # 按 |ΔL| 中位数挑 r_c
ppgen report -i output.json plot.json --pole     # 按极点是否重合 + 极点 ΔE 挑 r_c
```

两个输入都用 `-i` 给（顺序无所谓，**按内容识别**：`output.json` 里有 `tasks`、
`plot.json` 里有 `mode`/`figure`）：

* **`output.json`** —— `ppgen run` 写的运行记录：有哪些扫描点、在哪个目录、
  哪个轨道、`r_c` 多少、跑成没跑成；
* **`plot.json`** —— 分析模板：比哪两个文件（`ld1.dlog`/`ld1ps.dlog`）、
  能量窗口（`xlim`）、极点容差与开关（`analysis`）。

它对**每个扫描目录用同一套设置重算一遍**，再按轨道挑出排序统计量最小的 `r_c`，
写进 `output.json` 旁边的 **`report.md`**（`-o` 改目录，`-n` 只打印不落盘），
同时画一张 **`report_<统计量>.png`**（`--no-decision-figure` 可关掉）。

顺带**每个算成功的扫描点**也会在它自己的目录里产出（`--no-figures` 可关掉）：

| 产物 | 说明 |
|---|---|
| `logfunction.png` | 与 `ppgen plot` 画的一模一样（极点竖线/零点刻度、标 ΔE、图例带 `min/med/avg/max`） |
| `parser.json` | 机器可读的逐通道结果（极点/零点的极值与是否重合、`pairs`、`|ΔL|` 统计） |
| `plot.json` | 目录里**没有**才补一份分析模板 —— 已存在的不覆盖（可能是你手改过的） |

于是像 `runs/5F/rc_0p9` 这样的目录，跑完 `report` 就自带图与分析结果，不必再手写
`plot.json` 或逐个 `ppgen plot`。某个点画图失败只记一条提醒，不影响整份报告。

### 选择依据图：`report_<统计量>.png`

`report.md` 旁边那张图回答的是**"凭什么选这个 `r_c`"**：每个轨道一张子图，
**横轴都是各个扫描点的 `r_c`**，纵轴是**决定排序的那个量**，选中的点用星号 +
竖虚线标出：

| 命令 | 纵轴（主轴） | 右侧淡灰（对照） |
|---|---|---|
| `--pole` | **极点 ΔE 合计**（`pairs[].delta` 之和），虚线＝`pole_threshold` | `\|ΔL\|` 中位数（ΔE 一样时的次级判据） |
| `--median` / `--average` / `--maximum` / `--minimum` | **`\|ΔL\|` 的那个统计量**（对数轴） | 极点 ΔE（虚线＝`pole_threshold`） |

几个读图的细节：

* **实心点＝极点重合、空心点＝不重合** —— 所以 `--pole` 时第一眼就能看出哪些点
  在容差以内（`6S` 的 1.2/1.3/1.7/1.8 是实心，`5F` 全是空心）；
* **`∞`**：那个通道的极点配不上（窗口里没有 AE 极点可配），标在子图顶上；
* 跑失败的点不画线，只把 `r_c` 记在子图角落（`failed r_c: 0.6, 0.7, 1`）；
* `--median` 那张图上能直接看到"中位数最小的那个点，极点却已经出了容差线" ——
  这正是 `--pole` 会给出不同答案的原因。
（这些产物**不写进 `report.md`** —— 需要时看目录就行；报告里只列"画图失败的点"。）

```
最佳截断半径:
  5F   l=3  r_c = 1.3 bohr   |ΔL| 中位数 = 0.018965   rc_1p3
  6S   l=0  r_c = 1.2 bohr   |ΔL| 中位数 = 4.8892e-05   rc_1p2
  6P   l=1  r_c = 1.8 bohr   |ΔL| 中位数 = 0.0063025   rc_1p8
  6D   l=2  r_c = 2.4 bohr   |ΔL| 中位数 = 0.0028996   rc_2p4
  7S   l=0  r_c = 1.4 bohr   |ΔL| 中位数 = 6.2735e-05   rc_1p4
可用扫描点 : 45/55
```

* **排序统计量**：`--median`（默认）/ `--average` / `--maximum` / `--minimum`，
  指 `|L_PS − L_AE|` 的哪一个统计量最小者胜。
* **`--pole`：按"极点对不对得上"挑 r_c**（不看 `|ΔL|` 谁小）。次序是：

  1. **极点是否重合**（`parser.json` 里 `poles.coincide`）—— 重合的排前面；
  2. 都重合（或都不重合）时比 **极点 ΔE 合计** —— 也就是 `poles.pairs[].delta`
     之和（如 5F 的 l=3 是 0.32、6P 的 l=1 是 0.14），**越小越好**；窗口里
     压根没有 AE 极点可配的那些对（`delta: null`）没有 ΔE 可比，只能算"结构上
     就不一致"，记 **∞**，排在同一组的最后；
  3. ΔE 也一样才比 **`|ΔL|` 的中位数**（次级判据，固定用中位数）。

  也就是说"极点对得上优先，对得上之后再比谁对得准，全都一样才看曲线差多少"。
  `--scope` 一样管用：`scanned` 只看被扫轨道那个通道的极点，`all` 把所有通道
  的 ΔE 加起来。报告里会把**逐通道**的极点 ΔE 都列出来（加粗＝未重合），
  一眼看出是哪个通道对不上：

  ```
  最佳截断半径:
    5F   l=3  r_c = 1.2 bohr   极点不重合  ΔE 合计 0.32  |ΔL| 中位 = 0.019056   rc_1p2
    6S   l=0  r_c = 1.2 bohr   极点重合  ΔE 合计 0  |ΔL| 中位 = 4.8892e-05   rc_1p2
    6P   l=1  r_c = 1.8 bohr   极点不重合  ΔE 合计 0.14  |ΔL| 中位 = 0.0063025   rc_1p8
    6D   l=2  r_c = 2.4 bohr   极点重合  ΔE 合计 0.02  |ΔL| 中位 = 0.0028996   rc_2p4
    7S   l=0  r_c = 1.4 bohr   极点重合  ΔE 合计 0  |ΔL| 中位 = 6.2735e-05   rc_1p4
  ```

  实测这条判据会**改答案**：5F 的 `--median` 选 `r_c = 1.3`（中位数 0.018965
  最小），但 `--pole` 选 `r_c = 1.2` —— 因为 1.3 处 l=3 的极点 ΔE 是 0.34 Ry，
  而 1.2 处是 0.32 Ry；ΔE 相同的几个点（0.9 / 1.1 / 1.2）里才比中位数。
* **排序范围 `--scope`**（这个很关键）：一次只扫一个轨道，所以**其它**通道在
  同一次扫描里保持不变，池化进来会把差别稀释掉 —— 实测 5F 那一轮池化中位数在
  各点间只差第 4 位有效数字（0.0043246 vs 0.0043247），而**被扫通道**（l=3）的
  中位数从 0.0210 变到 0.0190，能分出好坏。所以默认 `--scope scanned`：
  用**被扫轨道对应的通道**排序（轨道标签末尾字母 → `l`：s=0 p=1 d=2 f=3…）；
  `--scope all` 则按所有通道池化。
* `report.md` 里除了逐点明细表（含每个通道的统计量与失败原因），还有
  **「最佳参数下的极点与零点」** 一节：对每个轨道选出的那个 `r_c`，逐通道列出
  极点（`L` 发散处）与零点（`L=0`）的**极值位置**（`min`/`max` 与全部位置）和
  **是否重合**（个数相同且每个都在 `pole_threshold` 内配上），并给出逐对差值。
  跳过的目录（跑失败、缺 `ld1.dlog`/`ld1ps.dlog`）单列一节说明原因。

## `configts`：逐个测试组态测同一份赝势（迁移性，论文表 3/表 5）

```json
{"workflow": "configts",
 "parameters": {"input": "rc_scan", "configts": "configts"},
 "qe_path": "/…/bin"}
```

```yaml
# parameters/configts.yml —— 只放"要测哪些组态"
configts:
  - '5f3 6s2 6p6 6d1 7s2'      # 参考组态（与生成用的 config 一致）
  - '5f3 6s2 6p6 6d0 7s2'
  - '5f3 6s2 6p6 6d2 7s1'
```

```yaml
# parameters/input/test_configts.yml —— 数值都在这里（数据的唯一来源）
input:
  iswitch: 2                   # 测试已生成的 PP（不是生成）
  # …config/dft/rel/nld/rlderiv/eminld/emaxld/deld 照写
test:
  nconf: 1
  'configts(1)': '5f3 6s2 6p6 6d1 7s2'
  file_pseudo: '/abs/path/U.pbe-n-nc.UPF'   # 要测的那一份（固定不变）
```

**两种模式**（由参数文件说了算）：

1. **单次运行（推荐、最简单）**：参数文件里直接写 `nconf = N` 与
   `configts(1..N)` → `ppgen run -i input.json` **只跑一次 ld1.x**（5 个组态约 5 s），
   QE 自己逐组态打印 `dEtot_ae` / `dEtot_ps` / `Delta E`（并写 `ld1<n>*.dlog`）。
   辅助组 `parameters.configts` 可以不写（写了也会被忽略）。
2. **逐组态扫描**：参数文件里只写一个 `configts(1)`（`nconf = 1`），把清单放到
   辅助组 `parameters.configts` 里 → 每个组态一个任务，落在 `runs/<组态>/`，
   每个目录里的 `U.in` **只换 `&test` 的 `configts(1)`**。

两种模式下 `iswitch`、`file_pseudo`、`nconf`、有没有 `&inputp`/卡片，一律照参数
文件来（工作流不替你做主）。两种模式算出来的 δ 完全一致（实测差 <1 µRy）。
* **参数文件必须是 `iswitch=2`**（`iswitch=2` 配卡片会被拦下）：`iswitch=3` 是
  "生成 PP + 测试"，它内存里的 PP 按生成组态做非屏蔽，测试组态一换 PS 总能就
  对不上（实测同一组态：`iswitch=3` 每组态单独跑给 −183.196 Ry，同一次 `nconf>1`
  或 `iswitch=2` + 固定 `file_pseudo` 给 −181.60 Ry，差 **1.6 Ry**）⇒ ΔE_PS 与 δ
  会完全失真。写成 `iswitch=3` 时 `prepare/run` 会**提醒**这一点（不改数据）。
* `ppgen report -i output.json`：**只给一个输入**（不需要 plot.json），产出
  `report.md`：迁移性表（每个组态一行 `E_AE / E_PS / ΔE_AE / ΔE_PS / δ`，
  δ = ΔE_AE − ΔE_PS 与 QE 的 `dEtot` 同号约定）+ 逐组态的 AE/PS 能级表。
  提醒会自动带出：价电子数与参考不同的组态、`|ΔE_AE|` 大得离谱（AE 落到别的
  解上）、参考组态不在任务列表里（那时改用第一个成功的组态当基准）。
* `--only` 可按序号、目录标签或组态串挑几个；`--rc` 在这个工作流里没有意义
  （报错拦下）。

## `rc_scan`：逐轨道扫 rcut

`input.json` 用的是**两个**参数组，一个是 ``ld1.x`` 的 namelist，另一个是扫描范围：

```json
{
  "workflow": "rc_scan",
  "parameters": {"input": "rc_scan", "rc_scan_range": "rc_scan_range"},
  "qe_path": "/home/liguozhou/oncv/work/qe/q-e-qe-7.4/bin"
}
```

* `parameters/input/rc_scan.yml`：`&input` / `&inputp` / **通道卡片** / `&test`，
  其中 `orbitals` 是逐轨道的 8 列卡片（`label n l occ ener rcut rcutus jj`），
  `n_orbitals` 是张数（对不上会报错）。
* `parameters/rc_scan_range.yml`：给偏移，以及可选的 `update`（见下）。

```yaml
rc_range: [-0.5, 0.5, 0.1]     # [起点, 终点, 步长]，含两端
update: true                   # 扫完一个轨道就把它的 rcut 更新成最佳值
```

**扫描规则**（一次只动一个轨道）：

* **扫描中心就是 `orbitals` 卡片里那个轨道的 rcut**（5F → 1.0641 bohr），
  不做任何四舍五入 —— 所以 `input` 参数是唯一的真值来源，不需要再另写一份"中心值"。
* 对轨道 X 的每个偏移 `卡片 rcut + offset`（5F → 0.5641、0.6641 … 1.5641），生成一份输入，
  落在 `runs/<轨道>/rc_<值>/U.in`，例如 `runs/5F/rc_1p0641/U.in`（**绝对值**不是偏移；
  小数位写足，名字就是那个点的真实 rcut）。
* 该轨道的 `rcut` **与** `rcutus` 一起改成同一个值（``ld1.x`` 的 ``read_psconfig``
  要求 `rcut <= rcutus`）；**其余轨道保持** `input` 参数里的精确值（6S 仍是 1.3527），
  这样一次只动一个变量。
* 每个点一个目录是必须的：所有点的 `file_pseudopw` 都叫 `U.pbe-n-nc.UPF`，
  挤在一个目录里会互相覆盖。

### 扫完自动更新 rcut（`update`）

`rc_scan_range.yml` 里写 **`update: true`** 时，`ppgen run -i input.json` 变成
**逐轨道串行**：

1. 按当前参数生成这个轨道的所有扫描点并逐个跑；
2. 拿跑出来的 `ld1.dlog` / `ld1ps.dlog` 挑最佳点 —— 判据与 `ppgen report --pole`
   **完全一致**：**极点是否重合 → 极点 ΔE 合计 → `\|ΔL\|` 中位数**（只看这个轨道
   对应的那个通道，即 report 的 `--scope scanned`）；
3. 把这个 rcut **就地写回提供 `orbitals` 卡片的那份参数**
   （`parameters/input/rc_scan.yml`），再用新值生成并扫描**下一个**轨道。

所以一次只动一个变量、而且**优化完的轨道从此固定在最佳值上**：6S 扫到 1.2 之后，
后面 7S/6P/6D/5F 的输入里 6S 就是 1.2。跑完一轮，那份参数里的 rcut 就都是当前最优值 ——
再跑一轮 `rc_scan` 就是以它们为中心重新细扫（迭代收敛）。

写回只动**那一行**里的两个数（`rcut`/`rcutus`），缩进、引号、行尾注释、其它列、
文件其余部分一律原样保留（这些 yml 是手写的，整份重排会把注释排版弄丢）。

```
更新设置   : 极点是否重合 → 极点 ΔE 合计 → |ΔL| 的中位数（极点容差 0.15 Ry）；能量窗口 &input 的 eminld/emaxld；逐轨道串行：扫完一个轨道就更新它的 rcut，再用新值扫下一个
…
更新 r_c   : 6S → 1.2 bohr（原 1.3）
           : 最佳点 .../runs/6S/rc_1p2，极点重合，ΔE 合计 0 Ry，|ΔL| 中位 4.8892e-05
           : 已写回 /home/liguozhou/oncv/uranium-ppgen/rc_scan/parameters/input/rc_scan.yml（后面的轨道用这个新值生成）
```

* **判据可以细调**（不写就是 `true` 那套）：

  ```yaml
  update: true                                    # 等于 {"criterion": "pole"}
  update: {"criterion": "pole", "pole_threshold": 0.02}
  update: {"criterion": "median"}                 # 只看 |ΔL| 中位数
  update: {"criterion": "pole", "window": [-5, 5]}  # 默认取 &input 的 eminld/emaxld
  ```

  * `criterion`：`pole`（默认）/ `minimum` / `median` / `average` / `maximum`；
  * `pole_threshold`（默认 `0.15` Ry，与 `report` 里那个默认值一致）决定"算不算重合"；
  * `window` 不写就用 `&input` 的 `eminld` / `emaxld`（与 `plot.json` 的 `xlim` 是同一个窗口）。
* **需要 `ld1.dlog` / `ld1ps.dlog`**（即 `iswitch=3` + `rlderiv`，这份参数就是）。
  缺这两个文件的点（跑失败、或者被 `--rc` 过滤掉）不参与比较；**一个都挑不出来就不更新**，
  并在报告与提醒里说明原因 —— 不会瞎猜一个值写回去。
* 更新记录写进 **`output.json` 的 `updates`**（新值 / 原值 / 判据 / 极点是否重合 / ΔE /
  最佳点目录 / 写回的文件），方便回头查"这个 rcut 是哪次扫描定的"。
* 干跑（`-n`）不跑也不更新，报告里会给一条提醒；直接跑（不给 `-i`）与它无关。

`prepare` 的报告：

```
扫描       : 5 个轨道 × 11 个点 = 55 个扫描点（中心取 orbitals 卡片里的 rcut + rc_range 偏移…）
           : 5F  rc = 0.5641、0.6641、0.7641、0.8641、0.9641、1.0641、1.1641、1.2641、1.3641、1.4641、1.5641   ->  runs/5F/rc_*/
产物目录   : /home/liguozhou/oncv/uranium-ppgen/rc_scan/runs
已生成     : 55 个 U.in（如 .../rc_scan/runs/6P/rc_1p6133/U.in）
```

`run` 逐个点跑（ld1.x 是单线程的，所以顺序执行），每个点一行，末尾汇总；
**中途失败不中断** —— 扫描的意义就是把点跑完再挑：

```
目标       : 55 个扫描点
逐个结果   :
  [  1/55] 6S rc=0.8527  成功  1.6 s  U.pbe-n-nc.UPF、ld1.dlog、ld1.pwe、ld1.test、ld1.wfc、ld1ps.dlog、…
  [  2/55] 6S rc=0.9527  成功  1.6 s  U.pbe-n-nc.UPF、ld1.dlog、…
  …
汇总       : 成功 55/55，总耗时 84.2 s
```

退出码：全部成功 0；**只有一个目标时**原样透传 `ld1.x` 的退出码，多目标时有失败就是 1
（失败清单在 `失败` 那行，具体错误行也一并打出来）。超时**只算那一个点失败**并
记进记录，不会把整批打断。

跑完会在 `input.json` 旁边写一份 **`output.json`**（干跑不写），逐个任务记下
**exit code / exit reason / 运行目录**，外加复现需要的输入文件名、耗时与产物清单：

```json
{
  "input": "/home/liguozhou/oncv/uranium-ppgen/rc_scan/input.json",
  "workflow": "rc_scan",
  "run_dir": "/home/liguozhou/oncv/uranium-ppgen/rc_scan",
  "qe_path": "/home/liguozhou/oncv/work/qe/q-e-qe-7.4/bin",
  "ld1": "/home/liguozhou/oncv/work/qe/q-e-qe-7.4/bin/ld1.x",
  "started_at": "2026-09-24T14:29:56+08:00",
  "finished_at": "2026-09-24T14:30:00+08:00",
  "summary": {"total": 3, "succeeded": 2, "failed": 1, "seconds": 3.975, "exit_code": 1},
  "tasks": [
    {
      "label": "5F rc=0.6", "orbital": "5F", "rcut": 0.6,
      "directory": ".../runs/5F/rc_0p6",
      "input": "U.in",
      "exit_code": 1,
      "exit_reason": "Error in routine run_pseudo (1): / Errors in PS-KS equation",
      "seconds": 0.858,
      "artifacts": ["CRASH", "input_tmp.in", "ld1.dlog", "ld1.wfc", "vx.pot"]
    },
    { "label": "5F rc=0.9", "...": "...", "exit_code": 0, "exit_reason": "正常结束" }
  ],
  "updates": [
    {
      "orbital": "6S", "rcut": 1.2, "previous": 1.3,
      "criterion": "pole", "pole_threshold": 0.15,
      "coincide": true, "delta": 0.0, "stat": 4.8892e-05,
      "file": ".../parameters/input/rc_scan.yml",
      "written": true, "directory": ".../runs/6S/rc_1p2", "note": ""
    }
  ]
}
```

`updates` 只有开了 `update: true` 才有内容：扫完每个轨道挑出的最佳 rcut、更新前后的值、
判据与判据里的数、写回哪份参数、最佳点目录（见「扫完自动更新 rcut」一节）。

`exit_reason` 的几种取值：`正常结束`（退出码 0）、QE 的 `Error in routine …` 那几行
（退出码非 0）、`跑了 N 秒还没结束，已被终止`（超时，`exit_code` 为 `null`）、
`找不到输入文件 U.in`（那个点缺输入、被跳过）。

**有些点会失败，这是扫描要发现的东西。** 扫描跑小了 rcut，QE 的
`run_pseudo` 会直接报

```
Error in routine run_pseudo (1):
Errors in PS-KS equation
```

（该点没有 UPF；AE 那半截是跑完的）。这不是工具出错：`ppgen run` 会照样把点全部
跑完，把失败清单与错误行打出来，并让退出码非 0。想重跑失败的点就 `--only`/`--rc`
挑出来重跑（实测同样的输入失败是可复现的，不是随机）。

批量时每个点的 stderr 只写进各自的 `U.err`，**不逐条回显**到终端（55 个点会有
55 段 OpenMPI/`MPI_ABORT` 噪声）；只有单个目标时才照旧回显。

扫描有筛选开关（`prepare` 与 `run` 都认，可逗号分隔、可合用）：

```bash
ppgen run -i .../rc_scan/input.json --only 5F                    # 只跑 5F 那 11 个点
ppgen run -i .../rc_scan/input.json --rc 1.0641                  # 只跑 rcut=1.0641 的那个点
ppgen prepare -i .../rc_scan/input.json --only 5F --rc 1.0641    # 只生成那一个点
```

`--rc` 取的是**真实 rcut**（卡片中心 + 偏移），不是偏移量；不确定有哪些值就先跑 `--dry-run`，报告里会列出每个轨道扫到的全部 rcut。

## `rlderiv`：扫 `&input` 的对数导数参考半径

`rlderiv` 决定"极点在哪个能量上出现"（极点 = 波函数在 `r=rlderiv` 处有过零），
所以换它可以把极点挪出你关心的能量段 —— 而且**它不改赝势**（见下）。它不是卡片上
的一列，而是 `&input` 里的一个标量，用的是和 `rcloc` 同一套 namelist 扫描机制：

```json
{
  "workflow": "rlderiv",
  "parameters": {"input": "energy_scan", "rlderiv": "rlderiv"},
  "qe_path": "/home/liguozhou/oncv/work/qe/q-e-qe-7.4/bin"
}
```

* `parameters.input`：**要测的那份 PP 设置**（`&input` / `&inputp` / 通道卡片 /
  `&test`，与本工作流无关的段照旧），里面的 `rlderiv` 就是扫描中心；
* `parameters.rlderiv`：范围与可选的 `update`：

```json
{"rlderiv": [4.0, 8.0, 1.0]}
```

**范围怎么解释**（`rcloc` 也适用）：

| 写法 | 含义 |
|---|---|
| `[4.0, 8.0, 1.0]` | 两端都为正、**又把中心（`&input` 里现在的 `rlderiv`）包在里面** → 按**绝对值**取：4、5、6、7、8 |
| `[-1.0, 3.0, 1.0]` | 没把中心包住（或有非正数）→ 按**相对中心的偏移**取（5.0 时是 4、5、6、7、8）|
| `{"offset": [1.0, 3.0, 1.0]}` | 明确按偏移 |
| `{"absolute": [6.0, 8.0, 1.0]}` | 明确按绝对值 |

每个值一个目录 `runs/rlderiv/rlderiv_<值>/U.in`，里面**只有 `&input` 的 `rlderiv`
被改成那个值**；`ppgen run -i input.json` 逐个跑完：

```
扫描       : &input 的 rlderiv × 5 个点 = 5 个扫描点（&input 的 rlderiv，范围写在 ['rlderiv'] 里，每个值一个目录）
           : rlderiv  = 4、5、6、7、8   ->  runs/rlderiv/rlderiv_*/
```

跑完在 run 目录里 `ppgen report -i output.json plot_dlog.json`：**每个点目录都会
生成自己的 `logfunction.png` 与 `parser.json`**（用的是你那份 `plot_dlog.json`
的窗口/容差/线型），并写出 `report.md` —— 每个 `rlderiv` 一行，**每个 l 一列给出
`|ΔL|` 的中位数**，末尾还有被选中那个 `rlderiv` 的极点/零点对照：

```
| rlderiv (a.u.) | 状态 | l=0 中位数 | l=1 中位数 | l=2 中位数 | l=3 中位数 | 池化中位数 |
| 4.0 | ok | 0.00021245 | 8.8633e-05 | 4.5817e-05 | 1.7052e-05 | 9.035e-05 |
| 5.0 | ok | 0.00011071 | 2.6424e-05 | 1.6769e-05 | 7.3193e-06 | 2.9083e-05 |
| 6.0 | ok | 8.6703e-05 | 5.6425e-05 | 6.108e-06 | 1.6319e-06 | 1.7644e-05 |
| 7.0 **←** | ok | 5.2452e-05 | 1.8828e-05 | 8.26e-06 | 9.9612e-06 | 1.3269e-05 |
| 8.0 | ok | 5.1653e-05 | 2.5881e-05 | 2.63e-05 | 2.3941e-05 | 2.7342e-05 |
```

因为 `rlderiv` 影响整份赝势（所有通道），排序与报告会自动按**所有通道池化**
（等价于 `--scope all`），并在报告里说明一句。

`update: true` 时把挑出来的最佳 `rlderiv` 写回 `parameters.input` 那份文件的
`&input`（只换那一行的数值，注释与别的键不动）。

> `rlderiv` 与 `rcut` / `ener` / `rcloc` 不同：**它不改变赝势本身**（同一套卡片、
> rlderiv=5 与 6 生成的 UPF 数值块逐字节相同，`ld1.test` 也相同），只改"对数
> 导数在哪儿看"。所以挑 `rlderiv` 是先从"哪套参数自己最一致"里挑一套**好看的
> 视角**，PP 的物理性质不会因此变化。

## `rcloc_scan`：扫 `&inputp` 的 `rcloc`（局域势的芯半径）

`rcloc` 是**局域势**的芯半径（`lloc=-2` 时局域势就是 l=2 那个通道的）：它不属于
任何一张通道卡片，而是 `&inputp` 里的一个标量，所以这个工作流的扫描轴不一样：

```json
{
  "workflow": "rcloc_scan",
  "parameters": {"input": "rcloc_scan", "rcloc_scan_range": "rcloc_scan_range"},
  "qe_path": "/home/liguozhou/oncv/work/qe/q-e-qe-7.4/bin"
}
```

```yaml
rcloc_range: [-0.5, 0.5, 0.25]   # [起点, 终点, 步长]，含两端
update: {"criterion": "pole", "pole_threshold": 0.02}
```

**扫描规则**：

* **扫描中心就是 `&inputp` 里现在的 `rcloc`**（没写就报错 —— 与其猜 `ld1.x` 的
  默认值，不如让你显式写一个中心）；`rcloc_range` 是**相对它的偏移**，与
  `rc_range` / `energy_range` 同一个约定；
* 每个点落在 `runs/rcloc/rcloc_<值>/U.in`，例如 `runs/rcloc/rcloc_2p5/U.in`；
* 每个点只改 `&inputp` 的 `rcloc`，`lloc`、`rclow`（如果写了）与别的段一律不动；
* 它影响**整份赝势**（所有通道），没有"被扫的那个轨道/卡片"这一说 —— 所以判据
  与报告都按**所有通道**来（极点全重合 → 全部极点 ΔE 合计 → 全部通道池化的
  `|ΔL|` 中位数），`ppgen report` 会自动把 `--scope` 当成 `all` 并给一条说明。

```
扫描       : &inputp 的 rcloc × 5 个点 = 5 个扫描点（中心取 &inputp 里 rcloc 现在值 + rcloc_range 偏移，每个值一个目录）
           : rcloc  = 2、2.25、2.5、2.75、3   ->  runs/rcloc/rcloc_*/
```

### 扫完自动更新 rcloc（`update`）

与其他两个扫描工作流同一套判据（逐组串行在这里只有一组：`rcloc` 本身）：

```
更新 rcloc: rcloc → 2.2 bohr（原 2.5）
          : 最佳点 .../runs/rcloc/rcloc_2p2，极点重合，ΔE 合计 0 Ry，|ΔL| 中位 4.8892e-05
          : 已写回 .../parameters/input/rcloc_scan.yml（后面的组用这个新值生成）
```

写回只换 `&inputp` 里 `rcloc:` 那一行的**那个数**：缩进、行尾注释、同段其它键、
文件其余部分一律原样保留；那个键不在就插到该段的第一行。

筛选：`--rc` 在这类工作流里就是"要哪些 rcloc 值"（`--rc 2.5` 只跑那一个点），
`--only` 用不上（只有一组）。

## `energy_scan`：给 `occ=0` 的额外投影算子扫 ener

多投影子（同一个 `l` 挂两张卡片）是"额外的那个投影算子用哪个参考能量"的问题：
第一张卡片的 `occ` 是价电子占据数、`ener` 写 `0.00`（`ld1.x` 取单电子能量），
**第二张卡片的 `occ` 必须是 0**、`ener` 必须自己给 —— 那个 `ener` 就是
`energy_scan` 要找的东西。

```json
{
  "workflow": "energy_scan",
  "parameters": {"input": "energy_scan", "energy_scan_range": "energy_scan_range"},
  "qe_path": "/home/liguozhou/oncv/work/qe/q-e-qe-7.4/bin"
}
```

* `parameters/input/energy_scan.yml`：与 `rc_scan` 同一份结构的参数
  （`&input` / `&inputp` / 通道卡片 / `&test`），卡片里至少要有一张 `occ == 0`
  且 `ener != 0` 的；
* `parameters/energy_scan_range.yml`：给偏移，以及可选的 `update`：

```yaml
energy_range: [-0.2, 0.2, 0.05]   # [起点, 终点, 步长]，含两端
update: true                      # 扫完一张卡片就把它的 ener 更新成最佳值
```

**扫描规则**（一次只动一张卡片）：

* 扫的对象是**所有 `occupancy == 0` 的卡片**（`occ` 是卡片第 4 列）——
  占据态卡片是"被描述的那个态"，不扫；
* **扫描中心就是那张卡片里现在的 `ener`**，`energy_range` 是**相对它的偏移**
  （与 `rc_range` 同一个约定）：5F 那张写 `-0.60`、范围 `[-0.2, 0.2, 0.05]`
  就扫 `-0.80、-0.75 … -0.40`；
* 每个点落在 `runs/<卡片键>/ener_<值>/U.in`，例如 `runs/5F-2/ener_m0p6/U.in`。
  **卡片键**：标签在整套卡片里唯一时就是标签（`5F`），同一个标签有多张
  （多投影子）时带序号 —— 上面那份参数里两张 `5F` 卡片，被扫的那张叫 `5F-2`；
* 只改**那一张**卡片的 `ener`；另一张 `occ=0` 卡片、以及所有占据态卡片
  （`ener`、`rcut`）保持参数里的精确值 —— 一次只动一个变量；
* `ener` 写 `0.0` 的 `occ=0` 卡片会被跳过并提醒（那种卡片 `ld1.x` 会用单电子
  能量，扫它没有意义）。

```
扫描       : 2 张 occ=0 卡片 × 9 个点 = 18 个扫描点（中心取每张 occ=0 卡片现在的 ener + energy_range 偏移，每张卡片一个目录）
           : 6D-2  ener = -0.25、-0.2、-0.15、-0.1、-0.05、0、0.05、0.1、0.15   ->  runs/6D-2/ener_*/
           : 5F-2  ener = -0.8、-0.75、-0.7、-0.65、-0.6、-0.55、-0.5、-0.45、-0.4   ->  runs/5F-2/ener_*/
产物目录   : .../energy_scan/runs
```

### 扫完自动更新 ener（`update`）

与 `rc_scan` 的 `update` 完全同一套（判据、容差、窗口、逐组串行都一样），
差别只有"改哪一列、写回哪一张"：

1. 生成并跑完**这一张卡片**的所有点；
2. 按 **极点是否重合 → 极点 ΔE 合计 → `|ΔL|` 中位数** 挑最佳点（通道取卡片上写的
   `l`，不是从标签猜的，所以 `5F-2` 这种键也能判对）；
3. 把这个 `ener` 就地写回**那张卡片**（按卡片序号定位，同一标签有多张也不会改错行），
   再用新值生成并扫描**下一张**卡片。

```yaml
update: true                                       # 等于 {"criterion": "pole"}
update: {"criterion": "pole", "pole_threshold": 0.02, "window": [-3, 3]}
```

```
更新 ener : 5F-2 → -0.60 Ry（原 -0.55）
          : 最佳点 .../runs/5F-2/ener_m0p6，极点重合，ΔE 合计 0 Ry，|ΔL| 中位 4.8892e-05
          : 已写回 .../parameters/input/energy_scan.yml（后面的组用这个新值生成）
```

> **窗口别开太宽。** 不写 `window` 时用的是 `&input` 的 `eminld`/`emaxld`；如果那
> 是 `-5…10`（15 Ry），高能端一堆跟赝势转移性无关的极点会把"极点重合"这一级顶掉，
> 排序落到 ΔE 合计上，挑出来的点可能毫无意义（实测 6P 曾被挑到
> `ΔE 3.88 Ry / |ΔL| 1.42` 的点上）。所以 `update` 开着而窗口宽过 6 Ry 时
> `run` 会给一条提醒，建议写 `window: [-3, 3]`。

`output.json` 里每个任务写的是 **`"scanned": "ener"`** 与 **`"ener": <值>`**
（`rc_scan` 是 `"scanned": "rcut"` 与 `"rcut"`，老记录没有 `scanned` 时按 `rcut` 读），
`updates` 里也是同一个键名：

```json
{"label": "5F-2 ener=-0.6", "orbital": "5F-2", "scanned": "ener", "ener": -0.6,
 "directory": ".../runs/5F-2/ener_m0p6", "exit_code": 0, "exit_reason": "正常结束"}
```

`ppgen report` 的表头、单位与"选择依据图"的横轴都跟着这个轴走
（`ener (Ry)` / `why this ener: pole ΔE vs ener`），不会写成 `r_c (bohr)`。

筛选开关：`--only` 用**卡片键**（`--only 5F-2`，也可以写标签一次性选同标签的所有
可扫卡片），`--ener` 用**真实 ener 值**（`--ener=-0.6`）：

```bash
ppgen run -i .../energy_scan/input.json --only 5F-2            # 只扫 5F 那张
ppgen run -i .../energy_scan/input.json --only 5F-2 --ener=-0.6  # 只跑那一个点
```

## `run` 的几条约定

* **给 `-i` 就会重新生成输入。** 那一套的语义是"按 `input.json` 生成文件并运行"，
  每次都覆盖同名 `.in`；要"原样跑手改过的那份"，就进那个目录直接用
  `ppgen run`（或 `ppgen run -i U.in`）—— 两种流程各管一件事。
* **cwd 是 run 目录。** `ld1.x` 的产物（`ld1.wfc`、`vx.pot`、`CRASH`…）落在
  `input.json` 旁边，不污染你当前所在的目录。
* **stdout → `U.out`，stderr → `U.err`**（并回显到终端）。OpenMPI 那句
  `librdmacm.so.1` 噪声会照样出现在终端，但也在 `U.err` 里留了底。
* **失败要看得见。** 退出码非 0 时，直接把 `U.out` 里 `Error in routine …`
  那几行与末尾若干行摆出来，并列出本次新出现的文件；退出码原样返回给 shell。
* **坐标范围**：报告里每根轴都标了是"给定"还是"自动按数据"（`x -5–5 Ry（给定），
  y 0.0003074–1.069e+04（对数轴）（自动按数据）`），所以删掉 `xlim`/`ylim` 之后
  仍然知道实际画的是哪一段。
* **产物清单是算出来的**：跑前跑后对比目录，列出新出现/被改写的文件及其
  「覆盖 / 新建」，所以不用猜"这次到底产出了什么、有没有覆盖上次的"。

## 三个输入是怎么对上的

```
run 目录/input.json                    pp-gen/src/pp_gen/parameters/
────────────────────                   ─────────────────────────────
workflow   : "ld1.x"       ──────┐
parameters : {"input":        │      （工作流登记表：ae_radial 只要 "input" 组）
              "ae_radial"} ───┼───>  input/ae_radial.yml ──┐
qe_path    : "xx"             │                           │
                              └─> pp_gen/workflows.py      └─> &input 段
                                                              ↓
                                             run 目录/U.in（<元素>.in，无注释）
```

* **`workflow`**：在 `pp_gen/workflows.py` 里登记的工作流名。现在只有两个：
  **`ld1.x`**（一次一份输入，只要 `input` 组）与 **`rc_scan`**（扫一批，要
  `input` + `rc_scan_range` 两组）。名字写错、或多配了工作流不用的组，都会直接报错
  —— 免得渲染出 `ld1.x` 不该读的段。
* **`parameters`**：`组名 -> YAML 名`。`{"input": "ae_radial"}` 就是
  `parameters/input/ae_radial.yml`（先找 `parameters/<组名>/<名>.yml`，再找
  `parameters/<名>.yml`）。YAML 的顶层键是 namelist 段名，段里的键值就是变量。
* **文件名**：由 `&input` 的 `atom` 定（只写了 `zed` 时由原子序数反推，pslibrary
  那批输入就是这么写的），所以 `atom: 'U'` → `U.in`。

## `qe_path`

`qe_path` 写 QE 的 **bin 目录**（推荐，`ld1.x` / `pw.x` 都在里面），也可以直接
写成 `ld1.x` 这个可执行文件本身。

`input.json` 里是占位符（`"xx"`、空、`TODO`…）时，`prepare` 会补上默认值
**并写回 `input.json`**（其它键原样保留）；已经写了真路径就尊重它，不动。

默认值 / 覆盖方式：

| | 值 |
|---|---|
| 默认 QE 路径 | `/home/liguozhou/oncv/work/qe/q-e-qe-7.4/bin`（`pp_gen/config.py` 里的 `DEFAULT_QE_PATH`） |
| 环境变量覆盖 | `PPGEN_QE_PATH`（QE 路径）、`PPGEN_PARAMETERS`（参数 YAML 根目录） |

本机默认路径下的 `ld1.x` 已确认存在且可执行；找不到时 `prepare` 只警告、不中断
（生成输入文件本身不需要 QE）。

## 安装与测试

```bash
python3 -m pip install -e /home/liguozhou/oncv/pp-gen     # 装出 ppgen 命令
python3 -m pip install -e '.[plot]'                       # 画图还需要 matplotlib
cd /home/liguozhou/oncv/pp-gen && python3 -m pytest -q    # 376 项
python3 -m pytest -q -m "not slow"                        # 373 项，不碰真 QE
```

没装也能跑：`PYTHONPATH=src python3 -m pp_gen prepare -i .../input.json`。

`run` 的测试大多用一个**假的 `ld1.x`**（读 stdin、写 stdout、在 cwd 落一个
文件的 bash 脚本），这样退出码、stderr、cwd 都能精确断言；`plot` 的测试直接
检查 matplotlib 的 `Axes` 状态（xlim、轴标题、线数、图例、线型/颜色），不去
比对 PNG 像素。标了 `slow` 的两项才真调本机 QE。

## 验证

用生成出来的 `U.in` 真跑了本机的 `ld1.x`（QE 7.4），实际就是
`ppgen run -i .../ae_radial/input.json`。

**与参考 run 逐字节一致**：当 `U.in` 与参考文件 `U.ae.iswitch1.in` 用同一套
网格时（`xmin=-8.0` / `dx=0.008` / `rmax=100`，即 ld1.x 对 Z≥19 的自动值），
产物 `ld1.wfc`、`vx.pot` 与
`ppgen/script/runs/ae_radial/` 里那份**逐字节相同**，`U.out` 只差时间戳、
可用内存和 title 那一行文字。

> `ae_radial.yml` 后来显式写了 `dx=0.005`（细网格），所以现在再跑出来的
> `ld1.wfc` 是 551902 B / 3427 点，而不是 344856 B / 2141 点。
> **这是网格变了，不是"没覆盖"** —— 覆盖与否看 `run` 报告里产物的
> 「覆盖 / 新建」标注，以及文件的 mtime。

失败路径也拿真 QE 验过：把 `U.in` 里塞一个不存在的变量名，`ld1.x` 退出码 1，
`ppgen run` 报出 `Error in routine ld1_readin (5010): reading input namelist`、
`U.out` 末尾，并把 QE 留下的 `CRASH`、`input_tmp.in` 一并列了出来。

`plot` 的结果与你早先那份分析逐一对照过（`ld1.wfc` -> 峰值表）：

* 七个态的最大 `|rR|` 与 `ae_radial_max_per_orbital.csv` **完全一致**
  （4F 1.7585、5D 1.2444、5F 0.8692、6S 1.0017、6P 0.8810、6D 0.5929、7S 0.5720），
  说明读文件与单位换算没错；
* 峰值位置（Å）与之相差 ≤ 0.004 Å（差别来自 `dx=0.005` 与参考 `dx=0.008`
  的网格不同）；
* 出图后再逐像素核对：七个态用的七种线色都真的画在画布上（`ae_radial.png`
  1260×864 = 7.0×4.8 in @180 dpi），没有哪条曲线因为范围设错而消失。
* **极值表与参考文件逐行对照**：在参考 run 那份 `ld1.wfc` 上，本工具算出的
  28 行极值与 `ae_radial_extrema.csv` **行集合完全相同**，`r_bohr`/`r_angstrom`/
  `rR` 的最大偏差 4e-11（就是 10 位小数的输出精度），`n_nodes` 与 `max/min`
  类型也逐个相同 —— 说明节点/极值算法与之前那份分析是一套。换成你的
  `dx=0.005` 网格后，行数、`n_nodes`、类型仍然一一对应，位置差 ≤ 7e-3 bohr。
* **`rc_scan` 真跑过一整轮**：55 个扫描点各出一份 `U.pbe-n-nc.UPF`
  （以及 `ld1.dlog`/`ld1.pwe`/`ld1.wfc`/`ld1ps.*`/`vx.pot`），日志以
  `End of Bessel function test` 收尾。卡片渲染还对着 QE 的
  `el_config.f90` 的 `read_psconfig` 核过：`rel<2` 只读 7 列
  （`rcut` 后紧跟 `rcutus`，且要求 `rcut <= rcutus`），第 8 列是可选
  `isws`/`jj`，多写的列会被 list-directed 读法在记录末尾丢掉 —— 所以本工具
  按 8 列写是安全的（与参考的 `U.gen.in` 一致），并且在渲染阶段就替 QE 把
  "标签与 l 不符、`n < l+1`、占据数超界、`rcut > rcutus`"这几类错拦住。
* **对数导数图的极点诊断与参考脚本逐通道一致**：把
  `rc_scan_logderiv.py` 的 `pole_report` 原样加载进来，对同一个点的 4 个通道比对
  `n_ae/n_ps/matched/shifted/extra`，三种 `r_c`（0.8、1.1、1.4）**全部相同**。
* **`rcloc_scan` 拿真 QE 跑过一遍**（5 个点扫 `&inputp` 的 `rcloc`：2.0、2.25、2.5、
  2.75、3.0，每个 1.5 s，5/5 成功）：`runs/rcloc/rcloc_2p5/U.in` 里只有 `rcloc`
  变了（`lloc`、卡片、`&test` 都没动）；判据按**全部通道**算，`update` 挑出
  **2.25**（极点重合、ΔE 合计 0.06 Ry、池化 `|ΔL|` 中位 1.4e-3）并写回
  `&inputp` 那一行；`report.md` 的表头是 `最佳 rcloc (bohr)`、通道列写 `l=全部`、
  并自动说明"按所有通道池化"。
* **`|ΔE|` 那一级修过一个洞**：`pairs` 是按 **PS 侧**的极点配的，AE 侧多出来的
  极点以前根本不出现在里面 —— 于是"PS 在窗口里连极点都没有"的点 `ΔE 合计` 算成
  `0`（看着最好），排序会挑中它。现在**两边个数不一样就记 `inf`**（报告里显示
  "配不上"）。上面那次 `rcloc=2.5 / 3.0` 就是这个情形：`|ΔL|` 中位最小
  （1.2e-3 / 4.6e-3）但 AE/PS 极点个数对不上，`--pole` 不再把它们排到前面。
* **`energy_scan` 拿真 QE 跑过一遍**（9 个点扫 `5F` 那张 `occ=0` 卡片的 ener，
  `-0.80…-0.40`，每个点 1.9 s，9/9 成功）：生成的 `U.in` 只改了那一张卡片的
  `ener`（第 5 张占据的 `5F` 仍是 `0.00`、另一张 `occ=0` 的 `6D` 仍是 `-0.05`）；
  `update: true` 把挑出来的最佳值写回了**第 7 张**卡片（同标签的另一张没动）；
  `output.json` 记的是 `"scanned": "ener"` / `"ener"`，`report.md` 的表头是
  `最佳 ener (Ry)`、结论行是 `ener = -0.60 Ry`，选择依据图的横轴也是 `ener (Ry)`。
* 对数横轴（`"xscale": "log"`）单独核过：`ax.get_xscale() == "log"`、
  主刻度落在 `0.01 / 0.1 / 1 / 10` 这四个十进制位上、纵轴仍是线性 `[-2, 2]`；
  曲线横向铺满 85% 的画布宽度，且这段宽度里**每一列**都有曲线经过
  （说明不是只有少数几个点被画上去了）。

## 目录

```
pp-gen/
├── pyproject.toml              # console script: ppgen = pp_gen.cli:main
├── src/pp_gen/
│   ├── cli.py                  # 命令行：ppgen prepare / run / plot
│   ├── prepare.py              # 用例：input.json -> <元素>.in
│   ├── runner.py               # 用例：跑 ld1.x，报退出码/耗时/产物
│   ├── plotter.py              # 用例：plot.json + ld1.wfc / dlog -> PNG
│   ├── dlog.py                 # ld1.dlog 读入、极点判据、AE/PS 三分类
│   ├── plot_spec.py            # plot.json 的读入/校验/推断（ylabel -> 量）
│   ├── figure.py               # 只管画：配色、线型、图例、落盘
│   ├── extrema.py              # 节点/极值 + 极值表 CSV（排序与格式）
│   ├── scan.py                 # rc 扫描：扫描点算术 + 把 rcut 打进卡片
│   ├── update.py               # 扫完一个轨道就挑最佳 rcut 并写回参数（update: true）
│   ├── report.py               # 用例：扫描结果 -> report.md（最佳 r_c + 极点情况）
│   ├── wfc.py                  # ld1.wfc 的读入与派生量（rR / R / 概率密度）
│   ├── run_input.py            # input.json 的读入/校验/写回（保序、留未知键）
│   ├── parameters.py           # 参数 YAML 的定位、读入、合并
│   ├── ld1_input.py            # 渲染 ld1.x 输入（无注释）
│   ├── workflows.py            # 工作流登记表（加新工作流只改这里）
│   ├── elements.py             # zed -> 元素符号（只为定文件名）
│   ├── config.py               # QE 默认路径 / 参数目录 / 超时 / 玻尔半径
│   └── parameters/             # 参数 YAML：input/、rc_scan_range.yml …
└── test/                       # pytest：渲染、run、plot、report、出错路径、CLI
```
