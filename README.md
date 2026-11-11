# MVkSP-SDP: Minimum Variance *k*-Subset Selection

**MVkSP-SDP** is branch-and-cut algorithm for the minimum variance *k*-subset problem (MVkSP). Given $n$ points in Euclidean space and a cardinality $k$, the problem selects exactly $k$ points that minimize their average squared distance from their centroid.

This repository contains the C++ and MATLAB source code associated with the paper ["Optimality Guarantees for Minimum Variance Subset Selection via Semidefinite Programming"](http://arxiv.org/abs/2611.13220) by Antonio M. Sudoso and Daniel Aloise.

## Repository structure

```text
.
├── cpp/
│   ├── main.cpp                  # Command-line entry point and data reader
│   ├── sdp_branch_and_bound.*    # Branch-and-cut algorithm
│   ├── ThreadPool.*              # Parallel node processing
│   ├── JobQueue.*                # Node queue and visiting strategies
│   ├── matlab_util.*             # MATLAB Engine interface
│   ├── config.txt                # Runtime configuration
│   ├── Makefile                  # C++ build instructions
├── matlab/
│   ├── solve_MVkSP_sdpnalplus.m  # SDP cutting-plane procedure
│   ├── build_MVkSP_*.m           # Original and facially reduced SDP models
│   ├── separate_RLT_triangle*.m  # Valid-inequality separation
│   ├── safe_bound_*.m            # Certified lower-bound post-processing
│   └── vns_*.m                   # Primal heuristics
├── paper.pdf
└── README.md
```

## Requirements

- A Linux C++ toolchain with GNU Make
- [Armadillo](https://arma.sourceforge.net/)
- MATLAB with the [MATLAB Engine API for C++](https://www.mathworks.com/help/matlab/calling-matlab-engine-from-cpp-programs.html)
- [SDPNAL+](https://blog.nus.edu.sg/mattohkc/softwares/sdpnalplus/)
- [Gurobi](https://www.gurobi.com/) with its MATLAB interface and a valid license

The implementation reported in the paper was tested on Ubuntu 22.04.3 LTS with MATLAB R2022b and Gurobi 13.0.

## Installation

1. Clone the repository:

   ```bash
   git clone https://github.com/antoniosudoso/MVkSP-SDP.git
   cd MVkSP-SDP
   ```

2. Install the C++ build tools and Armadillo on Ubuntu or Debian:

   ```bash
   sudo apt-get update
   sudo apt-get install build-essential libarmadillo-dev
   ```

3. Install MATLAB, Gurobi, and SDPNAL+.

4. Set `MATLAB_PATH` and `GUROBI_PATH` in `cpp/Makefile`:

   ```makefile
   MATLAB_PATH = /path/to/MATLAB/R2022b
   GUROBI_PATH = /path/to/gurobi1300
   ```

5. Set the following absolute paths in `cpp/config.txt`:

   ```text
   MATLAB_SOURCE_FOLDER = /path/to/MVkSP-SDP/matlab
   GUROBI_FOLDER = /path/to/gurobi1300
   SDPNAL_FOLDER = /path/to/SDPNAL+
   ```

   These directories and their subdirectories are added automatically to each MATLAB session started by the solver.

6. Compile the executable:

   ```bash
   cd cpp
   make
   ```

   This creates the executable `cpp/bb`. 

## Configuration

Runtime parameters are defined in `cpp/config.txt`.

### Branch-and-bound parameters

| Parameter | Description |
|---|---|
| `BB_TOL` | Relative optimality gap tolerance |
| `BB_THREADS` | Number of parallel worker threads |
| `BB_MAX_NODES` | Maximum number of processed nodes |
| `BB_TIME_LIMIT` | Time limit in seconds |
| `BB_VISITING_STRATEGY` | Node strategy: best-first (`0`), depth-first (`1`), or breadth-first (`2`) |

### MATLAB and solver parameters

| Parameter | Description |
|---|---|
| `MATLAB_SESSION_THREADS_ROOT` | MATLAB computational threads used at the root node |
| `MATLAB_SESSION_THREADS_CHILD` | MATLAB computational threads used at the child node |
| `MATLAB_SOURCE_FOLDER` | Absolute path to this repository's `matlab` directory |
| `GUROBI_FOLDER` | Absolute path to the Gurobi installation |
| `SDPNAL_FOLDER` | Absolute path to the SDPNAL+ installation |
| `SDP_FR` | Use the facially reduced SDP formulation (`1`) or the original formulation (`0`) |
| `SDP_SOLVER_TOL` | SDPNAL+ relative KKT-residual tolerance |
| `SDP_SOLVER_VERBOSE` | Disable (`0`) or enable (`1`) SDPNAL+ output |

### Cutting-plane parameters

| Parameter | Description |
|---|---|
| `CP_MAX_ITER` | Maximum number of cutting-plane iterations at a node |
| `CP_TOL` | Convergence tolerance used by the cutting-plane loop |
| `CP_MAX_INEQ` | Maximum number of candidate inequalities examined per iteration |
| `CP_NEW_INEQ` | Maximum number of violated inequalities added per iteration |
| `CP_EPS_INEQ` | Violation tolerance for separated inequalities |
| `CP_EPS_ACTIVE` | Tolerance used to retain active inequalities |

## Usage

From the `cpp` directory, run:

```bash
./bb <CONFIG_PATH> <DATA_PATH> <K> <LOG_PATH>
```

For example:

```bash
mkdir -p logs_BB
./bb ./config.txt /path/to/dataset.txt 10 ./logs_BB/dataset_k10.txt
```

The arguments are:

- `CONFIG_PATH`: path to the runtime configuration file;
- `DATA_PATH`: path to the input data matrix;
- `K`: number of points to select, with $1 < k < n$; and
- `LOG_PATH`: path of the log file to create.

To run several values of $k$, edit `INSTANCE` and the list of cardinalities in `cpp/run.sh`, ensure the data and output directories exist, and execute:

```bash
./run.sh
```

## Input format

Benchmark datasets are not included in this repository. Supply each instance as a whitespace-delimited text file. Its first line contains the number of points $n$ and the dimension $d$, followed by one point per line:

```text
n d
p_11 p_12 ... p_1d
p_21 p_22 ... p_2d
...
p_n1 p_n2 ... p_nd
```

## Output

The solver writes node-level progress and a final summary to the requested log file. The progress table contains:

| Column | Description |
|---|---|
| `N` | Number of unfixed points in the current node relaxation |
| `PAR ID` | Parent-node identifier |
| `ID` | Current-node identifier |
| `TYPE` | Root node (`-1`) or child created by fixing the branching variable to `0` or `1` |
| `BR IDX` | Zero-based index selected for the next branching decision; `-1` if no branching is needed |
| `PAR LB` | Parent-node lower bound |
| `LB` | Current-node valid lower bound |
| `ITER` | Number of cutting-plane iterations |
| `INEQ` | Number of retained valid inequalities |
| `TIME (s)` | Node solution time in seconds |
| `UB` | Best feasible objective value found at the node |
| `GUB` | Global incumbent objective value |
| `NGAP (%)` | Current-node relative gap |
| `GGAP (%)` | Global relative gap |
| `QUEUE` | Number of open nodes |

The final summary reports the root statistics, best objective value, final optimality gap, number of processed nodes, total running time, and the indices of the selected points.

## Citation

If you use this code, please cite:

> A. M. Sudoso and D. Aloise, "Optimality Guarantees for Minimum Variance Subset Selection via Semidefinite Programming" arXiv, 2026.

```bibtex
@article{sudoso2026optimality,
  title={Optimality Guarantees for Minimum Variance Subset Selection via Semidefinite Programming},
  author={Sudoso, Antonio M. and Aloise, Daniel},
  journal={arXiv preprint arXiv:2611.13220},
  year={2026},
  doi={10.48550/arXiv.2611.13220}
}
```
