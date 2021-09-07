# ina-model

*TODO: Write short description*

---
## Usage

### Activate environment
```shell
# in bash shell
cd PATH/TO/ina-model
julia

# in julia shell
# enter the Pkg REPL by pressing ] from the Julia REPL
(@v1.6) pkg> activate .
```
---
## Keep in mind
- **do not change** existing notebooks, **create new** based on old
- **clear** all unnecessary **cells outputs** in notebooks
---
## Troubleshooting

### Arpack

Problem:
```shell
ERROR: Error building `Arpack`
```
Solution:
```shell
sudo apt-get install libopenblas-dev
sudo ln -s /usr/lib/libopenblas.so /usr/lib/libopenblas64_.so.0
```
Refs:
https://discourse.julialang.org/t/resolved-arpack-deps-jl-linux-build-issue/30923
https://github.com/JuliaLinearAlgebra/Arpack.jl/issues/67#issuecomment-498692209
