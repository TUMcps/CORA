These is a test benchmark category. The folder contains .onnx and .vnnlib files used for the category. The instance.csv containts the full list of benchmark instances, one per line: onnx_file,vnn_lib_file,timeout_secs
 
The test properties correspond to ACAS Xu networks (from "Reluplex: An efficient SMT solver for verifying deep neural networks"), networks 1-6 (unsat) and 1-7 (sat). The properties correspond to property 3:

Description: If the intruder is directly ahead and is moving towards the
ownship, the score for COC will not be minimal.

Input constraints: 1500 <= rho <= 1800, −0.06 <= theta <= 0.06, psi >= 3.10, vown >= 980, vint >= 960.

Computing input bounds logic (python3 code):

```
init_lb = [1500, -0.06, 3.1, 980, 960]
init_ub = [1800, 0.06, 3.1415926535, 1200, 1200]

means_for_scaling = [19791.091, 0.0, 0.0, 650.0, 600.0, 7.5188840201005975]
range_for_scaling = [60261.0, 6.28318530718, 6.28318530718, 1100.0, 1200.0]

for i in range(len(init_lb)):
	print(f"\n; Unscaled Input {i}: {init_lb[i], init_ub[i]}")
	lb = (init_lb[i] - means_for_scaling[i]) / range_for_scaling[i]
	ub = (init_ub[i] - means_for_scaling[i]) / range_for_scaling[i]
	print(f"(assert (<= X_{i} {ub}))")
	print(f"(assert (>= X_{i} {lb}))")
```


Result:

```
; Unscaled Input 0: (1500, 1800)
(assert (<= X_0 -0.29855281193475053))
(assert (>= X_0 -0.30353115613746867))

; Unscaled Input 1: (-0.06, 0.06)
(assert (<= X_1 0.009549296585513092))
(assert (>= X_1 -0.009549296585513092))

; Unscaled Input 2: (3.1, 3.1415926535)
(assert (<= X_2 0.49999999998567607))
(assert (>= X_2 0.4933803235848431))

; Unscaled Input 3: (980, 1200)
(assert (<= X_3 0.5))
(assert (>= X_3 0.3))

; Unscaled Input 4: (960, 1200)
(assert (<= X_4 0.5))
(assert (>= X_4 0.3))
```

Desired output property: the score for COC (the first output) is not the minimal score.
