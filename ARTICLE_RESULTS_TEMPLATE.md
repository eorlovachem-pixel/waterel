# Section draft template — cannot be treated as completed science before computation

## Computational protocol
The starting arrangements of n water molecules around each monatomic cation were obtained from normalized sphere-packing directions, with initial ion–oxygen distances of 2.42 Å (Li+) or 3.04 Å (K+) as specified by the author-provided input generation script. In each fixed-stoichiometry batch, 100 randomly oriented candidate clusters were subjected to constrained orientational optimization in the GWH workflow using GFN2-xTB with ALPB(water). All successfully converged constrained configurations were subsequently optimized without the initial geometric restrictions at the same electronic-structure/solvation level. Importantly, no geometries were discarded on the basis of CoCoS equivalence or a relative-energy window before free optimization. After unconstrained relaxation, conformer equivalence was evaluated with CoCoS, and representative configurations were ranked by final relative energy within each stoichiometry. The local solvent environment and water orientations were then quantified with MDAtomLocalEnvironment. Exact executable versions, scripts, input parameters and the convergence statuses of all calculations are archived.

**Check on completion:** verify with actual programs and logs, especially author-confirmed ZYZ modification and CoCoS version, before putting this wording in a manuscript. *100* is the planned target, not achieved data yet.

## Results — TO BE FILLED WITH COMPUTED DATA ONLY
- Cluster: [Li+(H2O)7 or other case]
- Successful constrained calculations: [number/100]
- Successful full optimization: [number/100]
- Unique conformers after full optimization: [number; method parameters]
- Lowest electronic energy at GFN2-xTB/ALPB(water): [Eh]
- Range of relative energies (within n): [kcal/mol]
- Lowest-energy representative: [input ID, XYZ filename]
- Real CoCoS sensitivity tests: [count and survivor identity by settings]
- Distances and local water-orientation metrics: [values and method]
- Figure/Table IDs: [IDs]

## Cautions
This exploratory search does not demonstrate global-optimum completeness. Clustering depends on the CoCoS equivalence criterion and alignment strategy. Raw energy minima do not constitute standard hydration free energies; local water descriptors are not necessarily integer first-shell coordination numbers.
