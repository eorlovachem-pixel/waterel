# Scope and claims (snapshot, 30 October 2026)

The project intends hydrated clusters generated and characterized by GWH, xTB/ALPB, CoCoS, and MDAtomLocalEnvironment. **A case not visible in Drive is NOT proof it is unclaimed or has no calculations elsewhere.** Get the author's assignment first.

- **Rb+ (H2O)4–14:** 11 stoichiometries with source data folders; Rb explicitly claimed and calculations advanced. Do not duplicate.
- **Li+ (H2O)4–6:** directories with GWH output visible; Li n4–6 not in active package.
- **K+ (H2O)4:** GWH output on Drive; **K n5:** folder exists but was empty on this listing, so held because K was claimed by a participant. Do not assume availability.
- **Na+ (H2O)4–14:** no Na results directory visible, BUT participant said they would run sodium next day. Treat as verbally claimed; do not duplicate without confirmation.
- **Li+ (H2O)7–14 and K+ (H2O)6–14:** no active folders visible or claims in supplied chat, provisionally prepared in this kit. Explicit allocation pending.
- **F-:** original proposal considered this alternative. Author has not finalized inclusion. No F packing files.
- **Cl- and Br- (H2O)4–14:** technical packings added 2026-11-01 from mendeleev 1.3.0 ionic radii (X− VI + O²− VI, pm/100). Direction sets for n=6–14 copy author K unit vectors; n=4 tetrahedron; n=5 trigonal bipyramid. `CASE_STATUS` rows are **NOT_ASSIGNED**. Production GWH `--Opt` / 100-start xTB is **not** launched. F remains absent. Do not treat this prep as Denis assignment.

**Operational gate (2026-11-01):** `CASE_STATUS.csv` is `ON_HOLD` for `K_06` and `NOT_ASSIGNED` for Li_07–14, K_07–14, Cl_04–14, and Br_04–14. Production chemistry (`workflow.py gwh` … `dedup`, `RUN_ONE_APPROVED_CASE.sh`) refuses both `ON_HOLD` and `NOT_ASSIGNED` unless `DENIS_ASSIGNED=yes`. `ALLOW_FULL_BATCH=yes` is not a bypass. No Na, Rb, or F rows: Na verbally claimed; Rb already in progress elsewhere; F inclusion not prepared. `DENIS_ASSIGNED=yes` was **not** set. `workflow.py process-existing` only post-processes the existing K₆ pilot. See `ASSIGNMENT_AUDIT.md` and `MESSAGE_TO_DENIS.md`.

**All states reflect accessible Google Drive and supplied chat at snapshot time, not an authoritative real-time task register.**
