# Scope and claims (snapshot, 30 October 2026)

The project intends hydrated clusters generated and characterized by GWH, xTB/ALPB, CoCoS, and MDAtomLocalEnvironment. **A case not visible in Drive is NOT proof it is unclaimed or has no calculations elsewhere.** Get the author's assignment first.

- **Rb+ (H2O)4–14:** 11 stoichiometries with source data folders; Rb explicitly claimed and calculations advanced. Do not duplicate.
- **Li+ (H2O)4–6:** directories with GWH output visible; Li n4–6 not in active package.
- **K+ (H2O)4:** GWH output on Drive; **K n5:** folder exists but was empty on this listing, so held because K was claimed by a participant. Do not assume availability.
- **Na+ (H2O)4–14:** no Na results directory visible, BUT participant said they would run sodium next day. Treat as verbally claimed; do not duplicate without confirmation.
- **Li+ (H2O)7–14 and K+ (H2O)6–14:** no active folders visible or claims in supplied chat, provisionally prepared in this kit. Explicit allocation pending.
- **F-, Cl-, Br-:** original proposal considered these alternatives. Author has not finalized inclusion. User-provided vecfle.zip contains no anion inputs; anion-water geometry/orientation and charge -1 require separate validated protocol. No fabricated anion inputs included.

**Operational gate (2026-10-31):** `CASE_STATUS.csv` column `assignment_status` is `ON_HOLD` for `K_06` and `NOT_ASSIGNED` for the other Li/K rows in this kit. Production chemistry (`workflow.py gwh` … `dedup`, `RUN_ONE_APPROVED_CASE.sh`) refuses `ON_HOLD` unless `DENIS_ASSIGNED=yes`. `ALLOW_FULL_BATCH=yes` is not a bypass. No Na or Rb rows were added: those ions are supported in code (charge +1) but packing files and assignments are not in this kit (Na verbally claimed; Rb already in progress elsewhere). F/Cl/Br are supported in code (charge −1) with **no** CASE_STATUS rows and **no** packing files. Re-audit of kit docs and GWH vendor git on 2026-10-31 found **no confirmed Denis assignment**; statuses were left unchanged. `workflow.py process-existing` only post-processes the existing K₆ pilot. See `ASSIGNMENT_AUDIT.md` and `MESSAGE_TO_DENIS.md`.

**All states reflect accessible Google Drive and supplied chat at snapshot time, not an authoritative real-time task register.**
