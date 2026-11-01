# Assignment and Euler audit (2026-10-31)

This is a record of what is **in this kit**. It is not an assignment. `ON_HOLD` / `NOT_ASSIGNED` were **not** changed. `DENIS_ASSIGNED=yes` was **not** set. No production series was started.

## Assignments

| Source | Finding |
| --- | --- |
| `CASE_STATUS.csv` | `K_06` = `ON_HOLD`. Li_07–14, K_04–05, K_07–14, Cl_04–14, Br_04–14 = `NOT_ASSIGNED`. No `ASSIGNED` rows. No Na/Rb/F rows. |
| `CASE_STATUS.json` | Regenerated from the CSV so both files match. CSV remains the source of truth. |
| `ON_HOLD_AND_SCOPE.md` | Drive/chat snapshot 2026-10-30: Li 7–14 and K 6–14 “explicit allocation pending”. Rb claimed. Na verbally claimed. |
| `CHECKS_PERFORMED.md` | “Соавторство и распределение задач не согласованы”. |
| Kit git | This directory is **not** a git repository. No commit assigns K_06 to this participant. |
| GWH vendor git | HEAD `7411811374b6422a15f23d6d86a47f2145f8ac47` (Denis Tikhonov, 2026-10-30) is the Euler-string fix, not a case assignment. |

**Conclusion:** there is no confirmed Denis assignment of K_06 (or of Li_07–14 / K_07–14) in the files or vendor history available here.

Na⁺, Rb⁺, and F⁻ stay out of `CASE_STATUS.csv`. Cl/Br and K_04/K_05 rows are technical prep (`NOT_ASSIGNED`). Production launch refuses `NOT_ASSIGNED` as well as `ON_HOLD`.

## Euler: SciPy `zyz` vs `ZYZ`

Author code (`RotateMolecule`, `gloriouswaterhole.py` line 94):

```python
M = spst.Rotation.from_euler('zyz', [phi, theta, chi], degrees=True).as_matrix()
```

SciPy `Rotation.from_euler` docs: uppercase `ZYZ` = **intrinsic** (body-fixed); lowercase `zyz` = **extrinsic** (space-fixed). The two matrices for the same angles are different; `zyz(a,b,c)` equals `ZYZ(c,b,a)`.

| Claim | Status |
| --- | --- |
| Author wrote `'zyz'` on purpose in commit 7411811 | Fact |
| SciPy therefore applies **extrinsic** ZYZ | Fact (docs + matrix check in `tests/test_euler_zyz.py`) |
| `'ZYZ'` would be SciPy-intrinsic ZYZ | Fact |
| Whether Denis *meant* Wikipedia-style intrinsic ZYZ | **Open** — ask him |

**Recommendation (do not apply until Denis answers):** keep `from_euler('zyz')` as the production convention. Reasons: it is the author commit; changing to `'ZYZ'` would change every constrained GWH geometry; Li_07 already used old `'zyx'` and must not be mixed with a third convention; the kit previously mislabelled lowercase as “intrinsic”, but the source was never `'ZYZ'`.

Draft questions: `MESSAGE_TO_DENIS.md`.
