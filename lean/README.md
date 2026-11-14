# lean: kernel-checked algebraic steps

- `K3Chain.lean` (13 statements): step S1 of Proposition C0, the Route 2 trace-to-eigenvalue bounds and the symmetric
  sandwich, end to end; two statements are refutations kept as controls (a sign-flipped S1 and an exact form of a
  referee's square-root estimate are false).
- `K3Chain2.lean` (5 statements): the sandwich bounds and the trace lower bounds in the form the chain uses.

Lean 4 with Mathlib (`lean-toolchain`, `lakefile.toml`, `lake-manifest.json`). Check with `lake env lean K3Chain.lean`
inside a Lake project with Mathlib. Both files pass with standard axioms only (propext, Classical.choice, Quot.sound)
and no sorry. Lean checks the statements as written; that they match the steps of the chain is our reading.
