#!/usr/bin/env python3
"""Do the Tier-1 checks have teeth? Seeded-mutation test for audit.py.

    python3 -I provenance/test_checks.py .

Why this file exists. `audit.py` discarded `main()`'s return value: C1-C9
printed their findings and the exit code came only from the C10/C11 delegates.
A reviewer seeded eight simultaneous violations into a copy of the corpus -- a
disconnected counter-ion, a stale Molecular.Weight, ratios summing to 109,
Experiment_value='nan', a wrong Model, a wrong Route, a blanked IL_name, an
exact duplicate row -- and all eight printed while the command exited 0. Of 43
mutations, 30 appeared in the report and only 5 made anything exit non-zero.
A check that prints into a report a human skims is not a check (F-012).

So: mutate ONE field in ONE row of a COPY, run the real audit, and require a
non-zero exit. A control run with no mutation must still pass, or the suite is
only measuring that the audit is broken in general.

This file is the ratchet. Each new check gets a mutation here, and each
mutation that exits 0 is a coverage gap stated as a failing test rather than as
a note someone might read. Its first run found two: gating only
`blank_Publication_PMID` let a blanked `IL_name` through, and C8 checked 7 of
the 14 published descriptor columns, so a LogP of 42 passed. Both fixed; both
mutations kept, so they cannot come back.

NOT a substitute for the baseline file. This proves the gate fires; audit_baseline.json
says what it fires ON. The two fail in opposite directions: a baseline with no
teeth passes everything, and teeth with no baseline fail every clean run.

Copies the repo to a temp directory and mutates that. It never touches data/.
"""
import csv, os, shutil, subprocess, sys, tempfile

REPO = os.path.abspath(sys.argv[1] if len(sys.argv) > 1
                       else os.path.join(os.path.dirname(__file__), ".."))


def _epimerise(r):
    """Invert the steroid C3 stereocentre, reproducing F-013 exactly.

    Done with RDKit rather than a string substitution: the canonical SMILES
    writes that centre differently depending on the molecule, so a literal
    `[C@@H](OC` -> `[C@H](OC` swap silently matched nothing and the "mutation"
    left the row untouched. A seeded-mutation test whose mutation does not
    apply reports the CHECK as broken, which is the worst possible direction
    for this suite to fail in -- so this one verifies that it changed the
    string, and the row is skipped (leaving the control's clean result) only if
    there is no steroid ester to invert.
    """
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
    m = Chem.MolFromSmiles(r["IL_SMILES"])
    if m is None:
        return
    ri = m.GetRingInfo()
    for a in m.GetAtoms():
        if (a.GetSymbol() != "C" or not a.IsInRing()
                or ri.NumAtomRings(a.GetIdx()) != 1
                or a.GetChiralTag() == Chem.ChiralType.CHI_UNSPECIFIED):
            continue
        # C3: a ring carbon whose exocyclic oxygen is an ester oxygen
        for nb in a.GetNeighbors():
            if nb.GetSymbol() != "O" or nb.IsInRing() or nb.GetDegree() != 2:
                continue
            a.SetChiralTag(
                Chem.ChiralType.CHI_TETRAHEDRAL_CCW
                if a.GetChiralTag() == Chem.ChiralType.CHI_TETRAHEDRAL_CW
                else Chem.ChiralType.CHI_TETRAHEDRAL_CW)
            out = Chem.MolToSmiles(m)
            if out != r["IL_SMILES"]:
                r["IL_SMILES"] = out
                return
            a.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CW
                           if a.GetChiralTag() == Chem.ChiralType.CHI_TETRAHEDRAL_CCW
                           else Chem.ChiralType.CHI_TETRAHEDRAL_CCW)
MUT = {
 "C1 disconnected fragment":   ("vitro", lambda r: r.update(IL_SMILES=r["IL_SMILES"]+".[Cl-]")),
 "C1 radical electron":        ("vitro", lambda r: r.update(IL_SMILES="CCCCCCCC[CH2]")),
 "C1 non-canonical SMILES":    ("vitro", lambda r: r.update(IL_SMILES="C(C)(C)CCCCCCCCCCCN(C)C")),
 "C5 ratios sum to 109":       ("vitro", lambda r: r.update(IL_molratio=str(float(r["IL_molratio"] or 0)+9))),
 "C5 negative ratio":          ("vitro", lambda r: r.update(HL_molratio="-5")),
 "C6 Experiment_value nan":    ("vitro", lambda r: r.update(Experiment_value="nan")),
 "C6 |z| = 40":                ("vitro", lambda r: r.update(Experiment_value="40.0")),
 "C7 wrong Model":             ("vitro", lambda r: r.update(Model="in_vivo")),
 "C7 wrong Route":             ("vitro", lambda r: r.update(Route_of_administration="intravenous")),
 "C7 blanked IL_name":         ("vitro", lambda r: r.update(IL_name="")),
 "C7 blanked PMID (+1)":       ("vitro", lambda r: r.update(Publication_PMID="")),
 "C8 stale Molecular.Weight":  ("vitro", lambda r: r.update(**{"Molecular.Weight":"999.9"})),
 "C8 stale Rings":             ("vitro", lambda r: r.update(Rings="99")),
 "C13 bad confidence value":   ("vitro", lambda r: r.update(confidence="med")),
 "C13 bad value_is_bound":     ("vitro", lambda r: r.update(value_is_bound="floored")),
 "C13 bad norm unit":          ("vitro", lambda r: r.update(value_normalization_unit="per_plate")),
 "C5 vivo ratios (known=3)":   ("vivo",  lambda r: r.update(IL_molratio=str(float(r["IL_molratio"] or 0)+9))),
 "C8 vivo stale LogP":         ("vivo",  lambda r: r.update(LogP="42")),
 "C8 stale TPSA":              ("vitro", lambda r: r.update(**{"Topological.Polar.Surface.Area":"7"})),
 "C8 stale Molar.Refractivity": ("vitro", lambda r: r.update(**{"Molar.Refractivity":"3"})),
 "C8 stale Fraction.sp3":      ("vitro", lambda r: r.update(**{"Fraction.sp3.Carbons":"0.1"})),
 "C8 stale Aromatic.Rings":    ("vitro", lambda r: r.update(**{"Aromatic.Rings":"5"})),
 "C8 stale HB Acceptors":      ("vitro", lambda r: r.update(**{"Hydrogen.Bond.Acceptors":"77"})),
 "C7 blanked Model_type":      ("vitro", lambda r: r.update(Model_type="")),
 "C7 blanked Cargo":           ("vitro", lambda r: r.update(Cargo="")),
 "C7 blanked Experiment_method": ("vitro", lambda r: r.update(Experiment_method="")),
 # Needs a row that actually HAS a steroid ester, so it carries a selector:
 # the default middle row is RG_2023, with no sterol, and a mutation that
 # cannot land makes a working check look broken.
 "C14 epicholesterol (F-013)": ("vivo", _epimerise,
                                lambda r: r["Experiment_ID"] == "ZT_2026"),
 "C9 stereo stripped":         ("vitro", lambda r: r.update(
     IL_SMILES="CCCCC/C=C\\C/C=C\\CCCCCCCCC(CCCCCCCC/C=C\\C/C=C\\CCCCC)OC(=O)CCCN(C)C"
               .replace("/", "").replace("\\", ""))),
}

def run(tmp):
    env = dict(os.environ)
    r = subprocess.run([sys.executable, "-I", os.path.join(tmp, "provenance", "audit.py")],
                       capture_output=True, text=True, env=env)
    return r.returncode, r.stdout

results = []
# baseline: unmutated copy must PASS
for label, spec in [("(control: no mutation)", ("vitro", lambda r: None))] + list(MUT.items()):
    corpus, mutate = spec[0], spec[1]
    select = spec[2] if len(spec) > 2 else None
    with tempfile.TemporaryDirectory() as tmp:
        for d in ("provenance", "data", "curation"):
            # Only what audit.py and its C10/C11 delegates actually read. The
            # figure renders, PDFs and cloned repos are hundreds of MB and the
            # copy, not the audit, dominated the runtime: 29 iterations took
            # over 40 minutes of wall time for 20 seconds of CPU.
            shutil.copytree(os.path.join(REPO, d), os.path.join(tmp, d),
                            ignore=shutil.ignore_patterns(
                                "repo", "__pycache__", "00_pages", "review",
                                "*.pdf", "*.png", "*.jpg", "*.jpeg", "*.zst",
                                "*.tar", "*.ipynb"))
        path = os.path.join(tmp, "data", f"data_{corpus}.csv")
        rows = list(csv.DictReader(open(path)))
        cols = list(rows[0])
        target = rows[len(rows) // 2]
        if select is not None:
            picked = [r for r in rows if select(r)]
            assert picked, f"{label}: no row matches the selector"
            target = picked[len(picked) // 2]
        before = dict(target)
        mutate(target)
        assert target != before or label.startswith("(control"), (
            f"{label}: the mutation did not change the row, so this run would "
            f"test nothing. Fix the mutation, not the check.")
        with open(path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=cols); w.writeheader(); w.writerows(rows)
        rc, out = run(tmp)
        expect_fail = label != "(control: no mutation)"
        ok = (rc != 0) == expect_fail
        gate = [l for l in out.splitlines() if l.strip().startswith("*")]
        results.append((ok, label, rc, gate[0].strip()[:96] if gate else ""))

bad = 0
for ok, label, rc, gate in results:
    bad += not ok
    print(f"  {'PASS' if ok else 'FAIL':4}  rc={rc}  {label:30} {gate}")
print(f"\n{len(results)-bad}/{len(results)} as expected" + ("" if not bad else f"  -- {bad} WRONG"))
sys.exit(1 if bad else 0)
