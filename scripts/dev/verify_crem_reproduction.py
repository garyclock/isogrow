#!/usr/bin/env python
"""Reproduce a reference CReM variant set (stage B) from an official SB export.

This is the check we used to validate the pipeline against a previous run of the
same notebook:

    python scripts/dev/verify_crem_reproduction.py \\
        --official /path/to/SwissBioisostere_fragment_output.csv \\
        --crem-db /path/to/chembl22_sa2.db \\
        --reference-dir /path/to/1_crem_sdf \\
        --candidates 2

``--reference-dir`` must contain the reference SDFs (``*/*.sdf``, named
``mol_*.sdf``).  The script groups them by the candidate's ring system, then
compares with the variant set generated here.  It needs an external reference
run, so it is not part of the automated test suite.
"""
from __future__ import annotations

import argparse
import glob
import os
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from rdkit import Chem  # noqa: E402

from isogrow.pipeline import crem_variants, inputs  # noqa: E402

DUMMY = re.compile(r"\[\*:\d+\]")


def load_reference(reference_dir: str) -> dict[str, str]:
    ref: dict[str, str] = {}
    for f in glob.glob(os.path.join(reference_dir, "*", "*.sdf")) or \
            glob.glob(os.path.join(reference_dir, "*.sdf")):
        mol = Chem.MolFromMolFile(f, removeHs=False)
        if mol is None:
            continue
        ref[Chem.MolToSmiles(Chem.RemoveHs(mol))] = f
    return ref


def ring_core(candidate: str):
    """H-capped candidate used as a substructure filter."""
    return Chem.MolFromSmiles(DUMMY.sub("[H]", candidate))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--official", required=True)
    ap.add_argument("--crem-db", required=True)
    ap.add_argument("--reference-dir", required=True)
    ap.add_argument("--candidates", type=int, default=2)
    ap.add_argument("--logp-range", type=float, nargs=2, default=(-10, 10))
    ap.add_argument("--tpsa-range", type=float, nargs=2, default=(-300, 300))
    ap.add_argument("--mw-range", type=float, nargs=2, default=(-200, 80))
    args = ap.parse_args()

    print("=== 1. filtered candidates (official CSV, notebook ranges) ===")
    query, df = inputs.get_bioisosteres(
        source="official_csv", csv_path=args.official,
        logp_range=args.logp_range, tpsa_range=args.tpsa_range, mw_range=args.mw_range,
    )
    cands = df["Candidate Fragments"].tolist()
    print(f"query={query}  filtered_candidates={len(cands)}  first3={cands[:3]}")

    print("=== 2. reference variants ===")
    ref = load_reference(args.reference_dir)
    print(f"reference molecules: {len(ref)}")
    if not ref:
        print("no reference SDFs found — check --reference-dir")
        return 1

    print("=== 3. variant-set comparison ===")
    for cand in cands[: args.candidates]:
        lead, protected = crem_variants.replace_stars_with_carbons(cand)
        ours = crem_variants.crem_mutate(
            lead, protected, args.crem_db, radius=2, min_size=0, max_size=0,
            min_inc=-2, max_inc=2, replace_cycles=False, protect_H=True,
            max_replacements=None,
        )
        our_set = set(ours)
        core = ring_core(cand)
        theirs = {s for s in ref
                  if (Chem.MolFromSmiles(s) is not None
                      and Chem.MolFromSmiles(s).HasSubstructMatch(core))}
        inter = our_set & theirs
        print(f"\ncandidate: {cand}")
        print(f"  lead: {lead}  protected_ids={protected}")
        print(f"  ours: {len(our_set)}   reference: {len(theirs)}   overlap: {len(inter)}")
        if theirs:
            print(f"  recall: {len(inter)}/{len(theirs)} = {len(inter) / len(theirs):.1%}")
        print(f"  missing sample: {sorted(theirs - our_set)[:8]}")
        print(f"  extra   sample: {sorted(our_set - theirs)[:8]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
