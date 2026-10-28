"""Command line interface: ``isogrow <command>``.

Commands
--------
``doctor``     check the environment (chemistry stack, mmpdb, docking engine, data)
``fetch``      download ChEMBL / Vina artifacts into the data directory
``build``      build an mmpdb database from chemreps (optionally activity-aware)
``search``     fragment → replacements, SwissBioisostere-compatible output
``run``        the full A→G pipeline from a JSON config
``benchmark``  compare the local engine against an official export
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


# --------------------------------------------------------------------------- #
# doctor
# --------------------------------------------------------------------------- #
def cmd_doctor(args: argparse.Namespace) -> int:
    from isogrow import paths

    print(f"python          : {sys.version.split()[0]}  ({sys.platform})")
    checks = []
    for mod in ("rdkit", "pandas", "numpy", "meeko", "crem", "mmpdblib"):
        try:
            m = __import__(mod)
            checks.append((mod, getattr(m, "__version__", "ok")))
        except Exception as exc:  # noqa: BLE001 - report, never crash doctor
            checks.append((mod, f"MISSING ({type(exc).__name__})"))
    for name, ver in checks:
        print(f"{name:<15} : {ver}")

    for label, resolver in (("mmpdb", paths.mmpdb_exe), ("vina", paths.vina_exe)):
        try:
            print(f"{label:<15} : {resolver()}")
        except FileNotFoundError as exc:
            print(f"{label:<15} : not found — {exc}")

    data = paths.data_dir()
    print(f"ISOGROW_DATA   : {data} ({'exists' if data.is_dir() else 'missing'})")
    if data.is_dir():
        dbs = sorted(data.glob("mmpdb/*.mmpdb"))
        print(f"databases      : {len(dbs)}")
        for db in dbs[:10]:
            print(f"  {db.name:<40} {db.stat().st_size / 1e6:8.1f} MB")
    ok = not any("MISSING" in v for _, v in checks)
    return 0 if ok else 1


# --------------------------------------------------------------------------- #
# fetch
# --------------------------------------------------------------------------- #
def cmd_fetch(args: argparse.Namespace) -> int:
    from isogrow import fetch, paths

    out = Path(args.out_dir) if args.out_dir else paths.data_dir(create=True)
    if args.what == "chemreps":
        print(fetch.fetch_chemreps(args.chembl_version, out / "chembl"))
    elif args.what == "sqlite":
        print(fetch.fetch_sqlite(args.chembl_version, out / "chembl"))
    elif args.what == "activities":
        meta = fetch.save_activity_sample(args.target, args.std_type, args.max_rows,
                                          out / "chembl_sample")
        print(json.dumps(meta, ensure_ascii=False, indent=2))
    elif args.what == "vina":
        print(fetch.fetch_vina(args.vina_version, out / "bin"))
    return 0


# --------------------------------------------------------------------------- #
# build
# --------------------------------------------------------------------------- #
def cmd_build(args: argparse.Namespace) -> int:
    from isogrow import paths
    from isogrow.bioisostere import chembl

    out = Path(args.out) if args.out else paths.data_dir(create=True) / "mmpdb" / "chembl.mmpdb"
    prefix = str(out)[:-6] if str(out).endswith(".mmpdb") else str(out)
    exe = paths.mmpdb_exe(args.mmpdb_exe)

    if args.activity_table:
        db = chembl.build_from_activity(args.activity_table, args.chemreps, prefix,
                                        jobs=args.jobs, num_cuts=args.num_cuts, mmpdb_exe=exe)
    elif args.smi:
        db = chembl.build_from_smiles(args.smi, prefix, limit=args.limit, jobs=args.jobs,
                                      num_cuts=args.num_cuts, max_heavies=args.max_heavies,
                                      mmpdb_exe=exe)
    else:
        db = chembl.build_from_chemreps(args.chemreps, prefix, limit=args.limit,
                                        jobs=args.jobs, num_cuts=args.num_cuts,
                                        max_heavies=args.max_heavies, mmpdb_exe=exe)
    print(f"\nDatabase ready: {db}\nUse it with:\n  isogrow search --db {db} --query '[*:1]c1cc2ccccc2o1'")
    return 0


# --------------------------------------------------------------------------- #
# search
# --------------------------------------------------------------------------- #
def cmd_search(args: argparse.Namespace) -> int:
    import pandas as pd

    from isogrow.bioisostere import engine, plot
    from isogrow.bioisostere import io as sbio

    if args.activity:
        from isogrow.bioisostere.activity import add_similarity, query_activity

        df = query_activity(args.db, args.query)
        if df.empty:
            print("no replacements found")
            return 1
        df = add_similarity(df, args.query)
        show = [c for c in ("Candidate Fragments", "Frequency", "# Better", "# Similar",
                            "# Worse", "Δ pchembl", "success_score", "similarity")
                if c in df.columns]
    else:
        df = engine.search(args.query, args.db, topn=args.topn, add_l2=args.l2,
                           logp_range=args.logp_range, tpsa_range=args.tpsa_range,
                           mw_range=args.mw_range)
        if df.empty:
            print("no replacements found")
            return 1
        show = list(df.columns)

    with pd.option_context("display.width", 200, "display.max_columns", 30):
        print(df[show].head(args.topn or 20).to_string(index=False))
    print(f"\n{len(df)} replacements")

    if args.out:
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        if args.activity:
            df.to_csv(args.out, index=False)
        else:
            sbio.write_sb_csv(df, args.out, args.query)
        print(f"written -> {args.out}")

    if args.plot:
        Path(args.plot).parent.mkdir(parents=True, exist_ok=True)
        plot.plot_property_space(df, args.plot, title=f"{args.query} property space")
        print(f"plot -> {args.plot}")
    return 0


# --------------------------------------------------------------------------- #
# run / benchmark
# --------------------------------------------------------------------------- #
def cmd_run(args: argparse.Namespace) -> int:
    from isogrow import config as cfgmod
    from isogrow.pipeline.orchestrator import run

    overrides = {
        "source": args.source,
        "out": args.out,
        "max_candidates": args.max_candidates,
        "max_variants": args.max_variants,
        "variants_per_candidate": args.variants_per_candidate,
        "docking_backend": args.docking_backend,
    }
    cfg = cfgmod.load_config(args.config, overrides=overrides)
    res = run(cfg)
    print(f"query={res['query']} candidates={res['n_candidates']} "
          f"variants={res['n_variants']} docked={res['n_docked']}")
    ranked = res["ranked"]
    if len(ranked):
        print(ranked.head(15).to_string(index=False))
    print("out:", cfg["out"])
    return 0


def cmd_benchmark(args: argparse.Namespace) -> int:
    from isogrow.bioisostere.benchmark import benchmark

    report = benchmark(args.official, args.db, args.out, topn=args.topn)
    print(json.dumps(report, ensure_ascii=False, indent=2))
    return 0


# --------------------------------------------------------------------------- #
def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="isogrow", description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--version", action="version", version="isogrow 0.1.0")
    sub = p.add_subparsers(dest="command", required=True)

    d = sub.add_parser("doctor", help="check environment and data")
    d.set_defaults(func=cmd_doctor)

    f = sub.add_parser("fetch", help="download datasets / binaries")
    f.add_argument("what", choices=["chemreps", "sqlite", "activities", "vina"])
    f.add_argument("--out-dir", default=None, help="default: $ISOGROW_DATA")
    f.add_argument("--chembl-version", default="37")
    f.add_argument("--vina-version", default="1.2.5")
    f.add_argument("--target", default="CHEMBL203")
    f.add_argument("--std-type", default="IC50")
    f.add_argument("--max-rows", type=int, default=5000)
    f.set_defaults(func=cmd_fetch)

    b = sub.add_parser("build", help="build an mmpdb database")
    src = b.add_mutually_exclusive_group(required=True)
    src.add_argument("--chemreps", help="ChEMBL chemreps dump (.txt or .txt.gz)")
    src.add_argument("--smi", help="plain SMILES<TAB>id file (no downloads needed)")
    b.add_argument("--out", default=None, help="output .mmpdb path")
    b.add_argument("--limit", type=int, default=200_000, help="0 = all compounds")
    b.add_argument("--jobs", type=int, default=8)
    b.add_argument("--num-cuts", type=int, default=3,
                   help="bonds cut at most; SwissBioisostere uses 3")
    b.add_argument("--max-heavies", type=int, default=60)
    b.add_argument("--activity-table", default=None,
                   help="CSV with molecule_chembl_id,pchembl to add ΔpChEMBL")
    b.add_argument("--mmpdb-exe", default=None)
    b.set_defaults(func=cmd_build)

    s = sub.add_parser("search", help="fragment -> replacements")
    s.add_argument("--db", required=True)
    s.add_argument("--query", required=True, help="e.g. '[*:1]c1cc2ccccc2o1'")
    s.add_argument("--topn", type=int, default=None)
    s.add_argument("--activity", action="store_true",
                   help="include #Better/#Similar/#Worse and success score")
    s.add_argument("--l2", action="store_true", help="add ΔHAC/ΔRings/... columns")
    s.add_argument("--logp-range", type=float, nargs=2, default=None)
    s.add_argument("--tpsa-range", type=float, nargs=2, default=None)
    s.add_argument("--mw-range", type=float, nargs=2, default=None)
    s.add_argument("--out", default=None, help="write CSV (SwissBioisostere layout)")
    s.add_argument("--plot", default=None, help="write an HTML property-space plot")
    s.set_defaults(func=cmd_search)

    r = sub.add_parser("run", help="run the full pipeline from a JSON config")
    r.add_argument("--config", required=True)
    r.add_argument("--source", choices=["local", "official_csv"], default=None)
    r.add_argument("--out", default=None)
    r.add_argument("--max-candidates", type=int, default=None)
    r.add_argument("--max-variants", type=int, default=None)
    r.add_argument("--variants-per-candidate", type=int, default=None)
    r.add_argument("--docking-backend", default=None)
    r.set_defaults(func=cmd_run)

    m = sub.add_parser("benchmark", help="compare against an official export")
    m.add_argument("--official", required=True)
    m.add_argument("--db", required=True)
    m.add_argument("--out", required=True)
    m.add_argument("--topn", type=int, default=50)
    m.set_defaults(func=cmd_benchmark)
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
