"""Stage 2 orchestrator: bioisostere input -> ranked docking results."""
from __future__ import annotations

import math
import os

from . import align3d, crem_variants, dock, embed3d, inputs, pdb2pdbqt, rank_final


def run(config: dict) -> dict:
    out = config["out"]
    os.makedirs(out, exist_ok=True)

    # 1. bioisostere input (official CSV or local L1 DB)
    source = config.get("source", "official_csv")
    query, df = inputs.get_bioisosteres(
        source=source,
        query=config.get("query"),
        csv_path=config.get("csv_path"),
        mmpdb_path=config.get("mmpdb_path"),
        topn=config.get("input_topn"),
        logp_range=config.get("logp_range"),
        tpsa_range=config.get("tpsa_range"),
        mw_range=config.get("mw_range"),
    )
    df.to_csv(os.path.join(out, "1_bioisosteres.csv"), index=False)
    cands = df["Candidate Fragments"].tolist()
    if config.get("max_candidates"):
        cands = cands[: config["max_candidates"]]

    # 2. CReM variants
    #    SwissBioisostere ranks candidates and takes a top-k; we mirror that by
    #    allocating the variant budget *per candidate* (not a global slice, which
    #    would starve all but the first candidate).
    variants_all = crem_variants.expand_candidates(cands, config["crem_db"], **config.get("crem", {}))
    per_cand = config.get("variants_per_candidate")
    if per_cand is None and config.get("max_variants"):
        per_cand = max(1, math.ceil(config["max_variants"] / max(1, len(cands))))
    if per_cand:
        variants = []
        for cand in cands:
            vs = [v for v in variants_all if v["candidate"] == cand]
            variants.extend(vs[:per_cand])
    else:
        variants = variants_all
    smis = [v["variant"] for v in variants]
    cand_of = [v["candidate"] for v in variants]

    # 3. 3D embed
    sdf = embed3d.batch_embed(smis, os.path.join(out, "2_3d"))

    # 4. align to reference ligand (backend: "mcs" rigid, or "kcombu" flexible)
    aligned = align3d.batch_align(
        sdf,
        config["ref_ligand"],
        os.path.join(out, "3_aligned"),
        backend=config.get("align_backend", align3d.DEFAULT_BACKEND),
        **config.get("kcombu", {}),
    )

    # 5. PDBQT
    pdbqts = pdb2pdbqt.batch_pdbqt(aligned, os.path.join(out, "4_pdbqt"))

    # 6. docking
    center, size = align3d.box_from_ligand(config["ref_ligand"], padding=config.get("box_padding", 8.0))
    backend = dock.get_backend(config.get("docking_backend", "vina"))
    dock_dir = os.path.join(out, "5_dock")
    os.makedirs(dock_dir, exist_ok=True)
    recs = []
    for i, p in enumerate(pdbqts):
        r = backend(
            config["receptor"], p, center, size,
            out_pdbqt=os.path.join(dock_dir, f"dock_{i:06d}.pdbqt"),
            **config.get("dock", {}),
        )
        r["variant"] = smis[i] if i < len(smis) else ""
        r["candidate"] = cand_of[i] if i < len(cand_of) else ""
        r["center"] = center
        r["size"] = size
        recs.append(r)

    # 7. rank
    ranked = rank_final.rank_results(recs)
    rank_final.write_csv(ranked, os.path.join(out, "6_ranked.csv"))
    rank_final.write_html(ranked, os.path.join(out, "6_ranked.html"))

    return {
        "query": query,
        "n_candidates": len(cands),
        "n_variants": len(variants),
        "n_docked": len(pdbqts),
        "ranked": ranked,
    }
