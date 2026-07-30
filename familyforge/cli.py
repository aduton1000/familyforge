#!/usr/bin/env python3
"""familyforge CLI — interpret a CallForge annotated joint callset for one family."""
from __future__ import annotations
import argparse
import datetime
import json
import os
import shutil
import subprocess
import sys

from .inputs import load_config, load_pedigree, load_panel, load_cnv_calls, load_extra_scores
from .vcf_io import load_variants
from .inheritance import evaluate_variant, evaluate_compound_het, relatedness_sex_qc
from .acmg import tier_candidate
from .cnv import cnv_candidates
from .report import build_context, render


def _normalise_vcf(vcf, outdir):
    """Split multiallelics so genotype interpretation is unambiguous. Uses bcftools if
    present; otherwise returns the input unchanged (multiallelic sites are skipped)."""
    if not shutil.which("bcftools"):
        sys.stderr.write("[familyforge] WARN: bcftools not on PATH — skipping normalisation; "
                         "multiallelic sites will be ignored.\n")
        return vcf
    os.makedirs(outdir, exist_ok=True)
    out = os.path.join(outdir, "normalised.vcf.gz")
    subprocess.run(["bcftools", "norm", "-m-any", "-Oz", "-o", out, vcf], check=True)
    subprocess.run(["bcftools", "index", "-t", "-f", out], check=True)
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(prog="familyforge", description=__doc__)
    ap.add_argument("--vcf", required=True, help="CallForge annotated joint VCF.gz")
    ap.add_argument("--ped", required=True, help="PED pedigree (6-col)")
    ap.add_argument("--panel", required=True, help="gene panel TSV (gene,tier,inheritance,mechanism,...)")
    ap.add_argument("--cnv", help="CallForge cnv_calls.tsv (optional)")
    ap.add_argument("--config", help="override config YAML (optional)")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--study-name", default="Family study")
    ap.add_argument("--mz-twins", action="append", default=[],
                    help="comma-separated MZ-twin sample_ids, repeatable (e.g. TWIN_1,TWIN_2)")
    ap.add_argument("--provenance", help="CallForge provenance.json for the methods section")
    ap.add_argument("--extra-scores", help="TSV of finalist in-silico scores "
                    "(key,cadd_phred,revel,alphamissense_score,alphamissense_class); "
                    "key='chrom:pos:ref>alt'. Use when genome-wide dbNSFP is unavailable.")
    ap.add_argument("--skip-norm", action="store_true")
    a = ap.parse_args(argv)

    cfg = load_config(a.config)
    mz = [g.split(",") for g in a.mz_twins] if a.mz_twins else cfg.get("mz_twin_groups", [])
    cfg["mz_twin_groups"] = [[x.strip() for x in g] for g in mz]

    vcf = a.vcf if a.skip_norm else _normalise_vcf(a.vcf, os.path.join(a.outdir, "work"))

    panel = load_panel(a.panel)
    panel_genes = set(panel)

    # discover which pedigree members are actually in the VCF
    from cyvcf2 import VCF
    samples = list(VCF(vcf).samples)
    fam = load_pedigree(a.ped, sequenced_ids=set(samples), mz_twin_groups=cfg["mz_twin_groups"])
    sex_by_id = {m.iid: m.sex for m in fam.members.values()}

    variants, vcf_samples = load_variants(
        vcf, panel_genes, cfg["af_fields"], cfg["min_gq"], cfg["min_dp"],
        cfg["min_alt_ab"], cfg["max_ref_ab_for_hom"], sex_by_id)
    sys.stderr.write(f"[familyforge] loaded {len(variants)} panel variants across "
                     f"{len(vcf_samples)} samples\n")

    # attach externally-supplied finalist in-silico scores (AlphaMissense/REVEL/CADD)
    extra = load_extra_scores(a.extra_scores)
    if extra:
        nmatched = 0
        for v in variants:
            s = extra.get(v.key)
            if s:
                v.cadd, v.revel, v.am_score, v.am_class = s["cadd"], s["revel"], s["am_score"], s["am_class"]
                nmatched += 1
        sys.stderr.write(f"[familyforge] applied extra in-silico scores to {nmatched} variant(s)\n")

    cands = []
    for v in variants:
        cands += evaluate_variant(v, fam, panel)
    cands += evaluate_compound_het(variants, fam, panel)
    cnv_rows = load_cnv_calls(a.cnv)
    cands += cnv_candidates(cnv_rows, panel, fam, cfg)

    cands = [tier_candidate(c, panel, cfg) if c.variant is not None else c for c in cands]

    qc = relatedness_sex_qc(variants, vcf_samples, fam)

    prov = {}
    if a.provenance and os.path.isfile(a.provenance):
        with open(a.provenance) as fh:
            prov = json.load(fh)
    meta = {"study": a.study_name,
            "date": datetime.date.today().isoformat(),
            "vcf": os.path.abspath(a.vcf), "panel": os.path.abspath(a.panel),
            "n_panel_genes": len(panel), "n_cnv_rows": len(cnv_rows),
            "provenance": prov,
            "af_dominant_max": cfg["af_dominant_max"], "af_recessive_max": cfg["af_recessive_max"]}

    ctx = build_context(cands, qc, fam, panel, meta)
    written = render(ctx, a.outdir)
    # machine-readable candidate dump
    with open(os.path.join(a.outdir, "candidates.json"), "w") as fh:
        json.dump({"context_summary": {"n_candidates": ctx["n_candidates"],
                                       "tierA": len(ctx["tierA"]), "tierB": len(ctx["tierB"]),
                                       "tierC": len(ctx["tierC"])},
                   "rows": ctx["rows"], "qc": qc}, fh, indent=2, default=str)
    sys.stderr.write(f"[familyforge] candidates: {len(ctx['tierA'])}A / {len(ctx['tierB'])}B / "
                     f"{len(ctx['tierC'])}C\n[familyforge] wrote:\n  " + "\n  ".join(written) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
