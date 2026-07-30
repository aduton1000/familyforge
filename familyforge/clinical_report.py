#!/usr/bin/env python3
"""Reusable professional clinical-report engine.

Renders a letterhead clinical genomics report (HTML; render to PDF with headless Chrome)
from: a reusable lab-branding YAML, a per-report metadata YAML, a familyforge
candidates.json, and the PED. Disease-agnostic — all content is data-driven.

    python -m familyforge.clinical_report \
        --lab-config lab_branding.yaml --report-config report_meta.yaml \
        --candidates results/reports/candidates.json --ped config/family.ped \
        --out results/reports/clinical_report.html
"""
from __future__ import annotations
import argparse, base64, datetime, json, os, sys, yaml
from jinja2 import Environment, FileSystemLoader, select_autoescape

_SEX = {"M": "Male", "F": "Female", "U": "—", "1": "Male", "2": "Female", "0": "—"}


def _ped_members(ped_path):
    out = []
    with open(ped_path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.split()
            if len(f) >= 6:
                out.append({"iid": f[1], "sex": _SEX.get(f[4], "—"),
                            "sequenced": f[3] != "0" or f[2] != "0" or True})
    return out


def build_context(lab_cfg, rep_cfg, candidates, ped_path):
    banner_b64 = base64.b64encode(open(lab_cfg["logo_banner"], "rb").read()).decode()
    lab_logo_b64 = base64.b64encode(open(lab_cfg["logo_lab"], "rb").read()).decode()

    # family: order + labels from report_meta.samples; fall back to PED.
    # "Sequenced" is derived from the actual analysed sample set (candidates QC), so a
    # relationship-only member (e.g. an untested father) is correctly marked not sequenced.
    ped = {m["iid"]: m for m in _ped_members(ped_path)}
    seq_ids = set((candidates.get("qc", {}).get("sex", {})).keys())
    fam = []
    for iid, info in rep_cfg.get("samples", {}).items():
        pm = ped.get(iid, {})
        seq = "Yes" if iid in seq_ids else "No — not available"
        fam.append({"iid": iid, "relationship": info.get("relationship", ""),
                    "sex": pm.get("sex", info.get("sex", "—")),
                    "status": info.get("status", ""),
                    "sequenced": info.get("sequenced", seq)})

    rows = candidates.get("rows", [])
    findings = [r for r in rows if r.get("tier") in ("A", "B")]
    top = findings[0] if findings else None
    gm = rep_cfg.get("gene_management", {})
    top_mgmt = gm.get(top["gene"]) if top else None

    return {
        "lab": {**lab_cfg, "banner_b64": banner_b64, "lab_logo_b64": lab_logo_b64},
        "rep": rep_cfg,
        "date": datetime.date.today().strftime("%d %B %Y"),
        "family": fam,
        "findings": findings, "top": top, "top_mgmt": top_mgmt,
        "qc": candidates.get("qc", {}),
        "notes": rep_cfg.get("notes", []),
    }


def render(ctx, out_html):
    here = os.path.dirname(os.path.abspath(__file__))
    env = Environment(loader=FileSystemLoader(os.path.join(here, "templates")),
                      autoescape=select_autoescape(["html"]), trim_blocks=True, lstrip_blocks=True)
    os.makedirs(os.path.dirname(os.path.abspath(out_html)), exist_ok=True)
    outs = []
    with open(out_html, "w") as fh:
        fh.write(env.get_template("clinical_letterhead.html.j2").render(**ctx))
    outs.append(out_html)
    # Markdown twin (content mirror; markdown carries no page layout / banner / repeating header)
    if os.path.isfile(os.path.join(here, "templates", "clinical_letterhead.md.j2")):
        out_md = os.path.splitext(out_html)[0] + ".md"
        with open(out_md, "w") as fh:
            fh.write(env.get_template("clinical_letterhead.md.j2").render(**ctx))
        outs.append(out_md)
    return outs


def main(argv=None):
    ap = argparse.ArgumentParser(prog="familyforge-clinical-report", description=__doc__)
    ap.add_argument("--lab-config", required=True)
    ap.add_argument("--report-config", required=True)
    ap.add_argument("--candidates", required=True, help="familyforge candidates.json")
    ap.add_argument("--ped", required=True)
    ap.add_argument("--out", required=True, help="output HTML path")
    a = ap.parse_args(argv)
    lab = yaml.safe_load(open(a.lab_config))
    rep = yaml.safe_load(open(a.report_config))
    cand = json.load(open(a.candidates))
    ctx = build_context(lab, rep, cand, a.ped)
    outs = render(ctx, a.out)
    for o in outs:
        sys.stderr.write(f"[clinical_report] wrote {o}\n")
    sys.stderr.write(f"[clinical_report] render PDF: chrome --headless=new --no-pdf-header-footer "
                     f"--print-to-pdf=<out.pdf> file://{os.path.abspath(outs[0])}\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
