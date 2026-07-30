"""Build the report context and render the clinician + publication reports
(HTML + Markdown) from the ranked candidates and QC."""
from __future__ import annotations
import os
from jinja2 import Environment, FileSystemLoader, select_autoescape
from . import __version__

_TIER_RANK = {"A": 0, "B": 1, "C": 2, "": 3}
_CONF_RANK = {"high": 0, "moderate": 1, "low": 2}


def rank_candidates(cands, panel):
    def key(c):
        gt = int(panel.get(c.gene, {}).get("tier", "3") or 3)
        return (_TIER_RANK.get(c.tier, 3), _CONF_RANK.get(c.confidence, 3), gt)
    return sorted(cands, key=key)


def _var_row(c):
    v = c.variant
    base = {"gene": c.gene, "model": c.model, "tier": c.tier, "confidence": c.confidence,
            "zygosity": c.zygosity, "segregation": c.segregation,
            "inheritance_note": c.inheritance_note, "tier_reason": c.tier_reason,
            "flags": c.flags}
    if v is None:  # CNV
        base.update({"locus": (c.cnv or {}).get("chrom", "") + ":" + str((c.cnv or {}).get("start", "")),
                     "consequence": "deletion", "hgvsc": ".", "hgvsp": ".",
                     "clinvar": ".", "gnomad": ".", "gnomad_afr": ".", "phylop": ".",
                     "sift": ".", "polyphen": ".", "cadd": ".", "revel": ".", "alphamissense": ".",
                     "mane_refseq": ".", "rsid": ".", "filter": "."})
        return base
    # clinical-standard HGVS: prefer the MANE Select RefSeq accession; strip the protein-ID prefix
    c_part = v.gene.hgvsc.split(":")[-1] if ":" in v.gene.hgvsc else v.gene.hgvsc
    p_part = v.gene.hgvsp.split(":")[-1] if ":" in v.gene.hgvsp else v.gene.hgvsp
    mane = v.gene.mane_select if v.gene.mane_select not in (".", "") else None
    hgvsc_disp = f"{mane}:{c_part}" if (mane and c_part not in (".", "")) else v.gene.hgvsc
    base.update({
        "locus": f"{v.chrom}:{v.pos} {v.ref}>{v.alt}", "rsid": v.vid, "filter": v.filter,
        "consequence": v.gene.consequence, "hgvsc": hgvsc_disp, "hgvsp": p_part,
        "mane_refseq": (mane or "."),
        "clinvar": v.clinvar,
        "gnomad": ("absent" if v.max_af is None else f"{v.max_af:.2g}"),
        "gnomad_afr": next((f"{v.af_detail[k]:.2g}" for k in ("gnomADe_AFR_AF", "gnomADg_AFR_AF")
                            if v.af_detail.get(k) is not None), "."),
        "phylop": ("." if v.phylop is None else f"{v.phylop:.1f}"),
        "sift": v.gene.sift, "polyphen": v.gene.polyphen,
        "cadd": ("." if v.cadd is None else f"{v.cadd:.1f}"),
        "revel": ("." if v.revel is None else f"{v.revel:.3f}"),
        "alphamissense": ("." if v.am_score is None else f"{v.am_score:.3f} ({v.am_class})"),
    })
    if c.partner is not None:
        p = c.partner
        base["partner"] = f"{p.chrom}:{p.pos} {p.ref}>{p.alt} ({p.gene.hgvsp or p.gene.hgvsc})"
    return base


def classify(tier, clinvar, gnomad):
    """Map a finding to an ACMG-style classification (class code, label) for the
    colour convention: P=pathogenic, LP=likely pathogenic, VUS, LB=likely benign, B=benign."""
    cv = (clinvar or "").lower()
    try:
        af = float(gnomad) if gnomad not in ("absent", ".", "", "None", None) else 0.0
    except (TypeError, ValueError):
        af = 0.0
    if "conflicting" in cv:
        return "VUS", "Conflicting (VUS)"
    if "pathogenic" in cv and "benign" not in cv:
        return ("LP", "Likely pathogenic") if "likely" in cv else ("P", "Pathogenic")
    if "benign" in cv:
        return ("LB", "Likely benign") if "likely" in cv else ("B", "Benign")
    if tier == "A":
        return "LP", "Likely pathogenic"
    if af > 0.05:                       # common with no pathogenic ClinVar
        return "LB", "Likely benign"
    if tier == "B":
        return "VUS", "Uncertain significance"
    return "VUS", "Uncertain significance"


def build_context(cands, qc, fam, panel, meta):
    ranked = rank_candidates(cands, panel)
    rows = [_var_row(c) for c in ranked]
    for r in rows:
        r["acmg_class"], r["acmg_label"] = classify(r["tier"], r.get("clinvar", "."), r.get("gnomad", "."))
    tiersA = [r for r in rows if r["tier"] == "A"]
    tiersB = [r for r in rows if r["tier"] == "B"]
    tiersC = [r for r in rows if r["tier"] == "C"]
    fam_rows = [{"iid": m.iid, "role": ("affected" if m.affected else "unaffected"),
                 "sex": m.sex, "sequenced": m.sequenced,
                 "father": m.father, "mother": m.mother}
                for m in fam.members.values()]
    return {
        "version": __version__, "meta": meta,
        "family": fam_rows,
        "mz_groups": fam.mz_twin_groups,
        "n_candidates": len(rows), "rows": rows,
        "tierA": tiersA, "tierB": tiersB, "tierC": tiersC,
        "qc": qc,
        "panel": panel,
        "top": rows[0] if rows else None,
    }


def render(ctx, outdir):
    here = os.path.dirname(os.path.abspath(__file__))
    env = Environment(loader=FileSystemLoader(os.path.join(here, "templates")),
                      autoescape=select_autoescape(["html"]), trim_blocks=True, lstrip_blocks=True)
    os.makedirs(outdir, exist_ok=True)
    written = []
    for name, out in [("clinician.html.j2", "clinician_report.html"),
                      ("clinician.md.j2", "clinician_report.md"),
                      ("publication.html.j2", "publication_report.html"),
                      ("publication.md.j2", "publication_report.md")]:
        tpl = env.get_template(name)
        path = os.path.join(outdir, out)
        with open(path, "w") as fh:
            fh.write(tpl.render(**ctx))
        written.append(path)
    return written
