"""ACMG-style tiering of candidate variants (a pragmatic subset, not a full ACMG
classifier). Combines ClinVar, population frequency, predicted consequence, in-silico
damage, conservation, gene disease-mechanism, and segregation into a tier A/B/C with
an evidence trail for the publication report.
"""
from __future__ import annotations
from .model import Candidate


def _clinvar_class(clnsig: str):
    s = (clnsig or "").lower()
    if not s or s == ".":
        return "none"
    if "conflicting" in s:
        return "conflicting"
    if "pathogenic" in s and "benign" not in s:
        return "pathogenic"        # covers Pathogenic + Likely_pathogenic
    if "benign" in s:
        return "benign"            # covers Benign + Likely_benign
    if "uncertain" in s or "vus" in s:
        return "vus"
    return "other"


def _sift_del(sift: str):
    return "deleterious" in (sift or "").lower()


def _polyphen_dmg(pp: str):
    return "damaging" in (pp or "").lower()  # probably/possibly_damaging


def _is_stat3_hotspot(cand: Candidate, cfg: dict):
    if cand.gene != "STAT3":
        return False
    pp = cand.variant.gene.protein_pos
    if pp is None:
        return False
    if pp in set(cfg.get("stat3_hotspot_residues", [])):
        return True
    dbd = cfg.get("stat3_dbd_range", [321, 465])
    sh2 = cfg.get("stat3_sh2_range", [585, 688])
    return (dbd[0] <= pp <= dbd[1]) or (sh2[0] <= pp <= sh2[1])


def tier_candidate(cand: Candidate, panel: dict, cfg: dict) -> Candidate:
    v = cand.variant
    gene = cand.gene
    cons = v.gene.consequence
    lof = cons in set(cfg["lof_consequences"])
    missense = cons in set(cfg["missense_consequences"])
    splice = cons in set(cfg.get("splice_region", [])) or cons in ("splice_donor_variant", "splice_acceptor_variant")
    cvc = _clinvar_class(v.clinvar)
    recessive = cand.model in ("ar_homozygous", "ar_compound_het", "x_linked_recessive")
    af_max = cfg["af_recessive_max"] if recessive else cfg["af_dominant_max"]
    af = v.max_af
    rare = (af is None) or (af <= af_max)
    lof_gene = gene in set(cfg.get("lof_mechanism_genes", []))
    stat3_hot = _is_stat3_hotspot(cand, cfg)
    phylop_cons = (v.phylop is not None and v.phylop >= cfg["phylop_conserved"])
    # finalist in-silico (from --extra-scores): AlphaMissense / REVEL / CADD
    am_lp = (v.am_class or "").lower().startswith("likely_path") or (v.am_class or "").upper() == "P"
    revel_strong = (v.revel is not None and v.revel >= cfg.get("revel_strong", 0.7))
    cadd_high = (v.cadd is not None and v.cadd >= cfg.get("cadd_high", 25.0))
    strong_insilico = am_lp or revel_strong or cadd_high
    insilico_dmg = (_sift_del(v.gene.sift) or _polyphen_dmg(v.gene.polyphen) or phylop_cons
                    or strong_insilico or (v.revel is not None and v.revel >= cfg.get("revel_supportive", 0.5)))
    tier12 = str(panel.get(gene, {}).get("tier", "3")) in ("1", "2")

    ev = []
    if cvc == "pathogenic": ev.append("ClinVar pathogenic/likely-pathogenic")
    if cvc == "benign": ev.append("ClinVar benign/likely-benign")
    if cvc == "conflicting": ev.append("ClinVar conflicting")
    if lof: ev.append(f"predicted LoF ({cons})")
    if splice and not lof: ev.append(f"splice-region ({cons})")
    if missense: ev.append("missense")
    if stat3_hot: ev.append(f"STAT3 dominant-negative hotspot domain (aa {v.gene.protein_pos})")
    if rare: ev.append("rare/absent in gnomAD" + (f" (max AF {af:.2g})" if af is not None else " (absent)"))
    else: ev.append(f"too common for a {'recessive' if recessive else 'dominant'} cause (AF {af:.2g})")
    if _sift_del(v.gene.sift): ev.append(f"SIFT {v.gene.sift}")
    if _polyphen_dmg(v.gene.polyphen): ev.append(f"PolyPhen {v.gene.polyphen}")
    if phylop_cons: ev.append(f"conserved (PhyloP {v.phylop:.1f})")
    if v.am_score is not None: ev.append(f"AlphaMissense {v.am_score:.3f} ({v.am_class})")
    if v.revel is not None: ev.append(f"REVEL {v.revel:.3f}")
    if v.cadd is not None: ev.append(f"CADD {v.cadd:.1f}")
    ev.append(f"segregation consistent ({cand.model})")

    # ---- tier decision ----
    tier, reason = "C", "in-panel variant with limited supporting evidence"
    if not rare or cvc == "benign":
        tier = "C"
        reason = "downgraded: " + ("common in gnomAD" if not rare else "ClinVar benign")
    elif cvc == "pathogenic":
        tier, reason = "A", "ClinVar pathogenic/likely-pathogenic, rare, segregates"
    elif lof and lof_gene and recessive:
        tier, reason = "A", "biallelic predicted-LoF in a LoF-mechanism gene, rare, segregates"
    elif stat3_hot and (missense or lof) and cand.model == "de_novo_or_dominant":
        tier, reason = "A", "STAT3 dominant-negative hotspot missense, de-novo/dominant pattern"
    elif missense and strong_insilico and tier12 and cand.model in (
            "de_novo_or_dominant", "ar_homozygous", "ar_compound_het", "x_linked_recessive"):
        tier, reason = "A", ("likely-pathogenic candidate: rare damaging missense in a Tier-1/2 disease gene "
                             "with matching-inheritance segregation and strong in-silico consensus "
                             "(AlphaMissense/REVEL/CADD); pending functional + parental confirmation")
    elif lof and (recessive or gene not in set(cfg.get("lof_mechanism_genes", [])) is False):
        tier, reason = "B", "predicted-LoF in a panel gene, rare, segregates (mechanism/inheritance to confirm)"
    elif missense and insilico_dmg:
        tier, reason = "B", "rare damaging missense in a panel gene with consistent segregation"
    else:
        tier, reason = "C", "rare in-panel variant, weaker in-silico/segregation support"

    if cvc == "conflicting" and tier == "A":
        tier, reason = "B", reason + " (ClinVar conflicting — downgraded)"
    if cand.confidence == "low" and tier == "A":
        tier = "B"
        reason += " (capped: low genotype confidence / MZ-twin discordance)"

    cand.tier = tier
    cand.tier_reason = reason
    cand.flags = cand.flags + [f"evidence: {', '.join(ev)}"]
    return cand
