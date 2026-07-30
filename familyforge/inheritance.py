"""Pedigree-aware Mendelian inheritance models + segregation + MZ-twin concordance,
and genotype-based relatedness / sex QC. Father-not-sequenced is handled by inferring
the paternal allele and labelling it unconfirmed.
"""
from __future__ import annotations
from collections import defaultdict
from .model import Family, Variant, Candidate

CARRY = {"het", "hom_alt", "hemi_alt"}
HOM = {"hom_alt", "hemi_alt"}


def _g(v: Variant, sid: str):
    return v.geno.get(sid, {"zyg": "nocall", "ok": False})


def _carrier(v, sid):
    g = _g(v, sid); return g["zyg"] in CARRY and g["ok"]


def _hom(v, sid):
    g = _g(v, sid); return g["zyg"] in HOM and g["ok"]


def _homref(v, sid):
    return _g(v, sid)["zyg"] == "hom_ref"


def _het(v, sid):
    g = _g(v, sid); return g["zyg"] == "het" and g["ok"]


def _twin_concordant(v, fam: Family):
    for grp in fam.mz_twin_groups:
        zs = {_g(v, s)["zyg"] for s in grp if s in fam.members}
        if len(zs) > 1:
            return False
    return True


def _seg_string(v, fam: Family):
    parts = []
    for iid in fam.sequenced:
        g = _g(v, iid)
        tag = {"hom_ref": "ref", "het": "het", "hom_alt": "hom", "hemi_alt": "hemi",
               "nocall": "nocall"}.get(g["zyg"], g["zyg"])
        aff = "aff" if fam.members[iid].affected else "unaff"
        parts.append(f"{iid}({aff})={tag}")
    return "; ".join(parts)


def evaluate_variant(v: Variant, fam: Family, panel: dict) -> list[Candidate]:
    """Single-variant models (dominant/de novo, AR-hom, X-linked)."""
    gene = v.gene.symbol
    pinfo = panel.get(gene, {})
    inh = pinfo.get("inheritance", set())
    aff = fam.affected
    unaff = fam.unaffected
    cands = []
    seg = _seg_string(v, fam)
    twin_ok = _twin_concordant(v, fam)
    twin_flag = [] if twin_ok else ["MZ-twin discordant genotype — likely technical artifact"]

    all_aff_carry = aff and all(_carrier(v, a) for a in aff)
    all_aff_hom = aff and all(_hom(v, a) for a in aff)

    # ---- de novo / dominant (absent in all sequenced unaffected) ----
    if all_aff_carry and all(_homref(v, u) for u in unaff) and ("AD" in inh or not inh):
        conf = "high" if twin_ok else "low"
        cands.append(Candidate(
            variant=v, gene=gene, model="de_novo_or_dominant", tier="", tier_reason="",
            segregation=seg, zygosity="heterozygous (affected)",
            inheritance_note=("Present in both affected, absent in sequenced mother and sister. "
                              "Father not sequenced → de novo OR paternally inherited (unconfirmed)."),
            confidence=conf, flags=list(twin_flag)))

    # ---- dominant, carrier unaffected mother (reduced penetrance) ----
    elif all_aff_carry and any(_carrier(v, u) for u in unaff) and "AD" in inh:
        carriers = [u for u in unaff if _carrier(v, u)]
        cands.append(Candidate(
            variant=v, gene=gene, model="de_novo_or_dominant", tier="", tier_reason="",
            segregation=seg, zygosity="heterozygous (affected)",
            inheritance_note=(f"Dominant allele also carried by unaffected {', '.join(carriers)} → "
                              "would require reduced penetrance / variable expressivity."),
            confidence="low", flags=twin_flag + ["carried by an unaffected relative (penetrance?)"]))

    # ---- autosomal recessive, homozygous ----
    if all_aff_hom and v.chrom.replace("chr", "") not in ("X", "Y"):
        sister_hom = any(_hom(v, u) for u in unaff)
        mother_ok = True
        mo = fam.mothers_of_affected()
        # obligate-carrier mother expected het (or hom would make her affected)
        moms_seq = [m for m in mo if m in fam.sequenced]
        mom_states = {m: _g(v, m)["zyg"] for m in moms_seq}
        if any(st == "hom_alt" for st in mom_states.values()):
            mother_ok = False  # a hom unaffected mother argues against AR causal
        if not sister_hom and mother_ok and ("AR" in inh or not inh):
            conf = "high" if (twin_ok and all(m == "het" for m in mom_states.values())) else "moderate"
            note = "Homozygous in both affected; "
            note += (", ".join(f"{m} {s}" for m, s in mom_states.items()) or "mother n/a") + \
                    "; father not sequenced (obligate carrier)."
            cands.append(Candidate(
                variant=v, gene=gene, model="ar_homozygous", tier="", tier_reason="",
                segregation=seg, zygosity="homozygous (affected)",
                inheritance_note=note, confidence=conf, flags=list(twin_flag)))

    # ---- X-linked recessive ----
    if v.chrom.replace("chr", "") == "X" and aff:
        males_aff = [a for a in aff if fam.members[a].sex == "M"]
        if males_aff and all(_g(v, a)["zyg"] in HOM for a in males_aff):
            sister_hom = any(_hom(v, u) for u in unaff if fam.members[u].sex == "F")
            mo = [m for m in fam.mothers_of_affected() if m in fam.sequenced]
            mom_carrier = all(_g(v, m)["zyg"] in ("het",) for m in mo) if mo else True
            if not sister_hom and ("XLR" in inh or "XL" in inh or not inh):
                cands.append(Candidate(
                    variant=v, gene=gene, model="x_linked_recessive", tier="", tier_reason="",
                    segregation=seg, zygosity="hemizygous (affected males)",
                    inheritance_note=("Hemizygous in affected son(s); "
                                      + ("carrier mother; " if mom_carrier else "mother genotype check; ")
                                      + "sister not homozygous."),
                    confidence="high" if twin_ok else "low", flags=list(twin_flag)))
    return cands


def evaluate_compound_het(variants: list[Variant], fam: Family, panel: dict) -> list[Candidate]:
    """Gene-level AR compound-het: two variants, both het in all affected, split across
    parents (one maternal, one inferred paternal)."""
    by_gene = defaultdict(list)
    for v in variants:
        if v.chrom.replace("chr", "") in ("X", "Y"):
            continue
        by_gene[v.gene.symbol].append(v)
    aff = fam.affected
    unaff = fam.unaffected
    moms = [m for m in fam.mothers_of_affected() if m in fam.sequenced]
    out = []
    for gene, vs in by_gene.items():
        inh = panel.get(gene, {}).get("inheritance", set())
        if inh and "AR" not in inh:
            continue
        # variants het in ALL affected
        hetvars = [v for v in vs if aff and all(_het(v, a) for a in aff)]
        if len(hetvars) < 2:
            continue
        for i in range(len(hetvars)):
            for j in range(i + 1, len(hetvars)):
                v1, v2 = hetvars[i], hetvars[j]
                # opposite parental origin: one maternal (mother het), other not in mother
                m_v1 = all(_g(v1, m)["zyg"] in ("het", "hom_alt") for m in moms) if moms else None
                m_v2 = all(_g(v2, m)["zyg"] in ("het", "hom_alt") for m in moms) if moms else None
                mref_v1 = all(_homref(v1, m) for m in moms) if moms else None
                mref_v2 = all(_homref(v2, m) for m in moms) if moms else None
                phased_ok = moms and ((m_v1 and mref_v2) or (m_v2 and mref_v1))
                # no unaffected carries BOTH
                unaff_both = any(_carrier(v1, u) and _carrier(v2, u) for u in unaff)
                twin_ok = _twin_concordant(v1, fam) and _twin_concordant(v2, fam)
                if not unaff_both and (phased_ok or not moms):
                    matv, patv = (v1, v2) if m_v1 else (v2, v1)
                    conf = "high" if (phased_ok and twin_ok) else "moderate"
                    note = ("Two heterozygous variants in trans: "
                            + (f"{matv.gene.hgvsc} maternal, {patv.gene.hgvsc} paternal (inferred; "
                               "father not sequenced). " if phased_ok else
                               "parental phase unconfirmed (mother genotypes inconclusive). "))
                    out.append(Candidate(
                        variant=v1, partner=v2, gene=gene, model="ar_compound_het",
                        tier="", tier_reason="", segregation=_seg_string(v1, fam) + " || " + _seg_string(v2, fam),
                        zygosity="compound heterozygous (affected)", inheritance_note=note,
                        confidence=conf,
                        flags=([] if twin_ok else ["MZ-twin discordant — artifact?"])))
    return out


# ─────────────────────────── relatedness / sex QC ───────────────────────────
def relatedness_sex_qc(variants, samples, fam: Family):
    """Genotype-based cross-checks: MZ-twin concordance, parent-offspring sharing,
    and sex inference from chrX heterozygosity. Independent of somalier."""
    zyg = {s: [] for s in samples}          # autosomal zyg strings per sample (quality sites)
    x_het = {s: [0, 0] for s in samples}    # [het, called] on chrX non-PAR
    for v in variants:
        chrom = v.chrom.replace("chr", "")
        for s in samples:
            g = v.geno.get(s)
            if not g or not g["ok"] or g["zyg"] == "nocall":
                continue
            if chrom not in ("X", "Y"):
                zyg[s].append((v.key, g["zyg"]))
            if chrom == "X":
                x_het[s][1] += 1
                if g["zyg"] == "het":
                    x_het[s][0] += 1
    # pairwise genotype concordance over shared quality autosomal sites
    def concordance(a, b):
        da = dict(zyg[a]); db = dict(zyg[b])
        shared = set(da) & set(db)
        if not shared:
            return None, 0
        same = sum(1 for k in shared if da[k] == db[k])
        return same / len(shared), len(shared)

    pairs = {}
    sl = list(samples)
    for i in range(len(sl)):
        for j in range(i + 1, len(sl)):
            c, n = concordance(sl[i], sl[j])
            pairs[f"{sl[i]}~{sl[j]}"] = {"concordance": round(c, 4) if c is not None else None, "n_sites": n}
    sexcall = {}
    for s in samples:
        het, tot = x_het[s]
        rate = (het / tot) if tot else None
        sexcall[s] = {"x_het_rate": round(rate, 4) if rate is not None else None,
                      "x_sites": tot,
                      "inferred_sex": ("M" if (rate is not None and rate < 0.10)
                                       else "F" if rate is not None else "U"),
                      "declared_sex": fam.members.get(s).sex if s in fam.members else "U"}
    return {"pairwise": pairs, "sex": sexcall}
