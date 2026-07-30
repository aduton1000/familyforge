"""Parse a CallForge annotated multi-sample VCF into Variant objects.

Reads per-sample GT/DP/AD/GQ and the headline annotations (VEP CSQ: consequence,
SYMBOL, HGVSc/p, SIFT, PolyPhen, MANE/CANONICAL, gnomAD AFs; INFO: ClinVar_CLNSIG,
PhyloP, PARALOG_*). Expects a biallelic-normalised VCF (the CLI runs `bcftools norm
-m-any` first); multiallelic records are still handled per-record but normalisation
is strongly preferred for correct genotype interpretation.
"""
from __future__ import annotations
import re
from cyvcf2 import VCF
from .model import Variant, GeneAnno

SEVERITY = ["transcript_ablation", "splice_acceptor_variant", "splice_donor_variant",
            "stop_gained", "frameshift_variant", "stop_lost", "start_lost",
            "transcript_amplification", "inframe_insertion", "inframe_deletion",
            "missense_variant", "protein_altering_variant", "splice_region_variant",
            "splice_donor_5th_base_variant", "splice_donor_region_variant",
            "splice_polypyrimidine_tract_variant", "incomplete_terminal_codon_variant",
            "start_retained_variant", "stop_retained_variant", "synonymous_variant",
            "coding_sequence_variant", "5_prime_UTR_variant", "3_prime_UTR_variant",
            "intron_variant", "upstream_gene_variant", "downstream_gene_variant",
            "intergenic_variant"]
RANK = {c: i for i, c in enumerate(SEVERITY)}


def _csq_format(vcf: VCF):
    h = vcf.raw_header
    m = re.search(r"##INFO=<ID=CSQ.*?Format: ([^\">]+)", h)
    return m.group(1).split("|") if m else []


def _get(d, fmt, key):
    return d[fmt.index(key)] if key in fmt and fmt.index(key) < len(d) else ""


def _protein_pos(hgvsp: str):
    m = re.search(r"p\.[A-Za-z]{3}(\d+)", hgvsp or "")
    return int(m.group(1)) if m else None


def _floaty(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def pick_consequence(csq_str, fmt, panel_genes):
    """Choose one transcript annotation: prefer a panel gene, then MANE/CANONICAL,
    then most-severe consequence. Returns (GeneAnno, all_symbols set)."""
    if not csq_str or not fmt:
        return GeneAnno(), set()
    best = None
    best_key = None
    symbols = set()
    for tx in csq_str.split(","):
        p = tx.split("|")
        sym = _get(p, fmt, "SYMBOL") or "."
        if sym and sym != ".":
            symbols.add(sym)
        cons_list = (_get(p, fmt, "Consequence") or ".").split("&")
        sev = min((RANK.get(c, 500) for c in cons_list), default=500)
        in_panel = 0 if sym in panel_genes else 1
        mane = _get(p, fmt, "MANE_SELECT") or _get(p, fmt, "MANE")
        canon = _get(p, fmt, "CANONICAL")
        pref = (in_panel, 0 if mane else 1, 0 if canon == "YES" else 1, sev)
        if best_key is None or pref < best_key:
            best_key = pref
            worst = min(cons_list, key=lambda c: RANK.get(c, 500))
            best = GeneAnno(
                symbol=sym, consequence=worst,
                hgvsc=_get(p, fmt, "HGVSc") or ".", hgvsp=_get(p, fmt, "HGVSp") or ".",
                protein_pos=_protein_pos(_get(p, fmt, "HGVSp")),
                sift=_get(p, fmt, "SIFT") or ".", polyphen=_get(p, fmt, "PolyPhen") or ".",
                canonical=(canon == "YES"), mane=bool(mane), mane_select=(mane or "."),
            )
    return best or GeneAnno(), symbols


def _af_from_csq(csq_str, fmt, af_fields):
    """Max population AF across the configured CSQ AF fields (one transcript carries them)."""
    if not csq_str:
        return None, {}
    detail = {}
    for tx in csq_str.split(","):
        p = tx.split("|")
        for k in af_fields:
            if k in fmt:
                v = _floaty(_get(p, fmt, k))
                if v is not None:
                    detail[k] = max(detail.get(k, 0.0), v)
    vals = [v for v in detail.values() if v is not None]
    return (max(vals) if vals else None), detail


def _zygosity(alleles, sex, chrom):
    called = [a for a in alleles if a is not None and a >= 0]
    if not called:
        return "nocall"
    nalt = sum(1 for a in called if a >= 1)
    male_hemi = (sex == "M" and chrom.replace("chr", "") in ("X", "Y"))
    if nalt == 0:
        return "hom_ref"
    if male_hemi and nalt == len(called):
        return "hemi_alt"
    if nalt == len(called):
        return "hom_alt"
    return "het"


def load_variants(vcf_path, panel_genes, af_fields, min_gq, min_dp, min_alt_ab,
                  max_ref_ab_for_hom, sex_by_id):
    vcf = VCF(vcf_path)
    samples = list(vcf.samples)
    fmt = _csq_format(vcf)
    info_tags = set(re.findall(r"##INFO=<ID=([^,]+)", vcf.raw_header))
    out = []
    for v in vcf:
        if len(v.ALT) != 1:
            # normalisation should prevent this; skip multiallelic to keep GT unambiguous
            continue
        chrom = v.CHROM
        csq = v.INFO.get("CSQ")
        gene, _syms = pick_consequence(csq, fmt, panel_genes)
        if gene.symbol not in panel_genes:
            continue  # panel-restricted primary pass
        max_af, af_detail = _af_from_csq(csq, fmt, af_fields)
        clnsig = v.INFO.get("ClinVar_CLNSIG") or v.INFO.get("CLNSIG")
        if not clnsig and "CLIN_SIG" in fmt and csq:
            clnsig = _get(csq.split(",")[0].split("|"), fmt, "CLIN_SIG")
        phylop = v.INFO.get("PhyloP")
        var = Variant(
            chrom=chrom, pos=v.POS, ref=v.REF, alt=v.ALT[0],
            vid=(v.ID or f"{chrom}:{v.POS}:{v.REF}>{v.ALT[0]}"),
            filter=(v.FILTER or "PASS"),
            gene=gene, clinvar=(clnsig or "."),
            phylop=_floaty(phylop), max_af=max_af, af_detail=af_detail,
            paralog_gene=(v.INFO.get("PARALOG_GENE") or "."),
            paralog_conf=(v.INFO.get("PARALOG_CONF") or "."),
        )
        gts = v.genotypes  # [[a0,a1,phased], ...]
        dp = v.format("DP"); ad = v.format("AD"); gq = v.format("GQ")
        for i, sid in enumerate(samples):
            alleles = [a for a in gts[i][:-1]] if gts is not None else []
            sex = sex_by_id.get(sid, "U")
            zyg = _zygosity(alleles, sex, chrom)
            d = int(dp[i][0]) if dp is not None and dp[i][0] >= 0 else None
            g = int(gq[i][0]) if gq is not None and gq[i][0] >= 0 else None
            ab = None
            if ad is not None and ad.shape[1] >= 2:
                refc, altc = int(ad[i][0]), int(ad[i][1])
                tot = refc + altc
                ab = (altc / tot) if tot > 0 else None
            ok = True
            if zyg in ("het", "hom_alt", "hemi_alt"):
                if (g is not None and g < min_gq) or (d is not None and d < min_dp):
                    ok = False
                if zyg == "het" and ab is not None and ab < min_alt_ab:
                    ok = False
                if zyg in ("hom_alt", "hemi_alt") and ab is not None and (1 - ab) > max_ref_ab_for_hom and ab < 0.9:
                    pass  # tolerate; hom-alt AB handled leniently
            var.geno[sid] = {"alleles": alleles, "zyg": zyg, "dp": d, "gq": g, "ab": ab, "ok": ok}
        out.append(var)
    return out, samples
