"""CNV cross-check: surface deletions overlapping panel genes (esp. DOCK8, where large
deletions are a common AR-HIES mechanism that SNV calling misses). Schema-tolerant to
CallForge cnv_calls.tsv column-name variations.
"""
from __future__ import annotations
from .model import Candidate


def _find(row: dict, *names):
    low = {k.lower(): k for k in row}
    for n in names:
        if n in low:
            return row[low[n]]
    return ""


def _is_deletion(row: dict):
    call = _find(row, "call", "type", "cnv", "event", "status").lower()
    if any(t in call for t in ("del", "loss", "deletion")):
        return True
    cn = _find(row, "cn", "copy_number", "copynumber", "ncopies")
    try:
        return float(cn) < 1.5
    except (TypeError, ValueError):
        pass
    log2 = _find(row, "log2", "log2ratio", "logratio")
    try:
        return float(log2) < -0.4
    except (TypeError, ValueError):
        return False


def cnv_candidates(cnv_rows, panel: dict, fam, cfg: dict) -> list[Candidate]:
    priority = set(cfg.get("cnv_priority_genes", []))
    aff = set(fam.affected)
    out = []
    for r in cnv_rows:
        genes_field = _find(r, "gene", "genes", "symbol", "overlapping_genes")
        genes = {g.strip() for part in genes_field.replace(";", ",").split(",") for g in [part] if g.strip()}
        hit = genes & set(panel)
        if not hit or not _is_deletion(r):
            continue
        sample = _find(r, "sample", "sample_id", "sid")
        chrom = _find(r, "chrom", "chromosome", "chr")
        start = _find(r, "start", "begin")
        end = _find(r, "end", "stop")
        cn = _find(r, "cn", "copy_number", "log2", "log2ratio")
        callability = _find(r, "callability", "confidence", "label", "cnv_callable")
        for gene in hit:
            in_aff = (sample in aff) or (sample == "")
            flags = []
            if gene in priority:
                flags.append(f"{gene}: large deletions are a classic AR-HIES mechanism (SNV-blind)")
            if callability and "low" in callability.lower():
                flags.append(f"low-confidence CNV call ({callability})")
            tier = "A" if (gene in priority and in_aff) else "B"
            out.append(Candidate(
                variant=None, gene=gene, model="cnv_deletion",
                tier=tier, tier_reason=f"deletion overlapping {gene} in {sample or 'sample(s)'} "
                                       f"({'affected' if in_aff else 'unaffected'})",
                segregation=f"{sample}: deletion {chrom}:{start}-{end} (cn/log2={cn})",
                zygosity="deletion", inheritance_note="Confirm with read-depth/IGV and orthogonal assay (MLPA/aCGH).",
                confidence="moderate", flags=flags,
                cnv={"sample": sample, "chrom": chrom, "start": start, "end": end,
                     "cn": cn, "callability": callability, "gene": gene}))
    return out
