"""Loaders: pedigree (PED), gene panel (TSV), config (YAML), CNV calls (TSV)."""
from __future__ import annotations
import csv
import os
import yaml
from .model import Family, Member

_SEXMAP = {"1": "M", "2": "F", "0": "U", "m": "M", "f": "F",
           "male": "M", "female": "F", "u": "U", "unknown": "U"}


def load_config(path: str | None) -> dict:
    here = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(here, "config", "default.yaml")) as fh:
        cfg = yaml.safe_load(fh)
    if path:
        with open(path) as fh:
            user = yaml.safe_load(fh) or {}
        cfg.update(user)
    return cfg


def load_pedigree(path: str, sequenced_ids: set[str] | None = None,
                  mz_twin_groups: list[list[str]] | None = None) -> Family:
    """Parse a 6-column PED. A member is `sequenced` iff its iid is in the VCF
    (sequenced_ids); otherwise it is relationship-only (e.g. an unsequenced father)."""
    fam = Family(mz_twin_groups=mz_twin_groups or [])
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line.strip() or line.startswith("#"):
                continue
            f = line.split()
            if len(f) < 6:
                continue
            famid, iid, pat, mat, sex, pheno = f[:6]
            seq = True if sequenced_ids is None else (iid in sequenced_ids)
            fam.members[iid] = Member(
                fam=famid, iid=iid, father=pat, mother=mat,
                sex=_SEXMAP.get(sex.lower(), "U"),
                affected=(pheno.strip() == "2"),
                sequenced=seq,
            )
    return fam


def load_panel(path: str) -> dict[str, dict]:
    """HIES-style panel TSV -> {GENE: {tier,inheritance:set,mechanism:set,notes}}."""
    panel: dict[str, dict] = {}
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = [x.strip() for x in line.rstrip("\n").split("\t")]
            if f[0].lower() == "gene":
                continue
            gene = f[0]
            panel[gene] = {
                "tier": f[1] if len(f) > 1 else "3",
                "inheritance": {x.strip().upper() for x in (f[2] if len(f) > 2 else "").split(",") if x.strip()},
                "mechanism": {x.strip().upper() for x in (f[3] if len(f) > 3 else "").split(",") if x.strip()},
                "omim": f[4] if len(f) > 4 else "",
                "notes": f[5] if len(f) > 5 else "",
            }
    return panel


def load_extra_scores(path: str | None) -> dict[str, dict]:
    """Externally-retrieved finalist in-silico scores (genome-wide dbNSFP unavailable).
    TSV with header: key<TAB>cadd_phred<TAB>revel<TAB>alphamissense_score<TAB>alphamissense_class
    where key = 'chrom:pos:ref>alt' (matches Variant.key). Returns {key: {...}}."""
    if not path or not os.path.isfile(path):
        return {}
    out: dict[str, dict] = {}
    with open(path) as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        for r in rdr:
            key = (r.get("key") or "").strip()
            if not key:
                continue
            def num(x):
                try:
                    return float(x)
                except (TypeError, ValueError):
                    return None
            out[key] = {
                "cadd": num(r.get("cadd_phred")),
                "revel": num(r.get("revel")),
                "am_score": num(r.get("alphamissense_score")),
                "am_class": (r.get("alphamissense_class") or ".").strip() or ".",
            }
    return out


def load_cnv_calls(path: str | None) -> list[dict]:
    """Parse a CallForge cnv_calls.tsv (tolerant of column-name variations)."""
    if not path or not os.path.isfile(path):
        return []
    rows: list[dict] = []
    with open(path) as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        for r in rdr:
            rows.append({k: (v or "").strip() for k, v in r.items()})
    return rows
