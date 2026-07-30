"""Core data structures for familyforge."""
from __future__ import annotations
from dataclasses import dataclass, field
from typing import Optional


@dataclass
class Member:
    fam: str
    iid: str
    father: str          # "0" if founder/unknown
    mother: str          # "0" if founder/unknown
    sex: str             # "M" | "F" | "U"
    affected: bool
    sequenced: bool = True   # False for members present only to define relationships (e.g. father)


@dataclass
class Family:
    members: dict[str, Member] = field(default_factory=dict)
    mz_twin_groups: list[list[str]] = field(default_factory=list)

    @property
    def sequenced(self) -> list[str]:
        return [m.iid for m in self.members.values() if m.sequenced]

    @property
    def affected(self) -> list[str]:
        return [m.iid for m in self.members.values() if m.sequenced and m.affected]

    @property
    def unaffected(self) -> list[str]:
        return [m.iid for m in self.members.values() if m.sequenced and not m.affected]

    def mothers_of_affected(self) -> set[str]:
        out = set()
        for m in self.members.values():
            if m.sequenced and m.affected and m.mother in self.members:
                out.add(m.mother)
        return out


@dataclass
class GeneAnno:
    """The picked (most-severe, panel/MANE-preferred) VEP transcript consequence."""
    symbol: str = "."
    consequence: str = "."
    hgvsc: str = "."
    hgvsp: str = "."
    protein_pos: Optional[int] = None
    sift: str = "."          # e.g. deleterious(0.01)
    polyphen: str = "."      # e.g. probably_damaging(0.98)
    canonical: bool = False
    mane: bool = False
    mane_select: str = "."   # MANE Select RefSeq accession (e.g. NM_139276.3) for clinical HGVS


@dataclass
class Variant:
    chrom: str
    pos: int
    ref: str
    alt: str
    vid: str                       # rsID or chrom:pos:ref>alt
    filter: str
    gene: GeneAnno
    clinvar: str = "."             # ClinVar CLNSIG
    phylop: Optional[float] = None
    max_af: Optional[float] = None # population popmax across af_fields
    af_detail: dict = field(default_factory=dict)
    paralog_gene: str = "."
    paralog_conf: str = "."
    # externally-supplied finalist in-silico scores (e.g. via --extra-scores), keyed by variant
    cadd: Optional[float] = None
    revel: Optional[float] = None
    am_score: Optional[float] = None
    am_class: str = "."
    # per-sample genotype call: iid -> dict(gt=(a,b)|('.',), zyg, dp, gq, ab, ok)
    geno: dict = field(default_factory=dict)

    @property
    def key(self) -> str:
        return f"{self.chrom}:{self.pos}:{self.ref}>{self.alt}"


@dataclass
class Candidate:
    variant: Variant | None
    gene: str
    model: str                     # de_novo_or_dominant | ar_homozygous | ar_compound_het | x_linked_recessive | cnv_deletion
    tier: str                      # A | B | C
    tier_reason: str
    segregation: str               # human-readable segregation across the family
    zygosity: str
    inheritance_note: str
    confidence: str                # high | moderate | low
    flags: list[str] = field(default_factory=list)
    partner: Variant | None = None # for compound-het (the second variant)
    cnv: dict | None = None        # for CNV candidates
