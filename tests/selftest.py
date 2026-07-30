#!/usr/bin/env python3
"""Self-test: synthesize an annotated VCF with planted HIES variants across the four
family members and confirm familyforge recovers the expected models/tiers."""
import os, sys, tempfile
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from familyforge import cli

SAMPLES = ["MOTHER", "SIBLING", "TWIN_1", "TWIN_2"]
CSQFMT = "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature|BIOTYPE|HGVSc|HGVSp|CANONICAL|MANE_SELECT|SIFT|PolyPhen|gnomADe_AF|gnomADg_AF|MAX_AF"

HEADER = f"""##fileformat=VCFv4.2
##contig=<ID=chr17>
##contig=<ID=chr9>
##contig=<ID=chr6>
##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: {CSQFMT}">
##INFO=<ID=ClinVar_CLNSIG,Number=.,Type=String,Description="ClinVar">
##INFO=<ID=PhyloP,Number=1,Type=Float,Description="PhyloP">
##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="AD">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="DP">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="GQ">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t""" + "\t".join(SAMPLES)


def gt(g, alt=True):
    ad = {"0/0": "40,0", "0/1": "20,20", "1/1": "0,40"}[g]
    return f"{g}:{ad}:40:99"


def csq(allele, cons, sym, hgvsc, hgvsp, sift, pp, af):
    return (f"{allele}|{cons}|HIGH|{sym}|ENSG|ENST|protein_coding|{hgvsc}|{hgvsp}|YES|NM_1|"
            f"{sift}|{pp}|{af}|{af}|{af}")


def row(chrom, pos, ref, alt, info_csq, clnsig, phylop, gts):
    info = f"CSQ={info_csq};PhyloP={phylop}"
    if clnsig:
        info += f";ClinVar_CLNSIG={clnsig}"
    cells = "\t".join(gt(g) for g in gts)
    return f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t5000\tPASS\t{info}\tGT:AD:DP:GQ\t{cells}"


def build_vcf(path):
    lines = [HEADER]
    # 1) STAT3 de-novo hotspot missense R382W: affected het, parents/sib ref -> Tier A dominant
    lines.append(row("chr17", 42322000, "G", "A",
        csq("A", "missense_variant", "STAT3", "c.1144C>T", "p.Arg382Trp", "deleterious(0)", "probably_damaging(0.99)", ""),
        "Pathogenic", 5.2, ["0/0", "0/0", "0/1", "0/1"]))
    # 2) DOCK8 homozygous stop_gained: affected 1/1, mother+sister het -> Tier A AR-hom
    lines.append(row("chr9", 300000, "C", "T",
        csq("T", "stop_gained", "DOCK8", "c.100C>T", "p.Arg34Ter", "", "", ""),
        "", 3.1, ["0/1", "0/1", "1/1", "1/1"]))
    # 3a) PGM3 comp-het allele A (maternal): mother het, affected het, sister ref
    lines.append(row("chr6", 83150000, "C", "T",
        csq("T", "stop_gained", "PGM3", "c.200C>T", "p.Gln67Ter", "", "", ""),
        "", 4.0, ["0/1", "0/0", "0/1", "0/1"]))
    # 3b) PGM3 comp-het allele B (paternal, inferred): mother ref, affected het, sister ref
    lines.append(row("chr6", 83155000, "G", "GA",
        csq("GA", "frameshift_variant", "PGM3", "c.300dup", "p.Thr101fs", "", "", ""),
        "", 4.0, ["0/0", "0/0", "0/1", "0/1"]))
    # 4) common benign STAT3 variant -> must be filtered (Tier C / dropped)
    lines.append(row("chr17", 42330000, "A", "G",
        csq("G", "missense_variant", "STAT3", "c.50A>G", "p.Lys17Arg", "tolerated(0.4)", "benign(0.01)", "0.25"),
        "Benign", 0.1, ["0/1", "0/1", "0/1", "0/0"]))
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")


PED = """FAM01 FATHER 0 0 1 0
FAM01 MOTHER 0 0 2 1
FAM01 SIBLING FATHER MOTHER 2 1
FAM01 TWIN_1 FATHER MOTHER 1 2
FAM01 TWIN_2 FATHER MOTHER 1 2
"""

PANEL = """gene\ttier\tinheritance\tmechanism\tomim\tnotes
STAT3\t1\tAD\tDN\t147060\tDBD/SH2 hotspot incl R382
DOCK8\t1\tAR\tLOF\t243700\tbiallelic LOF, large deletions
PGM3\t2\tAR\tHYPO\t615816\thypomorphic/LOF
"""


def main():
    d = tempfile.mkdtemp(prefix="ff_selftest_")
    vcf = os.path.join(d, "annotated.vcf")
    ped = os.path.join(d, "family.ped")
    panel = os.path.join(d, "panel.tsv")
    build_vcf(vcf)
    open(ped, "w").write(PED)
    open(panel, "w").write(PANEL)
    out = os.path.join(d, "out")
    rc = cli.main(["--vcf", vcf, "--ped", ped, "--panel", panel, "--outdir", out,
                   "--study-name", "SELFTEST", "--mz-twins", "TWIN_1,TWIN_2", "--skip-norm"])
    import json
    res = json.load(open(os.path.join(out, "candidates.json")))
    genes_tierA = {r["gene"] for r in res["rows"] if r["tier"] == "A"}
    print("\n=== SELFTEST RESULT ===")
    print("Tier A genes:", genes_tierA)
    print("summary:", res["context_summary"])
    print("QC pairwise:", json.dumps(res["qc"]["pairwise"], indent=0))
    ok = True
    for g in ("STAT3", "DOCK8", "PGM3"):
        if g not in genes_tierA:
            print(f"  FAIL: expected {g} in Tier A"); ok = False
    # MZ twins TWIN_1~TWIN_2 concordance should be 1.0
    aly = res["qc"]["pairwise"].get("TWIN_1~TWIN_2", {}).get("concordance")
    if aly != 1.0:
        print(f"  WARN: TWIN_1~TWIN_2 concordance {aly} (expected 1.0)")
    reports = ["clinician_report.html", "clinician_report.md",
               "publication_report.html", "publication_report.md"]
    for r in reports:
        p = os.path.join(out, r)
        if not (os.path.isfile(p) and os.path.getsize(p) > 200):
            print(f"  FAIL: report missing/empty {r}"); ok = False
    print("OUTPUT DIR:", out)
    print("SELFTEST", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
