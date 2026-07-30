"""familyforge — pedigree-aware Mendelian + ACMG-style interpretation and dual
(clinician + publication) reporting on a CallForge annotated joint callset.

Standalone companion to CallForge (analysis) and CaptureForge (design). It consumes
a multi-sample annotated VCF + a PED pedigree + a gene panel and emits ranked,
segregation-checked candidate variants and two audience-specific reports.
Disease-agnostic: the panel and pedigree are inputs, never hardcoded.
"""
__version__ = "0.1.0"
