"""Bakta-name markers for the unknown-IE metals cargo class.

Literature systems (Bellanger 2017 + standard operons):
  Cad — CadA / Cd P-type ATPase; CadD
  Mer — MerA, MerB, MerC, MerE, MerP, MerT
  Ars — ArsA/B/C/D/I/N, ACR3, ArsM
  Cop — CopA/B/C/D/Z and Cu/Ag P-type ATPases
  Cus/Czc — RND cation-efflux of the same cargo type
  Chr — chromate efflux (ChrA)

Not cargo (same rule as C/Xre ≠ RM):
  ArsR / SmtB / MerR / CopG / CopY / CsoR / RicR
  generic HMA / cation transporter / multicopper / cytochrome / polymerase
"""
from __future__ import annotations

import re

# Hard rejects even if a loose token overlaps (arsenosugar, respiration Cu, regs).
EXCLUDE = re.compile(
    r"arsenosugar|"
    r"\baio[ab]\b|"
    r"chromosome|"
    r"cytochrome|"
    r"coproporphyrin|"
    r"\bcopg\b|"
    r"ribbon-helix-helix|"
    r"transcription|"
    r"regulator|"
    r"repressor|"
    r"sensor|"
    r"\bsmtb\b|"
    r"\barsr\b|"
    r"\bmerr\b|"
    r"\bcopy\b|"
    r"\bcutc\b|"
    r"pcu\(a\)c|"
    r"amine oxidase|"
    r"nitrite reductase|"
    r"superoxide dismutase|"
    r"polymerase|"
    r"hma domain|"
    r"heavy-metal-associated domain|"
    r"heavy metal-binding domain|"
    r"heavy metal-responsive|"
    r"heavy metal-dependent|"
    r"heavy metal response|"
    r"multicopper oxidase|"
    r"copper oxidase|"
    r"copper-binding protein|"
    r"blue \(type 1\) copper|"
    r"copper amine|"
    r"chrr-like|"
    r"ca\(2\+\) transporter|"
    r"cora family|"
    r"partition|"
    r"cation transporter$|"
    r"cation-transporting p-type|"
    r"cation efflux family|"
    r"cation efflux protein$",
    re.I,
)

# Structural / enzymatic names only. Regulators are not in this list.
INCLUDE = re.compile(
    r"cadmium|"
    r"\bcad[ad]\b|"
    r"cd\(?2\+?\)?-export|"
    r"lead,\s*cadmium|"
    r"mercur|"
    r"alkylmercury|"
    r"\bmer[abcept]\b|"
    r"arsenit|arsenat|arsenical|"
    r"arsenic resist|"
    r"\bars[abcdhin]\b|"
    r"\bacr3\b|"
    r"copper-transport|"
    r"copper-transloc|"
    r"copper-export|"
    r"copper/silver|"
    r"copper resist|"
    r"copper tolerance|"
    r"p-type cu|"
    r"\bcop[abcdz]\b|"
    r"\bcus[abcfs]\b|"
    r"\bczc[abcd]\b|"
    r"cobalt-zinc-cadmium|"
    r"chromate (transporter|resist|efflux)|"
    r"chromate resistance|"
    r"heavy metal translocating p-type|"
    r"heavy metal transporter|"
    r"cation efflux system protein (cus|czc)",
    re.I,
)

# Strong names win over EXCLUDE (CopD/cytochrome bifunctional, etc.).
STRONG = re.compile(
    r"cadmium-(transport|transloc|resist)|"
    r"putative cadmium-transport|"
    r"cd\(?2\+?\)?-export|"
    r"lead,\s*cadmium|"
    r"mercuric reductase|"
    r"mercury\(ii\) reductase|"
    r"\bmer[abcept]\b|"
    r"alkylmercury|"
    r"arsenate reductase|"
    r"arsenical|"
    r"\bars[abcdhin]\b|"
    r"\bacr3\b|"
    r"\bcop[abcdz]\b|"
    r"\bcus[abcfs]\b|"
    r"\bczc[abcd]\b|"
    r"copper-transport|"
    r"copper-transloc|"
    r"copper-export|"
    r"p-type cu|"
    r"chromate (transporter|resist|efflux)",
    re.I,
)

SUBSYS = [
    ("Cad", re.compile(
        r"cadmium|\bcad[ad]\b|cd\(?2\+?\)?-export|lead,\s*cadmium", re.I
    )),
    ("Mer", re.compile(r"mercur|alkylmercury|\bmer[abcept]\b", re.I)),
    ("Ars", re.compile(
        r"arsenit|arsenat|arsenical|arsenic resist|\bars[abcdhin]\b|\bacr3\b",
        re.I,
    )),
    ("Cop", re.compile(
        r"copper-transport|copper-transloc|copper-export|copper/silver|"
        r"copper resist|copper tolerance|p-type cu|\bcop[abcdz]\b",
        re.I,
    )),
    ("CusCzc", re.compile(
        r"\bcus[abcfs]\b|\bczc[abcd]\b|cobalt-zinc-cadmium|"
        r"cation efflux system protein (cus|czc)",
        re.I,
    )),
    ("Chr", re.compile(r"chromate (transporter|resist|efflux)|chromate resistance", re.I)),
    ("PtypeHMR", re.compile(r"heavy metal translocating p-type|heavy metal transporter", re.I)),
]


def is_metal_cargo(product: str) -> bool:
    p = product or ""
    if STRONG.search(p):
        return True
    if EXCLUDE.search(p):
        return False
    return bool(INCLUDE.search(p))


def metal_subsystems(product: str) -> list[str]:
    if not is_metal_cargo(product):
        return []
    hits = [name for name, rx in SUBSYS if rx.search(product)]
    return hits or ["other"]
