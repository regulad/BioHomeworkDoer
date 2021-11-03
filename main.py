"""
Copyright regulad 2021
"""

from typing import List, Dict
import argparse
from enum import Enum, auto


class NucleicAcidType(Enum):
    """Defines a type of Nucleic Acid."""

    DNA = auto()
    mRNA = auto()
    tRNA = auto()


class Protein(Enum):
    """Defines a protein."""

    MET = "Start"
    STOP = "Stop"
    TRP = "Trp"
    PHE = "Phe"
    LEU = "Leu"
    LLE = "Lle"
    VAL = "Val"
    SER = "Ser"
    PRO = "Pro"
    THR = "Thr"
    ALA = "Ala"
    TYR = "Tyr"
    HIS = "His"
    GLN = "Gln"
    ASN = "Asn"
    LYS = "Lys"
    ASP = "Asp"
    GLU = "Glu"
    CYS = "Cys"
    ARG = "Arg"
    GLY = "Gly"


mrna_protein_map: Dict[str, Protein] = {  # 3-Letter combos are special-cased
    "UUU": Protein.PHE,
    "UUC": Protein.PHE,
    "UU": Protein.LEU,
    "CU": Protein.LEU,
    "AUG": Protein.MET,
    "AU": Protein.LLE,
    "GU": Protein.VAL,
    "UC": Protein.SER,
    "CC": Protein.PRO,
    "AC": Protein.THR,
    "GC": Protein.ALA,
    "UAA": Protein.STOP,
    "UAG": Protein.STOP,
    "UA": Protein.TYR,
    "CAU": Protein.HIS,
    "CAC": Protein.HIS,
    "CA": Protein.GLN,
    "AAU": Protein.ASN,
    "AAC": Protein.ASN,
    "AA": Protein.LYS,
    "GAU": Protein.ASP,
    "GAC": Protein.ASP,
    "GA": Protein.GLU,
    "UGA": Protein.STOP,
    "UGG": Protein.TRP,
    "UG": Protein.CYS,
    "CG": Protein.ARG,
    "AGU": Protein.SER,
    "AGC": Protein.SER,
    "AG": Protein.ARG,
    "GG": Protein.GLY,
}

parser: argparse.ArgumentParser = argparse.ArgumentParser(description="Translates some DNA to/from mRNA and tRNA amino acids.")
parser.add_argument("input", help="The input nucleic acid base string.", type=str)
parser.add_argument("-i", dest="inputtype", help="Describes the type of the input. Options: DNA, mRNA.", default="DNA")
parser.add_argument("-o", dest="outputtype", help="Describes the type of the output. Options: DNA, mRNA, tRNA.", default="mRNA")
args: argparse.Namespace = parser.parse_args()

match args.inputtype.upper():
    case "DNA":
        input_nucleic_acid_type = NucleicAcidType.DNA
    case "MRNA":
        input_nucleic_acid_type = NucleicAcidType.mRNA
    case _:
        raise RuntimeError("Bad inputtype")

match args.outputtype.upper():
    case "DNA":
        output_nucleic_acid_type = NucleicAcidType.DNA
    case "MRNA":
        output_nucleic_acid_type = NucleicAcidType.mRNA
    case "TRNA":
        output_nucleic_acid_type = NucleicAcidType.tRNA
    case _:
        raise RuntimeError("Bad outputtype")

if input_nucleic_acid_type is NucleicAcidType.DNA:
    mrna_intermediate: str = args.input.upper().translate({65: 85, 84: 65, 71: 67, 67: 71})
else:
    mrna_intermediate: str = args.input.upper()

if output_nucleic_acid_type is NucleicAcidType.DNA:
    print(mrna_intermediate.translate({85: 65, 65: 84, 67: 71, 71: 67}))
elif output_nucleic_acid_type is NucleicAcidType.tRNA:
    new_codons: List[str] = [codon[:3] if len(codon) > 3 else codon for codon in mrna_intermediate.split()]
    for codon in new_codons:
        if len(codon) != 3:
            raise RuntimeError(f"codon of incorrect length: {codon}")
    proteins: List[Protein] = []
    for codon in new_codons:
        for known_codon, protein in mrna_protein_map.items():
            if codon.startswith(known_codon):
                proteins.append(protein)
                break
    print(" ".join([protein.value for protein in proteins]))
else:
    print(mrna_intermediate)
