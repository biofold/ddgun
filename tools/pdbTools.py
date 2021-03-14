#!/usr/bin/env python
'''
This package contains three classes RESIDUE,
PDBChain, and PDBProtein.
RESIDUE => class residue
PDBChain => only one chain is treated
PDBProtein => multiple PDBChains
This file also contains functions that read a pdb file
and return a PDBChain or PDBProtein.
'''


import sys
import string 
import numpy

TotAtomType=['C', 'CA', 'CB', 'CD', 'CD1', 'CD2', 'CE', 'CE1', 'CE2', 'CE3', 'CG', 'CG1', 'CG2', 'CH2', 'CZ', 'CZ2', 'CZ3', 'H', 'HA', 'HB', 'HD', 'HD1', 'HD2', 'HE', 'HE1', 'HE2', 'HE3', 'HG', 'HG1', 'HG2', 'HH', 'HH1', 'HH2', 'HZ', 'HZ2', 'HZ3', 'N', 'ND1', 'ND2', 'NE', 'NE1', 'NE2', 'NH1', 'NH2', 'NZ', 'O', 'OD1', 'OD2', 'OE1', 'OE2', 'OG', 'OG1', 'OH', 'OXT', 'SD', 'SG']

HeavyAtomType=['C', 'CA', 'CB', 'CD', 'CD1', 'CD2', 'CE', 'CE1', 'CE2', 'CE3', 'CG', 'CG1', 'CG2', 'CH2', 'CZ', 'CZ2', 'CZ3', 'N', 'ND1', 'ND2', 'NE', 'NE1', 'NE2', 'NH1', 'NH2', 'NZ', 'O', 'OD1', 'OD2', 'OE1', 'OE2', 'OG', 'OG1', 'OH', 'OXT', 'SD', 'SG']


# Dictionary that contains the relationship
# between A3 and A1 residue names
def _getResidueA1(res3):
   ''' 
   _getResidueA1(res3):
   returns single residue of 'res3' a 3-residue code
   '''
   resa3Toa1= {"00C": "C", "01W": "X", "02K": "A", "03Y": "C", "07O": "C",
      "08P": "C", "0A0": "D", "0A1": "Y", "0A2": "K", "0A8": "C",
      "0AA": "V", "0AB": "V", "0AC": "G", "0AD": "G", "0AF": "W",
      "0AG": "L", "0AH": "S", "0AK": "D", "0AM": "A", "0AP": "C",
      "0AU": "U", "0AV": "A", "0AZ": "P", "0BN": "F", "0C ": "C",
      "0CS": "A", "0DC": "C", "0DG": "G", "0DT": "T", "0FL": "A",
      "0G ": "G", "0NC": "A", "0SP": "A", "0U ": "U", "0YG": "YG",
      "10C": "C", "125": "U", "126": "U", "127": "U", "128": "N",
      "12A": "A", "143": "C", "175": "ASG", "193": "X", "1AP": "A",
      "1MA": "A", "1MG": "G", "1PA": "F", "1PI": "A", "1PR": "N",
      "1SC": "C", "1TQ": "W", "1TY": "Y", "1X6": "S", "200": "F",
      "23F": "F", "23S": "X", "26B": "T", "2AD": "X", "2AG": "A",
      "2AO": "X", "2AR": "A", "2AS": "X", "2AT": "T", "2AU": "U",
      "2BD": "I", "2BT": "T", "2BU": "A", "2CO": "C", "2DA": "A",
      "2DF": "N", "2DM": "N", "2DO": "X", "2DT": "T", "2EG": "G",
      "2FE": "N", "2FI": "N", "2FM": "M", "2GT": "T", "2HF": "H",
      "2LU": "L", "2MA": "A", "2MG": "G", "2ML": "L", "2MR": "R",
      "2MT": "P", "2MU": "U", "2NT": "T", "2OM": "U", "2OT": "T",
      "2PI": "X", "2PR": "G", "2SA": "N", "2SI": "X", "2ST": "T",
      "2TL": "T", "2TY": "Y", "2VA": "V", "2XA": "C", "32S": "X",
      "32T": "X", "3AH": "H", "3AR": "X", "3CF": "F", "3DA": "A",
      "3DR": "N", "3GA": "A", "3MD": "D", "3ME": "U", "3NF": "Y",
      "3QN": "K", "3TY": "X", "3XH": "G", "4AC": "N", "4BF": "Y",
      "4CF": "F", "4CY": "M", "4DP": "W", "4F3": "GYG", "4FB": "P",
      "4FW": "W", "4HT": "W", "4IN": "W", "4MF": "N", "4MM": "X",
      "4OC": "C", "4PC": "C", "4PD": "C", "4PE": "C", "4PH": "F",
      "4SC": "C", "4SU": "U", "4TA": "N", "4U7": "A", "56A": "H",
      "5AA": "A", "5AB": "A", "5AT": "T", "5BU": "U", "5CG": "G",
      "5CM": "C", "5CS": "C", "5FA": "A", "5FC": "C", "5FU": "U",
      "5HP": "E", "5HT": "T", "5HU": "U", "5IC": "C", "5IT": "T",
      "5IU": "U", "5MC": "C", "5MD": "N", "5MU": "U", "5NC": "C",
      "5PC": "C", "5PY": "T", "5SE": "U", "5ZA": "TWG", "64T": "T",
      "6CL": "K", "6CT": "T", "6CW": "W", "6HA": "A", "6HC": "C",
      "6HG": "G", "6HN": "K", "6HT": "T", "6IA": "A", "6MA": "A",
      "6MC": "A", "6MI": "N", "6MT": "A", "6MZ": "N", "6OG": "G",
      "70U": "U", "7DA": "A", "7GU": "G", "7JA": "I", "7MG": "G",
      "8AN": "A", "8FG": "G", "8MG": "G", "8OG": "G", "9NE": "E",
      "9NF": "F", "9NR": "R", "9NV": "V", "A  ": "A", "A1P": "N",
      "A23": "A", "A2L": "A", "A2M": "A", "A34": "A", "A35": "A",
      "A38": "A", "A39": "A", "A3A": "A", "A3P": "A", "A40": "A",
      "A43": "A", "A44": "A", "A47": "A", "A5L": "A", "A5M": "C",
      "A5N": "N", "A5O": "A", "A66": "X", "AA3": "A", "AA4": "A",
      "AAR": "R", "AB7": "X", "ABA": "A", "ABR": "A", "ABS": "A",
      "ABT": "N", "ACB": "D", "ACL": "R", "AD2": "A", "ADD": "X",
      "ADX": "N", "AEA": "X", "AEI": "D", "AET": "A", "AFA": "N",
      "AFF": "N", "AFG": "G", "AGM": "R", "AGT": "C", "AHB": "N",
      "AHH": "X", "AHO": "A", "AHP": "A", "AHS": "X", "AHT": "X",
      "AIB": "A", "AKL": "D", "AKZ": "D", "ALA": "A", "ALC": "A",
      "ALM": "A", "ALN": "A", "ALO": "T", "ALQ": "X", "ALS": "A",
      "ALT": "A", "ALV": "A", "ALY": "K", "AN8": "A", "AP7": "A",
      "APE": "X", "APH": "A", "API": "K", "APK": "K", "APM": "X",
      "APP": "X", "AR2": "R", "AR4": "E", "AR7": "R", "ARG": "R",
      "ARM": "R", "ARO": "R", "ARV": "X", "AS ": "A", "AS2": "D",
      "AS9": "X", "ASA": "D", "ASB": "D", "ASI": "D", "ASK": "D",
      "ASL": "D", "ASM": "X", "ASN": "N", "ASP": "D", "ASQ": "D",
      "ASU": "N", "ASX": "B", "ATD": "T", "ATL": "T", "ATM": "T",
      "AVC": "A", "AVN": "X", "AYA": "A", "AYG": "AYG", "AZK": "K",
      "AZS": "S", "AZY": "Y", "B1F": "F", "B1P": "N", "B2A": "A",
      "B2F": "F", "B2I": "I", "B2V": "V", "B3A": "A", "B3D": "D",
      "B3E": "E", "B3K": "K", "B3L": "X", "B3M": "X", "B3Q": "X",
      "B3S": "S", "B3T": "X", "B3U": "H", "B3X": "N", "B3Y": "Y",
      "BB6": "C", "BB7": "C", "BB8": "F", "BB9": "C", "BBC": "C",
      "BCS": "C", "BE2": "X", "BFD": "D", "BG1": "S", "BGM": "G",
      "BH2": "D", "BHD": "D", "BIF": "F", "BIL": "X", "BIU": "I",
      "BJH": "X", "BLE": "L", "BLY": "K", "BMP": "N", "BMT": "T",
      "BNN": "F", "BNO": "X", "BOE": "T", "BOR": "R", "BPE": "C",
      "BRU": "U", "BSE": "S", "BT5": "N", "BTA": "L", "BTC": "C",
      "BTR": "W", "BUC": "C", "BUG": "V", "BVP": "U", "BZG": "N",
      "C  ": "C", "C12": "TYG", "C1X": "K", "C25": "C", "C2L": "C",
      "C2S": "C", "C31": "C", "C32": "C", "C34": "C", "C36": "C",
      "C37": "C", "C38": "C", "C3Y": "C", "C42": "C", "C43": "C",
      "C45": "C", "C46": "C", "C49": "C", "C4R": "C", "C4S": "C",
      "C5C": "C", "C66": "X", "C6C": "C", "C99": "TFG", "CAF": "C",
      "CAL": "X", "CAR": "C", "CAS": "C", "CAV": "X", "CAY": "C",
      "CB2": "C", "CBR": "C", "CBV": "C", "CCC": "C", "CCL": "K",
      "CCS": "C", "CCY": "CYG", "CDE": "X", "CDV": "X", "CDW": "C",
      "CEA": "C", "CFL": "C", "CFY": "FCYG", "CG1": "G", "CGA": "E",
      "CGU": "E", "CH ": "C", "CH6": "MYG", "CH7": "KYG", "CHF": "X",
      "CHG": "X", "CHP": "G", "CHS": "X", "CIR": "R", "CJO": "GYG",
      "CLE": "L", "CLG": "K", "CLH": "K", "CLV": "AFG", "CM0": "N",
      "CME": "C", "CMH": "C", "CML": "C", "CMR": "C", "CMT": "C",
      "CNU": "U", "CP1": "C", "CPC": "X", "CPI": "X", "CQR": "GYG",
      "CR0": "TLG", "CR2": "GYG", "CR5": "G", "CR7": "KYG", "CR8": "HYG",
      "CRF": "TWG", "CRG": "THG", "CRK": "MYG", "CRO": "GYG", "CRQ": "QYG",
      "CRU": "EYG", "CRW": "ASG", "CRX": "ASG", "CS0": "C", "CS1": "C",
      "CS3": "C", "CS4": "C", "CS8": "N", "CSA": "C", "CSB": "C",
      "CSD": "C", "CSE": "C", "CSF": "C", "CSH": "SHG", "CSI": "G",
      "CSJ": "C", "CSL": "C", "CSO": "C", "CSP": "C", "CSR": "C",
      "CSS": "C", "CSU": "C", "CSW": "C", "CSX": "C", "CSY": "SYG",
      "CSZ": "C", "CTE": "W", "CTG": "T", "CTH": "T", "CUC": "X",
      "CWR": "S", "CXM": "M", "CY0": "C", "CY1": "C", "CY3": "C",
      "CY4": "C", "CYA": "C", "CYD": "C", "CYF": "C", "CYG": "C",
      "CYJ": "X", "CYM": "C", "CYQ": "C", "CYR": "C", "CYS": "C",
      "CZ2": "C", "CZO": "GYG", "CZZ": "C", "D11": "T", "D1P": "N",
      "D3 ": "N", "D33": "N", "D3P": "G", "D3T": "T", "D4M": "T",
      "D4P": "X", "DA ": "A", "DA2": "X", "DAB": "A", "DAH": "F",
      "DAL": "A", "DAR": "R", "DAS": "D", "DBB": "T", "DBM": "N",
      "DBS": "S", "DBU": "T", "DBY": "Y", "DBZ": "A", "DC ": "C",
      "DC2": "C", "DCG": "G", "DCI": "X", "DCL": "X", "DCT": "C",
      "DCY": "C", "DDE": "H", "DDG": "G", "DDN": "U", "DDX": "N",
      "DFC": "C", "DFG": "G", "DFI": "X", "DFO": "X", "DFT": "N",
      "DG ": "G", "DGH": "G", "DGI": "G", "DGL": "E", "DGN": "Q",
      "DHA": "S", "DHI": "H", "DHL": "X", "DHN": "V", "DHP": "X",
      "DHU": "U", "DHV": "V", "DI ": "I", "DIL": "I", "DIR": "R",
      "DIV": "V", "DLE": "L", "DLS": "K", "DLY": "K", "DM0": "K",
      "DMH": "N", "DMK": "D", "DMT": "X", "DN ": "N", "DNE": "L",
      "DNG": "L", "DNL": "K", "DNM": "L", "DNP": "A", "DNR": "C",
      "DNS": "K", "DOA": "X", "DOC": "C", "DOH": "D", "DON": "L",
      "DPB": "T", "DPH": "F", "DPL": "P", "DPP": "A", "DPQ": "Y",
      "DPR": "P", "DPY": "N", "DRM": "U", "DRP": "N", "DRT": "T",
      "DRZ": "N", "DSE": "S", "DSG": "N", "DSN": "S", "DSP": "D",
      "DT ": "T", "DTH": "T", "DTR": "W", "DTY": "Y", "DU ": "U",
      "DVA": "V", "DXD": "N", "DXN": "N", "DYG": "DYG", "DYS": "C",
      "DZM": "A", "E  ": "A", "E1X": "A", "ECC": "Q", "EDA": "A",
      "EFC": "C", "EHP": "F", "EIT": "T", "ENP": "N", "ESB": "Y",
      "ESC": "M", "EXB": "X", "EXY": "L", "EY5": "N", "EYS": "X",
      "F2F": "F", "FA2": "A", "FA5": "N", "FAG": "N", "FAI": "N",
      "FB5": "A", "FB6": "A", "FCL": "F", "FFD": "N", "FGA": "E",
      "FGL": "G", "FGP": "S", "FHL": "X", "FHO": "K", "FHU": "U",
      "FLA": "A", "FLE": "L", "FLT": "Y", "FME": "M", "FMG": "G",
      "FMU": "N", "FOE": "C", "FOX": "G", "FP9": "P", "FPA": "F",
      "FRD": "X", "FT6": "W", "FTR": "W", "FTY": "Y", "FVA": "V",
      "FZN": "K", "G  ": "G", "G25": "G", "G2L": "G", "G2S": "G",
      "G31": "G", "G32": "G", "G33": "G", "G36": "G", "G38": "G",
      "G42": "G", "G46": "G", "G47": "G", "G48": "G", "G49": "G",
      "G4P": "N", "G7M": "G", "GAO": "G", "GAU": "E", "GCK": "C",
      "GCM": "X", "GDP": "G", "GDR": "G", "GFL": "G", "GGL": "E",
      "GH3": "G", "GHG": "Q", "GHP": "G", "GL3": "G", "GLH": "Q",
      "GLJ": "E", "GLK": "E", "GLM": "X", "GLN": "Q", "GLQ": "E",
      "GLU": "E", "GLX": "Z", "GLY": "G", "GLZ": "G", "GMA": "E",
      "GMS": "G", "GMU": "U", "GN7": "G", "GND": "X", "GNE": "N",
      "GOM": "G", "GPL": "K", "GS ": "G", "GSC": "G", "GSR": "G",
      "GSS": "G", "GSU": "E", "GT9": "C", "GTP": "G", "GVL": "X",
      "GYC": "CYG", "GYS": "SYG", "H2U": "U", "H5M": "P", "HAC": "A",
      "HAR": "R", "HBN": "H", "HCS": "X", "HDP": "U", "HEU": "U",
      "HFA": "X", "HGL": "X", "HHI": "H", "HHK": "AK", "HIA": "H",
      "HIC": "H", "HIP": "H", "HIQ": "H", "HIS": "H", "HL2": "L",
      "HLU": "L", "HMR": "R", "HOL": "N", "HPC": "F", "HPE": "F",
      "HPH": "F", "HPQ": "F", "HQA": "A", "HRG": "R", "HRP": "W",
      "HS8": "H", "HS9": "H", "HSE": "S", "HSL": "S", "HSO": "H",
      "HTI": "C", "HTN": "N", "HTR": "W", "HV5": "A", "HVA": "V",
      "HY3": "P", "HYP": "P", "HZP": "P", "I  ": "I", "I2M": "I",
      "I58": "K", "I5C": "C", "IAM": "A", "IAR": "R", "IAS": "D",
      "IC ": "C", "IEL": "K", "IEY": "HYG", "IG ": "G", "IGL": "G",
      "IGU": "G", "IIC": "SHG", "IIL": "I", "ILE": "I", "ILG": "E",
      "ILX": "I", "IMC": "C", "IML": "I", "IOY": "F", "IPG": "G",
      "IPN": "N", "IRN": "N", "IT1": "K", "IU ": "U", "IYR": "Y",
      "IYT": "T", "IZO": "M", "JJJ": "C", "JJK": "C", "JJL": "C",
      "JW5": "N", "K1R": "C", "KAG": "G", "KCX": "K", "KGC": "K",
      "KNB": "A", "KOR": "M", "KPI": "K", "KST": "K", "KYQ": "K",
      "L2A": "X", "LA2": "K", "LAA": "D", "LAL": "A", "LBY": "K",
      "LC ": "C", "LCA": "A", "LCC": "N", "LCG": "G", "LCH": "N",
      "LCK": "K", "LCX": "K", "LDH": "K", "LED": "L", "LEF": "L",
      "LEH": "L", "LEI": "V", "LEM": "L", "LEN": "L", "LET": "X",
      "LEU": "L", "LEX": "L", "LG ": "G", "LGP": "G", "LHC": "X",
      "LHU": "U", "LKC": "N", "LLP": "K", "LLY": "K", "LME": "E",
      "LMF": "K", "LMQ": "Q", "LMS": "N", "LP6": "K", "LPD": "P",
      "LPG": "G", "LPL": "X", "LPS": "S", "LSO": "X", "LTA": "X",
      "LTR": "W", "LVG": "G", "LVN": "V", "LYF": "K", "LYK": "K",
      "LYM": "K", "LYN": "K", "LYR": "K", "LYS": "K", "LYX": "K",
      "LYZ": "K", "M0H": "C", "M1G": "G", "M2G": "G", "M2L": "K",
      "M2S": "M", "M30": "G", "M3L": "K", "M5M": "C", "MA ": "A",
      "MA6": "A", "MA7": "A", "MAA": "A", "MAD": "A", "MAI": "R",
      "MBQ": "Y", "MBZ": "N", "MC1": "S", "MCG": "X", "MCL": "K",
      "MCS": "C", "MCY": "C", "MD3": "C", "MD6": "G", "MDH": "X",
      "MDO": "ASG", "MDR": "N", "MEA": "F", "MED": "M", "MEG": "E",
      "MEN": "N", "MEP": "U", "MEQ": "Q", "MET": "M", "MEU": "G",
      "MF3": "X", "MFC": "GYG", "MG1": "G", "MGG": "R", "MGN": "Q",
      "MGQ": "A", "MGV": "G", "MGY": "G", "MHL": "L", "MHO": "M",
      "MHS": "H", "MIA": "A", "MIS": "S", "MK8": "L", "ML3": "K",
      "MLE": "L", "MLL": "L", "MLY": "K", "MLZ": "K", "MME": "M",
      "MMO": "R", "MMT": "T", "MND": "N", "MNL": "L", "MNU": "U",
      "MNV": "V", "MOD": "X", "MP8": "P", "MPH": "X", "MPJ": "X",
      "MPQ": "G", "MRG": "G", "MSA": "G", "MSE": "M", "MSL": "M",
      "MSO": "M", "MSP": "X", "MT2": "M", "MTR": "T", "MTU": "A",
      "MTY": "Y", "MVA": "V", "N  ": "N", "N10": "S", "N2C": "X",
      "N5I": "N", "N5M": "C", "N6G": "G", "N7P": "P", "NA8": "A",
      "NAL": "A", "NAM": "A", "NB8": "N", "NBQ": "Y", "NC1": "S",
      "NCB": "A", "NCX": "N", "NCY": "X", "NDF": "F", "NDN": "U",
      "NEM": "H", "NEP": "H", "NF2": "N", "NFA": "F", "NHL": "E",
      "NIT": "X", "NIY": "Y", "NLE": "L", "NLN": "L", "NLO": "L",
      "NLP": "L", "NLQ": "Q", "NMC": "G", "NMM": "R", "NMS": "T",
      "NMT": "T", "NNH": "R", "NP3": "N", "NPH": "C", "NPI": "A",
      "NRP": "LYG", "NRQ": "MYG", "NSK": "X", "NTY": "Y", "NVA": "V",
      "NYC": "TWG", "NYG": "NYG", "NYM": "N", "NYS": "C", "NZH": "H",
      "O12": "X", "O2C": "N", "O2G": "G", "OAD": "N", "OAS": "S",
      "OBF": "X", "OBS": "X", "OCS": "C", "OCY": "C", "ODP": "N",
      "OHI": "H", "OHS": "D", "OIC": "X", "OIP": "I", "OLE": "X",
      "OLT": "T", "OLZ": "S", "OMC": "C", "OMG": "G", "OMT": "M",
      "OMU": "U", "ONE": "U", "ONH": "A", "ONL": "X", "OPR": "R",
      "ORN": "A", "ORQ": "R", "OSE": "S", "OTB": "X", "OTH": "T",
      "OTY": "Y", "OXX": "D", "P  ": "G", "P1L": "C", "P1P": "N",
      "P2T": "T", "P2U": "U", "P2Y": "P", "P5P": "A", "PAQ": "Y",
      "PAS": "D", "PAT": "W", "PAU": "A", "PBB": "C", "PBF": "F",
      "PBT": "N", "PCA": "E", "PCC": "P", "PCE": "X", "PCS": "F",
      "PDL": "X", "PDU": "U", "PEC": "C", "PF5": "F", "PFF": "F",
      "PFX": "X", "PG1": "S", "PG7": "G", "PG9": "G", "PGL": "X",
      "PGN": "G", "PGP": "G", "PGY": "G", "PHA": "F", "PHD": "D",
      "PHE": "F", "PHI": "F", "PHL": "F", "PHM": "F", "PIA": "AYG",
      "PIV": "X", "PLE": "L", "PM3": "F", "PMT": "C", "POM": "P",
      "PPN": "F", "PPU": "A", "PPW": "G", "PQ1": "N", "PR3": "C",
      "PR5": "A", "PR9": "P", "PRN": "A", "PRO": "P", "PRS": "P",
      "PSA": "F", "PSH": "H", "PST": "T", "PSU": "U", "PSW": "C",
      "PTA": "X", "PTH": "Y", "PTM": "Y", "PTR": "Y", "PU ": "A",
      "PUY": "N", "PVH": "H", "PVL": "X", "PYA": "A", "PYO": "U",
      "PYX": "C", "PYY": "N", "QLG": "QLG", "QMM": "Q", "QPA": "C",
      "QPH": "F", "QUO": "G", "R  ": "A", "R1A": "C", "R4K": "W",
      "RC7": "HYG", "RE0": "W", "RE3": "W", "RIA": "A", "RMP": "A",
      "RON": "X", "RT ": "T", "RTP": "N", "S1H": "S", "S2C": "C",
      "S2D": "A", "S2M": "T", "S2P": "A", "S4A": "A", "S4C": "C",
      "S4G": "G", "S4U": "U", "S6G": "G", "SAC": "S", "SAH": "C",
      "SAR": "G", "SBL": "S", "SC ": "C", "SCH": "C", "SCS": "C",
      "SCY": "C", "SD2": "X", "SDG": "G", "SDP": "S", "SEB": "S",
      "SEC": "A", "SEG": "A", "SEL": "S", "SEM": "S", "SEN": "S",
      "SEP": "S", "SER": "S", "SET": "S", "SGB": "S", "SHC": "C",
      "SHP": "G", "SHR": "K", "SIB": "C", "SIC": "DC", "SLA": "P",
      "SLR": "P", "SLZ": "K", "SMC": "C", "SME": "M", "SMF": "F",
      "SMP": "A", "SMT": "T", "SNC": "C", "SNN": "N", "SOC": "C",
      "SOS": "N", "SOY": "S", "SPT": "T", "SRA": "A", "SSU": "U",
      "STY": "Y", "SUB": "X", "SUI": "DG", "SUN": "S", "SUR": "U",
      "SVA": "S", "SVV": "S", "SVW": "S", "SVX": "S", "SVY": "S",
      "SVZ": "X", "SWG": "SWG", "SYS": "C", "T  ": "T", "T11": "F",
      "T23": "T", "T2S": "T", "T2T": "N", "T31": "U", "T32": "T",
      "T36": "T", "T37": "T", "T38": "T", "T39": "T", "T3P": "T",
      "T41": "T", "T48": "T", "T49": "T", "T4S": "T", "T5O": "U",
      "T5S": "T", "T66": "X", "T6A": "A", "TA3": "T", "TA4": "X",
      "TAF": "T", "TAL": "N", "TAV": "D", "TBG": "V", "TBM": "T",
      "TC1": "C", "TCP": "T", "TCQ": "Y", "TCR": "W", "TCY": "A",
      "TDD": "L", "TDY": "T", "TFE": "T", "TFO": "A", "TFQ": "F",
      "TFT": "T", "TGP": "G", "TH6": "T", "THC": "T", "THO": "X",
      "THR": "T", "THX": "N", "THZ": "R", "TIH": "A", "TLB": "N",
      "TLC": "T", "TLN": "U", "TMB": "T", "TMD": "T", "TNB": "C",
      "TNR": "S", "TOX": "W", "TP1": "T", "TPC": "C", "TPG": "G",
      "TPH": "X", "TPL": "W", "TPO": "T", "TPQ": "Y", "TQI": "W",
      "TQQ": "W", "TRF": "W", "TRG": "K", "TRN": "W", "TRO": "W",
      "TRP": "W", "TRQ": "W", "TRW": "W", "TRX": "W", "TS ": "N",
      "TST": "X", "TT ": "N", "TTD": "T", "TTI": "U", "TTM": "T",
      "TTQ": "W", "TTS": "Y", "TY1": "Y", "TY2": "Y", "TY3": "Y",
      "TY5": "Y", "TYB": "Y", "TYI": "Y", "TYJ": "Y", "TYN": "Y",
      "TYO": "Y", "TYQ": "Y", "TYR": "Y", "TYS": "Y", "TYT": "Y",
      "TYU": "N", "TYW": "Y", "TYX": "X", "TYY": "Y", "TZB": "X",
      "TZO": "X", "U  ": "U", "U25": "U", "U2L": "U", "U2N": "U",
      "U2P": "U", "U31": "U", "U33": "U", "U34": "U", "U36": "U",
      "U37": "U", "U8U": "U", "UAR": "U", "UCL": "U", "UD5": "U",
      "UDP": "N", "UFP": "N", "UFR": "U", "UFT": "U", "UMA": "A",
      "UMP": "U", "UMS": "U", "UN1": "X", "UN2": "X", "UNK": "X",
      "UR3": "U", "URD": "U", "US1": "U", "US2": "U", "US3": "T",
      "US5": "U", "USM": "U", "VAD": "V", "VAF": "V", "VAL": "V",
      "VB1": "K", "VDL": "X", "VLL": "X", "VLM": "X", "VMS": "X",
      "VOL": "X", "WCR": "GYG", "X  ": "G", "X2W": "E", "X4A": "N",
      "X9Q": "AFG", "XAD": "A", "XAE": "N", "XAL": "A", "XAR": "N",
      "XCL": "C", "XCN": "C", "XCP": "X", "XCR": "C", "XCS": "N",
      "XCT": "C", "XCY": "C", "XGA": "N", "XGL": "G", "XGR": "G",
      "XGU": "G", "XPR": "P", "XSN": "N", "XTH": "T", "XTL": "T",
      "XTR": "T", "XTS": "G", "XTY": "N", "XUA": "A", "XUG": "G",
      "XX1": "K", "XXY": "THG", "XYG": "DYG", "Y  ": "A", "YCM": "C",
      "YG ": "G", "YOF": "Y", "YRR": "N", "YYG": "G", "Z  ": "C",
      "Z01": "A", "ZAD": "A", "ZAL": "A", "ZBC": "C", "ZBU": "U",
      "ZCL": "F", "ZCY": "C", "ZDU": "U", "ZFB": "X", "ZGU": "G",
      "ZHP": "N", "ZTH": "T", "ZU0": "T", "ZZJ": "A"}
   if res3 in resa3Toa1.keys():
      return (resa3Toa1[res3])
   else:
      return('X')
################################
#------------------------------#
def _listbsearch(vec,val,th,tol):
    ''' _listbsearch(vec,val,th,tol)
        return the list of neighbours of val in sorted vec    
    '''
    lenlist=len(vec)
    i=0
    j=lenlist
    while(i<j):
        m=abs((i+j)/2)
        if val > vec[m][0]: 
            i=m+1
        else:
            j=m
    # i points to the colsest value
    #if(vec[i][0] == val): # value located
    i=min(i,lenlist-1)
    lst=[vec[i][1]] #consider itself too
    # extend left
    j=i-1
    while (j >=0) and (abs(vec[j][0]-val)< th+tol):
        lst.insert(0,vec[j][1])
        j-=1
    j=i+1
    while (j <lenlist) and (abs(vec[j][0]-val)< th+tol):
        lst.append(vec[j][1])
        j+=1 
    return lst
#------------------------------#

def _neighborList(coords, threshold, tol=0, maincoords=None):
    ''' _neighborList(cooords, threshold, tol=0,maincoords=None)
        if maincoords==None => maincoords=cooords 
        => return for each residue in maincoords its neighbors in cooords 
                        
    '''
    if not maincoords:
       maincoords=coords
    protLen=len(maincoords[0])
    reshash={}
    for (val,i) in maincoords[0]:
        reshash[i]=[]
    for c in range(3):
        coords[c].sort()
        for (val,i) in maincoords[c]:
            reshash[i].append(_listbsearch(coords[c],val,threshold,tol)) 
    for i in reshash.keys(): # compute intersection
        x,y,z=reshash[i]
        reshash[i]=[]
        for e in x:
            if e in y and e in z:
               reshash[i].append(e) 
    return reshash    

#---------------------------------

##
# function that computes the outer product of 2 vectors
# and returns the vectorial product vector
def _vetP(a,b):
    ''' _vetP(a,b)
        returns the vectorial products of a and b
	a and b must be a 3-dim arrays
    '''
    assert(len(a)==len(b)==3)
    a=numpy.array(a)
    b=numpy.array(b)
    c=numpy.array([0.0]*3)
    c[0]=a[1]*b[2]-a[2]*b[1]
    c[1]=a[2]*b[0]-a[0]*b[2]
    c[2]=a[0]*b[1]-a[1]*b[0]
    return(c)

###########################
# computes the torsion angles
# around v2 and v3
# given the 4 points 
#
#  v1--v2
#        \
#         v3--v4
#
def _torsion(v1,v2,v3,v4):
    ''' _torsion(v1,v2,v3,v4)
        computes the torsion angles
	around v2 and v3
	given the 4 points
        v1--v2--v3--v4
    '''
    e=numpy.array([[0.0]*3]*3) # the versors
    v1=numpy.array(v1)
    v2=numpy.array(v2)
    v3=numpy.array(v3)
    v4=numpy.array(v4)
    e[0]=v2-v1  # versor e12
    e[1]=v3-v2  # versor e23
    e[2]=v4-v3  # versor e34
    for i in range(3): # normalize the 3 versors
        e[i]=e[i]/numpy.sqrt(numpy.dot(e[i],e[i])) 
    p1=_vetP(e[0],e[1]) # vectorial products
    p2=_vetP(e[1],e[2])
    seno=numpy.dot(e[0],p2)
    coseno=numpy.dot(p1,p2)
    return(numpy.arctan2(seno,coseno))
   
     

# Class Residue
#
'''
class RESIDUE is defined as
   Name => "residue name" (es. GLY ALA ...)
   AtomNumber => number of atoms composing the residue (eg 2)
   Coord => a list of alist o 3d coordinates for each atoms 
             (eg [[1.0, 2.0, 3.0],[2.33, 2.0, 3.0]])
   AtomNames => the list of the atom names (eg ['CA','CB'])
 
   There are 2 fnctions defined
   updateResidue(Coord,AtomNames,ResName=None) updating the
                 residue contents
   
   residueDistance(Residue,selfList=[],resList=[]) computing the
                 distance between another residue using all(default) or
		 lists of atoms for each residue


'''
class RESIDUE:
   ''' define the residue object
       RESIDUE is defined as:
       Name (string)
       NameOne (name one letter code)
       AtomNumber (number of atom composing the residues)
       Coord (array rows=AtomNumber coloum=the three coordinates)
   '''
   def __init__(self,Name,Coord,AtomNames,Others,Chain='',AltLoc=[]):
      assert(len(Coord)==len(AtomNames))
      self.Name=Name
      self.NameOne=_getResidueA1(Name)
      self.AtomNumber=len(AtomNames)
      self.Coord=Coord
      self.AtomNames=AtomNames
      self.Others=Others
      self.Chain=Chain
      self.AltLoc=AltLoc
   
   def updateResidue(self,Coord,AtomNames, Others, ResName=None, Chain='', AltLoc=[]):
      ''' update the residue content using a list of 
          coordinates and a list of the corresponding
	  atom names
      '''
      assert(len(Coord)==len(AtomNames))
      if(ResName!=None): 
          self.Name=ResName
	  self.NameOne=_getResidueA1(ResName)
      self.AtomNumber += len(AtomNames)
      self.Chain=Chain
      self.AltLoc=AltLoc
      for i in range(len(Coord)):
          self.Coord.append(Coord[i])
	  self.AtomNames.append(AtomNames[i])
	  self.Others.append(Others[i])

   
   def residueDistance(self,Residue,selfList=[],resList=[]):
       ''' compute the euclidean distance
           from residue and residue using minimal distance
	   among the atom in the residues
       '''
       if (selfList == []) :
           selfList=self.AtomNames # set the current list to all possible

       if (resList == []) :
           resList=Residue.AtomNames # set the current list to all possible
      
       minDistance=999999.9
       currentDistance=999999.9
       for i in selfList:
          if(i in self.AtomNames):     
             vi=numpy.array(self.Coord[self.AtomNames.index(i)],numpy.float)
	     for j in resList:
                if(j in Residue.AtomNames):
                   vj=numpy.array(Residue.Coord[Residue.AtomNames.index(j)],numpy.float)
                   numpy.add(vi,-vj,vj)
                   currentDistance= numpy.sqrt(numpy.dot(vj,vj))
	        if (currentDistance < minDistance):
		   minDistance=currentDistance
       return minDistance

#---- end class residue ----


# class PDBChain 
# this class concerns the pdb coordinate handling
class PDBChain:
   '''
   class PDBChain is defined as
       protLen => protein length
       __pdbPositions => label for eache residues in the
                      PDBChain that contains the 'number'
		      assigned to the residues in the pdb file
		      (it could be an alphanumeric es '123A')
       __Residues => list of residues in the chain
		 
      Constructor
      PDBChain(Positions,Residues)
         Positions => ['1','2',....,'123','123A',...]
	 Residues  => [Residue1, Residue2,.., Residue123,...]

   '''

   def __init__(self, Positions,Residues) :
       '''
          __init__(Positions,Residues)
          Positions => an array containing the lists  
              example of Positions => ['1','2',....,'123','123A',...]
          Residues => array of residues
	      example of Residues  => [Residue1, Residue2,.., Residue123,...]
       '''
       assert(len(Residues)==len(Positions))
       self.protLen=len(Positions)
       self.__pdbPositions=Positions
       self.__Residues=Residues
	  
   def updateChain(self,Position,Residue):
      ''' update the residue content using a 
          residue and a the corresponding
	  positions
	  if obj is a PDBChain object
	  obj.(Pos,Res)
	  add to the existing obj the residue Res with position Pos
      '''
      self.protLen=self.protLen+1
      self.__pdbPositions.append(Position)
      self.__Residues.append(Residue)
   
   def checkContinuity(self,CutOff=1.54,Atoms=['C','N'],minCutOff=0.5):
       ''' checkContinuity(CutOff=1.4,Atoms=['C','O'],minCutOff=0.5)
           tests if the bond distances between each
	   hypotetically bonded residues is below a given
	   threshold
       '''
       Errors=[]
       for i in range(self.protLen-1):
           d=self.distance(i,i+1,Atoms,Atoms)
	   if(d > CutOff or d < minCutOff):
              Errors.append([d,self.__pdbPositions[i],self.__pdbPositions[i+1]])
       return Errors
   
   def checkPDB(self):
       ''' check if there are holes or strange residues in the
           pdb enumerations, and print the lists
       '''
       Holes=[]
       Repited=[]
       Stranges=[]
       Errors={}
       Text='ABCDEFGHIJKLMNOPQRSTUVWXYZ'
       for i in range(self.protLen-1):
	  if(self.__pdbPositions[i][-1].upper() in Text):     
	     pi=string.atoi(self.__pdbPositions[i][0:-1])   
	  else: 
             pi=string.atoi(self.__pdbPositions[i])
	  if(self.__pdbPositions[i+1][-1].upper() in Text):     
             Stranges.append(self.__pdbPositions[i+1])		   
	     pii=string.atoi(self.__pdbPositions[i+1][0:-1])   
	  else: 
             pii=string.atoi(self.__pdbPositions[i+1])
          if(pi  == pii):
	     Repited.append([self.__pdbPositions[i],self.__pdbPositions[i+1]])	  
	  elif(pi +1 != pii):
             Holes.append([self.__pdbPositions[i],self.__pdbPositions[i+1]])
	  if(Holes !=[]):
             Errors['Holes']=Holes
	  if(Repited !=[]):
	     Errors['Repited']=Repited
	  if(Stranges!=[]):   
	     Errors['Stranges']=Stranges
       return Errors
   
   def getPDBNumbers(self):
       ''' return the list of the pdb numbers of each residue'''
       return self.__pdbPositions

   def getPDBNumber(self,i):
       ''' return the pdb number of the i-th residue'''
       return self.__pdbPositions[i]
   
   def getResidue(self,i):
       ''' return the i-th residue (starting from 0  '''
       assert(0<=i<self.protLen)
       return   self.__Residues[i]

   def TranslateChain(self,vector):
       for i in range(self.protLen):
          res=self.getResidue(i)
          for at in range(len(res.Coord)):
             for j in range(3):
                res.Coord[at][j]=res.Coord[at][j]-vector[j]
       return

   def RotateChain(self,matrix):
       for i in range(self.protLen):
          res=self.getResidue(i)
          for at in range(len(res.Coord)):
             ncoord=[0.0,0.0,0.0]
             for j in range(3):
                ncoord[j]=res.Coord[at][0]*matrix[j][0]+res.Coord[at][1]*matrix[j][1]+res.Coord[at][2]*matrix[j][2]
             res.Coord[at]=ncoord
       return
   
   def getPhiPsi(self,slice=[],inRadiants='Y'):
       ''' getPhiPsi(slice=[])
           return the list of the Phi and Psi
	   torsion angles for each residue in slice
	   if slice=[] returns all the residues in memory
	   if inRadiants != 'Y' returns degree
	   returns a list of len(slice)X2 ([phi,psi])
       '''
       if(slice==[]):
           slice=range(self.protLen)
       PhiPsi=[]
       for i in slice:
	   # Phi = C(i-1)-N(i)-CA(i)-C(i)
           if(i == 0):
	       phi=360.0 # default value
           else:
	       v1=self.getCoordResAtom(i-1,'C')
               v2=self.getCoordResAtom(i,'N')
	       v3=self.getCoordResAtom(i,'CA')
               v4=self.getCoordResAtom(i,'C')
	       phi=_torsion(v1,v2,v3,v4)
           # Psi = N(i)-CA(i)-C(i)-N(i+1)
           if(i == self.protLen-1):
               psi=360.0 # default value
	   else:
	       v1=self.getCoordResAtom(i,'N')
               v2=self.getCoordResAtom(i,'CA')
	       v3=self.getCoordResAtom(i,'C')
               v4=self.getCoordResAtom(i+1,'N')
	       psi=_torsion(v1,v2,v3,v4)
           PhiPsi.append((phi,psi)) 
       PhiPsi=numpy.array(PhiPsi,numpy.float)
       if(inRadiants != 'Y'):
           conversionFactor=180.0/numpy.pi
           PhiPsi=PhiPsi*conversionFactor
       return(PhiPsi.tolist())


   def getCoordResAtom(self,i,Atom):
       ''' getCoordResAtom(i,Atom) 
           returns the coordinates of the 'Atom' in the
	   residue 'i'
       '''
       assert(0<=i<self.protLen)
       if(Atom in self.__Residues[i].AtomNames):
	   coordPos=self.__Residues[i].AtomNames.index(Atom)    
           return   self.__Residues[i].Coord[coordPos]
       else:
           return None

   def getCoordResAtomByNumber(self,lab,Atom):
       ''' getCoordResAtomByNumber(lab,Atom) 
           returns the coordinates of the 'Atom' in the
           residue 'lab'
       '''
       i=self.__pdbPositions.index(lab)
       try:
          coordPos=self.__Residues[i].AtomNames.index(Atom)
          return self.__Residues[i].Coord[coordPos]
       except:
           return None
   
   def getResidueByNumber(self,i):
       ''' return the residue whose pdb number is i  '''
       if(i in self.__pdbPositions): 
           return   self.__Residues[self.__pdbPositions.index(i)]
       else:
           return None

   def distancePDBNumber(self,i,j,i_list=[],j_list=[]):
       ''' compute the eucliden distance between residues PDB positions i and j
           using i-list and j-list the two lists of atom type to consider
       '''
       Res_i=self.getResidueByNumber(i)
       Res_j=self.getResidueByNumber(j)
       return (Res_i.residueDistance(Res_j,i_list,j_list))
   
   def distance(self,i,j,i_list=[],j_list=[]):
       ''' compute the eucliden distance between residues i and j
           using i-list and j-list the two lists of atom type to consider
       '''
       Res_i=self.getResidue(i)
       Res_j=self.getResidue(j)
       return (Res_i.residueDistance(Res_j,i_list,j_list))

   def getChain(self):
       ''' getChain(self) -> Chain Label
       '''
       return self.getResidue(0).Chain

   def getSequence(self, slice=[]):
       ''' getSequence(self, slice=[]) -> string of protein one-lettercode sequence
           Starting from the PDBChain, retrive a string that is the chain sequence.
           It could also retrive a sequence slice. If slice is not provided (default 
           slice = []), the entire chain will be considered. 		
       '''
       if(slice == []): # if not provided use the whole Chain
		slice=range(self.protLen)
       seq=''
       for i in slice:
		if i in range(len(self.__Residues)):
			seq+= self.__Residues[i].NameOne
       return seq

   def getSequenceA3(self, slice=[]):
       ''' getSequenceA3(self, slice=[]) -> string of protein three-lettercode sequence
           Starting from the PDBChain, retrive a string that is the chain sequence.
           It could also retrive a sequence slice. If slice is not provided (default 
           slice = []), the entire chain will be considered. 		
       '''
       if(slice == []): # if not provided use the whole Chain
		slice=range(self.protLen)
       seq=[]
       for i in slice:
		if i in range(len(self.__Residues)):
			seq.append(self.__Residues[i].Name)
       return seq

   def neighborList(self, threshold, tol=0, pdblist=[]):
       ''' neighborList(pdbChain, threshold, tol=0, pdblist=[])
           pdbChain is an instance of PDBChain
           => return for each residue e list of neighbors 
       '''
       if not pdblist:
           pdblist=range(self.protLen)
       coords=[[],[],[]] # x=0, y=1, z=2
       for i in pdblist: 
           coord_i=self.getCoordResAtom(i,'CA')
           for c in range(3):
               coords[c].append((coord_i[c],i)) 
       return _neighborList(coords, threshold, tol)

   def neighborPDBList(self, threshold, tol=0, pdblist=[]):
       ''' neighborPDBList(pdbChain, threshold, tol=0, pdblist=[])
           pdbChain is an instance of PDBChain
	   pdblist are pdb labels
           => return for each residue e list of neighbors 
       '''
       if not pdblist:
           pdblist=self.__pdbPositions
       coords=[[],[],[]] # x=0, y=1, z=2
       for i in pdblist: 
           coord_i=self.getCoordResAtomByNumber(i,'CA')
           for c in range(3):
               coords[c].append((coord_i[c],i)) 
       return _neighborList(coords, threshold, tol)

   def neighbors(self,th_dist,num_nn,pdblist=[]):
       ''' neighbors(self,th_dist,num_nn,pdblist=[])
           computes the neighbors of each residue in pdblist
	   or in case pdblist=[] of the whole pdb
           The number of neighbors is num_nn
           The distance threshold cut-off is th_dist
           => returns the list [[neighbors of i]for each i in pdblist]
       '''
       if not pdblist:
           pdblist=self.__pdbPositions
       dist=[[100000 for j in range(num_nn)] for i in range(len(pdblist))]
       lst=[["--" for j in range(num_nn)] for i in range(len(pdblist))]
       neigh=self.neighborPDBList(th_dist,tol=0,pdblist=pdblist)
       i=0
       for lab1 in pdblist:
           for lab2 in neigh[lab1]: # 
               d=self.distancePDBNumber(lab1,lab2)
               if(d < th_dist):
                   k=0
                   while (k < num_nn) and (d > dist[i][k]):
                       k=k+1
                   if k != num_nn:
                       dist[i].insert(k,d)
                       dist[i].pop()
                       lst[i].insert(k,lab2)
                       lst[i].pop()
           i=i+1
       return(lst)


   def writePDB(self,fileName=None,slice=[],oldNumbers=None,chainName=' ',header='Y'):
       ''' writePDB(self,fileName=None,slice=[],oldNumbers=None,chainName=' ')
           write a pdb file coordinates using the existing
	   protein in memory. 
	   writePDB writes on <stdout> if a fileName is NOT provided (default)
	   or on the given fileName.
	   writePDB can also write only a slice (subset identified by the list
	   of the number of the residue we want (default=[] => all 
	   [0,1,..protLen])
	   oldNumbers='Y' => the old label numbers for each residue are used
	   chainName is set by default to ' ' 
       '''
       resCount=1  # residue conter
       resLable="" # residue label "number"
       atomCount=1 # atom counter
       if(slice == []): # if not provided use the whole protein
	    slice=range(self.protLen)
       strPrn="" # string to print
       if header=='Y' or header=='y':
          strPrn+= "%s %s\n" % ("HEADER ",fileName) # print header
          strPrn+= "%s %s\n" % ("COMPND ",fileName) # print header
          # print protein SEQRES fields
          k=1 # line counter for SEQRES field
          effProtLen=len(slice)  # effective protein length to save
          # print SEQRES fields
          for i in range(len(slice)):
             if(i%13 == 0): # add newline if the residue is %13
                if(i !=0) : strPrn+="\n" # and it is not == 0
                strPrn+= "%s%4d %1.1s%5d  " % ("SEQRES",k,chainName,effProtLen)
                k+=1
             strPrn+= "%3s " % self.getResidue(slice[i]).Name	
          strPrn+= "\n"
       # print ATOM fields
       for i in slice: # for each residue in the list
           currentResidue=self.getResidue(i)
           # if it is set to oldNumbers use them
           if(oldNumbers !=  None):
	       resLabel=self.__pdbPositions[i]
	   else :
               resLabel=str(resCount)
	   k=0
           if resLabel[-1] in '1234567890': resLabel=resLabel+' '
           for j in range(len(currentResidue.AtomNames)): # for each atom
              a=currentResidue.AtomNames[j]
              altloc=currentResidue.AltLoc[j]
              #strPrn+="ATOM%7d  %-4s%-4s%1.1s%4s    " % 
              strPrn+="ATOM  %5d  %-3s%1s%3s %1.1s%5s   " % \
	            (atomCount,a,altloc,currentResidue.Name,chainName,resLabel)
	      strPrn+="%8.3f%8.3f%8.3f" % tuple(currentResidue.Coord[k])
	      strPrn+="%6s%6s" % tuple(currentResidue.Others[k])+14*' '+'\n'
	      #strPrn+='  0.00  0.00'+14*' '+'\n'
              k+=1
	      atomCount+=1
	   resCount+=1
       strPrn+='TER\n'
       if(fileName == None):
           print strPrn
       else:
           try: 
	      f=open(fileName, 'w')
	   except:
              print "Can't open file ",fileName	   
              sys.exit()
           f.write(strPrn)
	   f.close()
	   
	   
#--- end class PDBChain

# class PDBProtein 
class PDBProtein:
   '''
   class PDBProtein
       seqNum => Number of sequences
       __chainLabels => label for each PDBchain in the PDBChain that contains
			the 'Label' (if exist), assigned to the  Chain in the 
		        pdb file (it could be a letter or a number,  for over 
			24 sequences proteins)
       __Chains => list of PDBChains instances in the protein
		 
      Constructor
      PDBProtein(Chains)
	 Chains =>  [PDBChains1,PDBChains2,...,PDBChainsN]
   '''

   def __init__(self, Chains) :
	'''
	  __init__(Chains)
	  Chains => array of Chains
	      example of Chains => [PDBChains1,PDBChains2,...,PDBChainsN]
	'''
	self.__Chains=Chains
	self.__chainLabels=self.getChains()
	self.seqNum=len(self.__chainLabels)
	self.resNum=0
	for PCh in self.__Chains:
		self.resNum+=PCh.protLen

   def updateProtein(self, Chain):
	''' update the Chains content using a PDbChain (the corresponding Label
	    is retrived by the PDB.Chain.getChain() function. 
	    if obj is a PDBProtein object obj.(Lab,Ch) add to the existing obj 
	    the PDBChain Chain with label Label
	'''
	self.seqNum=self.seqNum+1
	self.__chainLabels.append(Chain.getChain())
	self.__Chains.append(Chain)
	self.resNum=self.resNum+Chain.protLen

   def getPDBChain(self,i):
	'''Return i-th PDBChain instance'''
	if i in range(len(self.__Chains)):
		return self.__Chains[i]
	else: return None
	
   def getPDBChainByChain(self,ch):
	'''Return the PDBChain instance whose Chain label is i '''
	if ch in self.__chainLabels:
		i = self.__chainLabels.index(ch)
		return self.__Chains[i]
	else: return None
	
   def getChains(self):
	''' getChains() -> Return a list of the Protein Chains.
	'''
	chains=[]
	for ch in self.__Chains:
		chains.append(ch.getChain())
        return chains

   def getSequences(self):
	''' getSequences() -> Return a list of the Protein Sequences.'''
	seqs=[]
	for ch in self.__Chains:
		seqs.append(ch.getSequence())
	return seqs


   def getPDBNumbers(self, Chain=None):
       '''getPDBNumbers(self, Chain=None) -> List of PDB residue Labels
	  return the list of the pdb numbers of each residue.
	  If Chain is 'Y',  it  return  the  PDB Number  and its chain 
	  (example ['1 A'] for the residue 1 of the chain A)
       '''
       PDBNum = []
       for PCh in self.__Chains:
		PDBNumCh=[]
		PDBNumCh.extend(PCh.getPDBNumbers())
		if Chain != None:
			Ch = PCh.getResidue(0).Chain
			for i in range(len(PDBNumCh)):
				PDBNumCh[i]+=" "+Ch
				PDBNumCh[i]=PDBNumCh[i].strip()
		PDBNum.extend(PDBNumCh)
       return PDBNum
   
   def getResidueByNumber(self, PDBNum, Chain=None):
	'''getResidueByNumber(self, PDBNum,Chain=None)-> PDB residue instance
	   Return the Residue instance whose pdb number is PDBNum.
	   For multiple chains objects, there may be not unique identifiers
	   and so a Chain retrivial is possible. When Chain is not specified
	   the function returns the first Residue whose pdb Number is PDBNum.
	'''
	if Chain:
		if PDBNum in self.getPDBChainByChain(Chain).getPDBNumbers():
		      return self.getPDBChainByChain(Chain).getResidueByNumber(PDBNum)
		else: return None
	else:
		for PCh in self.__Chains:
			if PDBNum in PCh.getPDBNumbers():
				return PCh.getResidueByNumber(PDBNum)
			else:   return None

   def distancePDBNumber(self,i,j,chain_i,chain_j,i_list=[],j_list=[]):
       ''' compute the eucliden distance between residues PDB positions i and j
           of chain_i and chain_j
           using i-list and j-list the two lists of atom type to consider
       '''
       Res_i=self.getPDBChainByChain(chain_i).getResidueByNumber(i)
       if not Res_i:
           print "ERROR None",chain_i,i
       Res_j=self.getPDBChainByChain(chain_j).getResidueByNumber(j)
       if not Res_j:
           print "ERROR None",chain_j,j
       return (Res_i.residueDistance(Res_j,i_list,j_list))
   
   def distance(self,i,j,chain_i,chain_j,i_list=[],j_list=[]):
       ''' compute the eucliden distance between residues i and j
           of chain_i and chain_j
           using i-list and j-list the two lists of atom type to consider
       '''
       Res_i=self.getPDBChainByChain(chain_i).getResidue(i)
       Res_j=self.getPDBChainByChain(chain_j).getResidue(j)
       return (Res_i.residueDistance(Res_j,i_list,j_list))

   def patches(self,chain1,chain2,th_dist,list1=[],list2=[]):
       ''' patches(self,chain1,chain2,th_dist,list1=[],list2=[])
           returns 
               the dictionary k:1 or k:0 for k in reslist1 and
               the dictionary k:1 or k:0 for k in reslist2
               which are found in contact (1) on not (0) with some residues in the
               opposite chain
       '''
       if not list1:
           list1=self.getPDBChainByChain(chain1).getPDBNumbers()
       if not list2:
           list2=self.getPDBChainByChain(chain2).getPDBNumbers()

       coords1=[[],[],[]] # x=0, y=1, z=2
       coords2=[[],[],[]] # x=0, y=1, z=2
       for i in list1: 
           coord_i=self.getPDBChainByChain(chain1).getCoordResAtomByNumber(i,'CA')
           for c in range(3):
               coords1[c].append((coord_i[c],i)) 
       for i in list2: 
           coord_i=self.getPDBChainByChain(chain2).getCoordResAtomByNumber(i,'CA')
           for c in range(3):
               coords2[c].append((coord_i[c],i)) 
       neighbours=_neighborList(coords2, th_dist, tol=0, maincoords=coords1)
       patch1={}
       patch2={}
       for name1 in list1:
           patch1[name1]=0
       for name2 in list2:
           patch2[name2]=0
       for name1 in list1:
           for name2 in neighbours[name1]:
#           for name2 in list2:
               d=self.distancePDBNumber(name1,name2,chain1,chain2)
               if d < th_dist:              
                   patch1[name1]=1
                   patch2[name2]=1
       return(patch1,patch2) 

   def writePDBProtein(self,fileName=None,slice=[],oldNumbers=None,header='Y'):
      	''' writePDBProtein(self,fileName=None,slice=[],oldNumbers=None)
            write a pdb file coordinates using the existing
	    protein in memory. 
	    writePDB writes on <stdout> if a fileName is NOT provided (default)
	    or on the given fileName.
	    writePDB can also write only a slice (subset identified by the list
	    of the Chain Labels we want to write (default=[] => all Chains)
	    oldNumbers='Y' => the old label numbers for each residue are used
       	'''
	resCount=1  # residue counter
       	resLable="" # residue label "number"
       	atomCount=1 # atom counter
       	if not(slice): # if not provided use the whole protein
		slice=self.__chainLabels
      	strPrn="" # string to print
        if header=='Y' or header=='y':
           strPrn+= "%s %s\n" % ("HEADER ",fileName) # print header
           strPrn+= "%s %s\n" % ("COMPND ",fileName) # print header
	Seqres=""
	Atom=""
	# control chain lengths
	effProtLen=[]
	seqs=self.getSequences()
	for seq in seqs:
		effProtLen.append(len(seq))
	for j in range(self.seqNum):
	    ch = self.__chainLabels[j]
	    if ch in slice:
	        # print protein SEQRES fields
	        k=1 # line counter for SEQRES fields
	        # print SEQRES fields
	        for i in range(effProtLen[j]):
		   currentResidue=self.__Chains[j].getResidue(i)
	           if(i%13 == 0): # add newline if the residue is %13
		       if(i != 0) : Seqres+="\n" # and it is not == 0
	               Seqres+= "%s%4d %1.1s%5d  " % ("SEQRES",k,ch,effProtLen[j])
	               k+=1
	           Seqres+= "%3s " % currentResidue.Name	
	           # print ATOM fields
	           # if it is set to oldNumbers use them
	           if(oldNumbers):
		       resLabel=self.__Chains[j].getPDBNumbers()[i]
		   else:
	               resLabel=str(resCount)
                   if resLabel[-1] in '1234567890': resLabel=resLabel+' '
		   m=0
	           for l in range(len(currentResidue.AtomNames)): # for each atom
                      a=currentResidue.AtomNames[l]
                      altloc=currentResidue.AltLoc[l]
                      Atom+="ATOM  %5d  %-3s%1s%3s %1.1s%5s   " % \
                          (atomCount,a,altloc,currentResidue.Name,ch,resLabel)
	              #Atom+="ATOM%7d  %-4s%-4s%1.1s%4s    " % \
		      #      (atomCount,a,currentResidue.Name,ch,resLabel)
		      Atom+="%8.3f%8.3f%8.3f" % tuple(currentResidue.Coord[m])
		      Atom+="%6s%6s" % tuple(currentResidue.Others[m])+14*' '+'\n'
	              m+=1
		      atomCount+=1
		   resCount+=1 
		Seqres+="\n"
                Atom+='TER\n'
	if header=='Y' or header=='y':
           strPrn+=Seqres+Atom
        else:
           strPrn+=Atom
       	if(fileName == None):
           	print strPrn
       	else:
           	try: 
	      		f=open(fileName, 'w')
	   	except:
              		print "Can't open file ",fileName	   
              		sys.exit()
           	f.write(strPrn)
	   	f.close()
	   



#------------
# read a pdb file and build a PDBProtein object
# 
# taking into account the chain and a list of atoms
# to consider
# return only their  coordinates and positions
def readPDBProtein(file,chain='_',AtomType=['C','N','O','CA','CB'],HetAtm='N'):
   ''' read a pdb file and return
       the postions for the chain 
       for each atom type in the list AtomType
       If the chain == '_' means don't care 
       return a list PDBChain object
   '''
   
#       while not eof
#         se trovi Atom
#	    se trovi chain e atomtype compatibile
#	       se il residuo e' nuovo
#	         appendi quello vecchio al pdb
#		 inizializza il nuovo
#	       altrimenti
#	         appendi atomo e coordinate al vecchio
   
   lines = open(file).readlines()
   ChainRead=0      #  this means that the chain has not been read
   atom_coord=[]    #  positions
   atom_type=[]     #  type of atom
   atom_altloc=[]   #  alternative location
   atom_others=[]   #  occupancy and T-Factor
   res_name="NONE"
   old_res_name="NONE"   
   oldNum=-99999 
   old_chain='!' 
   rchain='!' 
   keepchain='!'
   pdbProt=PDBProtein([])     # init pdbProtein
   pdb=PDBChain([],[])        # init pdbChain
   already_done = []
   ter=0
   if HetAtm!='Y' and HetAtm!='y':
	atype=['ATOM','TER']
   else:
        atype=['ATOM','HETATM','TER']
   for line in lines:
      record=string.strip(line[:6])
      if ( record in atype ):
	if record=='TER':
		ter =1
	else:
		rchainold=rchain
        	rchain=line[21]  # current chain
	if rchainold!=rchain and rchain!=' ': ter = 0
	if keepchain != '!': rchainold = keepchain # for already_done chains problems
	if ( (rchain in chain) or (rchain in chain[0]) or (chain[0] == '_') ) and ter==0:
	  ChainRead=1
          atom_name=string.strip(line[13:16]) # current atom name
	  if(atom_name in AtomType):
            res_name=string.strip(line[17:20]) # residue
            altloc=line[16]       # alternative location
	    resNum=string.strip(line[22:27])  # current residue number
            x=string.atof(line[30:38]) # current x
            y=string.atof(line[38:46]) # current y 
            z=string.atof(line[46:54]) # current z
	    try:
            	occ=string.strip(line[54:60])
	    except:
		occ='  0.00'
	    try:
	    	tfac=string.strip(line[60:66])
	    except:
		tfac='  0.00'
	    if (resNum+" "+rchain+" "+atom_name).strip() not in already_done:
	      if(oldNum != (resNum+rchain).strip()):  # residue read != old number
                 if(oldNum != -99999): # this is for the first residue and first atom
                    residue.updateResidue(atom_coord,atom_type,atom_others,None,rchainold,atom_altloc)
	            pdb.updateChain(resNumold,residue)
                    atom_coord=[]    #  positions
                    atom_type=[]     #  type of atom
		    atom_altloc=[]   #  alternative location
                    atom_others=[]   #  occupancy and T-Factor
		 if (old_chain != rchain) and (old_chain != '!'): # take in account also the first chain
	    	    pdbProt.updateProtein(pdb)
		    ProtChains=pdbProt.getChains()
		    if rchain in ProtChains: # reload the correct chain to upgrade
		      pdb=pdbProt._PDBProtein__Chains[ProtChains.index(rchain)]
		      pdblen=pdb.protLen	
		      pdbProt._PDBProtein__Chains.pop(ProtChains.index(rchain))
		      pdbProt._PDBProtein__chainLabels.pop(ProtChains.index(rchain))
		      pdbProt.seqNum-=1
		      pdbProt.resNum-=pdblen
		    else:
	    	      pdb=PDBChain([],[]) # reinit the new chain
		 old_chain = rchain
	         residue=RESIDUE(res_name,[],[],[]) # init the residue
                 atom_coord.append([x,y,z])
                 atom_type.append(atom_name)
                 atom_altloc.append(altloc)   
		 atom_others.append([occ,tfac])
              else:   
	         atom_coord.append([x,y,z])
	         atom_type.append(atom_name)
                 atom_altloc.append(altloc)
		 atom_others.append([occ,tfac])
	      oldNum=(resNum+rchain).strip()
	      resNumold=resNum
	      already_done.append((resNum+" "+rchain+" "+atom_name).strip())
	      keepchain='!'
	    else:
	      keepchain = rchainold
   if keepchain != '!': rchain = keepchain # for already_done chains problems
   if(atom_type !=[]) and (ChainRead == 1):
      residue.updateResidue(atom_coord,atom_type,atom_others,None,rchain,atom_altloc) # for the last atom type(important: use rchain and not rchainold)
      pdb.updateChain(resNumold,residue) # for the last residue
      pdbProt.updateProtein(pdb) # for the last chain
   return pdbProt 
###################################################


#------------
# read a pdb file
# 
# taking into account the chain and a list of atoms
# to consider
# return only their  coordinates and positions
def readPDB(file,chain='_',AtomType=['C','N','O','CA','CB'],HetAtm='Y'):
   ''' read a pdb file and return
       the postions for the chain 
       for each atom type in the list AtomType
       If the chain == '_' means don't care 
       return a list PDBChain object
   '''
   
#       while not eof
#         se trovi Atom
#	    se trovi chain e atomtype compatibile
#	       se il residuo e' nuovo
#	         appendi quello vecchio al pdb
#		 inizializza il nuovo
#	       altrimenti
#	         appendi atomo e coordinate al vecchio
   
   lines = open(file).readlines()
   ChainRead=0      #  this means that the chain it has not been read
   atom_coord=[]    #  positions
   atom_type=[]     #  type of atom
   atom_others=[]   #  occupancy and T-Factor
   atom_altloc=[]   #  alternative location
   res_name="NONE"
   old_res_name="NONE"   
   oldNum=-99999      
   pdb=PDBChain([],[])        # init pdbchain
   if HetAtm!='Y' and HetAtm!='y':
	atype=['ATOM','TER']
   else:
        atype=['ATOM','HETATM','TER']
   for line in lines:
      record=string.strip(line[:6])
      if(ChainRead and record =='TER'):
         break
      if (record in atype):
         rchain=line[21]  # current chain
         atom_name=string.strip(line[13:16]) # current atom name
	 if((rchain==chain or chain =='_')and(atom_name in AtomType)):
	    ChainRead=1 # the chain it has beeen read
            altloc=line[16]       # alternative location
            res_name=string.strip(line[17:20]) # residue
	    resNum=string.strip(line[22:27])  # current residue number
            x=string.atof(line[30:38]) # current x
            y=string.atof(line[38:46]) # current y 
            z=string.atof(line[46:54]) # current z
	    try:
	    	occ=string.strip(line[54:60])
	    except:
		occ='  0.00'
	    try:
		tfac=string.strip(line[60:66])
	    except:
		tfac='  0.00'
	    if(oldNum != resNum):  # residue read != old number
               if(oldNum != -99999): # this is the first residue
                  residue.updateResidue(atom_coord,atom_type,atom_others,None,rchain,atom_altloc)   
	          pdb.updateChain(oldNum,residue)
                  atom_coord=[]    #  positions
                  atom_type=[]     #  type of atom
		  atom_others=[]   #  occupancy and T-Factor
                  atom_altloc=[]   #  alternative location
	       residue=RESIDUE(res_name,[],[],[]) # init the residue
               atom_coord.append([x,y,z])
               atom_type.append(atom_name)
	       atom_others.append([occ,tfac])
               atom_altloc.append(altloc)
            else:   
	       atom_coord.append([x,y,z])
	       atom_type.append(atom_name)
	       atom_others.append([occ,tfac])
               atom_altloc.append(altloc)
	    oldNum=resNum
   if(atom_type !=[]):
      residue.updateResidue(atom_coord,atom_type,atom_others,None,rchain,atom_altloc)
      pdb.updateChain(oldNum,residue)
   return pdb 
###################################################




if (__name__ == '__main__'):
   if(len(sys.argv) == 1):
      print "syntax prog inputfile [chain]"
#       z=readPDBProtein('/home/emidio/skype/p3.pdb','_',TotAtomType)
#       z.writePDBProtein('p4.pdb',[],'Y''Y')
#       z=readPDB('/home/emidio/skype/pdb1aac.ent','_',TotAtomType)
#       z.writePDB('p1.pdb',[],'Y','Z','N')
#       z.writePDB('p2.pdb',[],'Y','Y','N')
      sys.exit()
   else:
      try:
         if(len(sys.argv)>=3):
            Protein = readPDB(sys.argv[1],sys.argv[2],['CA','O','C','N'])
         else:	
            Protein = readPDBCA(sys.argv[1])
      except IOError, (errno, strerror):
         print "I/O error(%s): %s" % (errno, strerror)

      print Protein.protLen
      print Protein.checkContinuity()


