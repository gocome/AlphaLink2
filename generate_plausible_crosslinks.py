from Bio.PDB import PDBParser
import numpy as np


pdb_file = "./test_data/9DS3.pdb"
parser = PDBParser()
structure = parser.get_structure("1CDL", pdb_file)

# Map 3-letter to 1-letter amino acid codes
aa_map = {
    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E',
    'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N',
    'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S',
    'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
}

# Get Cα or Nζ (Lys sidechain) coordinates
def get_lys_nz(residue):
    if residue.get_resname() == "LYS" and "NZ" in residue:
        return residue["NZ"].get_coord()
    return None

# get chain sequence and Lys coordinates
def get_chain_info(chain):
    seq = []
    lys_coords = []
    lys_resids = []
    resid_start = None
    for residue in chain:
        if residue.get_id()[0] == ' ':  # Only standard amino acids
            resname = residue.get_resname().strip()
            if resname in aa_map:
                seq.append(aa_map[resname])
                if resid_start is None:
                    resid_start = residue.id[1]

            # Get Lys Nζ coordinates
            nz = get_lys_nz(residue)
            if nz is not None:
                lys_coords.append(nz)
                lys_resids.append(residue.id[1] - resid_start + 1)  # 1-based indexing
    return ''.join(seq), np.array(lys_coords), lys_resids


chain_ids = [chain.id for model in structure for chain in model]
print(f"Chains in structure: {chain_ids}")

chain_infos = {}
for model in structure:
    for chain in model:
        for chain_id in chain_ids:
            if chain.id == chain_id:
                chain_info = get_chain_info(chain)
                chain_infos[chain.id] = chain_info

for chain_id in chain_ids:
    print(f" {chain_id}: {chain_infos[chain_id][0]}")

# assume two unique protein chains
selected_chains = input("Select two chains (e.g., AC) to generate plausible crosslinks: ")

chain_ids = [chain_id for chain_id in chain_infos if chain_id in selected_chains]
if len(chain_ids) != 2:
    raise ValueError("Please select exactly two valid chains.")


def get_cross_linking(p50_coords, p50_resids, p65_coords, p65_resids):
    max_dist = 30.0  # Å (typical for DSS/BS3)
    cross_links = []
    for i, c1 in enumerate(p50_coords):
        for j, c2 in enumerate(p65_coords):
            dist = np.linalg.norm(c1 - c2)
            if dist <= max_dist:
                cross_links.append((p50_resids[i], p65_resids[j], dist))
    return cross_links


# Compute inter-chain Lys-Lys distances
cross_links = get_cross_linking(
    chain_infos[chain_ids[0]][1], chain_infos[chain_ids[0]][2],
    chain_infos[chain_ids[1]][1], chain_infos[chain_ids[1]][2]
)

# Map chain IDs to 'A' and 'B'
id_map = {chain_ids[0]: 'A', chain_ids[1]: 'B'}

# Write to FASTA
output_fasta = pdb_file.replace(".pdb", ".fasta")
with open(output_fasta, 'w') as f:
    for chain_id in chain_ids:
        f.write(f">{id_map[chain_id]}\n")
        sequence, _, _ = chain_infos[chain_id]
        f.write(sequence + "\n")

# Write cross-links to file
output_xl = pdb_file.replace(".pdb", ".csv")
with open(output_xl, 'w') as f:
    for cross_link in cross_links:
        f.write(f"{cross_link[0]} {id_map[chain_ids[0]]} {cross_link[1]} {id_map[chain_ids[1]]} 0.2\n")
