import os
import sys
import re
import argparse
import pandas as pd
from datetime import datetime
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import RDLogger
import warnings

RDLogger.DisableLog('rdApp.*')
warnings.filterwarnings("ignore")

def clean_control_chars(s):
    """Remove ANSI control characters, such as ^[[200~ and ^[[201~."""
    return re.sub(r'\x1b\[[0-9;]*[~A-Za-z]', '', s)

def count_heavy_atoms(molecule):
    """Count the number of heavy (non-hydrogen) atoms in a molecule."""
    return sum(1 for atom in molecule.GetAtoms() if atom.GetAtomicNum() != 1)

def get_heavy_atom_difference(smarts_reaction):
    """Parse a SMARTS reaction and calculate the difference in heavy atom counts between reactant and product."""
    parts = smarts_reaction.strip().split('>>')
    if len(parts) != 2:
        raise ValueError("Invalid SMARTS reaction format")
    before = Chem.MolFromSmarts(parts[0].strip())
    after = Chem.MolFromSmarts(parts[1].strip())
    if before is None or after is None:
        raise ValueError("Invalid SMARTS patterns in reaction")
    return count_heavy_atoms(after) - count_heavy_atoms(before)

def main():
    parser = argparse.ArgumentParser(description='Apply reactions to a SMILES molecule.')
    parser.add_argument('--smiles', type=str, help='Input SMILES string')
    args = parser.parse_args()

    # Method 1: Input from command-line arguments
    if args.smiles:
        input_smiles = clean_control_chars(args.smiles)
    else:
        # Method 2: Interactive input
        input_smiles = clean_control_chars(input("Enter SMILES: "))

    # Set file paths
    reaction_file_path = r'./Reaction_Fsp3-rich.csv'
    output_folder = datetime.now().strftime("%Y%m%d_%H%M%S")
    os.makedirs(output_folder, exist_ok=True)

    print("Loading reaction data...")
    df_reactions = pd.read_csv(reaction_file_path, encoding='latin-1')

    # Process the input molecule
    print(f"Processing molecule: {input_smiles}")
    mol = Chem.MolFromSmiles(input_smiles)
    if mol is None:
        print(f"Invalid SMILES: {input_smiles}")
        sys.exit(1)

    Chem.SanitizeMol(mol)
    input_heavy_atoms = count_heavy_atoms(mol)
    input_canonical = Chem.MolToSmiles(mol)
    seen_products = set()
    results = []

    print("Applying reactions and generating products...")
    for idx, row in df_reactions.iterrows():
        smarts_reaction = row['smarts_reaction']
        try:
            reaction = AllChem.ReactionFromSmarts(smarts_reaction)
            reaction.Initialize()
            products = reaction.RunReactants((mol,))
        except Exception:
            print("Skipping invalid reaction")
            continue

        for prod_set in products:
            for prod in prod_set:
                try:
                    prod_smiles = Chem.MolToSmiles(prod, isomericSmiles=True)
                    if not prod_smiles:
                        continue

                    try:
                        Chem.Kekulize(prod, clearAromaticFlags=True)
                        Chem.SanitizeMol(prod)
                    except Exception:
                        continue
                    prod_heavy_atoms = count_heavy_atoms(prod)
                    diff = get_heavy_atom_difference(smarts_reaction)
                    if abs(prod_heavy_atoms - input_heavy_atoms) <= abs(diff) and prod_smiles not in seen_products:
                        seen_products.add(prod_smiles)
                        results.append({
                            'Original_Molecule_SMILES': input_canonical,
                            'Reaction_SMARTS': smarts_reaction,
                            'Product_SMILES': prod_smiles
                        })
                except Exception:
                    continue

    if results:
        print("Saving generated products...")
        df_results = pd.DataFrame(results).drop_duplicates()
        result_csv = os.path.join(output_folder, 'reaction_products.csv')
        df_results.to_csv(result_csv, index=False)
        print(f"Product results saved to: {result_csv}")
    else:
        print("No products generated.")

if __name__ == '__main__':
    main()
