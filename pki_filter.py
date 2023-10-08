import csv
from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, Lipinski, rdMolDescriptors
from tqdm import tqdm
import concurrent.futures

def process_molecule(data):
    molregno, smiles = data
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            MW = Descriptors.MolWt(mol)
            CLogP = Crippen.MolLogP(mol)
            HBD = Lipinski.NumHDonors(mol)
            HBA = Lipinski.NumHAcceptors(mol)
            TPSA = rdMolDescriptors.CalcTPSA(mol)
            NRB = Lipinski.NumRotatableBonds(mol)
            NAR = rdMolDescriptors.CalcNumAromaticRings(mol)
            NCA = rdMolDescriptors.CalcNumRings(mol)

            if (314 <= MW <= 613) and (0.7 <= CLogP <= 6.3) and (0 <= HBD <= 4) and (3 <= HBA <= 10) and \
               (55 <= TPSA <= 138) and (1 <= NRB <= 11) and (1 <= NAR <= 5) and (0 <= NCA <= 2):
                return (True, molregno, smiles)
        return (False, None, None)
    except Exception as e:
        print(f"Error processing {smiles}: {e}")
        return (False, None, None)

if __name__ == "__main__":
    input_file = "filtered_chembl.tsv"
    output_file = "filtered_out.tsv"

    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        reader = csv.reader(infile, delimiter='\t')
        writer = csv.writer(outfile, delimiter='\t')
        
        # Write header
        writer.writerow(["molregno", "canonical_smiles"])
        
        # Skip header
        next(reader)
        
        data_list = [(row[0], row[1]) for row in reader]
        
        with concurrent.futures.ProcessPoolExecutor(max_workers=4) as executor:
            for is_valid, molregno, smiles in tqdm(executor.map(process_molecule, data_list), total=len(data_list)):
                if is_valid:
                    writer.writerow([molregno, smiles])
