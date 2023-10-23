import csv
import sys
from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, Lipinski, rdMolDescriptors
from tqdm import tqdm

def find_smiles_column(reader):
    for row in reader:
        for idx, field in enumerate(row):
            if Chem.MolFromSmiles(field) is not None:
                return idx
    return None

def process_molecule(data):
    _, smiles = data  # Ignorando molregno
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
            
            # Log das propriedades computadas
            print(f"SMILES: {smiles}, MW: {MW}, CLogP: {CLogP}, HBD: {HBD}, HBA: {HBA}, TPSA: {TPSA}, NRB: {NRB}, NAR: {NAR}, NCA: {NCA}")
            
            # Verificar as condições de filtragem
            if (314 <= MW <= 613) and (0.7 <= CLogP <= 6.3) and (0 <= HBD <= 4) and (3 <= HBA <= 10) and \
               (55 <= TPSA <= 138) and (1 <= NRB <= 11) and (1 <= NAR <= 5) and (0 <= NCA <= 2):
                return (True, smiles, MW, CLogP, HBD, HBA, TPSA, NRB, NAR, NCA)
        return (False, None, None, None, None, None, None, None, None, None)
    except Exception as e:
        print(f"Erro ao processar {smiles}: {e}")
        return (False, None, None, None, None, None, None, None, None, None)

if __name__ == "__main__":
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        reader = csv.reader(infile, delimiter='\t')
        writer = csv.writer(outfile, delimiter='\t')

        smiles_column = find_smiles_column(reader)
        if smiles_column is None:
            print("Não foi possível identificar a coluna SMILES.")
            sys.exit(1)

        infile.seek(0)  # Voltar ao início do arquivo
        writer.writerow(["canonical_smiles", "MW", "CLogP", "HBD", "HBA", "TPSA", "NRB", "NAR", "NCA"])

        for row in tqdm(reader):
            if len(row) > smiles_column:
                data = (None, row[smiles_column])
                is_valid, smiles, MW, CLogP, HBD, HBA, TPSA, NRB, NAR, NCA = process_molecule(data)
                if is_valid:
                    writer.writerow([smiles, MW, CLogP, HBD, HBA, TPSA, NRB, NAR, NCA])
