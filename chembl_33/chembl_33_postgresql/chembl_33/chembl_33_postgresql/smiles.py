import csv
from rdkit import Chem

# Nome dos arquivos
input_file = "pkidb.tsv"
output_file = "valid_smiles_only.tsv"

# Lê o arquivo TSV e escreve apenas SMILES válidos em um novo arquivo
with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
    reader = csv.DictReader(infile, delimiter='\t')
    writer = csv.writer(outfile, delimiter='\t')
    
    # Escreve os SMILES no arquivo de saída se forem válidos
    for row in reader:
        smiles = row['Canonical_Smiles']
        if not smiles:  # Se o SMILES estiver vazio ou None, continue para a próxima linha
            continue
        try:
            molecule = Chem.MolFromSmiles(smiles)
            if molecule:  # A conversão foi bem-sucedida
                writer.writerow([smiles])
        except:
            # Se houver um erro ao tentar converter o SMILES, apenas ignore e continue
            continue

print(f"SMILES válidos foram salvos em {output_file}")
