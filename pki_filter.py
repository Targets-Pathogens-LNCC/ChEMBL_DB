import csv
import sys
from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, Lipinski, rdMolDescriptors
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor

def find_smiles_column(reader):
    for row in reader:
        for idx, field in enumerate(row):
            if Chem.MolFromSmiles(field) is not None:
                return idx
    return None

def worker(args):
    row, smiles_column = args
    if len(row) > smiles_column:
        data = (None, row[smiles_column])
        return process_molecule(data)
    return (False, None, None, None, None, None, None, None, None, None)



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
    with open(input_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
        reader = csv.reader(infile, delimiter='\t')
        writer = csv.writer(outfile, delimiter='\t')

        smiles_column = find_smiles_column(reader)
        if smiles_column is None:
            print("Não foi possível identificar a coluna SMILES.")
            sys.exit(1)

        infile.seek(0)  # Voltar ao início do arquivo
        writer.writerow(["canonical_smiles", "MW", "CLogP", "HBD", "HBA", "TPSA", "NRB", "NAR", "NCA"])

        rows = list(reader)
        with ProcessPoolExecutor() as executor:
            # Cria um iterável de argumentos para passar para a função worker.
            # Cada item no iterável é uma tupla (row, smiles_column).
            args_iter = ((row, smiles_column) for row in rows)
            
            results = list(tqdm(executor.map(worker, args_iter), total=len(rows)))

            
            for is_valid, smiles, MW, CLogP, HBD, HBA, TPSA, NRB, NAR, NCA in results:
                if is_valid:
                    writer.writerow([smiles, MW, CLogP, HBD, HBA, TPSA, NRB, NAR, NCA])

### PLOT
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Carregar os conjuntos de dados
pkidb_data = pd.read_csv('../PKIDB/pkidb_2023-06-30.tsv', sep='\t')
out_pkidb_data = pd.read_csv('out_pkidb.tsv', sep='\t')

# Renomear as colunas do conjunto de dados out_pkidb para corresponder às colunas do conjunto de dados pkidb
out_pkidb_data = out_pkidb_data.rename(columns={"CLogP": "LogP"})

# Configurar o estilo dos gráficos
sns.set(style="whitegrid")

# Lista de descritores para análise
descritores = ['MW', 'LogP', 'HBD', 'HBA', 'TPSA', 'NRB']

# Criar uma figura e um conjunto de subplots
fig, axs = plt.subplots(nrows=3, ncols=2, figsize=(15, 12))

# Ajustar o espaço entre os plots
plt.tight_layout(h_pad=3)

# Plotar os histogramas para cada descritor
for i, descritor in enumerate(descritores):
    # Calcular a posição do subplot
    row = i // 2
    col = i % 2
    
    # Plotar o histograma para o conjunto de dados 'pkidb'
    sns.histplot(pkidb_data[descritor], bins=30, ax=axs[row, col], kde=False, color="blue", label="pkidb")
    
    # Plotar o histograma para o conjunto de dados 'out_pkidb'
    sns.histplot(out_pkidb_data[descritor], bins=30, ax=axs[row, col], kde=False, color="orange", label="out_pkidb")
    
    # Configurar o título e os rótulos
    axs[row, col].set_title(f'Distribuição de {descritor}')
    axs[row, col].set_xlabel(descritor)
    axs[row, col].set_ylabel('Frequência')

    # Adicionar a legenda
    axs[row, col].legend()

# Definir o caminho do arquivo de saída
output_directory = "histogramas"
os.makedirs(output_directory, exist_ok=True)
output_path = os.path.join(output_directory, "descritores_histogramas.png")

# Salvar a figura completa
plt.savefig(output_path)

# Mostrar o caminho do arquivo salvo
print("Histograma salvo em:", output_path)
