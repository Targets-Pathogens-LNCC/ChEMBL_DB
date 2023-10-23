import csv
import sys
from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, Lipinski, rdMolDescriptors
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

class Descritores:
    def __init__(self, input_file, output_file):
        self.input_file = input_file
        self.output_file = output_file

    @staticmethod
    def find_smiles_column(reader):
        for row in reader:
            for idx, field in enumerate(row):
                if Chem.MolFromSmiles(field) is not None:
                    return idx
        return None

    @staticmethod
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

    def worker(self, args):
        row, smiles_column = args
        if len(row) > smiles_column:
            data = (None, row[smiles_column])
            return self.process_molecule(data)
        return (False, None, None, None, None, None, None, None, None, None)

    def run(self):
        with open(self.input_file, 'r') as infile, open(self.output_file, 'w', newline='') as outfile:
            reader = csv.reader(infile, delimiter='\t')
            writer = csv.writer(outfile, delimiter='\t')

            smiles_column = self.find_smiles_column(reader)
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
                
                results = list(tqdm(executor.map(self.worker, args_iter), total=len(rows)))

                
                for is_valid, smiles, MW, CLogP, HBD, HBA, TPSA, NRB, NAR, NCA in results:
                    if is_valid:
                        writer.writerow([smiles, MW, CLogP, HBD, HBA, TPSA, NRB, NAR, NCA])

class VisualizarDescritores:
    def __init__(self, pkidb_file, out_pkidb_file, output_directory):
        self.pkidb_file = pkidb_file
        self.out_pkidb_file = out_pkidb_file
        self.output_directory = output_directory

    def run(self):
        # Carregar os conjuntos de dados
        pkidb_data = pd.read_csv(self.pkidb_file, sep='\t')
        out_pkidb_data = pd.read_csv(self.out_pkidb_file, sep='\t')

        # Renomear as colunas do conjunto de dados out_pkidb para corresponder às colunas do conjunto de dados pkidb
        out_pkidb_data = out_pkidb_data.rename(columns={"CLogP": "LogP"})

        # Configurar o estilo dos gráficos
        sns.set(style="whitegrid")

        # Lista de descritores para análise
        descritores = ['MW', 'LogP', 'HBD', 'HBA', 'TPSA', 'NRB']

        # Criar uma figura e um conjunto de subplots
        fig, axs = plt.subplots(nrows=3, ncols=2, figsize=(15, 12))

        # Ajustar o espaço entre os plots
        plt.tight_layout(h_pad=3, w_pad=2, rect=[0.05, 0, 1, 1])

        # Adicionar a legenda do eixo Y no meio do gráfico
        fig.text(0.04, 0.5, 'Frequency', va='center', rotation='vertical', fontsize=12)

        # Plotar os histogramas para cada descritor
        for i, descritor in enumerate(descritores):
            # Calcular a posição do subplot
            row = i // 2
            col = i % 2
            
            # Plotar o histograma para o conjunto de dados 'pkidb'
            sns.histplot(pkidb_data[descritor], bins=30, ax=axs[row, col], kde=False, color="blue", label="pkidb" if i==0 else "")
            
            # Plotar o histograma para o conjunto de dados 'out_pkidb'
            sns.histplot(out_pkidb_data[descritor], bins=30, ax=axs[row, col], kde=False, color="orange", label="out_pkidb" if i==0 else "")
            
            # Configurar o título e os rótulos
            axs[row, col].set_title(f'Distribution of {descritor}')
            axs[row, col].set_xlabel(descritor)
            axs[row, col].set(ylabel=None)

        # Adicionar a legenda fora dos plots
        fig.legend(loc="upper right", bbox_to_anchor=(1,1), bbox_transform=axs[0, 1].transAxes)

        # Definir o caminho do arquivo de saída
        os.makedirs(self.output_directory, exist_ok=True)
        output_path = os.path.join(self.output_directory, "descritores_histogramas.png")

        # Salvar a figura completa
        plt.savefig(output_path)
        plt.plot()

        # Mostrar o caminho do arquivo salvo
        print("Histograma salvo em:", output_path)

class QuantificarDescritores:
    def __init__(self, pkidb_file, output_file):
        self.pkidb_file = pkidb_file
        self.output_file = output_file

    def run(self):
        # Carregar o conjunto de dados
        pkidb_data = pd.read_csv(self.pkidb_file, sep='\t')

        # Lista de descritores para análise
        descritores = ['MW', 'LogP', 'HBD', 'HBA', 'TPSA', 'NRB']

        # Inicializar um dicionário para armazenar os resultados
        quantificacoes = {}

        # Calcular os valores de fronteira para cada descritor
        for descritor in descritores:
            minimo = np.min(pkidb_data[descritor])
            percentil_25 = np.percentile(pkidb_data[descritor], 25)
            mediana = np.median(pkidb_data[descritor])
            media = np.mean(pkidb_data[descritor])
            percentil_75 = np.percentile(pkidb_data[descritor], 75)
            maximo = np.max(pkidb_data[descritor])
            
            quantificacoes[descritor] = [minimo, percentil_25, mediana, media, percentil_75, maximo]

        # Salvar os resultados em um arquivo .tsv
        quantificacoes_df = pd.DataFrame.from_dict(quantificacoes, orient='index', columns=['Min', '25th Percentile', 'Median', 'Mean', '75th Percentile', 'Max'])
        quantificacoes_df.to_csv(self.output_file, sep='\t')

        print(f'Os valores de fronteira foram salvos em: {self.output_file}')

if __name__ == "__main__":
    descritores = Descritores('../PKIDB/pkidb_2023-06-30.tsv', 'out_pkidb222.tsv')
    descritores.run()
    
    visualizar_descritores = VisualizarDescritores('../PKIDB/pkidb_2023-06-30.tsv', 'out_pkidb.tsv', 'output_directory')
    visualizar_descritores.run()

    quantificar_descritores = QuantificarDescritores('../PKIDB/pkidb_2023-06-30.tsv', 'output_quantificacoes.tsv')
    quantificar_descritores.run()
