# ChEMBL_DB_exploring

Explore o banco de dados ChEMBL com comandos SQL úteis.

## Instalação, Configuração e Comandos SQL

```bash
# Instalar PostgreSQL
brew install postgresql

# Iniciar PostgreSQL
brew services start postgresql

# Iniciar o PostgreSQL e criar o banco de dados `chembl_23`
psql -U seu_usuario
create database chembl_23;
\q

# Restaurar Banco de Dados
pg_restore -U sulfierry -d chembl_23 --no-owner -n public ./chembl_23_postgresql.dmp
```

## Selecionar todos os compostos sem aplicar nenhum filtro e criar tabela persistente 'compounds_all'

```bash
CREATE TABLE public.compounds_all AS
SELECT molregno, canonical_smiles
FROM public.compound_structures;

-- Salvar esta seleção em .tsv
COPY (SELECT molregno, canonical_smiles FROM public.compound_structures) TO '/Users/sulfierry/Desktop/thil/chemblDB/chembl_33/chembl_33_molecules.tsv' WITH (FORMAT 'csv', DELIMITER E'\t', HEADER);

-- para visualizar a tabela recém criada
SELECT * FROM public.compounds_all;

-- Criar uma tabela persistente filtered_chembl_33
CREATE TABLE public.filtered_chembl_33_IC50 AS (
    SELECT DISTINCT cs.molregno, cs.canonical_smiles
    FROM public.compound_records AS cr
    JOIN public.compound_structures AS cs ON cr.molregno = cs.molregno
    JOIN public.activities AS act ON cr.molregno = act.molregno
    WHERE (act.standard_type = 'IC50' OR act.standard_type = 'Ki' OR act.standard_type = 'Kd') 
    AND act.standard_value > 6
    AND cs.canonical_smiles IS NOT NULL AND cs.canonical_smiles != ''
);

-- Exportar os dados para um arquivo .tsv
COPY public.filtered_chembl_33_IC50 TO '/Users/sulfierry/Desktop/thil/chemblDB/chembl_33/filtered_chembl_33_IC50.tsv' WITH (FORMAT 'csv', DELIMITER E'\t', HEADER);

```
# Cria a tabela 'kinase_ligand' no schema 'public' de 'chembl_33'

```bash
CREATE TABLE public.kinase_ligand AS
SELECT
-- Seleciona o nome preferencial do alvo (target) e o nomeia como 'kinase_name' 
    t.pref_name AS kinase_name,            			         
-- Conta os ligantes distintos associados ao alvo e nomeia a contagem como 'number_of_ligands'
    COUNT(DISTINCT act.molregno) AS number_of_ligands         
FROM
-- Tabela de atividades, que inclui informações sobre os ligantes e suas atividades 
    activities act                       				      
JOIN
-- Junta com a tabela 'assays' para obter informações sobre o ensaio em que o ligante foi testado 
    assays ass ON act.assay_id = ass.assay_id                
JOIN
-- Junta com a tabela 'target_dictionary' para obter informações sobre o alvo (target) 
    target_dictionary t ON ass.tid = t.tid                                      
WHERE
-- Filtra para considerar apenas alvos que são proteínas individuais 
    t.target_type = 'SINGLE PROTEIN' AND
-- Filtra para considerar apenas alvos com nomes que contêm a palavra 'kinase'                              
    t.pref_name LIKE '%kinase%'                                               
GROUP BY
-- Agrupa por nome preferencial do alvo para obter uma contagem única de ligantes por alvo 
    t.pref_name;                                                                         
```
