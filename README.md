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

```bash
-- Selecionar todos os compostos sem aplicar nenhum filtro e criar tabela persistente 'compounds_all'
CREATE TABLE public.compounds_all AS
SELECT molregno, canonical_smiles
FROM public.compound_structures;

-- Salvar seleção em .tsv
COPY (SELECT molregno, canonical_smiles FROM public.compound_structures) TO '/Users/sulfierry/Desktop/thil/chemblDB/chembl_33/chembl_33_molecules.tsv' WITH (FORMAT 'csv', DELIMITER E'\t', HEADER);

... (insira aqui os outros comandos SQL em sequência)

-- Remover Base de Dados `chembl_23`
DROP DATABASE chembl_23;

-- Excluir tabela (exemplo ‘kinase_all’)
DROP TABLE public.kinase_all;

-- Examinando o esquema da tabela ‘kinase_ligand'
\d public.kinase_ligand;
