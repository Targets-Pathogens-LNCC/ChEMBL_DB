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
-- Iniciar o PostgreSQL e criar o banco de dados `chembl_23`
psql -U seu_usuario
create database chembl_23;
\q

-- Restaurar Banco de Dados
pg_restore -U sulfierry -d chembl_23 --no-owner -n public ./chembl_23_postgresql.dmp

-- Listar Todas as Tabelas
\dt OR \dt public.*

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
