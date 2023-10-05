# ChEMBL_DB_exploring

# instalar
brew install postgresql

# iniciar postgresql
brew services start postgresql

# criar bd
psql -U seu_usuario
create database chembl_23;
\q

# restaurar banco de dados
pg_restore -U sulfierry -d chembl_23 --no-owner -n public ./chembl_23_postgresql.dmp

# listar todas as tabelas
\dt
# utilize este para o chembl_23
\dt public.*

-- Selecionar todos os compostos sem aplicar nenhum filtro e criar tabela persistente 'compounds_all'
CREATE TABLE public.compounds_all AS
SELECT molregno, canonical_smiles
FROM public.compound_structures;

-- Salvar esta seleção em .tsv
COPY (SELECT molregno, canonical_smiles FROM public.compound_structures) TO '/caminho/para/o/arquivo/chembl_23_molecules.tsv' WITH (FORMAT 'csv', DELIMITER E'\t', HEADER);


# para visualizar a tabela recem criada
SELECT * FROM public.compounds_all;


-- Criar uma tabela persistente filtered_chembl_23
CREATE TABLE public.filtered_chembl_23 AS (
    SELECT DISTINCT cs.molregno, cs.canonical_smiles
    FROM public.compound_records AS cr
    JOIN public.compound_structures AS cs ON cr.molregno = cs.molregno
    JOIN public.activities AS act ON cr.molregno = act.molregno
    WHERE (act.standard_type = 'IC50' OR act.standard_type = 'Ki' OR act.standard_type = 'Kd') 
    AND act.standard_value > 6
    AND cs.canonical_smiles IS NOT NULL AND cs.canonical_smiles != ''
);

-- Exportar os dados para um arquivo .tsv
COPY public.filtered_chembl_23 TO 'filtered_chembl_23.tsv' WITH (FORMAT 'csv', DELIMITER E'\t', HEADER);

— Excluir tabela, por exemplo ‘kinase_all’
DROP TABLE public.kinase_all;

— Examinando o esquema da tabela ‘kinase_ligand
\d public.kinase_ligand;

— Para obter e salvar todas as cinases únicas (sem redundância) presentes na tabela kinase_ligand:
COPY (SELECT DISTINCT kinase_name FROM public.kinase_ligand) TO '/Users/sulfierry/Desktop/thil/chemblDB/chembl_23/chembl_23_postgresql/kinase_all_chembl_23.tsv' WITH CSV HEADER DELIMITER E'\t';

— Para obter e salvar a quantidade de ligantes associados a cada cinase: 
COPY (SELECT kinase_name, SUM(number_of_ligands) FROM public.kinase_ligand GROUP BY kinase_name ORDER BY kinase_name) TO '/Users/sulfierry/Desktop/thil/chemblDB/chembl_23/chembl_23_postgresql/kinase_ligand_chembl_23.tsv' WITH CSV HEADER DELIMITER E'\t';

— Para obter e salvar o top 10 de cinases com o maior número de ligantes: 
COPY (SELECT kinase_name, SUM(number_of_ligands) as total_ligands FROM public.kinase_ligand GROUP BY kinase_name ORDER BY total_ligands DESC LIMIT 10) TO '/Users/sulfierry/Desktop/thil/chemblDB/chembl_23/chembl_23_postgresql/kinase_ligand_top10_chembl_23.tsv' WITH CSV HEADER DELIMITER E'\t';

