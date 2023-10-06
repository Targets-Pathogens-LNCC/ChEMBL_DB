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



## Selecionar todos os compostos sem aplicar nenhum filtro e criar tabela persistente 'compounds_all'
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


-- Cria a tabela 'kinase_ligand' no schema 'public' de 'chembl_33'
CREATE TABLE public.kinase_ligand AS
SELECT 
    t.pref_name AS kinase_name,            			          -- Seleciona o nome preferencial do alvo (target) e o nomeia como 'kinase_name'
    COUNT(DISTINCT act.molregno) AS number_of_ligands  -- Conta os ligantes distintos associados ao alvo e nomeia a contagem como 'number_of_ligands'
FROM 
    activities act                       				           -- Tabela de atividades, que inclui informações sobre os ligantes e suas atividades
JOIN 
    assays ass ON act.assay_id = ass.assay_id                        -- Junta com a tabela 'assays' para obter informações sobre o ensaio em que o ligante foi testado
JOIN 
    target_dictionary t ON ass.tid = t.tid                                      -- Junta com a tabela 'target_dictionary' para obter informações sobre o alvo (target)
WHERE 
    t.target_type = 'SINGLE PROTEIN' AND                              -- Filtra para considerar apenas alvos que são proteínas individuais
    t.pref_name LIKE '%kinase%'                                               -- Filtra para considerar apenas alvos com nomes que contêm a palavra 'kinase'
GROUP BY 
    t.pref_name;                                                                          -- Agrupa por nome preferencial do alvo para obter uma contagem única de ligantes por alvo


-- Exportar os dados para um arquivo .tsv
COPY public.filtered_chembl_33_IC50 TO '/Users/sulfierry/Desktop/thil/chemblDB/chembl_33/filtered_chembl_33_IC50.tsv' WITH (FORMAT 'csv', DELIMITER E'\t', HEADER);


— Para obter e salvar todas as cinases únicas (sem redundância) presentes na tabela kinase_ligand:
COPY (SELECT DISTINCT kinase_name FROM public.kinase_ligand) TO '/Users/sulfierry/Desktop/thil/chemblDB/chembl_33/kinase_all_chembl_33.tsv' WITH CSV HEADER DELIMITER E'\t';

— Para obter e salvar a quantidade de ligantes associados a cada cinase: 
-- Inicia uma operação de cópia de dados para um arquivo
COPY (
    -- Utiliza um CTE (Common Table Expression) para calcular o total de ligantes
    WITH TotalLigands AS (
        -- Seleciona a soma total de ligantes de todas as kinases
        SELECT SUM(number_of_ligands) AS total_ligands_of_kinases
        FROM public.kinase_ligand
    )
    -- Seleciona os dados da tabela kinase_ligand
    SELECT 
        k.kinase_name, 
        -- Calcula a soma de ligantes para cada kinase
        SUM(k.number_of_ligands) AS ligands_per_kinase,
        -- Verifica se é a primeira linha do resultado
        CASE
            WHEN ROW_NUMBER() OVER(ORDER BY SUM(k.number_of_ligands) DESC, k.kinase_name) = 1 THEN t.total_ligands_of_kinases
            -- Se não for a primeira linha, coloca NULL
            ELSE NULL
        END AS total_ligands_of_kinases
    FROM 
        -- Junta os dados da tabela kinase_ligand com o total de ligantes
        public.kinase_ligand k
    CROSS JOIN TotalLigands t
    -- Agrupa os resultados por nome da kinase e total de ligantes
    GROUP BY 
        k.kinase_name, t.total_ligands_of_kinases 
    -- Ordena os resultados pela soma de ligantes em ordem decrescente e, em seguida, pelo nome da kinase
    ORDER BY 
        ligands_per_kinase DESC, k.kinase_name
)
-- Especifica o local e o formato do arquivo para onde os dados serão copiados
TO '/Users/sulfierry/Desktop/thil/chemblDB/chembl_33/kinase_ligand_chembl_33.tsv' WITH CSV HEADER DELIMITER E'\t';

— Para obter e salvar o top 10 de cinases com o maior número de ligantes: 
COPY (SELECT kinase_name, SUM(number_of_ligands) as total_ligands FROM public.kinase_ligand GROUP BY kinase_name ORDER BY total_ligands DESC LIMIT 10) TO '/Users/sulfierry/Desktop/thil/chemblDB/chembl_33/kinase_ligand_top10_chembl_33.tsv' WITH CSV HEADER DELIMITER E'\t';


# remover base de dados chembl_23
DROP DATABASE chembl_23;

# Excluir tabela, por exemplo ‘kinase_all’
DROP TABLE public.kinase_all;

# Examinando o esquema da tabela ‘kinase_ligand
\d public.kinase_ligand;

