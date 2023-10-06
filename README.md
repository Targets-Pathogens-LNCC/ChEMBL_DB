# ChEMBL_DB_exploring

Explore the ChEMBL ([Chembl_33](https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/)) database using SQL.

The following example is applied to the identification of ligands related to protein kinases.


## Installation, Configuration, and SQL Commands

```bash
# Install PostgreSQL
brew install postgresql

# Start PostgreSQL
brew services start postgresql

# Start PostgreSQL and create the `chembl_23` database
psql -U your_user
create database chembl_23;
\q

# Restore Database
pg_restore -U sulfierry -d chembl_23 --no-owner -n public ./chembl_23_postgresql.dmp
```

## Select all compounds without any filter and create persistent table 'compounds_all'

```bash
CREATE TABLE public.compounds_all AS
SELECT molregno, canonical_smiles
FROM public.compound_structures;

-- Save this selection in .tsv
COPY (SELECT molregno, canonical_smiles FROM public.compound_structures) TO '/Users/sulfierry/Desktop/thil/chemblDB/chembl_33/chembl_33_molecules.tsv' WITH (FORMAT 'csv', DELIMITER E'\t', HEADER);

-- to view the newly created table
SELECT * FROM public.compounds_all;

-- Create a persistent table filtered_chembl_33
CREATE TABLE public.filtered_chembl_33_IC50 AS (
    SELECT DISTINCT cs.molregno, cs.canonical_smiles
    FROM public.compound_records AS cr
    JOIN public.compound_structures AS cs ON cr.molregno = cs.molregno
    JOIN public.activities AS act ON cr.molregno = act.molregno
    WHERE (act.standard_type = 'IC50' OR act.standard_type = 'Ki' OR act.standard_type = 'Kd') 
    AND act.standard_value > 6
    AND cs.canonical_smiles IS NOT NULL AND cs.canonical_smiles != ''
);

-- Export the data to a .tsv file
COPY public.filtered_chembl_33_IC50 TO '/Users/sulfierry/Desktop/thil/chemblDB/chembl_33/filtered_chembl_33_IC50.tsv' WITH (FORMAT 'csv', DELIMITER E'\t', HEADER);

```
## Create the 'kinase_ligand' table in the 'public' schema of 'chembl_33'

```bash
CREATE TABLE public.kinase_ligand AS
SELECT
-- Select the targets preferred name and label it as 'kinase_name'
    t.pref_name AS kinase_name,                                             
-- Count distinct ligands associated with the target and label the count as 'number_of_ligands'
    COUNT(DISTINCT act.molregno) AS number_of_ligands                     
FROM
-- Activities table, which includes information about ligands and their activities
    activities act                                                          
JOIN
-- Join with 'assays' table to get information about the assay in which the ligand was tested
    assays ass ON act.assay_id = ass.assay_id                            
JOIN
-- Join with 'target_dictionary' table to get information about the target
    target_dictionary t ON ass.tid = t.tid                              
WHERE
-- Filter to consider only targets that are individual proteins
    t.target_type = 'SINGLE PROTEIN' AND
-- Filter to consider only targets with names containing the word 'kinase'
    t.pref_name LIKE '%kinase%'                                       
GROUP BY
-- Group by the targets preferred name to get a unique count of ligands per target
    t.pref_name;

-- Export the data to a .tsv file
COPY public.filtered_chembl_33_IC50 TO '/Users/sulfierry/Desktop/thil/chemblDB/chembl_33/filtered_chembl_33_IC50.tsv' WITH (FORMAT 'csv', DELIMITER E'\t', HEADER);

-- To retrieve and save all unique kinases (non-redundant) present in the kinase_ligand table:
COPY (SELECT DISTINCT kinase_name FROM public.kinase_ligand) TO '/Users/sulfierry/Desktop/thil/chemblDB/chembl_33/kinase_all_chembl_33.tsv' WITH CSV HEADER DELIMITER E'\t';

                                                                  
```
## To retrieve and save the number of ligands associated with each kinase:

```bash
-- Start a data copy operation to a file
COPY (
    -- Use a CTE (Common Table Expression) to calculate the total ligands
    WITH TotalLigands AS (
        -- Select the total sum of ligands for all kinases
        SELECT SUM(number_of_ligands) AS total_ligands_of_kinases
        FROM public.kinase_ligand
    )
    -- Select data from the kinase_ligand table
    SELECT 
        k.kinase_name, 
        -- Calculate the sum of ligands for each kinase
        SUM(k.number_of_ligands) AS ligands_per_kinase,
        -- Check if its the first row of the result
        CASE
            WHEN ROW_NUMBER() OVER(ORDER BY SUM(k.number_of_ligands) DESC, k.kinase_name) = 1 THEN t.total_ligands_of_kinases
            -- If it's not the first row, set to NULL
            ELSE NULL
        END AS total_ligands_of_kinases
    FROM 
        -- Join the data from the kinase_ligand table with the total ligands
        public.kinase_ligand k
    CROSS JOIN TotalLigands t
    -- Group the results by kinase name and total ligands
    GROUP BY 
        k.kinase_name, t.total_ligands_of_kinases 
    -- Order the results by the sum of ligands in descending order and then by kinase name
    ORDER BY 
        ligands_per_kinase DESC, k.kinase_name
)
-- Specify the location and format of the file to which the data will be copied
TO '/Users/sulfierry/Desktop/thil/chemblDB/chembl_33/kinase_ligand_chembl_33.tsv' WITH CSV HEADER DELIMITER E'\t';

-- To retrieve and save the top 10 kinases with the highest number of ligands:
COPY (SELECT kinase_name, SUM(number_of_ligands) as total_ligands FROM public.kinase_ligand GROUP BY kinase_name ORDER BY total_ligands DESC LIMIT 10) TO '/Users/sulfierry/Desktop/thil/chemblDB/chembl_33/kinase_ligand_top10_chembl_33.tsv' WITH CSV HEADER DELIMITER E'\t';

```
## SQL Utils
```bash
-- remove the chembl_23 database
DROP DATABASE chembl_23;

-- Drop table, for example, ‘kinase_all’
DROP TABLE public.kinase_all;

-- Inspecting the schema of the ‘kinase_ligand table
\d public.kinase_ligand;
```


