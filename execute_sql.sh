#!/bin/bash

# Define path variables
PATH_CHEMBL_33_MOLECULES="/Users/sulfierry/Desktop/thil/chemblDB/chembl_33/chembl_33_molecules.tsv"
PATH_FILTERED_CHEMBL_33_IC50="/Users/sulfierry/Desktop/thil/chemblDB/chembl_33/filtered_chembl_33_IC50.tsv"
PATH_KINASE_ALL_CHEMBL_33="/Users/sulfierry/Desktop/thil/chemblDB/chembl_33/kinase_all_chembl_33.tsv"
PATH_KINASE_LIGAND_CHEMBL_33="/Users/sulfierry/Desktop/thil/chemblDB/chembl_33/kinase_ligand_chembl_33.tsv"
PATH_KINASE_LIGAND_TOP10_CHEMBL_33="/Users/sulfierry/Desktop/thil/chemblDB/chembl_33/kinase_ligand_top10_chembl_33.tsv"

# Execute SQL commands
psql -U sulfierry <<EOF

-- Create compounds_all table
CREATE TABLE public.compounds_all AS
SELECT molregno, canonical_smiles
FROM public.compound_structures;

-- Save this selection in .tsv
COPY (SELECT molregno, canonical_smiles FROM public.compound_structures) TO '$PATH_CHEMBL_33_MOLECULES' WITH (FORMAT 'csv', DELIMITER E'\t', HEADER);

-- View the newly created table
SELECT * FROM public.compounds_all;

-- Create filtered_chembl_33_IC50 table
CREATE TABLE public.filtered_chembl_33_IC50 AS (
    SELECT DISTINCT cs.molregno, cs.canonical_smiles
    FROM public.compound_records AS cr
    JOIN public.compound_structures AS cs ON cr.molregno = cs.molregno
    JOIN public.activities AS act ON cr.molregno = act.molregno
    WHERE (act.standard_type = 'IC50' OR act.standard_type = 'Ki' OR act.standard_type = 'Kd') 
    AND act.standard_value > 6
    AND cs.canonical_smiles IS NOT NULL AND cs.canonical_smiles != ''
);

-- Export the data to .tsv file
COPY public.filtered_chembl_33_IC50 TO '$PATH_FILTERED_CHEMBL_33_IC50' WITH (FORMAT 'csv', DELIMITER E'\t', HEADER);

-- Create kinase_ligand table
CREATE TABLE public.kinase_ligand AS
SELECT
    t.pref_name AS kinase_name,                                             
    COUNT(DISTINCT act.molregno) AS number_of_ligands                     
FROM
    activities act                                                          
JOIN
    assays ass ON act.assay_id = ass.assay_id                            
JOIN
    target_dictionary t ON ass.tid = t.tid                              
WHERE
    t.target_type = 'SINGLE PROTEIN' AND
    t.pref_name LIKE '%kinase%'                                       
GROUP BY
    t.pref_name;

-- Export the data to a .tsv file
COPY public.filtered_chembl_33_IC50 TO '$PATH_KINASE_ALL_CHEMBL_33' WITH (FORMAT 'csv', DELIMITER E'\t', HEADER);

-- Save all unique kinases
COPY (SELECT DISTINCT kinase_name FROM public.kinase_ligand) TO '$PATH_KINASE_LIGAND_CHEMBL_33' WITH CSV HEADER DELIMITER E'\t';

-- Start a data copy operation to a file
COPY (
    WITH TotalLigands AS (
        SELECT SUM(number_of_ligands) AS total_ligands_of_kinases
        FROM public.kinase_ligand
    )
    SELECT 
        k.kinase_name, 
        SUM(k.number_of_ligands) AS ligands_per_kinase,
        CASE
            WHEN ROW_NUMBER() OVER(ORDER BY SUM(k.number_of_ligands) DESC, k.kinase_name) = 1 THEN t.total_ligands_of_kinases
            ELSE NULL
        END AS total_ligands_of_kinases
    FROM 
        public.kinase_ligand k
    CROSS JOIN TotalLigands t
    GROUP BY 
        k.kinase_name, t.total_ligands_of_kinases 
    ORDER BY 
        ligands_per_kinase DESC, k.kinase_name
)
TO '$PATH_KINASE_LIGAND_TOP10_CHEMBL_33' WITH CSV HEADER DELIMITER E'\t';

EOF

echo "SQL commands executed successfully!"
