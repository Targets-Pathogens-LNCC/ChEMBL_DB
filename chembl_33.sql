
CREATE TABLE public.compounds_all AS
SELECT molregno, canonical_smiles
FROM public.compound_structures;

SELECT * FROM public.compounds_all;

CREATE TABLE public.filtered_chembl_33_IC50 AS (
    SELECT DISTINCT cs.molregno, cs.canonical_smiles
    FROM public.compound_records AS cr
    JOIN public.compound_structures AS cs ON cr.molregno = cs.molregno
    JOIN public.activities AS act ON cr.molregno = act.molregno
    WHERE (act.standard_type = 'IC50' OR act.standard_type = 'Ki' OR act.standard_type = 'Kd') 
    AND act.standard_value > 6
    AND cs.canonical_smiles IS NOT NULL AND cs.canonical_smiles != ''
);


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

COPY (SELECT molregno, canonical_smiles FROM public.compound_structures) TO '/Users/sulfierry/Desktop/thil/chemblDB/chembl_33/chembl_33_molecules.tsv' WITH (FORMAT 'csv', DELIMITER E'\t', HEADER);

COPY public.filtered_chembl_33_IC50 TO '/Users/sulfierry/Desktop/thil/chemblDB/chembl_33/filtered_chembl_33_IC50.tsv' WITH (FORMAT 'csv', DELIMITER E'\t', HEADER);

COPY (SELECT DISTINCT kinase_name FROM public.kinase_ligand) TO '/Users/sulfierry/Desktop/thil/chemblDB/chembl_33/kinase_all_chembl_33.tsv' WITH CSV HEADER DELIMITER E'\t';

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
TO '/Users/sulfierry/Desktop/thil/chemblDB/chembl_33/kinase_ligand_chembl_33.tsv' WITH CSV HEADER DELIMITER E'\t';

COPY (SELECT kinase_name, SUM(number_of_ligands) as total_ligands FROM public.kinase_ligand GROUP BY kinase_name ORDER BY total_ligands DESC LIMIT 10) TO '/Users/sulfierry/Desktop/thil/chemblDB/chembl_33/kinase_ligand_top10_chembl_33.tsv' WITH CSV HEADER DELIMITER E'\t';

