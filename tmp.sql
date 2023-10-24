CREATE TABLE public.compounds_all AS
SELECT molregno, canonical_smiles
FROM public.compound_structures;

-- to view the newly created table
SELECT * FROM public.compounds_all;

-- Save this selection in .tsv
COPY (SELECT molregno, canonical_smiles FROM public.compound_structures) TO '/Users/sulfierry/Desktop/thil/chemblDB/chembl_33/chembl_33_molecules.tsv' WITH (FORMAT 'csv', DELIMITER E'\t', HEADER);
