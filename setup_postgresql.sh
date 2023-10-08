#!/bin/bash

# Install PostgreSQL
brew install postgresql

# Start PostgreSQL
brew services start postgresql

# Give a short delay to ensure PostgreSQL has started
sleep 5

# Start PostgreSQL and create the `chembl_23` database
PSQL_USER="sulfierry" # Change this to your user 
psql -U $PSQL_USER <<EOF
create database chembl_33;
\q
EOF

# Restore Database
DUMP_FILE="./chembl_23_postgresql.dmp" # Change this to your dump file path if different
pg_restore -U $PSQL_USER -d chembl_33 --no-owner -n public $DUMP_FILE

echo "Database restoration complete!"
