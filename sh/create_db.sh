#!/bin/bash

path='taxdump/'
pypath='../python/'
db='ncbi.db'

cd $path
echo 'Downloading taxonomy files...'
wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/gi_taxid_nucl.dmp.gz
gzip -d gi_taxid_nucl.dmp.gz
wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.Z
tar -xvzf taxdump.tar.Z

echo 'Creating taxonomy database...'
sqlite3 $db < create_db.sql
echo 'Speeding up database...'
python3 ${pypath}pre_pyphy.py names.dmp nodes.dmp ncbi.db
sqlite3 $db < speed_db.sql
echo 'Testing database'
sqlite3 $db < test_db.sql

rm taxdump.tar.Z
rm *.dmp
rm *.txt
rm *.prt
