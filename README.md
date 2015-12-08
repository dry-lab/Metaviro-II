# Metaviro-II

Repository of the *Metaviro-II Machine Learning* project.

Its goal is to provide a complete pipeline for metagenomics classification:
- Creating and feeding a classifier
  - Downloading reference data from [NCBI](http://www.ncbi.nlm.nih.gov/) databases
  - Preprocessing and learning from downloaded data
- Running a classification daemon

## Features

- Metagenomics classification on *contigs*, using *k-mers frequencies*
- *Online* learning
- *Space efficient* and *updatable* models

## Prerequisite software

- [Vowpal Wabbit](https://github.com/JohnLangford/vowpal_wabbit/) machine learning system (package `vowpal-wabbit` available for some Debian-based distributions)
- Python 3 with module [Biopython](http://biopython.org/wiki/Main_Page) (package `python3-biopython` or `python-biopython` available for some distributions)
- SQLite3

## Installation

```
$ sh/create_db.sh
```

Will take some time. At least 25G of space needed.

Create a SQLite3 taxonomic database associating GIs and taxIDs. This database is constructed from [NCBI taxonomy database](ftp://ftp.ncbi.nih.gov/pub/taxonomy/)

## Begin learning

```
$ sh/init_vw.sh
$ python3 python/get_data.py
```

1. Start the classification daemon with an empty model.
2. Start the learning process (downloading, preprocessing, learning).

## Classification

The Vowpal Wabbit daemon is listening on port 1234.

Each instance to be classified has to be processed (k-merization, [Vowpal Wabbit formatting](https://github.com/JohnLangford/vowpal_wabbit/wiki/Input-format)) and then can be sent to the listening port.

Classification results will be sent back by the daemon.

