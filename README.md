# Metaviro-II
Repository of the *Metaviro-II Machine Learning* project

## Prerequisite software

- [Vowpal Wabbit](https://github.com/JohnLangford/vowpal_wabbit/) machine learning system (package `vowpal-wabbit` available for some Debian-based distributions)
- Python 3 with module [Biopython](http://biopython.org/wiki/Main_Page) (package `python3-biopython` or `python-biopython` available for some distributions)
- SQLite3

## Installation

```
$ sh/create_db.sh
```
Will take some time. At least 25G of space needed.

## Begin classification

```
$ sh/init_vw.sh
$ python3 python/get_data.py
```
