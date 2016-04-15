# encoding utf-8

from Bio import Entrez
from Bio import SeqIO

import time
import logging

from urllib.error import HTTPError, URLError
from http.client import HTTPException
import socket
socket.setdefaulttimeout(1800)

Entrez.email = 'marie.gasparoux@u-bordeaux.fr'

error_sleep_delay = 60
encountered_errors = (HTTPError, HTTPException, URLError,
                      ConnectionError, socket.timeout, TimeoutError)


def get_NCBI_slice_sequence(gi, start, stop):
    handle = Entrez.efetch(db='nuccore', id=str(gi),
                           seq_start=str(start), seq_stop=str(stop), rettype='fasta')
    records = list(SeqIO.parse(handle, "fasta"))
    handle.close()
    return str(records[0].seq)


def get_NCBI_gis_pos(ids_list, error_msg=False):
    handle = Entrez.esummary(db='nuccore', id=ids_list)
    records = list(Entrez.parse(handle))
    handle.close()
    return records


get_gis_pos = { 'name':'esummary', 'func':get_NCBI_gis_pos}
get_gis_slice_seq = { 'name':'efetch', 'func':get_NCBI_slice_sequence}


def exception_procedure(origin, e, delay):
    logger.error("%s: %s" % (origin, str(e)))
    logger.info("Sleeping %ds, then trying again" % (delay * error_sleep_delay))
    time.sleep(delay * error_sleep_delay)


def get(ncbi_action_dict, *args):
    result = False
    delay = 1
    while not result:
        try:
            records = ncbi_action_dict['func'](*args)
            result = True
        except encountered_errors as e:
            exception_procedure(ncbi_action_dict['name'], e, delay)
            delay *= 2
    if records == []:
        logger.info("%s: obtained empty record list, trying again"
                    % ncbi_action_dict['name'])
        return ncbi_action_dict['func'](*args)
    return records
