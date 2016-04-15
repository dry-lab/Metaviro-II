# encoding utf-8

import socket


def print_vw_features(g_label, sk_dict, g_cat=None):
    if g_cat != None:
        feats = "%d '%d |" % (g_cat, g_label)
    else:
        feats = "'%d |" % g_label
    for k, v in sorted(sk_dict.items()):
        feats += " %s:%.2f" % (k, v)
    return feats + "\n"


def send_to_daemon(port, vw_data):
    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    s.connect(('localhost', int(port)))
    s.sendall(vw_data.encode())
    s.shutdown(socket.SHUT_WR)
