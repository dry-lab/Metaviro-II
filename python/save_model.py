# encoding utf-8

import sys
import socket
from datetime import datetime
import pickle


def main(argv=None):
    now = datetime.now()
    timestr = now.strftime("%d%B_%H:%M")

    with open("pattern_port_dict.pickle", 'rb') as input:
        pattern_port_dict = pickle.load(input)
    with open("port_dest_dict.pickle", 'rb') as input:
        port_dest_dict = pickle.load(input)

    for pattern in pattern_port_dict.keys():
        for port in pattern_port_dict[pattern]:
            dest = port_dest_dict[port]
            send_to_daemon(model_list[0], "log/models/%s/%s_%s.mod" %(dest, pattern, timestr))
    return 0

if __name__ == '__main__':
    sys.exit(main())
