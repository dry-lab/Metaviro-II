# encoding utf-8

import sqlite3
from contextlib import closing

db = "taxdump/ncbi.db"
unknown = -1
no_rank = "no rank"

#pyphy.getTaxidByName("Bacteria",1)
def getTaxidByName(name,limit=1):
    with sqlite3.connect(db) as conn:
        with closing(conn.cursor()) as cursor:
            cursor.execute("SELECT taxid FROM tree WHERE name = ?", (name,))
            results = cursor.fetchall()

    taxid_list = [result[0] for result in results]
    if len(taxid_list) > 0:
        taxid_list.sort()
        return taxid_list[:limit]
    else:
        return [unknown]

#pyphy.getRankByTaxid("2")
def getRankByTaxid(taxid):
    with sqlite3.connect(db) as conn:
        with closing(conn.cursor()) as cursor:
            cursor.execute("SELECT rank FROM tree WHERE taxid = ?", (taxid,))
            result = cursor.fetchone()

    if not result:
        return no_rank
    return '{0}'.format(result[0])

#pyphy.getRankByName("Bacteria")
def getRankByName(name):
    try:
        return getRankByTaxid(getTaxidByName(name)[0])
    except:
        return no_rank

#pyphy.getNameByTaxid("2")
def getNameByTaxid(taxid):
    with sqlite3.connect(db) as conn:
        with closing(conn.cursor()) as cursor:
            cursor.execute("SELECT name FROM tree WHERE taxid = ?", (taxid,))
            result = cursor.fetchone()

    if not result:
        return "unknown"
    return '{0}'.format(result[0])

def getParentByTaxid(taxid):
    with sqlite3.connect(db) as conn:
        with closing(conn.cursor()) as cursor:
            cursor.execute("SELECT parent FROM tree WHERE taxid = ?", (taxid,))
            result = cursor.fetchone()

    if not result:
        return unknown
    return result[0]

#pyphy.getParentByName("Flavobacteriia")
def getParentByName(name):
    try:
        return getParentByTaxid(getTaxidByName(name)[0])
    except:
        return unknown

def getPathByTaxid(taxid):
    path = []
    current_id = int(taxid)
    path.append(current_id)
    
    while current_id != 1 and current_id != unknown:
        current_id = int(getParentByTaxid(current_id))
        path.append(current_id)
    
    return path[::-1]
    
def getTaxidByGi(gi):
    with sqlite3.connect(db) as conn:
        with closing(conn.cursor()) as cursor:
            cursor.execute("SELECT taxid FROM gi_taxid WHERE gi = ?", (gi,))
            result = cursor.fetchone()

    if not result:
        return unknown
    return result[0]
    
def getSonsByTaxid(taxid):
    with sqlite3.connect(db) as conn:
        with closing(conn.cursor()) as cursor:
            cursor.execute("SELECT taxid FROM tree WHERE parent = ?", (taxid,))
            result = [row[0] for row in cursor]

    return result

def getSonsByName(name):
    return getSonsByTaxid(getTaxidByName(name)[0])

def getGiByTaxid(taxid):
    conn = sqlite3.connect(db)
    command = "SELECT gi FROM gi_taxid WHERE taxid = '" + str(taxid) +  "';"
    cursor = conn.cursor()
    #print(command)
    try:
        result = [row[0] for row in cursor.execute(command)] 
        
        cursor.close()
        return result
    except Exception as e:
        print(e)

def getAllSonsByTaxid(taxid):
    from threading import Thread
    import queue
    in_queue = queue.Queue()
    out_queue = queue.Queue()

#what a thread is supposed to do
    def work():
        while True:
            sonId = in_queue.get()
#if here use getAllSonsByTaxid, it will be true recursive,
#but it will run over the thread limit set by the os by trying "Flavobacteriia" (6000+)
#error: can't start new thread
#it is more elegant here

            for s_s_id in getSonsByTaxid(sonId):
                #print("s_s_id%s" % str(s_s_id))
                out_queue.put(s_s_id)
                in_queue.put(s_s_id)
       
            in_queue.task_done()

#spawn 20 threads
    for i in range(20):
        
        t = Thread(target=work)
        t.daemon = True
        t.start()

#feed the queue
    for son in getSonsByTaxid(taxid):
        out_queue.put(son)
        in_queue.put(son)
    
    in_queue.join()   
    result = []
#reap the results
    while not out_queue.empty():
        result.append( out_queue.get())
    return result

def getAllGiByTaxid(taxid):
    from threading import Thread
    import queue
    in_queue = queue.Queue()
    out_queue = queue.Queue()

    allSons = getAllSonsByTaxid(taxid)

    def work():
        while True:
            sonId = in_queue.get()
            out_queue.put(getGiByTaxid(sonId))
            in_queue.task_done()

    for i in range(20):
        
        t = Thread(target=work)
        t.daemon = True
        t.start()
        
    in_queue.put(taxid)

    for son in allSons:
        in_queue.put(son)
    
    in_queue.join()        
    
    result = []
    while not out_queue.empty():
        result += out_queue.get()
    return result


