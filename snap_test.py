#/usr/bin/python
#
#

import snap

def tester():
    Graph = snap.(snap.PNGraph, 10, 100)
    print('generated')
    FOut = snap.TFOut("test.graph")
    print('file stream generated')
    Graph.Save(FOut)
    print('saved')
    FOut.Flush()

if __name__ == "__main__":
    tester()