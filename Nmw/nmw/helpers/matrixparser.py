#!/usr/bin/env python
import pickle
import os

class MatrixParser(object):
    def __init__(self, matrixf):
        self.matrix = self.readMatrix(matrixf)

    def readMatrix(self, matrixf):
        first = True
        aa1 = []
        aa2d = {}
        with open(matrixf, 'r') as inf:
            for l in inf:
                l = l.strip()
                if l.startswith("#"):
                    pass
                else:
                    if first:
                        aa1 = l.split(" ")
                        aa1 = [a for a in aa1 if a !=" " and a!="" ]
                        first = False
                    else:
                        tmp = l.split(" ")
                        tmp = [t for t in tmp if t !=" " and t !=""]
                        aa2d[tmp[0]] ={}
                        for i, a in enumerate(aa1[:-1]):
                            aa2d[l.split(" ")[0]] [a] = int(tmp[i+1])

        return(aa2d)

    def dumpMatrix(self, matrixp):
        with open(matrixp, 'wb') as out:
            pickle.dump(self.matrix, out)
    #def dumpMatrix(self, matrixf, matrixp):
    #    print("overload?")


class MatrixReader(object):
    def __init__(self, matrixp):
        self.matrix = pickle.load(open(matrixp,"rb"))


def main():
    sourcedir = "data"
    mp = MatrixParser("data/EBLOSUM62")
    print(mp.matrix)
    mp.dumpMatrix("databin/EBLOSUM62")
    mr = MatrixReader("databin/EBLOSUM62")
    print(mr.matrix)
    print("\n")
    for f in mr.matrix:
        for s in mr.matrix[f]:
            print(mr.matrix[f][s],f,s)
    for f in os.listdir(sourcedir):
        print(f)
        fn = os.path.basename(f)
        mp = MatrixParser(sourcedir+os.sep+f)
        mp.dumpMatrix(matrixp = "databin"+os.sep+fn)
if __name__ == "__main__":
     main()
