#!/usr/bin/env python

#########################
# last update:
# Tue 17 Jun 2014 06:25:26 PM CEST
# [JMass]
#########################

import pickle
import getopt
import sys
import helpers.matrixparser as matrixparser
import helpers.fastahelper as fastahelper
import os
import multiprocessing
import queue
import time

package_directory = os.path.dirname(os.path.abspath(__file__))


def usage():
    """Usage"""
    print("""
	nmw.py -m MAT -f FILE_A -g FILE_B -o OUTFILE

  general options:
    -m, --matrix=MAT    (see matrix options)
    -f, --fasta_a=FILE    multifasta file a
    -g, --fasta_b=FILE    multifasta file b
    -o, --outfile=FILE    store results in FILE
    -h, --help            prints this
    -n, --num_cores=NUM   number of processors to use [default 1]

  alignment options:
     -O, --gap_open=FLOAT           gap opening cost [default 10]
     -E, --gap_extend=FLOAT         gap extension cost [default 0.5]

  fasta options:
    -d, --delim=STRING               split fasta headers at STRING
    -a, --asID=INT                   use INTth part of fasta header as transcript-ID
                                     (default:0)
  matrix options:
    MAT can be any of   ["EBLOSUM62", "EBLOSUM30",
                         "EBLOSUM45", "EBLOSUM60",
                         "EBLOSUM65", "EBLOSUM80",
                         "EBLOSUMN", "EBLOSUM35",
                         "EBLOSUM50","EBLOSUM62",
                         "EBLOSUM70","EBLOSUM85",
                         "EBLOSUM40","EBLOSUM55",
                         "EBLOSUM62-12","EBLOSUM75",
                         "EBLOSUM90"]
    """)
    sys.exit(2)

class Sequence(object):
    def __init__(self, header, sequence):
        self.header = header
        self.sequence = sequence.strip()


class Score(object):
    """Sequence seqA, Sequence seqB"""
    def __init__(self, seqA, seqB, scoringMatrix, gapOpen=10, gapExtend=0.5):
        self.seqA = seqA
        self.seqB = seqB
        self.scoringMatrix = scoringMatrix
        self.gapOpen = gapOpen
        self.gapExtend = gapExtend
        self.score = self.nmw()


    def nmw(self):
        score = None
        tmp =  Matrix(len(self.seqA.sequence)+1, len(self.seqB.sequence)+1)
        tmpmat = tmp.matrix
        tmpmat[0][0] = 0
        #print(tmp.gap)

        """ needleman wunsch algo """
        for i,si in enumerate(self.seqA.sequence):
            for j,sj in enumerate(self.seqB.sequence):
                tmp.gap[i][j] = False
                if tmp.gap[i][j-1]:
                    left = tmpmat[i][j-1]-self.gapExtend
                else:
                    left = tmpmat[i][j-1]-self.gapOpen
                if tmp.gap[i-1][j]:
                    top = tmpmat[i-1][j]- self.gapExtend
                else:
                    top = tmpmat[i-1][j] - self.gapOpen
                diag = tmpmat[i-1][j-1]+int(self.scoringMatrix.matrix[si][sj])
                choice = [diag, left, top]
                tmpmat[i][j] = max(choice)
                if max(choice) in [left,top]:
                    tmp.gap[i][j] = True

        score = tmp.getMax()
        return (score)

    def __repr__(self):
        if not self.score:
            self.nmw()
        s = ""
        s += "{}\t{}\t{}".format(self.seqA.header.split(" ")[0], self.seqB.header.split(" ")[0],self.score)
        return(s)

class Matrix(object):
    def __init__(self, n, m):
        self.n = n
        self.m = m
        self.matrix = [[0 for x in range(0,m)] for y in range(0,n)]
        self.gap = [[False for x in range(0,m)] for y in range(0,n)]

    def __repr__(self):
        s = ""
        for i in range(0,self.n):
            s+=str(i)+" "
            for j in range(0,self.m):
                s+=" "+str(self.matrix[i][j])
            s+="\n"
        return(s)

    def getMax(self):
        mx = 0
        for i in range(self.n-2, self.n-1):
            for j in range(0,self.m):
                if self.matrix[i][j] > mx:
                    mx = self.matrix[i][j]
        return(mx)

class HelperMatrix(object):
    """Stores tuples (i,j-1), (i-1,j), (i-1,j-1) """
    def __init__(self):
        self.arr = []

class ScoringMatrix(object):
    """Get SubstitutionMatrices from databin"""
    avail = ["EBLOSUM62", "EBLOSUM30", "EBLOSUM45", "EBLOSUM60", "EBLOSUM65",
            "EBLOSUM80", "EBLOSUMN", "EBLOSUM35", "EBLOSUM50",
            "EBLOSUM62", "EBLOSUM70","EBLOSUM85", "EBLOSUM40",
            "EBLOSUM55", "EBLOSUM62-12","EBLOSUM75", "EBLOSUM90"]

    def __init__(self,name):
        name = name.upper()
        self.matrix = {}
        if name in ScoringMatrix.avail:
            self.matrix = pickle.load(open(package_directory+os.sep+"databin"+os.sep+name,'rb'))
        else :
            raise Error("no such thing {}".format(name))

def main():
    avail = ["EBLOSUM62", "EBLOSUM30", "EBLOSUM45", "EBLOSUM60", "EBLOSUM65",
            "EBLOSUM80", "EBLOSUMN", "EBLOSUM35", "EBLOSUM50",
            "EBLOSUM62", "EBLOSUM70","EBLOSUM85", "EBLOSUM40",
            "EBLOSUM55", "EBLOSUM62-12","EBLOSUM75", "EBLOSUM90"]


    matrix = None
    fasta_a = None
    fasta_b = None
    outfile = None
    num_cores = 1
    try:
        opts, args = getopt.gnu_getopt(sys.argv[1:], "m:f:g:o:n:h",["matrix=","fasta_a=","fasta_b=","outfile=","num_cores","help"])
    except getopt.GetoptError as err:
        sys.stderr.write(str(err))
        usage()
    for o, a in opts:
        if o in ("-m", "--matrix"):
            matrix = a
        elif o in ("-f", "--fasta_a"):
            fasta_a = a
        elif o in ("-g", "--fasta_b"):
            fasta_b = a
        elif o in ("-h", "--help"):
            usage()
        elif o in ("-o", "--outfile"):
            outfile = a
        elif o in ("-n", "--num_cores"):
            num_cores = int(a)
        else:
            assert False, "unhandled option"

    if not matrix:
        usage()
    if matrix not in avail:
        usage()
    if not fasta_a or not fasta_b:
        usage()
    if not outfile:
        outfile = "out"
    if num_cores < 2:
        seen = set()
        res = " "
        sm = ScoringMatrix(matrix)
        fpa = fastahelper.FastaParser().read(fasta_a, " ", 0)
        for ha,sa in fpa:
            seen.add(ha)
            a = Sequence(ha,sa)
            fpb = fastahelper.FastaParser().read(fasta_b, " ", 0)
            for hb,sb in fpb:
                if hb in seen:
                    continue
                b = Sequence(hb,sb)
                s = Score(seqA = a, seqB = b, scoringMatrix = sm)
                res +=str(s)+"\n"
        with open(outfile,'w')as out:
            out.write(res)
    else:
        global SEMAPHORE
        global LOCK
        SEMAPHORE = multiprocessing.BoundedSemaphore(num_cores)
        LOCK = multiprocessing.Lock()
        nmw_multi(matrix=matrix, fasta_a= fasta_a, fasta_b = fasta_b, outfile= outfile, num_cores = num_cores)



def nmw_multi(matrix, fasta_a, fasta_b, outfile, num_cores):
    seen = set()
    tasks = multiprocessing.JoinableQueue()
    sm = ScoringMatrix(matrix)
    fpa = fastahelper.FastaParser().read(fasta_a, " ", 0)
    for ha,sa in fpa:
        seen.add(ha)
        a = Sequence(ha,sa)
        fpb = fastahelper.FastaParser().read(fasta_b, " ", 0)
        for hb,sb in fpb:
            if hb in seen:
                continue
            b = Sequence(hb,sb)
            tasks.put(Task(seqA=a, seqB=b, scoringMatrix = sm))
    for i in range(num_cores):
        tasks.put(None)

    consumers = [ Consumer(tasks, outfile) for i in range(num_cores) ]
    for c in consumers:
        c.start()
    for c in consumers:
        c.join()


class Task(object):
    def __init__(self, seqA, seqB, scoringMatrix):
        self.seqA = seqA
        self.seqB = seqB
        self.scoringMatrix = scoringMatrix

    def call(self):
        #print(self,"called")
        return(Score(seqA = self.seqA, seqB = self.seqB, scoringMatrix=self.scoringMatrix))

class Consumer(multiprocessing.Process):
    def __init__(self,taskq, outfile):
        multiprocessing.Process.__init__(self)
        self.taskq = taskq
        self.outfile = outfile
        global SEMAPHORE
        global LOCK
    def run(self):
        SEMAPHORE.acquire()
        while True:
            t = self.taskq.get()
            if t is None:
                print("quit", self)
                break
            res = t.call()
            LOCK.acquire()
            with open(self.outfile,'a') as out:
                out.write(str(res)+"\n")
            LOCK.release()
        SEMAPHORE.release()





if __name__=="__main__":
    main()
