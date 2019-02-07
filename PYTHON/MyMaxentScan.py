#!/usr/bin/env python3
'''
Author: James Casaletto
Date: 28.Dec.2018
'''


from maxentpy import maxent
import sys


def main():

    mySpliceSiteScores = SpliceSiteScores(sys.argv[1])
    mySpliceSiteScores.findSpliceSites()


class SpliceSiteScores:

    def __init__(self, inputFile):
        self.fivePrimeSites = dict()
        self.threePrimeSites = dict()
        self.inputData = open(inputFile, 'r')
        self.fivePrimeSiteLength = 9
        self.threePrimeSiteLength = 23

    def findSpliceSites(self):
        # chr15:28235701-28235871,chr15:28235789-28235790,G,A,1921,682,1467,315,3223,1142,1276,367,attgggactgtgac...
        # first line is schema, every line after that is data
        i = 0
        for line in self.inputData:
            if i == 0:
                print(line.strip() + "," + "wt5start" + "," + "wt5sequence" + "," "wt5score"  +
                      "," "mu5start" + "," + "mu5sequence" + ","  + "mu5score" +
                      "," + "wt3start" + "," + "wt3sequence" + "," + "wt3score" +
                      "," + "mu3start" + "," + "mu3sequence" + "," + "mu3score")
            else:
                splitLine = line.strip().split(',')
                sequenceStart = int(splitLine[0].split(":")[1].split("-")[0])
                mutationStart = int(splitLine[1].split(":")[1].split("-")[0])
                mutationTuple = (splitLine[2], splitLine[3])
                matrix5 = maxent.load_matrix5()
                self.fivePrimeSites[splitLine[0]] = self.check5Prime(splitLine[12].strip(), self.fivePrimeSiteLength,
                                                                sequenceStart, mutationStart, mutationTuple, matrix5)
                matrix3 = maxent.load_matrix3()
                self.threePrimeSites[splitLine[0]] = self.check3Prime(splitLine[12].strip(), self.threePrimeSiteLength,
                                                                sequenceStart, mutationStart, mutationTuple, matrix3)

                print(line.strip() + "," +
                      str(self.fivePrimeSites[splitLine[0]][0]) + "," +
                      str(self.fivePrimeSites[splitLine[0]][1]) + "," +
                      str(self.fivePrimeSites[splitLine[0]][2]) + "," +
                      str(self.fivePrimeSites[splitLine[0]][3]) + "," +
                      str(self.fivePrimeSites[splitLine[0]][4]) + "," +
                      str(self.fivePrimeSites[splitLine[0]][5]) + "," +

                      str(self.threePrimeSites[splitLine[0]][0]) + "," +
                      str(self.threePrimeSites[splitLine[0]][1]) + "," +
                      str(self.threePrimeSites[splitLine[0]][2]) + "," +
                      str(self.threePrimeSites[splitLine[0]][3]) + "," +
                      str(self.threePrimeSites[splitLine[0]][4]) + "," +
                      str(self.threePrimeSites[splitLine[0]][5]))
            i += 1

    def check5Prime(self, sequence, length, sequenceStart, mutationStart, mutationTuple, matrix5):
        wtMaxScore = -99.0
        muMaxScore = -99.0
        wtMaxStart = 0
        muMaxStart = 0
        wtMaxSequence = None
        muMaxSequence = None
        mutationOffset = mutationStart - sequenceStart
        mutatedSequence = sequence[:mutationOffset] + mutationTuple[1].lower() + sequence[mutationOffset+1:]

        for i in range(length, 0, -1):
            start = mutationOffset - i
            end = start + length
            if (end - start) != length:
                continue
            wtSequence = sequence[start:end]
            muSequence = mutatedSequence[start:end].strip()
            try:
                wtSequenceScore = maxent.score5(wtSequence, matrix5)
                muSequenceScore = maxent.score5(muSequence, matrix5)
            except:
                sys.stderr.write("maxent failure")
                continue
            if(wtSequenceScore > wtMaxScore):
                wtMaxScore = wtSequenceScore
                wtMaxStart = start
                wtMaxSequence = wtSequence
            if (muSequenceScore > muMaxScore):
                muMaxScore = muSequenceScore
                muMaxStart = start
                muMaxSequence = muSequence
        return (wtMaxStart + sequenceStart, wtMaxSequence, wtMaxScore, muMaxStart + sequenceStart, muMaxSequence,
                muMaxScore)

    def check3Prime(self, sequence, length, sequenceStart, mutationStart, mutationTuple, matrix3):
        wtMaxScore = -99.0
        muMaxScore = -99.0
        wtMaxStart = 0
        muMaxStart = 0
        wtMaxSequence = None
        muMaxSequence = None
        mutationOffset = mutationStart - sequenceStart
        mutatedSequence = sequence[:mutationOffset] + mutationTuple[1].lower() + sequence[mutationOffset+1:]

        for i in range(length, 0, -1):
            start = mutationOffset - i
            end = start + length
            if (end - start) != length:
                continue
            wtSequence = sequence[start:end]
            muSequence = mutatedSequence[start:end].strip()
            try:
                wtSequenceScore = maxent.score3(wtSequence, matrix3)
                muSequenceScore = maxent.score3(muSequence, matrix3)
            except:
                sys.stderr.write("maxent failure")
                continue
            if(wtSequenceScore > wtMaxScore):
                wtMaxScore = wtSequenceScore
                wtMaxStart = start
                wtMaxSequence = wtSequence
            if (muSequenceScore > muMaxScore):
                muMaxScore = muSequenceScore
                muMaxStart = start
                muMaxSequence = muSequence
        return (wtMaxStart + sequenceStart, wtMaxSequence, wtMaxScore, muMaxStart + sequenceStart, muMaxSequence,
                muMaxScore)

    def printSitesAndScores(self):
        print("3 prime sequences, offsets, and scores: ")
        for locus in self.threePrimeSites:
            print(locus + ":" + str(self.threePrimeSites[locus]))
        print("5 prime sequences, offsets, and scores: ")
        for locus in self.fivePrimeSites:
            print(locus + ":" + str(self.fivePrimeSites[locus]))

if __name__ == "__main__":
    main()
