# Author: Conner Warnock
# This program takes an input txt file, parses it, generates a Huffman Code for it, and also computes
# a 2nd-Order Markov Chain based on the input
# October 17, 2020

import copy
import heapq
import binascii
import math

# Reads in file
FileData = []
with open("EE6743_grail_testfile.txt") as TestFile:
    for line in TestFile:
        FileData.append(line)


# The word doc lists a carriage return '\r', of which none were found in the text file
# Total of 29 characters including new line, space, and end of file
SourceAlphabet = ['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z',' ','\n','\x03']

# Gets frequencies/probabilities of characters in text
def getProbabilities(SourceAlphabet, FileData):
    Count = 0
    Sum = 0
    ProbList = []
    Temp = []
    for i in range(0, len(SourceAlphabet)):
        Temp.clear()
        for j in range(0, len(FileData)):
            for k in range(0, len(FileData[j])):
                if FileData[j][k] == SourceAlphabet[i]:
                    Count += 1
        Temp.append(Count)
        Temp.append(SourceAlphabet[i])
        ProbList.append(Temp[:])
        Sum = Sum + Count
        Count = 0
    for i in range(0, len(SourceAlphabet)):
        ProbList[i][0] = ProbList[i][0] / Sum

    return ProbList



# Creates binary Huffman code. Smallest probabilities get added together, and tree branch added
# 0 value assigned for values closer to end of ProbList, 1 assigned for values closer to beginning
def binaryHuffman(ProbList):
    # Heap contains weight, symbol, and space for Huffman code
    ProbHeap = [[weight, [character, '']] for weight, character in ProbList]
    # Min heap
    heapq.heapify(ProbHeap)
    while (len(ProbHeap)) > 1:
        # Takes lowest probabilities, pops from heap
        first = heapq.heappop(ProbHeap)
        second = heapq.heappop(ProbHeap)

        # Add 0 to Huffman code of first value
        for value in first[1:]:
            value[1] = '0' + value[1]
        # Add 1 to Huffman code of second value
        for value in second[1:]:
            value[1] = '1' + value[1]

        # Push added probabilities back into heap, with new Huffman code added
        heapq.heappush(ProbHeap, [first[0]+second[0]] + first[1:]+second[1:])

    return sorted(heapq.heappop(ProbHeap)[1:], key=lambda p: (len(p[-1]), p))


    return HuffmanCode

# Table showing character, probability, Huffman code representation, and code length
def createTable(ProbList, HuffmanCode):
    print("Character       |       Probability     |       Code      |   Length")
    print("____________________________________________________________________")

    for i in range(0, len(HuffmanCode)):
        character, code = HuffmanCode[i]
        prob = 0
        for i in range(0, len(ProbList)):
            prob2, character2 = ProbList[i]
            if character == character2:
                prob = prob2
                break
        line = [repr(character), prob, code, len(code)]
        print('{:<8}           {:<30} {:<15} {:<10}'.format(*line))

    return

# Computes entropy = p1*log2(1/p1)+p2*log2(1/p2)+...+pn*log2(1/pn)
def getEntropy(ProbList):
    Entropy = 0
    for i in range(0, len(ProbList)):
        Entropy = Entropy + (ProbList[i][0]*(math.log2(1/ProbList[i][0])))

    print("Entropy:", Entropy,"bits/source symbol")

    return Entropy

# Computes Lavg = p1*(length of Huffman Code symbol 1)+p2*(length of Huffman Code symbol 2)+...
def getLavg(ProbList, HuffmanCode):
    Lavg = 0
    for i in range(0, len(ProbList)):
        prob, character = ProbList[i]
        code = 2
        for j in range(0, len(HuffmanCode)):
            character2, code2 = HuffmanCode[j]
            if character == character2:
                code = code2
                break
        Lavg = Lavg + prob*len(code)

    print("Lavg:", Lavg, "bits/bit")

    return Lavg

# Gets array of first-order Markov probabilities, and conditional/source entropies
def getFirstMarkov(SourceAlphabet, FileData, ProbList):
    FirstMarkov = []
    Count = 0
    Sum = 0
    CharProb = [0]*len(SourceAlphabet)
    for i in range(0, len(SourceAlphabet)):
        CharProb = [0]*len(SourceAlphabet)
        for j in range(0, len(FileData)):
            for k in range(0, len(FileData[j])):
                if FileData[j][k] == SourceAlphabet[i]:
                    Count += 1
                    if SourceAlphabet[i] != '\x03':
                        if SourceAlphabet[i] != '\n':
                            NextLetter = FileData[j][k+1]
                        else:
                            NextLetter = FileData[j+1][0]
                        for m in range(0, len(SourceAlphabet)):
                            if NextLetter == SourceAlphabet[m]:
                                CharCount = CharProb[m]
                                CharCount = CharCount + 1
                                CharProb[m] = CharCount
        Sum = Sum + Count
        FirstMarkov.append(CharProb[:])
        for j in range(0, len(FirstMarkov[i])):
            FirstMarkov[i][j] = FirstMarkov[i][j] / Count
        Count = 0

    FirstMarkovSourceEntropy = 0
    FirstMarkovConditionalEntropies = []
    # Get conditional entropies
    for i in range(0, len(FirstMarkov)):
        ConditionalEntropy = 0
        for j in range(0, len(FirstMarkov[i])):
            if FirstMarkov[i][j] != 0:
                ConditionalEntropy = ConditionalEntropy + (FirstMarkov[i][j])*(math.log2(1/FirstMarkov[i][j]))
        FirstMarkovConditionalEntropies.append(ConditionalEntropy)
    # Get source entropy
    for i in range(0, len(FirstMarkovConditionalEntropies)):
        FirstMarkovSourceEntropy = FirstMarkovSourceEntropy + (ProbList[i][0])*(FirstMarkovConditionalEntropies[i])
    print("Source Entropy Based on First-Order Markov: ", FirstMarkovSourceEntropy)

    return FirstMarkov, FirstMarkovSourceEntropy

# Gets array of second-order Markov probabilities, and conditional/source entropies
def getSecondMarkov(SourceAlphabet, FileData, ProbList):
    SecondMarkov = []
    TwoCharProb = []
    TwoCharProbRow = []
    for i in range(0, len(SourceAlphabet)):
        TwoCharProbRow.clear()
        for j in range(0, len(SourceAlphabet)):
            TwoCharProbRow.append(0)
        TwoCharProb.append(TwoCharProbRow[:])
    for i in range(0, len(SourceAlphabet)):
        TwoCharProb.clear()
        TwoCharProb = []
        for x in range(0, len(SourceAlphabet)):
            TwoCharProbRow.clear()
            for y in range(0, len(SourceAlphabet)):
                TwoCharProbRow.append(0)
            TwoCharProb.append(TwoCharProbRow[:])
        for j in range(0, len(FileData)):
            for k in range(0, len(FileData[j])):
                if FileData[j][k] == SourceAlphabet[i]:
                    if SourceAlphabet[i] != '\x03':
                        if SourceAlphabet[i] != '\n':
                            SecondChar = FileData[j][k+1]
                        else:
                            SecondChar = FileData[j+1][0]
                            if SecondChar != '\x03':
                                if SecondChar != '\n':
                                    NextLetter = FileData[j+1][1]
                                else:
                                    NextLetter = FileData[j+2][1]
                        for m in range(0, len(SourceAlphabet)):
                            if SecondChar == SourceAlphabet[m]:
                                # Find next letter
                                if SecondChar != '\x03':
                                    if SecondChar != '\n' and SourceAlphabet[i] != '\n':
                                        NextLetter = FileData[j][k+2]
                                    else:
                                        NextLetter = FileData[j + 1][0]
                                NextLetterIndex = 0
                                for n in range(0, len(SourceAlphabet)):
                                    if SourceAlphabet[n] == NextLetter:
                                        NextLetterIndex = n
                                CharCount = TwoCharProb[m][NextLetterIndex]
                                CharCount = CharCount + 1
                                TwoCharProb[m][NextLetterIndex] = CharCount
        SecondMarkov.append(TwoCharProb[:])
    # Find state probs
    CountList = copy.deepcopy(SecondMarkov)
    StateCountList = []
    StateCountListRow = []
    StateProbList = []
    StateProbListRow = []
    CountSum = 0
    StateCount = 0
    for i in range(0, len(CountList)):
        StateCountListRow.clear()
        for j in range(0, len(CountList)):
            StateCount = 0
            for k in range(0, len(CountList)):
                StateCount = StateCount + CountList[i][j][k]
            CountSum = CountSum + StateCount
            StateCountListRow.append(StateCount)
        StateCountList.append(StateCountListRow[:])
    for i in range(0, len(StateCountList)):
        StateProbListRow.clear()
        for j in range(0, len(StateCountList)):
            StateProbListRow.append(StateCountList[i][j] / CountSum)
        StateProbList.append(StateProbListRow[:])


    # Find probabilities from count list
    for i in range(0, len(SecondMarkov)):
        for j in range(0, len(SecondMarkov)):
            Sum = 0
            for k in range(0, len(SecondMarkov)):
                Sum = Sum + SecondMarkov[i][j][k]
            for k in range(0, len(SecondMarkov)):
                if SecondMarkov[i][j][k] != 0:
                    SecondMarkov[i][j][k] = SecondMarkov[i][j][k] / Sum

    # Find conditional entropies
    SecondMarkovSourceEntropy = 0
    SecondMarkovConditionalEntropies = []
    SecondMarkovConditionalEntropiesRow = []
    for i in range(0, len(FirstMarkov)):
        SecondMarkovConditionalEntropiesRow.clear()
        for j in range(0, len(FirstMarkov)):
            ConditionalEntropy = 0
            for k in range(0, len(FirstMarkov)):
                if SecondMarkov[i][j][k] != 0:
                    ConditionalEntropy = ConditionalEntropy + (SecondMarkov[i][j][k])*(math.log2(1/SecondMarkov[i][j][k]))
            SecondMarkovConditionalEntropiesRow.append(ConditionalEntropy)
        SecondMarkovConditionalEntropies.append(SecondMarkovConditionalEntropiesRow[:])
    # Get source entropy
    for i in range(0, len(SecondMarkovConditionalEntropies)):
        for j in range(0, len(SecondMarkovConditionalEntropies)):
            SecondMarkovSourceEntropy = SecondMarkovSourceEntropy + (StateProbList[i][j])*(SecondMarkovConditionalEntropies[i][j])
    print("Source Entropy Based on Second-Order Markov: ", SecondMarkovSourceEntropy)

    return SecondMarkov, SecondMarkovSourceEntropy

ProbList = getProbabilities(SourceAlphabet, FileData)
HuffmanCode = binaryHuffman(ProbList)
createTable(ProbList, HuffmanCode)
Entropy = getEntropy(ProbList)
Lavg = getLavg(ProbList, HuffmanCode)
FirstMarkov, FirstMarkovSourceEntropy = getFirstMarkov(SourceAlphabet, FileData, ProbList)
SecondMarkov, SecondMarkovSourceEntropy = getSecondMarkov(SourceAlphabet, FileData, ProbList)