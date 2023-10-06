#!/usr/bin/env python3
### 2022.10.27
### 2023.03.06 change rpLongestSoftClipped: add I D process, ignoring them causes failure of recheck because soft-clipped part is clipped by I or D which shouled have been transformed to S in replaced cigar, its potential bug about index movement is sovled too.
### 2023.04.11
### change lowFI to SomTD
### 2023.05.30
### add getDirection method for readGT
### 2023.08.02
### modify for the new data form of combine.pair.sam for the arbitrary mode and neural network mode
### 2023.09.03
### add the module for extracting read pairs in the arbitrary mode
### 2023.09.14
### add the module for parsing insertion information in the arbitrary mode
### 2023.09.20
### pseudoalignments are checked whether clipped after location parsing. By the way Genome alignments are checked at the very beginning
### when collecting reads while TE alignments just before location parsing.
### 2023.09.22
### add read pair accordance check for both clipped read pair
### 2023.09.27
### start deep mode

import re
from SomTD_GlobalParam import convertDict,converChrDict,convertTEDict,TElenDict



# ReadGTP functions
def revTuple(align,alignLen):
    """
    function for generating new cigartuples whose direction follows the original read
    """
    if align.cigartuples:
        if align.flag & 16 == 16:
            return list(reversed(align.cigartuples))
        else:
            return(align.cigartuples)
    else:
        return [(4,alignLen)] # 4 means soft-clipped which will be annotated as not aligned

def convertNuc(nuc):
    """
    function for converting a nucleotide string to an integer list, A:0, T:1, C:2, G:3, N:4
    """
    intList = []
    for i in nuc:
        intList.append(convertDict[i])
    return intList

def convertChr(chrName):
    """
    function for converting chr name to a integer list
    """
    if chrName:
        return [converChrDict[chrName]]
    else:
        return [0] # not aligned situation

def convertTE(TEname):
    """
    function for converting TE name to a integer list
    """
    if TEname:
        return [convertTEDict[TEname]]
    else:
        return [0] # not aligned situation

def convertP(XP):
    """
    function for converting XP information to a integer list
    """
    if XP != "NULL":
        if len(XP.split(";")) == 2:
            return[convertTEDict[XP.split(",")[0]]]
        else:
            return [-1]
    else:
        return [0] # not aligned situation

def locMapped(align,mode="origin"):
    """
    function for checking whether alignment is soft-clipped and where the mapped part locates, mode == origin for the original location
    and mode == forward for the location in chromsome/TE reference forward strand
    """
    if mode not in ["origin","forward"]:
        print("unexpected mode parameter")
        return None
    if align.reference_name == "polyA" and mode == "origin":
        if align.is_reverse:
            return 2 # in read tail
        else:
            return 1 # in read head
    cigarSimpString = ""
    if mode == "origin" and align.cigartuples:
        cigartuples = list(reversed(align.cigartuples)) if align.is_reverse else align.cigartuples
    elif not align.cigartuples:
        return None
    else:
        cigartuples = align.cigartuples
    for i in range(len(cigartuples)):
        if cigartuples[i][0] == 4 and cigartuples[i][1] >= 10: # soft-clipped part less than 10 bp in read head or tail is tolerated, make sure TE alignment also contains soft-clipped part
            cigarSimpString += "S"
        if cigartuples[i][0] == 5 and cigartuples[i][1] >= 10:
            cigarSimpString += "H"
        if cigartuples[i][0] == 0:
            cigarSimpString += "M"
    if "S" not in cigarSimpString and "H" not in cigarSimpString:
        return None # alignment not soft-clipped
    elif cigarSimpString.startswith("M"):
        return 1 # read head
    elif cigarSimpString.endswith("M"):
        return 2 # read tail
    else:
        return None # internal read

def returnMapped(loc):
    """
    function for mapping return mapped list
    """
    # if align.is_mapped:
    #     return True
    # else:
    #     return False
    if loc is not None:
        return True
    else:
        return False

def alignState(locList):
    """
    function for checking read pair alignments state
    """
    match list(map(returnMapped,locList)):
        case [True,False,True,False,False,False]:
            return 0
        case [True,True,True,False,False,False]:
            return 1
        case [True,False,True,True,False,False]:
            return 2
        case [True,True,True,True,False,False]:
            return 3
        case [True,False,False,False,True,False]:
            return 4
        case [True,True,False,False,True,False]:
            return 5
        case [True,False,False,False,True,True]:
            return 6
        case [True,True,False,False,True,True]:
            return 7
        case [True,False,True,False,True,False]:
            return 8
        case [True,True,True,False,True,False]:
            return 9
        case [True,False,True,True,True,False]:
            return 10
        case [True,False,True,False,True,True]:
            return 11
        case [True,True,True,True,True,False]:
            return 12
        case [True,True,True,False,True,True]:
            return 13
        case [True,False,True,True,True,True]:
            return 14
        case [True,True,True,True,True,True]:
            return 15
        case _:
            return None

def compareLoc(loc1,loc2):
    """
    function for comparing mapped locations of different references
    """
    if loc1 is not None and loc2 is not None and loc1 == loc2:
        return 1
    elif loc1 is not None and loc2 is not None and loc1 != loc2:
        return 0
    else:
        return None

def locState(locList):
    """
    function for giving relative location state of different references
    """
    locStateList=[] # [GPrivsTPri,GPrivsTSup,GPrivsPPri,GPrivsPSup,GSupvsTPri,GSupvsTSup,GSupvsPPri,GSupvsPSup,TPrivsPPri,TPrivsPSup,TSupvsPPri,TSupvsPSup]
    locStateList.extend([compareLoc(locList[0],locList[2]),compareLoc(locList[0],locList[3]),compareLoc(locList[0],locList[4]),compareLoc(locList[0],locList[5]),compareLoc(locList[1],locList[2]),compareLoc(locList[1],locList[3]),compareLoc(locList[1],locList[4]),compareLoc(locList[1],locList[5]),compareLoc(locList[2],locList[4]),compareLoc(locList[2],locList[5]),compareLoc(locList[3],locList[4]),compareLoc(locList[3],locList[5])])
    return locStateList

def TPState(alignT,alignP):
    """
    function for checking TE alignment and polyA alignment accordance
    """
    # check type and direction
    Ttype = alignT.reference_name
    Ptype = alignP.get_tag("XP").split(";")[0].split(",")[0] if len(alignP.get_tag("XP").split(";")) == 2 else "NULL"
    if Ttype == Ptype and alignT.flag & 16 == alignP.flag & 16:
        return True
    else:
        return False

def PState(alignP):
    """
    function for checking cutadapt pseudoalignment state
    """
    if len(alignP.get_tag("XP").split(";")) == 2:
        return True
    else:
        return False

# ReadGT functions
def TalignState(Talign,readLen):
    """
    function for checking the mapped portion of TE alignments of discordant reads
    """
    if Talign.reference_length/readLen >= 0.5:
        return True
    else:
        return False


# other functions

def string2tuple(cigarstring):
    """
    function for converting cigarstring to cigartuples
    """
    cigarstringSplit = re.split(r'([M|I|D|N|S|H|P|=|X])',cigarstring)[0:-1] # r for escape
    character2number = {"M":0,"I":1,"D":2,"N":3,"S":4,"H":5,"P":6,"=":7,"X":8}
    cigartuples = [(character2number[cigarstringSplit[i+1]],int(cigarstringSplit[i])) for i in range(len(cigarstringSplit)) if i%2==0]
    return cigartuples

def complseq(lt1): 
    """
    function for generating reverse complement
    """
    lt2 = lt1[::-1]
    lt3=''
    for x in lt2:
        if x == 'A': lt3 = lt3 + 'T'
        elif x == 'T': lt3 = lt3 + 'A'
        elif x == 'C': lt3 = lt3 + 'G'
        elif x == 'G': lt3 = lt3 + 'C'
        else: lt3 = lt3 + x
    return lt3

def checkClipped(cigartuples, cutoff):
    """
    function for checking whether there is one soft/hard-clipped part and whether its length is longer than cutoff 
    """
    if cigartuples is None:
        return None
    for i in cigartuples:
        if i[0] in [4,5] and i[1] >= cutoff:
            return True
    return False

def locateSoftClipped(cigartuples):
    """
    function for generating the location of the longest soft-clipped part
    """
    location_list=[]
    length_list=[]
    start_base=0
    end_base=-1
    for i in cigartuples:
        if i[0] in [2,3,5,6]:
            continue
        start_base = end_base + 1 # closed interval
        end_base = end_base + i[1] # closed interval
        if i[0] == 4:
            location_list.append([start_base,end_base])
            length_list.append(end_base - start_base + 1)
    max_value = max(length_list)
    return location_list[length_list.index(max_value)]

def mergeSame(cigartuples):
    """
    function for merging new adjacent same state bases
    """
    for i in range(len(cigartuples)-1,0,-1):
        if cigartuples[i][0] == cigartuples[i-1][0]:
            cigartuples[i-1] = (cigartuples[i][0],cigartuples[i][1] + cigartuples[i-1][1])
            del cigartuples[i]

def rpLongestSoftClipped(oriCigartuples,reCigartuples,oriFlag,reFlag): # attention to the direction of alignments! 
    """
    function for replacing overlapping alignments with realignment results
    """
    cigarlist = []
    if oriFlag & 16 != reFlag & 16:
        oriCigartuples.reverse()
    for i in range(len(oriCigartuples)): # replace M and deal with I from oriC
        if oriCigartuples[i][0] == 4:
            cigarlist.append(oriCigartuples[i][1])
        elif oriCigartuples[i][0] != 2:
            cigarlist.append(0)
        if oriCigartuples[i][0] == 0:
            oriCigartuples[i] = (4,oriCigartuples[i][1])
        if oriCigartuples[i][0] == 1:
            oriCigartuples[i] = (4,oriCigartuples[i][1])
    oriCigartuples = [x for x in oriCigartuples if x[0] != 2] # remove D from oriC
    softindex = cigarlist.index(max(cigarlist))
    for i in range(len(reCigartuples)): # replace ori S
        oriCigartuples.insert(softindex+i,reCigartuples[i])
    del oriCigartuples[i+softindex+1]
    mergeSame(oriCigartuples)
    return oriCigartuples

def getDirection(Gflag,Tflag):
    """
    function for determining inseriton direction
    """
    if Gflag & 16 == Tflag & 16:
        return "forward"
    else:
        return "reverse"

def reCheckSoftClipped(cigartuples,cutoff):
    """
    function for rechecking soft-clipped part number after replacing genome realignments, read with more than 1 soft-clipped part is removed 
    """
    Snum = 0
    Spass = 0
    for i in cigartuples:
        if i[0] == 4:
            Snum = Snum + 1
            if i[1] >= cutoff:
                Spass = Spass + 1
    if Spass == 1 and Snum == 1:
        return True

def checkClippedNum(cigartuples, cutoff):
    """
    function for checking the number of soft/hard-clipped parts whose its length is not shorted than cutoff 
    """
    clippedNum = 0
    if cigartuples is None:
        return None
    for i in cigartuples:
        if i[0] in [4,5] and i[1] >= cutoff:
            clippedNum += 1
    return clippedNum

class ReadGTP:
    """
    Class which combines Genome and TE primary/supplementary alignment and Cutadapt primary/supplementary pseudoalignment. A neural network is used for classification of read pairs. An alternative arbitrary classification depending on information of overlapping part is also available.
    """
    def __init__(self, alignList, alignLen,TElenDict,deepType="test"):
        self.alignList = alignList
        self.TElenDict = TElenDict
        self.deepType = deepType
        if len(self.alignList) != 12:
            print("length of alignment list is not 12")
            return None
        # generate basic information
        self.alignLen,self.Sequence1,self.Sequence2 = alignLen,alignList[0].get_forward_sequence(),alignList[6].get_forward_sequence()
        self.Phred1 = list(reversed(alignList[0].query_qualities.tolist())) if alignList[0].flag & 16 == 16 else alignList[0].query_qualities.tolist()
        self.Phred2 = list(reversed(alignList[6].query_qualities.tolist())) if alignList[0].flag & 16 == 16 else alignList[6].query_qualities.tolist()
        self.GPriQuality1,self.GPriQuality2,self.GSupQuality1,self.GSupQuality2,self.TPriQuality1,self.TPriQuality2,self.TSupQuality1,self.TSupQuality2 = alignList[0].mapping_quality,alignList[6].mapping_quality,alignList[1].mapping_quality,alignList[7].mapping_quality,alignList[2].mapping_quality,alignList[8].mapping_quality,alignList[3].mapping_quality,alignList[9].mapping_quality
        # indicate read pair type
        if self.deepType == "train":
            if self.alignList[0].get_tag("XG") == "CLIP" and (self.alignList[2].get_tag("XT") == "CLIP" or self.alignList[4].is_mapped) and self.alignList[6].get_tag("XG") == "CLIP" and (self.alignList[8].get_tag("XT") == "CLIP" or self.alignList[10].is_mapped):
                self.alignPairMode = "both"
            elif self.alignList[0].get_tag("XG") == "CLIP" and (self.alignList[2].get_tag("XT") == "CLIP" or self.alignList[4].is_mapped):
                self.alignPairMode = "single:1"
            elif self.alignList[6].get_tag("XG") == "CLIP" and (self.alignList[8].get_tag("XT") == "CLIP" or self.alignList[10].is_mapped):
                self.alignPairMode = "single:2"
            else:
                self.alignPairMode = "none"

    def aligns2array(self):
        """
        method for constructing a 2-dimension list output from 12 alignments of 1 read pair, for both neural network mode and arbitrary mode
        """
        self.alignArray=[] # output initiation
        # channel 1: nucleotide sequence
        # channel 2: sequencing quality
        # channel 3: primary genome alignment mapping quality
        # channel 4: supplementary genome alignment mapping quality
        # channel 5: primary TE alignment mapping quality
        # channel 6: supplementary TE alignment mapping quality
        self.alignArray.extend([convertNuc(self.Sequence1),convertNuc(self.Sequence2),self.Phred1,self.Phred2,[self.GPriQuality1] * self.alignLen,[self.GPriQuality2] * self.alignLen,[self.GSupQuality1] * self.alignLen,[self.GSupQuality2] * self.alignLen,[self.TPriQuality1] * self.alignLen,[self.TPriQuality2] * self.alignLen,[self.TSupQuality1] * self.alignLen,[self.TSupQuality2] * self.alignLen])
        # channel 7-8: primary/supplementary genome alignment and direction
        # 0:not aligned, 1: forward, 2:reverse
        for i in [0,6,1,7]:
            cigarNumList = []
            for j in revTuple(self.alignList[i],self.alignLen):
                if j[0] == 0:
                    if self.alignList[i].flag & 16 == 16:
                        cigarNumList.extend([2]*j[1])
                    else:
                        cigarNumList.extend([1]*j[1])
                elif j[0] in [1,4,5]:
                    cigarNumList.extend([0]*j[1])
            self.alignArray.append(cigarNumList)
        # channel 9-10: primary/supplementary genome alignment chr
        self.alignArray.extend([convertChr(self.alignList[0].reference_name) * self.alignLen,convertChr(self.alignList[6].reference_name) * self.alignLen,convertChr(self.alignList[1].reference_name) * self.alignLen,convertChr(self.alignList[7].reference_name) * self.alignLen])
        # channel 11-12: primary/supplementary TE alignment with location
        for i in [2,8,3,9]:
            cigarNumList = []
            startLoc=self.alignList[i].reference_start + 1 # 0based -> 1based
            for j in revTuple(self.alignList[i],self.alignLen):
                if j[0] == 0:
                    if self.alignList[i].flag & 16 == 16: # reverse direction
                        cigarNumList.extend(reversed(range(startLoc,startLoc+j[1])))
                    else: # forward reverse
                        cigarNumList.extend(range(startLoc,startLoc+j[1]))
                elif j[0] in [1,4,5]:
                    cigarNumList.extend([0]*j[1])
            self.alignArray.append(cigarNumList)
        # channel 13-14 primary/supplementary TE alignment type
        self.alignArray.extend([convertTE(self.alignList[2].reference_name) * self.alignLen,convertTE(self.alignList[8].reference_name) * self.alignLen,convertTE(self.alignList[3].reference_name) * self.alignLen,convertTE(self.alignList[9].reference_name) * self.alignLen])
        # channel 15-16 primary/supplementary cutadapt pseudoalignment with speculated location
        # -1: polyA
        # loc: speculated location
        # 0: other unmapped
        for i in [4,10,5,11]:
            cigarNumList = []
            tailLoc = self.alignList[i].get_tag("XP") # XP information helps discriminating aligned and not aligned reads
            if len(tailLoc.split(";")) == 2:
                TEtype = tailLoc.split(",")[0]
                TEpartLen = int(tailLoc.split(";")[0].split(",")[3]) - int(tailLoc.split(";")[0].split(",")[2]) + 1
                endLoc = TElenDict[TEtype]
                startLoc = endLoc - TEpartLen + 1
                if self.alignList[i].flag & 16 == 16:
                    markTail = 2 # parameter for indicating TE tail part location in read tail
                    locationRan = reversed(range(startLoc,endLoc+1))
                else:
                    markTail = 1 # parameter for indicating TE tail part location in read head
                    locationRan = range(startLoc,endLoc+1)
                for j in revTuple(self.alignList[i],self.alignLen):
                    if j[0] == 0:
                        cigarNumList.extend([-1]*j[1])
                    elif j[0] in [1,5]:
                        cigarNumList.extend([0]*j[1])
                    elif j[0] == 4 and markTail == 1:
                        cigarNumList.extend(locationRan)
                        markTail = 0
                    elif j[0] == 4 and markTail == 0:
                        cigarNumList.extend([0]*j[1])
                    elif j[0] == 4 and markTail == 2:
                        cigarNumList.extend([0]*j[1])
                        markTail = 1
            elif len(tailLoc.split(";")) == 1:
                for j in revTuple(self.alignList[i],self.alignLen):
                    if j[0] == 0:
                        cigarNumList.extend([-1]*j[1])
                    elif j[0] in [1,4,5]:
                        cigarNumList.extend([0]*j[1])
            else:
                cigarNumList.extend([0]*j[1])
            self.alignArray.append(cigarNumList)
        # channel 17-18 primary/supplementary cutadapt pseudoalignment type
        self.alignArray.extend([convertP(self.alignList[4].get_tag("XP")) * self.alignLen,convertP(self.alignList[10].get_tag("XP")) * self.alignLen,convertP(self.alignList[5].get_tag("XP")) * self.alignLen,convertP(self.alignList[11].get_tag("XP")) * self.alignLen])


    def deepMode(self):
        """
        method for executing classification based on the deep learning model
        """

    def arbMode(self):
        """
        method for executing classification based on the arbitrary rules
        """
        self.arbList1,self.arbList2,self.arbPolyAPass1,self.arbPolyAPass2,self.arbDisc,self.arbOutputList = [False,"NULL",[[None],[None],[None]],"NULL","NULL",False],[False,"NULL",[[None],[None],[None]],"NULL","NULL",False],False,False,False,[] # prepare basic information
        self.alignsArbLocParse() # T clipped check is before loc parse while G clipped check is executed during collecting reads
        self.arbPass1,self.arbCode1,self.arbAlignMain1,self.arbMainType1,self.arbLocType1,self.arbTail1 = self.arbList1
        # deal with read2 information
        self.arbPass2,self.arbCode2 = self.arbList2[:2]
        self.arbMainType2,self.arbLocType2,self.arbTail2 = self.arbList2[3:]
        modifiedList = [] # change the ID of read 2 by adding 6
        for i in self.arbList2[2]:
            modifiedList.append(list(map(lambda x:x+6 if x is not None else None,i))) 
        self.arbAlignMain2 = modifiedList
        # other parse after relative location parsing
        self.alignsArbOverlapParse()
        self.alignsArbTEParse()
        self.alignsArbPolyACheck() # P clipped check is after loc parse
        self.alignsArbGmainACheck()
        self.alignsArbAccordParse()
        self.alignsArbDiscCheck()
        self.arbGenerateOutput()

    def alignsArbLocParse(self):
        """
        method for parsing whether read pairs are aligned and the state of relative location of alignments to different references of read pairs, arbitrary mode
        """
        # check genome and TE mapped part relative locations in the read, arbLocList = [GPriLoc1, GSupLoc1, TPriLoc1, TSupLoc1, PPriLoc1, PSupLoc1]
        self.arbLocList1,self.arbLocList2 = list(map(locMapped,self.alignList[0:6])),list(map(locMapped,self.alignList[6:12]))
        # exchange primary alignments corresponding supplenmentary alignments if the latter is more suitable, these alignments will NOT be reverted afterwards
        self.arbChangeList = [] # record the lines exchanged
        for i in [0,2]: # no need to change polyA pseudoalignment because if primary pseudoalignment exists, their alignState will be either 1 or 2 without None, and polyA check will be done in method self.arbPolyACheck
            if self.arbLocList1[i] is None and self.arbLocList1[i+1] is not None:
                self.alignList[i],self.alignList[i+1] = self.alignList[i+1],self.alignList[i]
                self.arbChangeList.append([i,i+1])
                self.arbLocList1 = list(map(locMapped,self.alignList[0:6]))
            if self.arbLocList2[i] is None and self.arbLocList2[i+1] is not None:
                self.alignList[i+6],self.alignList[i+7] = self.alignList[i+7],self.alignList[i+6]
                self.arbChangeList.append([i+6,i+7])
                self.arbLocList2 = list(map(locMapped,self.alignList[6:12]))
        # parse and phrase read pairs
        self.arbAlignState1, self.arbAlignState2, self.arbLocStateList1, self.arbLocStateList2 = alignState(self.arbLocList1), alignState(self.arbLocList2), locState(self.arbLocList1), locState(self.arbLocList2)
        match self.arbAlignState1: # parse read 1
            case 0:
                if self.arbLocStateList1[0] == 0:
                    self.arbList1 = [True,"0:0",[[0],[2],[None]],"CONCLUSIVE","GEN",False]
            case 1:
                if [self.arbLocStateList1[i] for i in [0,4]] == [0,1]:
                    self.arbList1 = [True,"1:0",[[0],[2],[None]],"CONCLUSIVE","GEN",False]
                elif [self.arbLocStateList1[i] for i in [0,4]] == [1,0]:
                    self.arbList1 = [True,"1:1",[[1],[2],[None]],"CONCLUSIVE","GEN",False]
            case 2:
                if [self.arbLocStateList1[i] for i in [0,1]] == [1,0]:
                    self.arbList1 = [True,"2:0",[[0],[3],[None]],"CONCLUSIVE","REFINS",False]
                elif [self.arbLocStateList1[i] for i in [0,1]] == [0,1]:
                    self.arbList1 = [True,"2:1",[[0],[2],[None]],"CONCLUSIVE","REFINS",False]
            case 3:
                if [self.arbLocStateList1[i] for i in [0,1,4,5]] == [1,0,0,1]:
                    self.arbList1 = [True,"3:0",[[0,1],[3,2],[None,None]],"ACCESSORY:2REF","REFINS",False]
                elif [self.arbLocStateList1[i] for i in [0,1,4,5]] == [0,1,1,0]:
                    self.arbList1 = [True,"3:1",[[0,1],[2,3],[None,None]],"ACCESSORY:2REF","REFINS",False]
            case 4:
                if self.arbLocStateList1[2] == 0:
                    if PState(self.alignList[4]):
                        self.arbList1 = [True,"4:0",[[0],[None],[4]],"CONCLUSIVE","GEN",True]
                    else:
                        self.arbList1 = [True,"4:0",[[0],[None],[4]],"ACCESSORY:POLYA","GEN",True]
            case 5:
                if [self.arbLocStateList1[i] for i in [2,6]] == [1,0]:
                    if PState(self.alignList[4]):
                        self.arbList1 = [True,"5:0",[[1],[None],[4]],"CONCLUSIVE","GEN",True]
                    else:
                        self.arbList1 = [True,"5:0",[[1],[None],[4]],"ACCESSORY:POLYA","GEN",True]
                elif [self.arbLocStateList1[i] for i in [2,6]] == [0,1]:
                    if PState(self.alignList[4]):
                        self.arbList1 = [True,"5:1",[[0],[None],[4]],"CONCLUSIVE","GEN",True]
                    else:
                        self.arbList1 = [True,"5:1",[[0],[None],[4]],"ACCESSORY:POLYA","GEN",True]
            case 6:
                if [self.arbLocStateList1[i] for i in [2,3]] == [1,0]:
                    if PState(self.alignList[4]) and PState(self.alignList[5]):
                        self.arbList1 = [True,"6:0",[[0],[None],[5]],"CONCLUSIVE","REFINS",True]
                    elif PState(self.alignList[4]) and not PState(self.alignList[5]):
                        self.arbList1 = [True,"6:0",[[0],[None],[5]],"ACCESSORY:POLYA","REFINS",True]
                elif [self.arbLocStateList1[i] for i in [2,3]] == [0,1]:
                    if PState(self.alignList[4]) and PState(self.alignList[5]):
                        self.arbList1 = [True,"6:1",[[0],[None],[4]],"CONCLUSIVE","REFINS",True]
                    elif not PState(self.alignList[4]) and PState(self.alignList[5]):
                        self.arbList1 = [True,"6:1",[[0],[None],[4]],"ACCESSORY:POLYA","REFINS",True]
            case 7:
                if [self.arbLocStateList1[i] for i in [2,3,6,7]] == [1,0,0,1]:
                    if PState(self.alignList[4]) and PState(self.alignList[5]):
                        self.arbList1 = [True,"7:0",[[0,1],[None,None],[5,4]],"ACCESSORY:2REF","REFINS",True]
                    elif PState(self.alignList[4]) and not PState(self.alignList[5]):
                        self.arbList1 = [True,"7:0",[[0],[None],[5]],"ACCESSORY:POLYA","REFINS",True]
                    elif not PState(self.alignList[4]) and PState(self.alignList[5]):
                        self.arbList1 = [True,"7:0",[[1],[None],[4]],"ACCESSORY:POLYA","REFINS",True]
                elif [self.arbLocStateList1[i] for i in [2,3,6,7]] == [0,1,1,0]:
                    if PState(self.alignList[4]) and PState(self.alignList[5]):
                        self.arbList1 = [True,"7:1",[[0,1],[None,None],[4,5]],"ACCESSORY:2REF","REFINS",True]
                    elif PState(self.alignList[4]) and not PState(self.alignList[5]):
                        self.arbList1 = [True,"7;1",[[1],[None],[5]],"ACCESSORY:POLYA","REFINS",True]
                    elif not PState(self.alignList[4]) and PState(self.alignList[5]):
                        self.arbList1 = [True,"7:1",[[0],[None],[4]],"ACCESSORY:POLYA","REFINS",True]
            case 8:
                if [self.arbLocStateList1[i] for i in [0,2,8]] == [0,0,1] and TPState(self.alignList[2],self.alignList[4]):
                    self.arbList1 = [True,"8:0",[[0],[2],[4]],"CONCLUSIVE","GEN",True]
                elif [self.arbLocStateList1[i] for i in [0,2,8]] == [0,1,0] and PState(self.alignList[4]):
                    self.arbList1 = [True,"8:1",[[0],[2],[None]],"CONCLUSIVE","REFINS",False]
                elif [self.arbLocStateList1[i] for i in [0,2,8]] == [1,0,0]:
                    if PState(self.alignList[4]):
                        self.arbList1 = [True,"8:2",[[0],[None],[4]],"CONCLUSIVE","REFINS",True]
                    else:
                        self.arbList1 = [True,"8:2",[[0],[None],[4]],"ACCESSORY:POLYA","REFINS",True]
            case 9:
                if [self.arbLocStateList1[i] for i in [0,2,4,6,8]] == [1,1,0,0,1] and TPState(self.alignList[2],self.alignList[4]):
                    self.arbList1 = [True,"9:0",[[1],[2],[4]],"CONCLUSIVE","GEN",True]
                elif [self.arbLocStateList1[i] for i in [0,2,4,6,8]] == [1,0,0,1,0]:
                    if PState(self.alignList[4]):
                        self.arbList1 = [True,"9:1",[[0,1],[None,2],[4,None]],"ACCESSORY:2REF","REFINS",False]
                    else:
                        self.arbList1 = [True,"9:1",[[0],[None],[4]],"ACCESSORY:POLYA","REFINS",True]
                elif [self.arbLocStateList1[i] for i in [0,2,4,6,8]] == [0,1,1,0,0]:
                    if PState(self.alignList[4]):
                        self.arbList1 = [True,"9:2",[[0,1],[2,None],[None,4]],"ACCESSORY:2REF","REFINS",False]
                    else:
                        self.arbList1 = [True,"9:2",[[1],[None],[4]],"ACCESSORY:POLYA","REFINS",True]
                elif [self.arbLocStateList1[i] for i in [0,2,4,6,8]] == [0,0,1,1,1] and TPState(self.alignList[2],self.alignList[4]):
                    self.arbList1 = [True,"9:3",[[0],[2],[4]],"CONCLUSIVE","GEN",True]
            case 10:
                if [self.arbLocStateList1[i] for i in [0,1,2,8,10]] == [1,0,1,1,0] and TPState(self.alignList[2],self.alignList[4]):
                    self.arbList1 = [True,"10:0",[[0],[3],[None]],"CONCLUSIVE","REFINS",False]
                elif [self.arbLocStateList1[i] for i in [0,1,2,8,10]] == [0,1,1,0,1] and TPState(self.alignList[3],self.alignList[4]):
                    self.arbList1 = [True,"10:1",[[0],[2],[None]],"CONCLUSIVE","REFINS",False]
                elif [self.arbLocStateList1[i] for i in [0,1,2,8,10]] == [1,0,0,0,1] and TPState(self.alignList[3],self.alignList[4]):
                    self.arbList1 = [True,"10:2",[[0],[3],[4]],"CONCLUSIVE","REFINS",True]
                elif [self.arbLocStateList1[i] for i in [0,1,2,8,10]] == [0,1,0,1,0] and TPState(self.alignList[2],self.alignList[4]):
                    self.arbList1 = [True,"10:3",[[0],[2],[4]],"CONCLUSIVE","REFINS",True]
            case 11:
                if [self.arbLocStateList1[i] for i in [0,2,3,8,9]] == [1,1,0,1,0] and TPState(self.alignList[2],self.alignList[4]):
                    if PState(self.alignList[5]):
                        self.arbList1 = [True,"11:0",[[0],[None],[5]],"CONCLUSIVE","REFINS",True]
                    else:
                        self.arbList1 = [True,"11:0",[[0],[None],[5]],"ACCESSORY:POLYA","REFINS",True]
                elif [self.arbLocStateList1[i] for i in [0,2,3,8,9]] == [1,0,1,0,1] and TPState(self.alignList[2],self.alignList[5]):
                    if PState(self.alignList[4]):
                        self.arbList1 = [True,"11:1",[[0],[None],[4]],"CONCLUSIVE","REFINS",True]
                    else:
                        self.arbList1 = [True,"11:1",[[0],[None],[4]],"ACCESSORY:POLYA","REFINS",True]
                elif [self.arbLocStateList1[i] for i in [0,2,3,8,9]] == [0,1,0,0,1] and TPState(self.alignList[2],self.alignList[5]) and PState(self.alignList[4]):
                    self.arbList1 = [True,"11:2",[[0],[2],[5]],"CONCLUSIVE","REFINS",True]
                elif [self.arbLocStateList1[i] for i in [0,2,3,8,9]] == [0,0,1,1,0] and TPState(self.alignList[2],self.alignList[4]) and PState(self.alignList[5]):
                    self.arbList1 = [True,"11:3",[[0],[2],[4]],"CONCLUSIVE","REFINS",True]
            case 12:
                if [self.arbLocStateList1[i] for i in [0,1,2,4,5,6,8,10]] == [1,0,1,0,1,0,1,0] and TPState(self.alignList[2],self.alignList[4]):
                    self.arbList1 = [True,"12:0",[[0,1],[3,2],[None,4]],"ACCESSORY:2REF","REFINS",True]
                elif [self.arbLocStateList1[i] for i in [0,1,2,4,5,6,8,10]] == [1,0,0,0,1,1,0,1] and TPState(self.alignList[3],self.alignList[4]):
                    self.arbList1 = [True,"12:1",[[0,1],[3,2],[4,None]],"ACCESSORY:2REF","REFINS",True]
                elif [self.arbLocStateList1[i] for i in [0,1,2,4,5,6,8,10]] == [0,1,1,1,0,0,0,1] and TPState(self.alignList[3],self.alignList[4]):
                    self.arbList1 = [True,"12:2",[[0,1],[2,3],[None,4]],"ACCESSORY:2REF","REFINS",True]
                elif [self.arbLocStateList1[i] for i in [0,1,2,4,5,6,8,10]] == [0,1,0,1,0,1,1,0] and TPState(self.alignList[2],self.alignList[4]):
                    self.arbList1 = [True,"12:3",[[0,1],[2,3],[4,None]],"ACCESSORY:2REF","REFINS",True]
            case 13:
                if [self.arbLocStateList1[i] for i in [0,2,3,4,6,7,8,9]] == [1,1,0,0,0,1,1,0] and TPState(self.alignList[2],self.alignList[4]):
                    if PState(self.alignList[5]):
                        self.arbList1 = [True,"13:0",[[0,1],[None,2],[5,4]],"ACCESSORY:2REF","REFINS",True]
                    else:
                        self.arbList1 = [True,"13:0",[[0],[None],[5]],"ACCESSORY:POLYA","REFINS",True]
                elif [self.arbLocStateList1[i] for i in [0,2,3,4,6,7,8,9]] == [0,1,0,1,0,1,0,1] and TPState(self.alignList[2],self.alignList[5]):
                    if PState(self.alignList[4]):
                        self.arbList1 = [True,"13:1",[[0,1],[2,None],[5,4]],"ACCESSORY:2REF","REFINS",True]
                    else:
                        self.arbList1 = [True,"13:1",[[1],[None],[4]],"ACCESSORY:POLYA","REFINS",True]
                elif [self.arbLocStateList1[i] for i in [0,2,3,4,6,7,8,9]] == [1,0,1,0,1,0,0,1] and TPState(self.alignList[2],self.alignList[5]):
                    if PState(self.alignList[4]):
                        self.arbList1 = [True,"13:2",[[0,1],[None,2],[4,5]],"ACCESSORY:2REF","REFINS",True]
                    else:
                        self.arbList1 = [True,"13:2",[[0],[None],[4]],"ACCESSORY:POLYA","REFINS",True]
                elif [self.arbLocStateList1[i] for i in [0,2,3,4,6,7,8,9]] == [0,0,1,1,1,0,1,0] and TPState(self.alignList[2],self.alignList[4]):
                    if PState(self.alignList[5]):
                        self.arbList1 = [True,"13:3",[[0,1],[2,None],[4,5]],"ACCESSORY:2REF","REFINS",True]
                    else:
                        self.arbList1 = [True,"13:3",[[1],[None],[5]],"ACCESSORY:POLYA","REFINS",True]
            case 14:
                if [self.arbLocStateList1[i] for i in [0,1,2,3,8,9,10,11]] == [1,0,1,0,1,0,0,1] and TPState(self.alignList[2],self.alignList[4]) and TPState(self.alignList[3],self.alignList[5]):
                    self.arbList1 = [True,"14:0",[[0],[3],[5]],"CONCLUSIVE","REFINS",True]
                elif [self.arbLocStateList1[i] for i in [0,1,2,3,8,9,10,11]] == [0,1,0,1,1,0,0,1] and TPState(self.alignList[2],self.alignList[4]) and TPState(self.alignList[3],self.alignList[5]):
                    self.arbList1 = [True,"14:1",[[0],[2],[4]],"CONCLUSIVE","REFINS",True]
                elif [self.arbLocStateList1[i] for i in [0,1,2,3,8,9,10,11]] == [1,0,0,1,0,1,1,0] and TPState(self.alignList[2],self.alignList[5]) and TPState(self.alignList[3],self.alignList[4]):
                    self.arbList1 = [True,"14:2",[[0],[3],[4]],"CONCLUSIVE","REFINS",True]
                elif [self.arbLocStateList1[i] for i in [0,1,2,3,8,9,10,11]] == [0,1,1,0,0,1,1,0] and TPState(self.alignList[2],self.alignList[5]) and TPState(self.alignList[3],self.alignList[4]):
                    self.arbList1 = [True,"14:3",[[0],[2],[5]],"CONCLUSIVE","REFINS",True]
            case 15:
                if [self.arbLocStateList1[i] for i in [0,1,2,3,4,5,6,7,8,9,10,11]] == [1,0,1,0,0,1,0,1,1,0,0,1] and TPState(self.alignList[2],self.alignList[4]) and TPState(self.alignList[3],self.alignList[5]):
                    self.arbList1 = [True,"15:0",[[0,1],[3,2],[5,4]],"ACCESSORY:2REF","REFINS",True]
                elif [self.arbLocStateList1[i] for i in [0,1,2,3,4,5,6,7,8,9,10,11]] == [1,0,0,1,0,1,1,0,0,1,1,0] and TPState(self.alignList[2],self.alignList[5]) and TPState(self.alignList[3],self.alignList[4]):
                    self.arbList1 = [True,"15:1",[[0,1],[3,2],[4,5]],"ACCESSORY:2REF","REFINS",True]
                elif [self.arbLocStateList1[i] for i in [0,1,2,3,4,5,6,7,8,9,10,11]] == [0,1,1,0,1,0,0,1,0,1,1,0] and TPState(self.alignList[2],self.alignList[5]) and TPState(self.alignList[3],self.alignList[4]):
                    self.arbList1 = [True,"15:2",[[0,1],[2,3],[5,4]],"ACCESSORY:2REF","REFINS",True]
                elif [self.arbLocStateList1[i] for i in [0,1,2,3,4,5,6,7,8,9,10,11]] == [0,1,0,1,1,0,1,0,1,0,0,1] and TPState(self.alignList[2],self.alignList[4]) and TPState(self.alignList[3],self.alignList[5]):
                    self.arbList1 = [True,"15:3",[[0,1],[2,3],[4,5]],"ACCESSORY:2REF","REFINS",True]
        match self.arbAlignState2: # parse read 2
            case 0:
                if self.arbLocStateList2[0] == 0:
                    self.arbList2 = [True,"0:0",[[0],[2],[None]],"CONCLUSIVE","GEN",False]
            case 1:
                if [self.arbLocStateList2[i] for i in [0,4]] == [0,1]:
                    self.arbList2 = [True,"1:0",[[0],[2],[None]],"CONCLUSIVE","GEN",False]
                elif [self.arbLocStateList2[i] for i in [0,4]] == [1,0]:
                    self.arbList2 = [True,"1:1",[[1],[2],[None]],"CONCLUSIVE","GEN",False]
            case 2:
                if [self.arbLocStateList2[i] for i in [0,1]] == [1,0]:
                    self.arbList2 = [True,"2:0",[[0],[3],[None]],"CONCLUSIVE","REFINS",False]
                elif [self.arbLocStateList2[i] for i in [0,1]] == [0,1]:
                    self.arbList2 = [True,"2:1",[[0],[2],[None]],"CONCLUSIVE","REFINS",False]
            case 3:
                if [self.arbLocStateList2[i] for i in [0,1,4,5]] == [1,0,0,1]:
                    self.arbList2 = [True,"3:0",[[0,1],[3,2],[None,None]],"ACCESSORY:2REF","REFINS",False]
                elif [self.arbLocStateList2[i] for i in [0,1,4,5]] == [0,1,1,0]:
                    self.arbList2 = [True,"3:1",[[0,1],[2,3],[None,None]],"ACCESSORY:2REF","REFINS",False]
            case 4:
                if self.arbLocStateList2[2] == 0:
                    if PState(self.alignList[4+6]):
                        self.arbList2 = [True,"4:0",[[0],[None],[4]],"CONCLUSIVE","GEN",True]
                    else:
                        self.arbList2 = [True,"4:0",[[0],[None],[4]],"ACCESSORY:POLYA","GEN",True]
            case 5:
                if [self.arbLocStateList2[i] for i in [2,6]] == [1,0]:
                    if PState(self.alignList[4+6]):
                        self.arbList2 = [True,"5:0",[[1],[None],[4]],"CONCLUSIVE","GEN",True]
                    else:
                        self.arbList2 = [True,"5:0",[[1],[None],[4]],"ACCESSORY:POLYA","GEN",True]
                elif [self.arbLocStateList2[i] for i in [2,6]] == [0,1]:
                    if PState(self.alignList[4+6]):
                        self.arbList2 = [True,"5:1",[[0],[None],[4]],"CONCLUSIVE","GEN",True]
                    else:
                        self.arbList2 = [True,"5:1",[[0],[None],[4]],"ACCESSORY:POLYA","GEN",True]
            case 6:
                if [self.arbLocStateList2[i] for i in [2,3]] == [1,0]:
                    if PState(self.alignList[4+6]) and PState(self.alignList[5+6]):
                        self.arbList2 = [True,"6:0",[[0],[None],[5]],"CONCLUSIVE","REFINS",True]
                    elif PState(self.alignList[4+6]) and not PState(self.alignList[5+6]):
                        self.arbList2 = [True,"6:0",[[0],[None],[5]],"ACCESSORY:POLYA","REFINS",True]
                elif [self.arbLocStateList2[i] for i in [2,3]] == [0,1]:
                    if PState(self.alignList[4+6]) and PState(self.alignList[5+6]):
                        self.arbList2 = [True,"6:1",[[0],[None],[4]],"CONCLUSIVE","REFINS",True]
                    elif not PState(self.alignList[4+6]) and PState(self.alignList[5+6]):
                        self.arbList2 = [True,"6:1",[[0],[None],[4]],"ACCESSORY:POLYA","REFINS",True]
            case 7:
                if [self.arbLocStateList2[i] for i in [2,3,6,7]] == [1,0,0,1]:
                    if PState(self.alignList[4+6]) and PState(self.alignList[5+6]):
                        self.arbList2 = [True,"7:0",[[0,1],[None,None],[5,4]],"ACCESSORY:2REF","REFINS",True]
                    elif PState(self.alignList[4+6]) and not PState(self.alignList[5+6]):
                        self.arbList2 = [True,"7:0",[[0],[None],[5]],"ACCESSORY:POLYA","REFINS",True]
                    elif not PState(self.alignList[4+6]) and PState(self.alignList[5+6]):
                        self.arbList2 = [True,"7:0",[[1],[None],[4]],"ACCESSORY:POLYA","REFINS",True]
                elif [self.arbLocStateList2[i] for i in [2,3,6,7]] == [0,1,1,0]:
                    if PState(self.alignList[4+6]) and PState(self.alignList[5+6]):
                        self.arbList2 = [True,"7:1",[[0,1],[None,None],[4,5]],"ACCESSORY:2REF","REFINS",True]
                    elif PState(self.alignList[4+6]) and not PState(self.alignList[5+6]):
                        self.arbList2 = [True,"7;1",[[1],[None],[5]],"ACCESSORY:POLYA","REFINS",True]
                    elif not PState(self.alignList[4+6]) and PState(self.alignList[5+6]):
                        self.arbList2 = [True,"7:1",[[0],[None],[4]],"ACCESSORY:POLYA","REFINS",True]
            case 8:
                if [self.arbLocStateList2[i] for i in [0,2,8]] == [0,0,1] and TPState(self.alignList[2+6],self.alignList[4+6]):
                    self.arbList2 = [True,"8:0",[[0],[2],[4]],"CONCLUSIVE","GEN",True]
                elif [self.arbLocStateList2[i] for i in [0,2,8]] == [0,1,0] and PState(self.alignList[4+6]):
                    self.arbList2 = [True,"8:1",[[0],[2],[None]],"CONCLUSIVE","REFINS",False]
                elif [self.arbLocStateList2[i] for i in [0,2,8]] == [1,0,0]:
                    if PState(self.alignList[4+6]):
                        self.arbList2 = [True,"8:2",[[0],[None],[4]],"CONCLUSIVE","REFINS",True]
                    else:
                        self.arbList2 = [True,"8:2",[[0],[None],[4]],"ACCESSORY:POLYA","REFINS",True]
            case 9:
                if [self.arbLocStateList2[i] for i in [0,2,4,6,8]] == [1,1,0,0,1] and TPState(self.alignList[2+6],self.alignList[4+6]):
                    self.arbList2 = [True,"9:0",[[1],[2],[4]],"CONCLUSIVE","GEN",True]
                elif [self.arbLocStateList2[i] for i in [0,2,4,6,8]] == [1,0,0,1,0]:
                    if PState(self.alignList[4+6]):
                        self.arbList2 = [True,"9:1",[[0,1],[None,2],[4,None]],"ACCESSORY:2REF","REFINS",False]
                    else:
                        self.arbList2 = [True,"9:1",[[0],[None],[4]],"ACCESSORY:POLYA","REFINS",True]
                elif [self.arbLocStateList2[i] for i in [0,2,4,6,8]] == [0,1,1,0,0]:
                    if PState(self.alignList[4+6]):
                        self.arbList2 = [True,"9:2",[[0,1],[2,None],[None,4]],"ACCESSORY:2REF","REFINS",False]
                    else:
                        self.arbList2 = [True,"9:2",[[1],[None],[4]],"ACCESSORY:POLYA","REFINS",True]
                elif [self.arbLocStateList2[i] for i in [0,2,4,6,8]] == [0,0,1,1,1] and TPState(self.alignList[2+6],self.alignList[4+6]):
                    self.arbList2 = [True,"9:3",[[0],[2],[4]],"CONCLUSIVE","GEN",True]
            case 10:
                if [self.arbLocStateList2[i] for i in [0,1,2,8,10]] == [1,0,1,1,0] and TPState(self.alignList[2+6],self.alignList[4+6]):
                    self.arbList2 = [True,"10:0",[[0],[3],[None]],"CONCLUSIVE","REFINS",False]
                elif [self.arbLocStateList2[i] for i in [0,1,2,8,10]] == [0,1,1,0,1] and TPState(self.alignList[3+6],self.alignList[4+6]):
                    self.arbList2 = [True,"10:1",[[0],[2],[None]],"CONCLUSIVE","REFINS",False]
                elif [self.arbLocStateList2[i] for i in [0,1,2,8,10]] == [1,0,0,0,1] and TPState(self.alignList[3+6],self.alignList[4+6]):
                    self.arbList2 = [True,"10:2",[[0],[3],[4]],"CONCLUSIVE","REFINS",True]
                elif [self.arbLocStateList2[i] for i in [0,1,2,8,10]] == [0,1,0,1,0] and TPState(self.alignList[2+6],self.alignList[4+6]):
                    self.arbList2 = [True,"10:3",[[0],[2],[4]],"CONCLUSIVE","REFINS",True]
            case 11:
                if [self.arbLocStateList2[i] for i in [0,2,3,8,9]] == [1,1,0,1,0] and TPState(self.alignList[2+6],self.alignList[4+6]):
                    if PState(self.alignList[5+6]):
                        self.arbList2 = [True,"11:0",[[0],[None],[5]],"CONCLUSIVE","REFINS",True]
                    else:
                        self.arbList2 = [True,"11:0",[[0],[None],[5]],"ACCESSORY:POLYA","REFINS",True]
                elif [self.arbLocStateList2[i] for i in [0,2,3,8,9]] == [1,0,1,0,1] and TPState(self.alignList[2+6],self.alignList[5+6]):
                    if PState(self.alignList[4+6]):
                        self.arbList2 = [True,"11:1",[[0],[None],[4]],"CONCLUSIVE","REFINS",True]
                    else:
                        self.arbList2 = [True,"11:1",[[0],[None],[4]],"ACCESSORY:POLYA","REFINS",True]
                elif [self.arbLocStateList2[i] for i in [0,2,3,8,9]] == [0,1,0,0,1] and TPState(self.alignList[2+6],self.alignList[5+6]) and PState(self.alignList[4+6]):
                    self.arbList2 = [True,"11:2",[[0],[2],[5]],"CONCLUSIVE","REFINS",True]
                elif [self.arbLocStateList2[i] for i in [0,2,3,8,9]] == [0,0,1,1,0] and TPState(self.alignList[2+6],self.alignList[4+6]) and PState(self.alignList[5+6]):
                    self.arbList2 = [True,"11:3",[[0],[2],[4]],"CONCLUSIVE","REFINS",True]
            case 12:
                if [self.arbLocStateList2[i] for i in [0,1,2,4,5,6,8,10]] == [1,0,1,0,1,0,1,0] and TPState(self.alignList[2+6],self.alignList[4+6]):
                    self.arbList2 = [True,"12:0",[[0,1],[3,2],[None,4]],"ACCESSORY:2REF","REFINS",True]
                elif [self.arbLocStateList2[i] for i in [0,1,2,4,5,6,8,10]] == [1,0,0,0,1,1,0,1] and TPState(self.alignList[3+6],self.alignList[4+6]):
                    self.arbList2 = [True,"12:1",[[0,1],[3,2],[4,None]],"ACCESSORY:2REF","REFINS",True]
                elif [self.arbLocStateList2[i] for i in [0,1,2,4,5,6,8,10]] == [0,1,1,1,0,0,0,1] and TPState(self.alignList[3+6],self.alignList[4+6]):
                    self.arbList2 = [True,"12:2",[[0,1],[2,3],[None,4]],"ACCESSORY:2REF","REFINS",True]
                elif [self.arbLocStateList2[i] for i in [0,1,2,4,5,6,8,10]] == [0,1,0,1,0,1,1,0] and TPState(self.alignList[2+6],self.alignList[4+6]):
                    self.arbList2 = [True,"12:3",[[0,1],[2,3],[4,None]],"ACCESSORY:2REF","REFINS",True]
            case 13:
                if [self.arbLocStateList2[i] for i in [0,2,3,4,6,7,8,9]] == [1,1,0,0,0,1,1,0] and TPState(self.alignList[2+6],self.alignList[4+6]):
                    if PState(self.alignList[5+6]):
                        self.arbList2 = [True,"13:0",[[0,1],[None,2],[5,4]],"ACCESSORY:2REF","REFINS",True]
                    else:
                        self.arbList2 = [True,"13:0",[[0],[None],[5]],"ACCESSORY:POLYA","REFINS",True]
                elif [self.arbLocStateList2[i] for i in [0,2,3,4,6,7,8,9]] == [0,1,0,1,0,1,0,1] and TPState(self.alignList[2+6],self.alignList[5+6]):
                    if PState(self.alignList[4+6]):
                        self.arbList2 = [True,"13:1",[[0,1],[2,None],[5,4]],"ACCESSORY:2REF","REFINS",True]
                    else:
                        self.arbList2 = [True,"13:1",[[1],[None],[4]],"ACCESSORY:POLYA","REFINS",True]
                elif [self.arbLocStateList2[i] for i in [0,2,3,4,6,7,8,9]] == [1,0,1,0,1,0,0,1] and TPState(self.alignList[2+6],self.alignList[5+6]):
                    if PState(self.alignList[4+6]):
                        self.arbList2 = [True,"13:2",[[0,1],[None,2],[4,5]],"ACCESSORY:2REF","REFINS",True]
                    else:
                        self.arbList2 = [True,"13:2",[[0],[None],[4]],"ACCESSORY:POLYA","REFINS",True]
                elif [self.arbLocStateList2[i] for i in [0,2,3,4,6,7,8,9]] == [0,0,1,1,1,0,1,0] and TPState(self.alignList[2+6],self.alignList[4+6]):
                    if PState(self.alignList[5+6]):
                        self.arbList2 = [True,"13:3",[[0,1],[2,None],[4,5]],"ACCESSORY:2REF","REFINS",True]
                    else:
                        self.arbList2 = [True,"13:3",[[1],[None],[5]],"ACCESSORY:POLYA","REFINS",True]
            case 14:
                if [self.arbLocStateList2[i] for i in [0,1,2,3,8,9,10,11]] == [1,0,1,0,1,0,0,1] and TPState(self.alignList[2+6],self.alignList[4+6]) and TPState(self.alignList[3+6],self.alignList[5+6]):
                    self.arbList2 = [True,"14:0",[[0],[3],[5]],"CONCLUSIVE","REFINS",True]
                elif [self.arbLocStateList2[i] for i in [0,1,2,3,8,9,10,11]] == [0,1,0,1,1,0,0,1] and TPState(self.alignList[2+6],self.alignList[4+6]) and TPState(self.alignList[3+6],self.alignList[5+6]):
                    self.arbList2 = [True,"14:1",[[0],[2],[4]],"CONCLUSIVE","REFINS",True]
                elif [self.arbLocStateList2[i] for i in [0,1,2,3,8,9,10,11]] == [1,0,0,1,0,1,1,0] and TPState(self.alignList[2+6],self.alignList[5+6]) and TPState(self.alignList[3+6],self.alignList[4+6]):
                    self.arbList2 = [True,"14:2",[[0],[3],[4]],"CONCLUSIVE","REFINS",True]
                elif [self.arbLocStateList2[i] for i in [0,1,2,3,8,9,10,11]] == [0,1,1,0,0,1,1,0] and TPState(self.alignList[2+6],self.alignList[5+6]) and TPState(self.alignList[3+6],self.alignList[4+6]):
                    self.arbList2 = [True,"14:3",[[0],[2],[5]],"CONCLUSIVE","REFINS",True]
            case 15:
                if [self.arbLocStateList2[i] for i in [0,1,2,3,4,5,6,7,8,9,10,11]] == [1,0,1,0,0,1,0,1,1,0,0,1] and TPState(self.alignList[2+6],self.alignList[4+6]) and TPState(self.alignList[3+6],self.alignList[5+6]):
                    self.arbList2 = [True,"15:0",[[0,1],[3,2],[5,4]],"ACCESSORY:2REF","REFINS",True]
                elif [self.arbLocStateList2[i] for i in [0,1,2,3,4,5,6,7,8,9,10,11]] == [1,0,0,1,0,1,1,0,0,1,1,0] and TPState(self.alignList[2+6],self.alignList[5+6]) and TPState(self.alignList[3+6],self.alignList[4+6]):
                    self.arbList2 = [True,"15:1",[[0,1],[3,2],[4,5]],"ACCESSORY:2REF","REFINS",True]
                elif [self.arbLocStateList2[i] for i in [0,1,2,3,4,5,6,7,8,9,10,11]] == [0,1,1,0,1,0,0,1,0,1,1,0] and TPState(self.alignList[2+6],self.alignList[5+6]) and TPState(self.alignList[3+6],self.alignList[4+6]):
                    self.arbList2 = [True,"15:2",[[0,1],[2,3],[5,4]],"ACCESSORY:2REF","REFINS",True]
                elif [self.arbLocStateList2[i] for i in [0,1,2,3,4,5,6,7,8,9,10,11]] == [0,1,0,1,1,0,1,0,1,0,0,1] and TPState(self.alignList[2+6],self.alignList[4+6]) and TPState(self.alignList[3+6],self.alignList[5+6]):
                    self.arbList2 = [True,"15:3",[[0,1],[2,3],[4,5]],"ACCESSORY:2REF","REFINS",True]

    def alignsArbOverlapParse(self):
        """
        method for parsing the state of overlapping proption of alignments to different references of read pairs. arbitrary mode
        """
        self.aligns2array()
        self.arrayIDcompList = [[12, 20], [12, 22], [12, 28], [12, 30], [14, 20], [14, 22], [14, 28], [14, 30], [20, 28], [20, 30], [22, 28], [22, 30], [13, 21], [13, 23], [13, 29], [13, 31], [15, 21], [15, 23], [15, 29], [15, 31], [21, 29], [21, 31], [23, 29], [23, 31]]
        # [GPrivsTPri, GPrivsTSup, GPrivsPPri, GPrivsPSup, GSupvsTPri, GSupvsTSup, GSupvsPPri, GSupvsPSup, TPrivsPPri, TPrivsPSup, TSupvsPPri, TSupvsPSup]
        self.arrayListIDcompDict={12:0,14:1,20:2,22:3,28:4,30:5,13:6,15:7,21:8,23:9,29:10,31:11}
        self.overlapList = []
        for calId1,calID2 in self.arrayIDcompList:
            if self.alignList[self.arrayListIDcompDict[calId1]].flag & 4 == 4 or self.alignList[self.arrayListIDcompDict[calID2]].flag & 4 == 4:
                self.overlapList.append([None,None])
            else:
                overlapLen = len([i for i in range(self.alignLen) if self.alignArray[calId1][i] != 0 and self.alignArray[calID2][i] != 0])
                ID1Len = len([i for i in range(self.alignLen) if self.alignArray[calId1][i] != 0])
                ID2Len = len([i for i in range(self.alignLen) if self.alignArray[calID2][i] != 0])
                self.overlapList.append([round(overlapLen/ID1Len,4),round(overlapLen/ID2Len,4)])

    def alignsArbTEParse(self):
        """
        method for checking whether read accords with TE head or tail etc.
        """
        self.headTailList1,self.headTailList2,self.arbTEdirecList1,self.arbTEdirecList2=[],[],[],[]
        if self.arbPass1:
            for i in range(len(self.arbAlignMain1[0])):
                GlineID = self.arbAlignMain1[0][i]
                TlineID = self.arbAlignMain1[1][i]
                Gloc = locMapped(self.alignList[GlineID],"forward")
                if TlineID is not None: # if corresponding T alignment is not None
                    self.arbTEdirecList1.append(getDirection(self.alignList[GlineID].flag,self.alignList[TlineID].flag))
                    if getDirection(self.alignList[GlineID].flag,self.alignList[TlineID].flag) =="reverse" and Gloc == 2 or getDirection(self.alignList[GlineID].flag,self.alignList[TlineID].flag)=="forward" and Gloc == 1: # check insertion direction and relative location of different read parts in main Genome alignment
                        if self.arbTail1:
                            self.headTailList1.append("headTailConflict")
                        elif self.alignList[TlineID].reference_start <= 10: # indicate head by alignment location
                            self.headTailList1.append("head")
                        else: # indicate headTruncation by alignment location
                            self.headTailList1.append("headTruncation")
                    elif getDirection(self.alignList[GlineID].flag,self.alignList[TlineID].flag) =="forward" and Gloc == 2 or getDirection(self.alignList[GlineID].flag,self.alignList[TlineID].flag)=="reverse" and Gloc == 1:
                        if self.TElenDict[self.alignList[TlineID].reference_name] - self.alignList[TlineID].reference_end <= 10:
                            self.headTailList1.append("tail")
                        else:
                            self.headTailList1.append("tailTruncation")
                    else:
                        self.headTailList1.append("unknown")
                else: # only PolyA pseudoalignment exist
                    self.headTailList1.append("tail")
                    PlineID = self.arbAlignMain1[2][i]
                    self.arbTEdirecList1.append(getDirection(self.alignList[GlineID].flag,self.alignList[PlineID].flag))
        if self.arbPass2:
            for i in range(len(self.arbAlignMain2[0])):
                GlineID = self.arbAlignMain2[0][i]
                TlineID = self.arbAlignMain2[1][i]
                Gloc = locMapped(self.alignList[GlineID],"forward")
                if TlineID is not None: # if corresponding T alignment is not None
                    self.arbTEdirecList2.append(getDirection(self.alignList[GlineID].flag,self.alignList[TlineID].flag))
                    if getDirection(self.alignList[GlineID].flag,self.alignList[TlineID].flag) =="reverse" and Gloc == 2 or getDirection(self.alignList[GlineID].flag,self.alignList[TlineID].flag)=="forward" and Gloc == 1: # check insertion direction and relative location of different read parts in main Genome alignment
                        if self.arbTail2:
                            self.headTailList2.append("headTailConflict")
                        elif self.alignList[TlineID].reference_start <= 10: # indicate head by alignment location
                            self.headTailList2.append("head")
                        else: # indicate headTruncation by alignment location
                            self.headTailList2.append("headTruncation")
                    elif getDirection(self.alignList[GlineID].flag,self.alignList[TlineID].flag) =="forward" and Gloc == 2 or getDirection(self.alignList[GlineID].flag,self.alignList[TlineID].flag)=="reverse" and Gloc == 1:
                        if self.TElenDict[self.alignList[TlineID].reference_name] - self.alignList[TlineID].reference_end <= 10:
                            self.headTailList2.append("tail")
                        else:
                            self.headTailList2.append("tailTruncation")
                    else:
                        self.headTailList2.append("unknown")
                else: # only PolyA pseudoalignment exist
                    self.headTailList2.append("tail")
                    PlineID = self.arbAlignMain2[2][i]
                    self.arbTEdirecList2.append(getDirection(self.alignList[GlineID].flag,self.alignList[PlineID].flag))

    def alignsArbPolyACheck(self):
        """
        method for checking whether pseudoalignment is clipped after exchanging locations
        """
        XPlist =[]
        if self.arbPass1:
            for j in self.arbAlignMain1[2]:
                if j is not None:
                    if self.alignList[j].reference_name == "polyA":
                        XPlist = self.alignList[j].get_tag("XP").split(";")
                        XPPlen = 0
                    for i in XPlist:
                        XPPlen += (int(i.split(",")[3]) - int(i.split(",")[2]) + 1)
                    XPGlen = self.alignLen - XPPlen
                    if XPGlen >= 10:
                        self.arbPolyAPass1,self.arbPass1 = True,True
                    else:
                        self.arbPolyAPass1,self.arbPass1 = False,False
        if self.arbPass2:
            for j in self.arbAlignMain2[2]:
                if j is not None:
                    if self.alignList[j].reference_name == "polyA":
                        XPlist = self.alignList[j].get_tag("XP").split(";")
                        XPPlen = 0
                    for i in XPlist:
                        XPPlen += (int(i.split(",")[3]) - int(i.split(",")[2]) + 1)
                    XPGlen = self.alignLen - XPPlen
                    if XPGlen >= 10:
                        self.arbPolyAPass2,self.arbPass2 = True,True
                    else:
                        self.arbPolyAPass2,self.arbPass2 = False,False

    def alignsArbAccordParse(self):
        """
        method for checking whether the state of two reads of the read pair is matched
        """
        if self.arbPass1 and self.arbPass2: # both mode
            if self.arbMainType1 == 'CONCLUSIVE' and self.arbMainType2 == "CONCLUSIVE" or self.arbMainType1 == 'ACCESSORY:POLYA' and self.arbMainType2 == "CONCLUSIVE" or self.arbMainType1 == 'CONCLUSIVE' and self.arbMainType2 == "ACCESSORY:POLYA":
                Galign1,Galign2 = self.alignList[self.arbAlignMain1[0][0]],self.alignList[self.arbAlignMain2[0][0]]
                if not (Galign1.reference_name == Galign2.reference_name and abs(Galign1.reference_start - Galign2.reference_start) <= 1000):
                    self.arbPass1,self.arbPass2,self.alignPairMode = False,False,"both:dif"
            elif self.arbMainType1 == 'CONCLUSIVE' and self.arbMainType2 == "ACCESSORY:2REF":
                for i in range(len(self.arbAlignMain2[0])):
                    Galign1,Galign2 = self.alignList[self.arbAlignMain1[0][0]],self.alignList[self.arbAlignMain2[0][i]]
                    if Galign1.reference_name == Galign2.reference_name and abs(Galign1.reference_start - Galign2.reference_start) <= 1000:
                        self.arbAlignMain2=[[self.arbAlignMain2[0][i]],[self.arbAlignMain2[1][i]],[self.arbAlignMain2[2][i]]]
                        break
                if len(self.arbAlignMain2[0]) == 2:
                    self.arbPass1,self.arbPass2,self.alignPairMode = False,False,"both:dif"
            elif self.arbMainType1 == 'ACCESSORY:2REF' and self.arbMainType2 == "CONCLUSIVE":
                for i in range(len(self.arbAlignMain1[0])):
                    Galign1,Galign2 = self.alignList[self.arbAlignMain1[0][i]],self.alignList[self.arbAlignMain2[0][0]]
                    if Galign1.reference_name == Galign2.reference_name and abs(Galign1.reference_start - Galign2.reference_start) <= 1000:
                        self.arbAlignMain1=[[self.arbAlignMain1[0][i]],[self.arbAlignMain1[1][i]],[self.arbAlignMain1[2][i]]]
                        break
                if len(self.arbAlignMain1[0]) == 2:
                    self.arbPass1,self.arbPass2,self.alignPairMode = False,False,"both:dif"

    def alignsArbGmainACheck(self):
        """
        method for checking whether G main record is consisted of polyA
        """
        self.GmainAratioList1 = []
        self.GmainAratioList2 = []
        if self.arbPass1:
            for i in self.arbAlignMain1[0]:
                GmappedSeq = self.alignList[i].query_alignment_sequence
                GmappedLen = len(GmappedSeq)
                ALen = len([j for j in GmappedSeq if j == "A" or j == "T"])
                self.GmainAratioList1.append(ALen/GmappedLen)
            for i in reversed(range(len(self.GmainAratioList1))):
                if self.GmainAratioList1[i] > 0.8:
                    del self.GmainAratioList1[i]
                    del self.arbAlignMain1[0][i],self.arbAlignMain1[1][i],self.arbAlignMain1[2][i]
            if len(self.arbAlignMain1[0]) == 0:
                self.arbPass1 = False
                self.GmainPolyA1 = True
        if self.arbPass2:
            for i in self.arbAlignMain2[0]:
                GmappedSeq = self.alignList[i].query_alignment_sequence
                GmappedLen = len(GmappedSeq)
                ALen = len([j for j in GmappedSeq if j == "A" or j == "T"])
                self.GmainAratioList2.append(ALen/GmappedLen)
            for i in reversed(range(len(self.GmainAratioList2))):
                if self.GmainAratioList2[i] > 0.8:
                    del self.GmainAratioList2[i]
                    del self.arbAlignMain2[0][i],self.arbAlignMain2[1][i],self.arbAlignMain2[2][i]
            if len(self.arbAlignMain2[0]) == 0:
                self.arbPass2 = False
                self.GmainPolyA2 = True

    def alignsArbDiscCheck(self):
        """
        method for checking alignments failing in clip-check to in discordant read pair standards after exchanging locations
        """
        if not self.arbPass1 and not self.arbPass2:
            if self.alignList[0].is_paired and not self.alignList[0].is_unmapped and not self.alignList[0].mate_is_unmapped and (self.alignList[0].reference_id != self.alignList[0].next_reference_id or abs(self.alignList[0].reference_start - self.alignList[0].next_reference_start) >= 400 + 3 * 100 or (self.alignList[0].is_forward == self.alignList[0].mate_is_forward and self.alignList[0].is_reverse == self.alignList[0].mate_is_reverse)):
                self.arbDisc = True # only output for further parsing in ReadGT

    def arbGenerateOutput(self):
        """
        method for generating the output list
        """
        if self.arbPass1 and self.arbPass2:
            self.arbPairMode = "both"
            for i in [1,2]:
                arbMainType = self.arbMainType1 if i == 1 else self.arbMainType2
                arbAlignMain = self.arbAlignMain1 if i == 1 else self.arbAlignMain2
                headTailList = self.headTailList1 if i == 1 else self.headTailList2
                arbTEdirecList = self.arbTEdirecList1 if i == 1 else self.arbTEdirecList2
                for j in range(len(arbAlignMain[0])):
                    Galign = self.alignList[arbAlignMain[0][j]]
                    TrefName,TlocStart,TlocEnd,Plen = "*","*","*","*"
                    if arbAlignMain[1][j] is not None:
                        Talign = self.alignList[arbAlignMain[1][j]]
                        TrefName,TlocStart,TlocEnd = Talign.reference_name,Talign.reference_start,Talign.reference_end
                    if arbAlignMain[2][j] is not None:
                        Palign = self.alignList[arbAlignMain[2][j]]
                        Plen = Palign.reference_length
                        PalignList = Palign.get_tag("XP").split(";")
                        if arbAlignMain[1][j] is None and len(PalignList) == 2:
                            TrefName = PalignList[0].split(",")[0]
                            TlocStart = TElenDict[PalignList[0].split(",")[0]] - (int(PalignList[0].split(",")[3]) - int(PalignList[0].split(",")[2]) + 1) # keep 0based leftmost
                            TlocEnd = TElenDict[PalignList[0].split(",")[0]] - 1
                    self.arbOutputList.append([Galign.reference_name,Galign.reference_start,Galign.reference_end,Galign.query_name + "/" + Galign.get_tag("XO"),str(j+1) + "/" + str(len(arbAlignMain[0])),arbMainType,TrefName,TlocStart,TlocEnd,headTailList[j],arbTEdirecList[j],Plen])
        elif self.arbPass1 and not self.arbPass2:
            self.arbPairMode = "single:1"
            for j in range(len(self.arbAlignMain1[0])):
                    Galign = self.alignList[self.arbAlignMain1[0][j]]
                    TrefName,TlocStart,TlocEnd,Plen = "*","*","*","*"
                    if self.arbAlignMain1[1][j] is not None:
                        Talign = self.alignList[self.arbAlignMain1[1][j]]
                        TrefName,TlocStart,TlocEnd = Talign.reference_name,Talign.reference_start,Talign.reference_end
                    if self.arbAlignMain1[2][j] is not None:
                        Palign = self.alignList[self.arbAlignMain1[2][j]]
                        Plen = Palign.reference_length
                        PalignList = Palign.get_tag("XP").split(";")
                        if self.arbAlignMain1[1][j] is None and len(PalignList) == 2:
                            TrefName = PalignList[0].split(",")[0]
                            TlocStart = TElenDict[PalignList[0].split(",")[0]] - (int(PalignList[0].split(",")[3]) - int(PalignList[0].split(",")[2]) + 1) # keep 0based leftmost
                            TlocEnd = TElenDict[PalignList[0].split(",")[0]] - 1
                    self.arbOutputList.append([Galign.reference_name,Galign.reference_start,Galign.reference_end,Galign.query_name + "/" + Galign.get_tag("XO"),str(j+1) + "/" + str(len(self.arbAlignMain1[0])),self.arbMainType1,TrefName,TlocStart,TlocEnd,self.headTailList1[j],self.arbTEdirecList1[j],Plen])
        elif not self.arbPass1 and self.arbPass2:
            self.arbPairMode = "single:2"
            for j in range(len(self.arbAlignMain2[0])):
                    Galign = self.alignList[self.arbAlignMain2[0][j]]
                    TrefName,TlocStart,TlocEnd,Plen = "*","*","*","*"
                    if self.arbAlignMain2[1][j] is not None:
                        Talign = self.alignList[self.arbAlignMain2[1][j]]
                        TrefName,TlocStart,TlocEnd = Talign.reference_name,Talign.reference_start,Talign.reference_end
                    if self.arbAlignMain2[2][j] is not None:
                        Palign = self.alignList[self.arbAlignMain2[2][j]]
                        Plen = Palign.reference_length
                        PalignList = Palign.get_tag("XP").split(";")
                        if self.arbAlignMain2[1][j] is None and len(PalignList) == 2:
                            TrefName = PalignList[0].split(",")[0]
                            TlocStart = TElenDict[PalignList[0].split(",")[0]] - (int(PalignList[0].split(",")[3]) - int(PalignList[0].split(",")[2]) + 1) # keep 0based leftmost
                            TlocEnd = TElenDict[PalignList[0].split(",")[0]] - 1
                    self.arbOutputList.append([Galign.reference_name,Galign.reference_start,Galign.reference_end,Galign.query_name + "/" + Galign.get_tag("XO"),str(j+1) + "/" + str(len(self.arbAlignMain2[0])),self.arbMainType2,TrefName,TlocStart,TlocEnd,self.headTailList2[j],self.arbTEdirecList2[j],Plen])
        elif self.arbDisc:
            self.arbPairMode = "disc"
        else:
            self.arbPairMode = "none"

    def alignNeuralNetworkParse(self):
        """
        method for parsing the type of read pair, neural network mode
        """


class ReadGT:
    """
    Class which combines Genome and TE primary/supplementary alignments of discordant read pairs.
    """
    def __init__(self, alignList, alignLen,TElenDict):
        self.alignList = alignList
        self.TElenDict = TElenDict
        if len(self.alignList) != 8:
            print("length of the alignment list is not 8")
            return None
        # generate basic information
        # indicate read pair type
        self.alignPairMode = "disc"

    def arbDiscParse(self):
        """
        method for parsing discordant read pairs
        """
        
    def generateOutput(self):
        """
        method for generating output
        """
