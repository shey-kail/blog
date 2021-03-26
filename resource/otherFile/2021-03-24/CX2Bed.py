#######
#######CX2Bed
#######

# give me two parameters，inputFilePath outputFilePath

import sys
from collections import OrderedDict



############################################
################# CONFIG ###################
############################################
STEP=50
BIN_SIZE=100
############################################
################# CONFIG ###################
############################################





class Bin :
    def __init__(self,size:int=BIN_SIZE,step:int=STEP,begin:int=0):
        self.beginSite=begin
        self.step=step
        self.endSite=self.beginSite+self.step
    def move(self):
        self.beginSite+=self.step
        self.endSite=self.beginSite+self.step
    def reset(self,size:int=BIN_SIZE,step:int=STEP,begin:int=0):
        self.beginSite=begin
        self.step=step
        self.endSite=self.beginSite+self.step

    def isInBin(self,site):
        if site>self.beginSite and site<=self.endSite :
            return "InBin"
        elif site>self.endSite:
            return "BehindBin"
        else :
            return "ForewardBin"




#### This class is used to create collection of all chromosomes' sites, C_counts and T_counts, and to simplify search those imformation
#### information can be searched by keys, and keys follow "chromosome,site" format, such as "chr1,1002"
class ChrDict:
    def __init__(self,CX_data:list):
        self.ChrDict=OrderedDict()
        self.lastKey="NULL"

        for i in range(len(CX_data)):
            Chr=CX_data[i].split("\t")[0]
            site=int(CX_data[i].split("\t")[1])
            C=int(CX_data[i].split("\t")[3])
            T=int(CX_data[i].split("\t")[4])
            context=CX_data[i].split("\t")[5]
            
            key=Chr+","+str(site)


            self.ChrDict[key]={
                    "site":site,
                    "C":C,
                    "T":T,
                    "context":context
            }

    # return a list that restores all keys
    def allKeys(self):
        return list(self.ChrDict.keys())

    # return a site's site number in a chromosome 
    def Site(self,key):
        return self.ChrDict[key]["site"]

    # return a site's C_count in a chromosome 
    def Ccount(self,key):
        return self.ChrDict[key]["C"]

    # return a site's T_count in a chromosome 
    def Tcount(self,key):
        return self.ChrDict[key]["T"]

    # return a site's context in a chromosome 
    def Context(self,key):
        return self.ChrDict[key]["context"]

    # return the key's chromosome
    def Chrom(self,key):
        return key.split(",")[0]

    # judge whether the chromosome is changed
    def chrChange(self,key):
        nowKey=key

        if self.lastKey=="NULL":
            self.lastKey=nowKey
            return False
        elif nowKey.split(",")[0]!=self.lastKey.split(",")[0]:
            self.lastKey=nowKey
            return True
        else:
            self.lastKey=nowKey
            return False


class resultNotebook :
    def __init__(self):
        self.resultList=[]


    def addResult(self,Chr,BIN:Bin,Counter):

        binResult=[Chr,
                str(BIN.beginSite), str(BIN.endSite), 
                str(Counter.C_CG), 
                str(Counter.T_CG), 
                str(Counter.C_CHG), 
                str(Counter.T_CHG), 
                str(Counter.C_CHH), 
                str(Counter.T_CHH)
        ]

        self.resultList.append(binResult)


    def writeToFile(self,outputFile):
        for i in range(len(self.resultList)):
            line="\t".join(self.resultList[i])+"\n"
            outputFile.write(line)


# this class is designed to simplify counting C and T
class counter :
    def __init__(self):
        self.C_CG=0
        self.T_CG=0
        self.C_CHG=0
        self.T_CHG=0
        self.C_CHH=0
        self.T_CHH=0

    def add(self,context:str,C:int,T:int):
        if (context=="CG"):
            self.C_CG+=C
            self.T_CG+=T
        elif (ChrDB.Context(key)=="CHG"):
            self.C_CHG+=C
            self.T_CHG+=T
        elif (ChrDB.Context(key)=="CHH"):
            self.C_CHH+=C
            self.T_CHH+=T

    def empty(self):
        self.C_CG=0
        self.T_CG=0
        self.C_CHG=0
        self.T_CHG=0
        self.C_CHH=0
        self.T_CHH=0






if __name__=='__main__':
    CX_File=open(sys.argv[1])
    outputFile=open(sys.argv[2],"w")


    CX_line=CX_File.readlines()

    result=resultNotebook()

    #ChrDB restores information of CX_report, and organize it into orderdDict format. the information should be searched by keys, and keys can be got by ChrDict.allKeys() function
    ChrDB=ChrDict(CX_data=CX_line)

    #theBin is an object that represents alignment regions
    theBin=Bin()

    #Counter is an object that is used to count C and T according to context 
    Counter=counter()


    #ergodic ChrDB to count
    for key in ChrDB.allKeys():

        # if chromosome changed, save the C's and T's counting result to outputResult
        if (ChrDB.chrChange(key)==True):
            result.addResult(ChrDB.Chrom(key),theBin,Counter)
        # and reset the bin and recount C and T
            theBin.reset()
            Counter.empty()

        # if the site is out of Bin, move the Bin until the site is in it
        while (theBin.isInBin(ChrDB.Site(key))!="InBin"):
            result.addResult(ChrDB.Chrom(key),theBin,Counter)
            theBin.move()
            Counter.empty()


        # when the site is in the Bin, start to count the C and T
        ##### the count of C and T were added according to the context
        Counter.add(context=ChrDB.Context(key),
                C=ChrDB.Ccount(key),
                T=ChrDB.Tcount(key))



        # To ensure the information in last region of last chromosome not being lost, C's and T's count must be saved when it turns to last ChrDB's key.
        if (key==ChrDB.allKeys()[len(ChrDB.allKeys())-1]):
            result.addResult(ChrDB.Chrom(key),theBin,Counter)
           


    result.writeToFile(outputFile)
    outputFile.close()
    CX_File.close()


            


