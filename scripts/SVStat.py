import math
import mmap
import matplotlib.pyplot as plt
import os
import pandas as pd


class SVStat:
    def __getABSPath(self,path)->str:
        """
        Converts a relative path to an absolute path
        :param path: input path to a file or folder
        :return: absolute path to the file or folder
        """
        return os.path.abspath(os.path.expanduser(path))

    @classmethod
    def extractAF(cls, filepath:str):
        """
        Extracts all allele frequencies form a given SVTYPE in .vcf file
        :param filepath: path to .vcf
        :return: dictionary with SVTYPE->AF:Count
        """

        filepath = cls.__getABSPath(cls, filepath)
        mutationDict = dict()
        # data extraction
        if not os.path.exists(filepath):
            print("File does not exist")
        else:
            with open(filepath, "r") as fp:
                mm = mmap.mmap(fp.fileno(), length=0, access=mmap.ACCESS_READ)
                for line in iter(mm.readline, b""):
                    if str(line).__contains__("AF="):
                        line = line.decode("utf-8")
                        i = line.find("SVTYPE=")
                        j = line[i:].find(";")
                        svtype = str(line[i + 7:i + j])
                        i = line.find("AF=")
                        j = line[i:].find("\t")
                        af = str(line[i + 3:i + j])
                        if svtype not in mutationDict.keys():
                            mutationDict[str(svtype)] = {}
                        if float(af) not in mutationDict[str(svtype)].keys():
                            mutationDict[str(svtype)][float(af)] = 0
                        mutationDict[str(svtype)][float(af)] += 1

        return mutationDict
    @classmethod
    def getAFDistribution(cls, filepath:str, outputpath:str):
        """
        Generates a visualisation of the allele frequency in vcf file.
        :param filepath: path to vcf file
        :param outputpath: output path of the visualisation (.png or .pdf)
        :return: figure, pandas dataframe
        """
        outputpath = cls.__getABSPath(cls, outputpath)
        if not os.path.exists(filepath):
            print("path not valid")
            outputpath = os.path.join(os.getcwd()+"out.pdf")
            print(outputpath)

        # plotting
        df = pd.DataFrame.from_dict(cls.extractAF(filepath)).sort_index().fillna(0)

        cols = 3
        rows = math.ceil(len(df.columns) / cols)

        plt.figure(figsize=(20, 15))
        fig = plt.figure(1)

        for i, (k, v) in enumerate(df.items(), 1):
            ax = fig.add_subplot(rows, cols, i, title=k)
            ax.plot(df[k])
            ax.autoscale()
            ax.set_xlabel("Allele frequencies")
            ax.set_ylabel("Count")

        plt.savefig(outputpath)
        return plt, df



if __name__ == '__main__':
    fig, df = SVStat.getAFDistribution(
        filepath=r"~/HG002.HiFi.GRCh37.SVLEN50.RE10.vcf",
        outputpath=r"~/out.pdf"
    )
    #x.show()
    #print(df)