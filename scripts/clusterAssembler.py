import os
import glob
import subprocess
import sys

class FlyeClusterAssembler:
   def __getFilesFromDir(self,dir: str) -> list:
       """
       Get all .fasta files from a directory
       :param dir: input dir with all cluster files
       :return: list of filepaths
       """
       files = []
       for file in glob.glob(dir + "/*.fasta"):
           files.append(file)
       return files

   @classmethod
   def MultipleClusterAssembler(cls,flyePath: str, inputDir: str, outputDir: str) -> None:
       """
       Run multiple cluster assemblies
       :param flyePath: path to flye
       :param inputDir: path to input dir with all files
       :param outputDir: output dir
       :return: Nothing
       """
       files = cls.__getFilesFromDir(cls,inputDir)
       for file in files:
           cls.SingleClusterAssembler( flyePath, file, outputDir)

   @classmethod
   def SingleClusterAssembler(cls,flyePath: str, inputFile: str, outputDir: str) -> None:
       """
       Run a single cluster assembly
       :param flyePath: path to flye
       :param inputFile: path to inputfile
       :param outputDir: output dir
       :return: None
       """
       # find clustername
       x = str(inputFile.split("\\")[-1])
       i = x.find("cluster")
       j = x[i + 7:].find(".")
       flyeOutputPath = os.path.join(outputDir, str(x[i + 7:i + 7 + j]))

       cmd = flyePath + " --pacbio-hifi " + x + " --out-dir "+flyeOutputPath
       cls.__runCOMMAND(cls,cmd)


   def __runCOMMAND(self, command)->None:
       """
       Executes a shell command
       :param command: command
       :return: None
       """
       print(command)
       sp = subprocess.Popen(command, shell=True)
       sp.wait()


if __name__ == '__main__':
    """
    for multiple clusters in a folder :     clusterAssembler.py multi path/flye inputDir/File outputDir
    for a single cluster:                   clusterAssembler.py path/flye inputDir/File outputDir
    """

    print(sys.argv)
    if len(sys.argv) > 3:
        flyePath = str(sys.argv[2])
        inputDir = str(sys.argv[3])
        outputDir = str(sys.argv[4])
        if str(sys.argv[1]) == "multi":
            FlyeClusterAssembler.MultipleClusterAssembler(flyePath, inputDir, outputDir)
        else:
            FlyeClusterAssembler.SingleClusterAssembler(flyePath, inputDir, outputDir)
    else:
        print("use cmd like this: clusterAssembler.py [nothing/multi] inputDir/File outputDir")