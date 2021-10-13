import os
import glob
import subprocess
import sys

class ClusterAssembler:
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
   def MultipleClusterAssembler(cls,assemblerPath: str, assemblerName: str, inputDir: str, outputDir: str) -> None:
       """
       Run multiple cluster assemblies
       :param assemblerPath: path to assembler
       :param inputDir: path to input dir with all files
       :param outputDir: output dir
       :return: Nothing
       """
       files = cls.__getFilesFromDir(cls,inputDir)
       for file in files:
           cls.SingleClusterAssembler( assemblerPath, assemblerName, file, outputDir)

   @classmethod
   def SingleClusterAssembler(cls,assemblerPath: str, assemblerName: str, inputFile: str, outputDir: str) -> None:
       """
       Run a single cluster assembly
       :param assemblerPath: path to assembler
       :param inputFile: path to inputfile
       :param outputDir: output dir
       :return: None
       """
       # find clustername
       x = str(inputFile.split("\\")[-1])
       i = x.find("cluster")
       j = x[i + 7:].find(".")
       assemblerOutputPath = os.path.join(outputDir, str(x[i + 7:i + 7 + j]))
       import pdb
       pdb.set_trace()
       if assemblerName == 'flye':
           cmd = assemblerPath + " --pacbio-hifi " + x + " --out-dir "+assemblerOutputPath
       elif assemblerName == 'spades':
           cmd = assemblerPath + " --only-assembler -s " + x + " -o "+assemblerOutputPath
       else:
           print('Unknown assembler. Exiting')
           sys.exit()

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
    for multiple clusters in a folder :     clusterAssembler.py multi assemblername path/assembler inputDir/File outputDir
    for a single cluster:                   clusterAssembler.py assemblername path/assembler inputDir/File outputDir
    """

    print(sys.argv)
    if len(sys.argv) > 4:
        assemblerName = str(sys.argv[2])
        assemblerPath = str(sys.argv[3])
        inputDir = str(sys.argv[4])
        outputDir = str(sys.argv[5])
        if str(sys.argv[1]) == "multi":
            ClusterAssembler.MultipleClusterAssembler(assemblerPath, assemblerName, inputDir, outputDir)
        else:
            ClusterAssembler.SingleClusterAssembler(assemblerPath, assemblerName, inputDir, outputDir)
    else:
        print("use cmd like this: clusterAssembler_general.py [nothing/multi] assemblername path/assembler inputDir/File outputDir")