import ROOT
import argparse


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Create list of files from grid')
    parser.add_argument('--input', default="/home/fchinu/Run3/Dplus_PbPb_QC/utils/AODsForDownload.txt", help='Input file name')
    parser.add_argument('--output', default="/home/fchinu/Run3/Dplus_PbPb_QC/utils/outfile.txt", help='Output file name')

    args = parser.parse_args()
    inputfilename = args.input
    outputfilename = args.output

    with open(inputfilename,"r") as f:
        InputFiles = f.readlines()
        InputFiles = [x.strip() for x in InputFiles]

    grid = ROOT.TGrid.Connect("alien://")
    fileList = []
    for file in InputFiles:
        gridresult = grid.Ls(file)
        fileinfolist = gridresult.GetFileInfoList()
        for i in range(fileinfolist.GetSize()):
            directoryname = gridresult.GetFileName(i)
            if "aod_collection" not in directoryname and "analysis" not in directoryname and "full_config" not in directoryname:
                fileList.append(file + '/' + directoryname)
    
    with open(outputfilename,"w") as f:
        for file in fileList:
            f.write(file + "\n")


