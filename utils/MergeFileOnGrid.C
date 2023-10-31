#include <fstream>
#include <string>
#include <iostream>
#include "TFileMerger.h"
#include "TGrid.h"
#include "TString.h"
#include "TTree.h"
#include "TFile.h"
#include "TList.h"
#include "TCollection.h"

void MergeFileOnGrid()
{
    //std::ifstream file("/home/fchinu/Run3/Ds_pp_13TeV/utils/fileList.txt");
    std::ifstream file("/home/fchinu/Run3/Ds_pp_13TeV/utils/fileList130748.txt");
    std::string str; 
    int nFilesMax = 60;

    // connect to grid
    if(!TGrid::Connect("alien://")) {
       printf("no grid connection available... Exiting!");
       return;}
    TFileMerger *mergeAOD = new TFileMerger();
    TFileMerger *mergeAnalysisResults = new TFileMerger();
    int count = 0;
    while (std::getline(file, str) && count < nFilesMax)
    {
        count++;
        std::cout << "Adding file: " << Form("alien://%s/AO2D.root",str.c_str()) << std::endl;
        mergeAOD->AddFile(Form("alien://%s/AO2D.root",str.c_str()));
        mergeAnalysisResults->AddFile(Form("alien://%s/AnalysisResults.root",str.c_str()));
    }
    //mergeAOD->OutputFile(Form("LHC22b1a_Train128762_AO2D.root"), "RECREATE");
    mergeAOD->OutputFile(Form("/home/fchinu/Run3/Ds_pp_13TeV/Datasets/Ds_pp_run3_ml/LHC22o_Train130748_AO2D.root"), "RECREATE");
    mergeAOD->Merge();
    mergeAnalysisResults->OutputFile(Form("/home/fchinu/Run3/Ds_pp_13TeV/Datasets/Ds_pp_run3_ml/LHC22o_Train130748_AnalysisResults.root"), "RECREATE");
    mergeAnalysisResults->Merge();
}

//void ReduceSize() {
//    std::vector<std::string> BranchesToRemove =  {"fIndexBCs", "fNumContrib", "fPosX", "fPosY", "fPosZ","fXSecondaryVertex", "fYSecondaryVertex", //"fZSecondaryVertex", "fErrorDecayLength","fErrorDecayLengthXY", "fRSecondaryVertex", "fImpactParameterNormalised0", "fPtProng0", "fPProng0", //"fImpactParameterNormalised1", "fPtProng1", "fPProng1", "fImpactParameterNormalised2", "fPtProng2", "fPProng2", "fPxProng0", "fPyProng0", "fPzProng0", //"fPxProng1", "fPyProng1", "fPzProng1", "fPxProng2", "fPyProng2", "fPzProng2", "fImpactParameter0", "fImpactParameter1", "fImpactParameter2", //"fErrorImpactParameter0", "fErrorImpactParameter1", "fErrorImpactParameter2", "fP", "fCt", "fEta", "fE"};
//
//    std::string inputfilename = "/home/fchinu/Run3/Ds_pp_13TeV/Datasets/Ds_pp_run3_ml/LHC22o_Train130748_AO2D.root";
//    std::string outputfilename = "/home/fchinu/Run3/Dplus_PbPb_QC/Dataset/Dplus_PbPb_run3_ml/Light_LHC22o_Train130748_AO2D.root";
//
//    TFile *input = TFile::Open(inputfilename.c_str(), "READ");
//    if (!input || input->IsZombie()) {
//        std::cerr << "Error: Could not open input file " << inputfilename << std::endl;
//        return;
//    }
//
//    TFile *output = TFile::Open(outputfilename.c_str(), "RECREATE");
//    if (!output || output->IsZombie()) {
//        std::cerr << "Error: Could not create output file " << outputfilename << std::endl;
//        input->Close();
//        delete input;
//        return;
//    }
//
//    TIter next(input->GetListOfKeys());
//    TKey *key;
//
//    std::vector<std::string> folders;
//    while ((key = (TKey*)next())) {
//        if (key->IsFolder() && std::string(key->GetName()) != "parentFiles") {
//            folders.push_back(std::string(key->GetName()));
//        }
//    }
//
//    TTree *inputtree;
//    TTree *outputtree;
//
//    // Loop through all TDirectoryFile objects in the input file
//    int count = 0;
//
//    for (auto folder : folders) {
//        count++;
//        std::cout << "Processing folder " << count << ": " << folder << std::endl;
//        // Get the TDirectoryFile
//        TDirectoryFile *directory = dynamic_cast<TDirectoryFile*>(input->Get(folder.c_str()));
//
//        // Check if the TDirectoryFile contains a TTree named "O2hfcanddpfull"
//        inputtree = dynamic_cast<TTree*>(directory->Get("O2hfcanddpfull"));
//        if (inputtree) {
//            std::cout << "Processing TTree " << inputtree->GetName() << std::endl;
//
//            // Create a new TTree with selected branches
//            outputtree = inputtree->CloneTree(0);
//            for (auto branch : BranchesToRemove) {
//                outputtree->SetBranchStatus(branch.c_str(), 0);
//            }
//
//            // Copy the selected branches to the new TTree
//            outputtree->CopyEntries(inputtree);
//
//            // Write the modified TTree to the output file
//            outputtree->SetDirectory(output);
//            outputtree->Write();
//
//            // Cleanup
//            delete outputtree;
//        }
//    }
//
//    // Close files and cleanup
//    output->Close();
//    input->Close();
//
//    delete output;
//    delete input;
//}

