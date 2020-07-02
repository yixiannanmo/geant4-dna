//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// This example is provided by the Geant4-DNA collaboration
// Any report or published results obtained using the Geant4-DNA software
// shall cite the following Geant4-DNA collaboration publication:
// Med. Phys. 37 (2010) 4692-4708
// Delage et al. PDB4DNA: implementation of DNA geometry from the Protein Data
  
//                  Bank (PDB) description for Geant4-DNA Monte-Carlo
//                  simulations (submitted to Comput. Phys. Commun.)
// The Geant4-DNA web site is available at http://geant4-dna.org
//
//
/// \file RunAction.cc
/// \brief Implementation of the RunAction class

#include "RunAction.hh"
#include "RunInitObserver.hh"

#include "Analysis.hh"
#include "G4Run.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction() : G4UserRunAction()
{
//在构造函数中创建analysisManager，建议（但没有必要）在RunAction的构造函数中创建anaysisManager，并在其析构函数中删除它。这保证了在多线程模式下的正确行为。
//10.6版本可以利用字符串来G4AnalysisManager* analysisManager = G4Analysis::ManagerInstance("root");
//AnalysisManager一次只能处理一个基本文件，要想多个相见B5例子以及他的run2.mac

  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetFirstHistoId(1);//设置直方图初始id为1

  // Creating histograms
  //创建直方图可以在构造函数中构造，在EndOfRunAction中填充直方图，见B4例子
  analysisManager->CreateH1("1",
                            "Energy deposit in the target (eV)",
                            1000,0.,1000.);
  analysisManager->CreateH1("2",
                            "Number of SSB",
                            10,0.,10.);
  analysisManager->CreateH1("3",
                            "Number of DSB",
                            10,0.,10.);

  // Open an output file
  //这段代码可以在BeginOfAction中创建文本，在EndOfRunAction中关闭
  G4String fileName = "pdb4dna_output";
  analysisManager->OpenFile(fileName);
  G4String extension = analysisManager->GetFileType();
  fileName = fileName + "." + extension;

  G4cout << "\n----> Histogram file is opened in " << fileName << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->CloseFile();

  delete G4AnalysisManager::Instance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run*)
{
  RunInitManager::Instance()->Initialize();
//可以在这里打开一个文本
//auto analysisManager = G4AnalysisManager::Instance();

  // Open an output file
  //analysisManager->OpenFile("pdb4dna");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run*)
{
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
//在一个run结束可以写入数据，并且关闭文件
  // save histograms 
  //auto analysisManager = G4AnalysisManager::Instance();
  //analysisManager->Write();
  //analysisManager->CloseFile();
  analysisManager->Write();
}

