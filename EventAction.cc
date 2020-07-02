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
/// \file EventAction.cc
/// \brief Implementation of the EventAction class

#include "EventAction.hh"

#include "Analysis.hh"
#include "EventActionMessenger.hh"
#include "G4Event.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

#include <algorithm>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction():G4UserEventAction()
{
  //default parameter values
  //
  fThresEdepForSSB=8.2*eV;
  fThresDistForDSB=10;
  fTotalEnergyDeposit=0;

  //create commands
  //
  fpEventMessenger = new EventActionMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{
  delete fpEventMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction( const G4Event*)
{
  // Initialization of parameters
  //
  fTotalEnergyDeposit=0.;
  fEdepStrand1.clear();
  fEdepStrand2.clear();//清空容器为空，后面的empty详见，std：：vector(网址www.cplusplus.com)
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction( const G4Event*)
{
  // At the end of an event, compute the number of strand breaks
  //定义一个有两个数据的数组，数组名为sb
  G4int sb[2] = {0,0};
  ComputeStrandBreaks(sb);
  // Fill histograms
  //这里见RunAction中定义的直方图（1：能量沉积图；2：SSB的个数；3：DSB的个数）
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
//如果能量沉积有了就是>0就记录，填充进直方图
  if ( fTotalEnergyDeposit>0. )
  {
    analysisManager->FillH1(1,fTotalEnergyDeposit);
  }
//如果数组sb中第一个数据>0，就填充为SSB中
  if ( sb[0]>0 )
  {
    analysisManager->FillH1(2,sb[0]);
  }
//如果第二个数据大于0，就填充在DSB中
  if ( sb[1]>0 )
  {
    analysisManager->FillH1(3,sb[1]);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//下面是一个函数计算锻炼的个数，参数为sb
void EventAction::ComputeStrandBreaks(G4int* sb)
{
  // sb quantities
  //
  G4int ssb1=0;
  G4int ssb2=0;
  G4int dsb=0;

  // nucleotide id and energy deposit for each strand
  G4int nucl1;
  G4int nucl2;
  G4double edep1;
  G4double edep2;

  //Read strand1
  //读取链1
  while ( !fEdepStrand1.empty() )
  {
//这里利用stl中的map容器存取了链1的能量沉积，对于map里利用迭代器fEdepStrand1.begin()，读取了他的key和value(first就是key这里是用的dna核苷酸之间的距离，second就是value指能量沉积)
    nucl1 = fEdepStrand1.begin()->first;
    edep1 = fEdepStrand1.begin()->second;
    fEdepStrand1.erase( fEdepStrand1.begin() );//容器的迭代器删除erase函数(表示删除这个元素，并指向下一个元素)

    // SSB in strand1
    //
    if ( edep1 >= fThresEdepForSSB/eV )
    {
      ssb1++;
    }

    // Look at strand2
    //
    if ( !fEdepStrand2.empty() )
    {
      do
      {

        nucl2 = fEdepStrand2.begin()->first;
        edep2 = fEdepStrand2.begin()->second;
        if ( edep2 >= fThresEdepForSSB/eV )
        {
          ssb2++;
        }
        fEdepStrand2.erase( fEdepStrand2.begin() );
      } while ( ((nucl1-nucl2)>fThresDistForDSB) && (!fEdepStrand2.empty()) );

      // no dsb
      //如果两条链的距离大于给定的阈值，第二条链上的能量沉积就是第二条链上的第一个数据，这个能量沉积大于给定的ssb的阈值，ssb2就减少一个
      if ( nucl2-nucl1 > fThresDistForDSB )
      {
        fEdepStrand2[nucl2]=edep2;
        if ( edep2 >= fThresEdepForSSB/eV )
        {
          ssb2--;
        }
      }

      // one dsb
      //如果两条链的距离小于给定的阈值，并且两条链能量都大于给定的阈值，就ssb都减少1，dsb+1
      if ( std::abs(nucl2-nucl1) <= fThresDistForDSB )
      {
        if ( ( edep2 >= fThresEdepForSSB/eV ) &&
            ( edep1 >= fThresEdepForSSB/eV ) )
        {
          ssb1--;
          ssb2--;
          dsb++;
        }
      }
    }
  }

  // End with not processed data
  //
  while ( !fEdepStrand1.empty() )
  {
    nucl1 = fEdepStrand1.begin()->first;
    edep1 = fEdepStrand1.begin()->second;
    if ( edep1 >= fThresEdepForSSB/eV )
    {
      ssb1++;
    }
    fEdepStrand1.erase( fEdepStrand1.begin() );
  }

  while ( !fEdepStrand2.empty() )
  {
    nucl2 = fEdepStrand2.begin()->first;
    edep2 = fEdepStrand2.begin()->second;
    if ( edep2 >= fThresEdepForSSB/eV )
    {
      ssb2++;
    }
    fEdepStrand2.erase( fEdepStrand2.begin() );
  }

//最后sb数组保存的数据
  sb[0]=ssb1+ssb2;
  sb[1]=dsb;
}
