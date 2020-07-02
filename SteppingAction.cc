#include "SteppingAction.hh"

#include "DetectorConstruction.hh"
#include "G4EventManager.hh"
#include "EventAction.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4RunManager.hh"
#include "G4SteppingManager.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction()
:G4UserSteppingAction(),RunInitObserver(),fpEventAction(0),fpDetector(0)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::Initialize()
{
  fpEventAction = (EventAction*) G4EventManager::GetEventManager()->
      GetUserEventAction();
  fpDetector = (DetectorConstruction*)G4RunManager::GetRunManager()->
      GetUserDetectorConstruction();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* theStep)
{
  if(theStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()!=
      "Transportation")
  {
    // Get position and edep of current step
    //获取当前的step的位置和能量沉积
    G4double x = theStep->GetPreStepPoint()->GetPosition().x()/nanometer;
    G4double y = theStep->GetPreStepPoint()->GetPosition().y()/nanometer;
    G4double z = theStep->GetPreStepPoint()->GetPosition().z()/nanometer;
    G4double edepStep = theStep->GetTotalEnergyDeposit()/eV;

    G4LogicalVolume* targetVolume =
        G4LogicalVolumeStore::GetInstance()->GetVolume("BoundingLV");
    G4LogicalVolume* theVolume =
        theStep->GetPreStepPoint()->GetPhysicalVolume()->GetLogicalVolume();
//沉积能量>0，并且逻辑体是靶区，就认为这个event沉积能量
    if ((edepStep > 0.) && (theVolume==targetVolume))
    {
      // Add edep to this event
      //
      fpEventAction->AddEdepEvent(edepStep);
//判断探测器指向了哪
      if (fpDetector->GetBarycenterList()==NULL)
      {
        G4cout << "Barycenter list is null!!!" << G4endl;
      }
      else
      {
        CheckAndProcessDNAHit(x,y,z,edepStep);
      }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool SteppingAction::CheckAndProcessDNAHit(G4double x,G4double y, G4double z,
    G4double edepStep)
{
  int numStrand=0;
  int numNucl=0;
  int intResidue=-1; // 0 for Phospat, 1 for Sugar, 2 for Base
  unsigned short int hit = (fpDetector->GetPDBlib()).ComputeMatchEdepDNA(
      fpDetector->GetBarycenterList(),
      fpDetector->GetMoleculeList(),
      x*10., y*10., z*10.,// x10 => angstrom<->nm
      numStrand, numNucl, intResidue);

  if (hit==1)
  {
    if ((intResidue==0)||(intResidue==1)) //Edep in Phosphate or Sugar
    {
      fpEventAction->AddEdepToNucleotide(numStrand,numNucl,edepStep);
      return true;
    }
    else
    {
      return false;
    }
  }
  else
  {
    return false;
  }
}
