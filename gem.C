#include <cstdlib>
#include <iostream>
#include <vector>
#include <numeric>

#include <TApplication.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <TH1F.h>

#include "Garfield/AvalancheMicroscopic.hh"
#include "Garfield/AvalancheMC.hh"
#include "Garfield/ComponentAnsys123.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/Random.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/ViewField.hh"
#include "Garfield/ViewFEMesh.hh"
#include "Garfield/ViewDrift.hh"
#include "Garfield/ComponentComsol.hh"
#include "Garfield/Random.hh"

using namespace Garfield;

int main(int argc, char * argv[]) {

  TApplication app("app", &argc, argv);

  // Load de field map
  ComponentComsol fm;
  fm.Initialise("MMTHGEM_mesh.mphtxt","MMTHGEM_dielectrics.dat","MMTHGEM_fields.txt","mm");
  fm.EnableMirrorPeriodicityX();
  fm.EnableMirrorPeriodicityY();
  fm.PrintRange();

  // Setup the gas.
  MediumMagboltz gas("cf4");
  gas.SetTemperature(293.15);
  gas.SetPressure(150.);
  gas.Initialise(true); 
  /*
  // Set the Penning transfer efficiency.
  constexpr double rPenning = 0.51;
  constexpr double lambdaPenning = 0.;
  gas.EnablePenningTransfer(rPenning, lambdaPenning, "ar");
  // Load the ion mobilities.
  gas.LoadIonMobility("IonMobility_Ar+_Ar.txt");
  */

  // Associate the gas with the corresponding field map material.
  fm.SetGas(&gas); 
  fm.PrintMaterials();
  // fm.Check();

  // Dimensions of the GEM [cm]
  constexpr double pitch = 0.1;

  // Create the sensor.
  Sensor sensor(&fm);
  sensor.SetArea(-4 * pitch, -4 * pitch, -2 * pitch,
                  4 * pitch,  4 * pitch,  2 * pitch);

  // Define the avalanche components                  
  AvalancheMicroscopic aval(&sensor);
  AvalancheMC drift(&sensor);
  drift.SetDistanceSteps(2.e-4);

  // Plotting options
  ViewField fieldView(&fm);
  ViewFEMesh meshView(&fm);
  ViewDrift driftView;
  constexpr bool plotField = true;
  constexpr bool plotVUV = true;
  constexpr bool plotDrift = true;
  if (plotDrift) {
    // Plot every tenth collision.
    aval.EnablePlotting(&driftView, 10);
    drift.EnablePlotting(&driftView);
  }
  TH1D hEn("hEn", "energy distribution", 1000, 0., 100.);
  aval.EnableElectronEnergyHistogramming(&hEn);
  
  // Count the total number of ions produced the back-flowing ions.
  unsigned int nTotal = 0;
  unsigned int nBF = 0;
  constexpr unsigned int nEvents = 5;
  std::vector<unsigned int> nElectrons;
  std::vector<unsigned int> nVUV;

  for (unsigned int i = 0; i < nEvents; ++i) { 
    std::cout << i << "/" << nEvents << "\n";
    // Randomize the initial position. 
    const double x0 = -0.25 * pitch + 0.5 * RndmUniform() * pitch;
    const double y0 = 0;
    const double z0 = 0.15; 
    const double t0 = 0.;
    const double e0 = 0.1;
    aval.AvalancheElectron(x0, y0, z0, t0, e0, 0., 0., 0.);
    int ne = 0, ni = 0;
    aval.GetAvalancheSize(ne, ni);

    for (const auto& electron : aval.GetElectrons()) {
      const auto& p0 = electron.path[0];
      drift.DriftIon(p0.x, p0.y, p0.z, p0.t);
      ++nTotal;
      const auto& endpoint = drift.GetIons().front().path.back();
      if (endpoint.z > 0.005) ++nBF;
    }
    nElectrons.push_back(ne);
    std::cout << "Number of electrons: " << ne << "\n";

    // Get the number of electron excitations.
    unsigned int nEl = 0;
    unsigned int nIon = 0;
    unsigned int nAtt = 0;
    unsigned int nInel = 0;
    unsigned int nExc = 0;
    unsigned int nSup = 0;
    gas.GetNumberOfElectronCollisions(nEl, nIon, nAtt, nInel, nExc, nSup);
    gas.ResetCollisionCounters();
    nVUV.push_back(nExc + ni); // Multiple ionizations??
  }

  std::cout << "Number of electrons (gas gain): " 
            << std::accumulate(nElectrons.begin(), nElectrons.end(), 0)/nEvents << "\n";
  std::cout << "Fraction of back-flowing ions: " 
            << double(nBF) / double(nTotal) << "\n";

  // Plot the field map.
  if (plotField) {
    TCanvas* cf = new TCanvas("cf", "", 600, 800);
    cf->SetLeftMargin(0.16);

    fieldView.SetCanvas(cf);
    fieldView.SetPlane(0, -1, 0, 0, 0, 0); //XZ plane
    fieldView.SetArea(-2 * pitch, -0.1 * pitch, 2 * pitch, 2 * pitch);
    fieldView.SetVoltageRange(-1500., 0.);
    fieldView.PlotContour();

    meshView.SetCanvas(cf);
    meshView.SetPlane(0, -1, 0, 0, 0, 0);
    meshView.SetArea(-2 * pitch, -0.1 * pitch, 2 * pitch, 2 * pitch); 
    meshView.SetFillMesh(true);
    meshView.SetColor(0, kGray);
    meshView.SetColor(1, kYellow + 3);
    meshView.SetColor(2, kRed);
    meshView.Plot(true);
  }

  // Plot the drift lines.
  if (plotDrift) {
    TCanvas* cd = new TCanvas("cd", "", 600, 800);
    constexpr bool plotMesh = true;
    constexpr bool twod = true;
    if (plotMesh) {
      meshView.SetCanvas(cd);
      meshView.SetComponent(&fm);
      meshView.SetPlane(0, -1, 0, 0, 0, 0);
      if (twod) {
        meshView.SetArea(-2 * pitch, -0.1 * pitch, 2 * pitch, 2 * pitch);
      } else {
        meshView.SetArea(-0.5 * pitch, -0.5 * pitch, -0.02, 0.5 * pitch, 0.5 * pitch, 0.02);
      }
      meshView.SetFillMesh(true);
      meshView.SetColor(0, kGray);
      meshView.SetColor(1, kYellow + 3);
      meshView.SetColor(2, kRed);
      meshView.EnableAxes();
      meshView.SetViewDrift(&driftView);
      const bool outline = twod ? false : true;
      meshView.Plot(twod, outline);
    } else {
      driftView.SetPlane(0, -1, 0, 0, 0, 0);
      driftView.SetArea(-2 * pitch, -0.02, 2 * pitch, 0.02);
      driftView.SetCanvas(cd);
      constexpr bool twod = true;
      driftView.Plot(twod);
    }
  }

  // Plot the electron and VUV distributions.
  if (plotVUV) {
    auto nMinVUV = *std::min_element(nVUV.begin(), nVUV.end());
    auto nMaxVUV = *std::max_element(nVUV.begin(), nVUV.end());
    TH1D hVUV("hVUV", "number of VUV photons", nMaxVUV - nMinVUV, nMinVUV, nMaxVUV);
    hVUV.StatOverflows(true);
    for (const auto& n : nVUV) {hVUV.Fill(n);}
    std::cout << "Mean number of VUV photons: " << hVUV.GetMean() << "\n";
    auto nMinElectrons = *std::min_element(nElectrons.begin(), nElectrons.end());
    auto nMaxElectrons = *std::max_element(nElectrons.begin(), nElectrons.end());
    TH1D hEl("hEn", "number of electrons", nMaxElectrons - nMinElectrons, nMinElectrons, nMaxElectrons);
    hEl.StatOverflows(true);
    for (const auto& n : nElectrons) {hEl.Fill(n);}
    std::cout << "Mean number of electron avalanches: " << hEl.GetMean() << "\n";

    TCanvas* cv = new TCanvas("cv", "", 600, 600);
    //cv->Divide(1, 2);
    //cv->cd(1);
    //cv->SetLogy();
    //hEn.SetTitle("Energy distribution");
    //hEn.Draw();
    //cv->cd(2);
    //cv->SetLogy();
    cv->cd();
    hVUV.SetTitle("Number of VUV photons");
    hVUV.Draw();
    hEl.Draw("same");
    cv->Update();
  }

  app.Run();
}
