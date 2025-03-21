#include "Pythia8/Pythia.h" // access to Pythia objects.
#include "Pythia8/HIInfo.h" // access to Pythia objects.
using namespace Pythia8;
// allow simplified notation.

int main(int argc, char *argv[])
{

  char *configfile = argv[1];
  char *outfile = argv[2];

  // --- Initialization ---
  Pythia pythia;
  // Define Pythia object.
  Event &event = pythia.event; // quick access to current event.

  // Read settings from file
  pythia.readFile(configfile);
  pythia.settings.parm("HISigmaDiffractive:PomFluxEpsilon", -0.04, true);

  // Define histograms, external links,
  // local variables etc. here. E.g.
  pythia.init();

  int nD0 = 0;
  int nCharm = 0;
  int nPions = 0;
  int nKaons = 0;
  int evGen = 0;

  FILE *outtxt = fopen(outfile, "w");

  int charmHadPDG[7] = {411, 421, 431, 4122, 4132, 4232, 4332};
  // Initialize
  // --- The event loop ---
  int nAbort = pythia.mode("Main:timesAllowErrors");
  int nEvent = pythia.mode("Main:numberOfEvents");

  int iAbort = 0;
  for (int iEvent = 0; iEvent < nEvent; iEvent++)
  {
    // Generate next event;
    // Produce the next event, returns true on success.
    if (!pythia.next())
    {
      if (++iAbort < nAbort)
        continue;
      cout << " Event generation aborted prematurely, owing to error!\n";
    }
    //    pythia.event.list();
    fprintf(outtxt, "Event %d\n", iEvent);
    int nCharged = 0;
    int nFinal = 0;
    int nChargedy05 = 0;
    int nChargedFV0 = 0;
    int nPhysPrimy05 = 0;
    int nPhysPrimFV0 = 0;
    printf("Event %d size %d\n", iEvent, pythia.event.size());
    for (int i = 0; i < pythia.event.size(); ++i)
    {
      int abspdg = abs(pythia.event[i].id());
      double eta = event[i].eta();
      if (pythia.event[i].isFinal())
      {
        ++nFinal;
        if (pythia.event[i].isCharged())
        {
          ++nCharged;
          if (abs(eta) < 0.8)
            ++nChargedy05;
          if (eta > 3.5 && eta < 4.9)
            ++nChargedFV0;
          if (eta > -3.3 && eta < -2.1)
            ++nChargedFV0;
        }
        if (abspdg == 211)
          ++nPions;
        else if (abspdg == 321)
          ++nKaons;
      }
      if (abspdg == 11 || abspdg == 13 || abspdg == 211 || abspdg == 321 || abspdg == 2212)
      {
        double rad = sqrt(pythia.event[i].xProd() * pythia.event[i].xProd() + pythia.event[i].yProd() * pythia.event[i].yProd());
        if (rad < 10.)
        {
          if (abs(eta) < 0.8)
            ++nPhysPrimy05;
          if (eta > 3.5 && eta < 4.9)
            ++nPhysPrimFV0;
          if (eta > -3.3 && eta < -2.1)
            ++nPhysPrimFV0;
        }
      }
      if (abspdg == 421)
        ++nD0;
      else if (abspdg == 4)
        ++nCharm;

      for (int j = 0; j < 7; j++)
      {
        if (abspdg == charmHadPDG[j])
        {
          int moth1 = pythia.event[i].mother1();
          int moth2 = pythia.event[i].mother2();
          int idMoth1 = pythia.event[moth1].id();
          int idMoth2 = pythia.event[moth2].id();
          int moth = moth1;
          if (abs(idMoth2) > abs(idMoth1))
            moth = moth2;
          bool fromB = false;
          while (moth > 0)
          {
            int idMoth = pythia.event[moth].id();
            if (abs(idMoth) == 5 || (abs(idMoth) > 500 && abs(idMoth) < 600) || (abs(idMoth) > 5000 && abs(idMoth) < 6000))
            {
              fromB = true;
              break;
            }
            moth1 = pythia.event[moth].mother1();
            moth2 = pythia.event[moth].mother2();
            idMoth1 = pythia.event[moth1].id();
            idMoth2 = pythia.event[moth2].id();
            moth = moth1;
            if (abs(idMoth2) > abs(idMoth1))
              moth = moth2;
          }
          fprintf(outtxt, "%d %d %f %f %f %f\n", pythia.event[i].id(), fromB, pythia.event[i].px(), pythia.event[i].py(), pythia.event[i].pz(), pythia.event[i].e());
        }
      }
    }
    fprintf(outtxt, "nCharged %d %d %d %d %d\n", nCharged, nChargedy05, nChargedFV0, nPhysPrimy05, nPhysPrimFV0);
    //    printf("nCharged %d %d %d %d %d\n",nCharged,nChargedy05,nChargedFV0,nPhysPrimy05,nPhysPrimFV0);
    ++evGen;
  } // End event loop.
  // --- Calculate final statistics ---
  pythia.stat();
  double xsec = pythia.info.sigmaGen();
  double nacc = pythia.info.nAccepted();
  double xsecHI = (pythia.info.hiInfo)->sigmaTot();
  double naccHI = (pythia.info.hiInfo)->nAccepted();
  fprintf(outtxt, "---------------------------\n");
  fprintf(outtxt, "N events = %d  Cross sectionHI = %f nAcceptedHI = %f Cross section = %f nAccepted = %f\n", evGen, xsecHI, naccHI, xsec, nacc);
  fclose(outtxt);

  cout << "number of pions = " << nPions << endl;
  cout << "number of kaons = " << nKaons << endl;
  cout << "number of D0 = " << nD0 << endl;
  cout << "number of c quarks = " << nCharm << endl;
  cout << "NAcc = " << nacc << endl;
  cout << "xSec = " << xsec << endl;

  return 0;
}
