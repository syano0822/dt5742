#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstdio>

#include "TFile.h"
#include "TDirectory.h"
#include "TGraph.h"

static const double TSAMPLE_NS   = 0.2;    // 5 GHz → 0.2 ns
static const int    NPEDESTAL    = 100;    
static const double PED_TARGET   = 3500.0; 

void convert2root(const char* inPattern = "wave_%d.dat",
		  int nch = 16,
		  const char* outFile = "wave_form.root")
{

  std::vector<std::ifstream> fins(nch);
  for (int ch = 0; ch < nch; ++ch) {
    char fname[512];

    if (ch == 3) {
      std::snprintf(fname, sizeof(fname), "TR_0_0.dat");
    } else {
      std::snprintf(fname, sizeof(fname), inPattern, ch);
    }

    fins[ch].open(fname, std::ios::binary);

    if (!fins[ch].is_open()) {
      std::cerr << "ERROR: cannot open " << fname << std::endl;
      return;
    }
    std::cout << "Opened " << fname << std::endl;
  }

  TFile* fout = TFile::Open(outFile, "RECREATE");
  if (!fout || fout->IsZombie()) {
    std::cerr << "ERROR: cannot create ROOT file " << outFile << std::endl;
    return;
  }

  int evt = 0;
  const uint32_t headerBytes = 8 * sizeof(uint32_t);
  
  while (true) {
    uint32_t header[8];
    if (!fins[0].read(reinterpret_cast<char*>(header), 8*sizeof(uint32_t))) {
      std::cout << "EOF reached at event " << evt << " (ch0)" << std::endl;
      break;
    }

    uint32_t eventSize0   = header[0];
    uint32_t boardId0     = header[1];
    uint32_t channel0     = header[3];
    uint32_t eventCnt0    = header[4];

    if (eventSize0 <= headerBytes) {
      std::cerr << "ERROR: strange EventSize=" << eventSize0
		<< " at event " << evt << " ch0" << std::endl;
      break;
    }

    uint32_t payloadBytes0 = eventSize0 - headerBytes;
    if (payloadBytes0 % 4 != 0) {
      std::cerr << "ERROR: payloadBytes0 " << payloadBytes0
		<< " not multiple of 4 at event " << evt << " ch0" << std::endl;
      break;
    }

    int nsamples = payloadBytes0 / 4;  // 4 byte / sample (float)

    std::vector<uint32_t> eventSize(nch);
    eventSize[0] = eventSize0;

    for (int ch = 1; ch < nch; ++ch) {
      if (!fins[ch].read(reinterpret_cast<char*>(header), 8*sizeof(uint32_t))) {
	std::cerr << "EOF or read error at event " << evt
		  << " (channel " << ch << ")" << std::endl;
	goto end;
      }
      eventSize[ch] = header[0];
      if (eventSize[ch] != eventSize0) {
	std::cerr << "WARNING: EventSize mismatch at event " << evt
		  << " ch" << ch << " (" << eventSize[ch]
		  << " vs " << eventSize0 << ")" << std::endl;
      }
    }

    std::vector<double> x(nsamples);
    for (int i = 0; i < nsamples; ++i) {
      x[i] = i * TSAMPLE_NS;
    }

    std::vector< std::vector<double> > y(nch, std::vector<double>(nsamples));

    for (int ch = 0; ch < nch; ++ch) {
      uint32_t evSize = eventSize[ch];
      uint32_t payloadBytes = evSize - headerBytes;
      int nsamp_ch = payloadBytes / 4;  // float 4byte/サンプル

      if (nsamp_ch != nsamples) {
	std::cerr << "WARNING: nsamples mismatch at event " << evt
		  << " ch" << ch << " (" << nsamp_ch
		  << " vs " << nsamples << ")" << std::endl;
      }

      std::vector<float> fbuf(nsamp_ch);
      fins[ch].read(reinterpret_cast<char*>(fbuf.data()),
		    nsamp_ch * sizeof(float));
      if (!fins[ch].good()) {
	std::cerr << "ERROR: failed to read data at event "
		  << evt << " ch" << ch << std::endl;
	goto end;
      }

      for (int i = 0; i < nsamp_ch; ++i) {
	y[ch][i] = static_cast<double>(fbuf[i]);
      }
    }
    
    char dname[64];
    std::snprintf(dname, sizeof(dname), "event_%06d", evt);
    TDirectory* evtDir = fout->mkdir(dname);
    evtDir->cd();
    
    for (int ch = 0; ch < nch; ++ch) {
      
      TGraph* g_raw = new TGraph(nsamples, x.data(), y[ch].data());
      char gname_raw[64];
      std::snprintf(gname_raw, sizeof(gname_raw), "ch%02d_raw", ch);
      char gtitle_raw[128];
      std::snprintf(gtitle_raw, sizeof(gtitle_raw),
		    "Event %d Channel %d (raw);time [ns];ADC", evt, ch);
      g_raw->SetName(gname_raw);
      g_raw->SetTitle(gtitle_raw);
      g_raw->Write();

      int nPed = (nsamples < NPEDESTAL) ? nsamples : NPEDESTAL;
      double ped = 0.0;
      for (int i = 0; i < nPed; ++i) {
	ped += y[ch][i];
      }
      ped /= (double)nPed;

      std::vector<double> y_ped(nsamples);
      for (int i = 0; i < nsamples; ++i) {
	y_ped[i] = y[ch][i] - ped + PED_TARGET;
      }

      TGraph* g_ped = new TGraph(nsamples, x.data(), y_ped.data());
      char gname_ped[64];
      std::snprintf(gname_ped, sizeof(gname_ped), "ch%02d_ped", ch);
      char gtitle_ped[128];
      std::snprintf(gtitle_ped, sizeof(gtitle_ped),
		    "Event %d Channel %d (ped->%.1f);time [ns];ADC", evt, ch, PED_TARGET);
      g_ped->SetName(gname_ped);
      g_ped->SetTitle(gtitle_ped);
      g_ped->Write();

      delete g_raw;
      delete g_ped;
    }

    fout->cd();
    ++evt;
  }

 end:
  std::cout << "Total events written: " << evt << std::endl;
  fout->Write();
  fout->Close();
  delete fout;

  for (int ch = 0; ch < nch; ++ch) {
    if (fins[ch].is_open()) fins[ch].close();
  }

  std::cout << "Done." << std::endl;
}
