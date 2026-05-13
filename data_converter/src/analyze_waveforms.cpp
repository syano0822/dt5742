#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <limits>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "utils/file_io.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TTree.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMarker.h"
#include "TLatex.h"
#include "TH2F.h"
#include "TF1.h"
#include "TMath.h"

#include "config/analysis_config.h"
#include "analysis/waveform_math.h"
#include "analysis/waveform_plotting.h"

using namespace std;

namespace {

std::string to6digits(int n) {
    std::ostringstream oss;
    oss << std::setw(6) << std::setfill('0') << n;
    return oss.str();
}

  
bool EnsureParentDirectory(const std::string &path) {
  size_t lastSlash = path.find_last_of('/');
  if (lastSlash != std::string::npos) {
    std::string dirPath = path.substr(0, lastSlash);
    if (!CreateDirectoryIfNeeded(dirPath)) {
      return false;
    }
  }
  return true;
}

enum class NsamplesPolicy { kStrict, kPad };

NsamplesPolicy ResolveNsamplesPolicy(const std::string &policyText) {
  std::string lowered = policyText;
  std::transform(lowered.begin(), lowered.end(), lowered.begin(),
                 [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
  if (lowered == "pad") {
    return NsamplesPolicy::kPad;
  }
  if (lowered != "strict") {
    std::cerr << "WARNING: unknown nsamples_policy '" << policyText
              << "', defaulting to 'strict'" << std::endl;
  }
  return NsamplesPolicy::kStrict;
}

// Helper to check if sensor should be displayed horizontally
bool IsSensorHorizontal(int sensorID, const AnalysisConfig& cfg) {
  // Find unique sensor IDs in current config and map to local index
  std::set<int> uniqueSensors(cfg.sensor_ids.begin(), cfg.sensor_ids.end());
  std::vector<int> sortedSensors(uniqueSensors.begin(), uniqueSensors.end());
  std::sort(sortedSensors.begin(), sortedSensors.end());
  
  // Find index of this sensorID in the sorted unique list
  auto it = std::find(sortedSensors.begin(), sortedSensors.end(), sensorID);
  if (it == sortedSensors.end()) {
    return false;  // Sensor not found
  }

  int localIndex = std::distance(sortedSensors.begin(), it);
  if (localIndex < 0 || localIndex >= static_cast<int>(cfg.sensor_orientations.size())) {
    return false;  // Out of bounds, default to vertical
  }

  return cfg.sensor_orientations[localIndex] == "horizontal";
}

// Lightweight linear calibration: mv = offset + slope * adc
// Matches the interface of TF1::Eval() used in QuickCheckHitMapDUT.C.
struct CalibPol1 {
  float offset = 0.f;
  float slope  = 0.f;
  bool  valid  = false;
  float Eval(float adc) const { return offset + slope * adc; }
};

// ADC-to-mV linear calibration table, indexed by [daq][ch] (both 0-based).
// Channels without a calibration entry have valid=false.
// DAQ 0: ch1-ch15 calibrated; ch0 has no calibration.
// DAQ 1: ch0-ch14 calibrated; ch15 has no calibration.
static const CalibPol1 g_calib_pol1[2][16] = {
  // DAQ 0
  {
    {0, 1, true},                          // ch0: no calibration
    {-3.176676f, 0.286681f, true},         // ch1
    {-3.555424f, 0.288232f, true},         // ch2
    {-3.115196f, 0.287901f, true},         // ch3
    {-2.451742f, 0.281876f, true},         // ch4
    {-2.674749f, 0.284403f, true},         // ch5
    {-2.813437f, 0.285116f, true},         // ch6
    {-2.851733f, 0.286678f, true},         // ch7
    {-2.480619f, 0.279976f, true},         // ch8
    {-2.652653f, 0.285591f, true},         // ch9
    {-1.574811f, 0.277765f, true},         // ch10
    {-1.371317f, 0.277248f, true},         // ch11
    {-1.307235f, 0.276787f, true},         // ch12
    {-1.884459f, 0.281404f, true},         // ch13
    {-1.992954f, 0.281581f, true},         // ch14
    {-1.935413f, 0.281919f, true},         // ch15
  },
  // DAQ 1
  {
    {-1.315517f, 0.304584f, true},         // ch0
    {-1.498006f, 0.312468f, true},         // ch1
    {-1.754055f, 0.305675f, true},         // ch2
    {-1.586169f, 0.306521f, true},         // ch3
    {-1.790790f, 0.306217f, true},         // ch4
    {-1.626740f, 0.307081f, true},         // ch5
    {-1.884997f, 0.309002f, true},         // ch6
    {-1.910199f, 0.309617f, true},         // ch7
    {-1.815818f, 0.306846f, true},         // ch8
    {-1.786031f, 0.306699f, true},         // ch9
    {-1.809334f, 0.307955f, true},         // ch10
    {-1.761285f, 0.308772f, true},         // ch11
    {-1.900867f, 0.306163f, true},         // ch12
    {-1.852977f, 0.308358f, true},         // ch13
    {-1.938787f, 0.309390f, true},         // ch14
    {0, 1, true},                          // ch15: no calibration
  },
};

// ── Fit-based feature extraction ─────────────────────────────────────────────
//
// Results use polarity-corrected space (peak is always positive) to stay
// consistent with the existing ampMax / peakTime output branches.
//
// DUT0-2:
//   1. pol2 fit in [peak_time ± 0.4 ns]  → ampMax_Fit, peakTime_Fit
//   2. erf fit in [t_5%, peak_time]       → timeCFD_Fit, riseTime_Fit, leadingEdge_Fit
//
// DUT3 (sensor_id == 3):
//   1. Data peak → peakTime_Fit (raw sample), amp_ref = data peak amplitude
//   2. erf fit in [t_5%, t_5% + 5 ns]    → ampMax_Fit (= erf |A|),
//                                            timeCFD_Fit, riseTime_Fit, leadingEdge_Fit
//
// All timing quantities are in ns.  Amplitude is in polarity-corrected ADC units.
// Failed / skipped values are set to kBad = -999.
struct FitFeatures {
    static constexpr float kBad = -999.f;
    float ampMax_Fit      = kBad;
    float peakTime_Fit    = kBad;
    float riseTime_Fit    = kBad;
    float leadingEdge_Fit = kBad;   // time at 10% of amplitude (LE discriminator)
    std::vector<float> timeCFD_Fit; // one entry per cfd_thresholds element
};

// Analytical least-squares pol2 fit (y = a0 + a1·x + a2·x²).
// Gaussian elimination with partial pivoting. Returns false if system is singular.
static bool FitPol2Analytical(const double* x, const double* y, int n,
                               double& a0, double& a1, double& a2) {
    double S0=n, Sx=0, Sx2=0, Sx3=0, Sx4=0, Sy=0, Sxy=0, Sx2y=0;
    for (int i = 0; i < n; ++i) {
        double xi=x[i], yi=y[i], xi2=xi*xi;
        Sx   += xi;  Sx2 += xi2;  Sx3 += xi2*xi;  Sx4 += xi2*xi2;
        Sy   += yi;  Sxy += xi*yi; Sx2y += xi2*yi;
    }
    double M[3][4] = {
        {S0,  Sx,  Sx2, Sy},
        {Sx,  Sx2, Sx3, Sxy},
        {Sx2, Sx3, Sx4, Sx2y}
    };
    for (int col = 0; col < 3; ++col) {
        int pivot = col;
        for (int r = col+1; r < 3; ++r)
            if (std::abs(M[r][col]) > std::abs(M[pivot][col])) pivot = r;
        if (pivot != col) std::swap(M[col], M[pivot]);
        if (std::abs(M[col][col]) < 1e-12) return false;
        double inv = 1.0 / M[col][col];
        for (int r = col+1; r < 3; ++r) {
            double f = M[r][col] * inv;
            for (int k = col; k <= 3; ++k) M[r][k] -= f * M[col][k];
        }
    }
    a2 = M[2][3] / M[2][2];
    a1 = (M[1][3] - M[1][2]*a2) / M[1][1];
    a0 = (M[0][3] - M[0][2]*a2 - M[0][1]*a1) / M[0][0];
    return true;
}

FitFeatures ComputeFitFeatures(
    const std::vector<float>& amp_raw,      // ped-subtracted ADC (from ch%02d_ped)
    const std::vector<float>& time,
    float  baseline,                         // already computed by AnalyzeWaveform
    int    polarity,                         // cfg.signal_polarity[ch]: +1 or -1
    int    sensor_id,                        // cfg.sensor_ids[ch]
    const std::vector<int>& cfd_thresholds, // cfg.cfd_thresholds (in percent)
    double rt_low  = 0.1,                   // cfg.rise_time_low  (fraction, e.g. 0.4)
    double rt_high = 0.9)                   // cfg.rise_time_high (fraction, e.g. 0.6)
{
    FitFeatures res;
    res.timeCFD_Fit.assign(cfd_thresholds.size(), FitFeatures::kBad);

    int n = static_cast<int>(std::min(amp_raw.size(), time.size()));
    if (n < 5) return res;

    // Build polarity-corrected waveform (peak is always positive)
    std::vector<double> gx(n), gy(n);
    for (int i = 0; i < n; ++i) {
        gx[i] = time[i];
        gy[i] = static_cast<double>((amp_raw[i] - baseline) * polarity);
    }

    // Peak of corrected signal (maximum value)
    int    peak_idx  = static_cast<int>(
                           std::max_element(gy.begin(), gy.end()) - gy.begin());
    double peak_amp  = gy[peak_idx];
    double peak_time = gx[peak_idx];

    if (peak_amp < 10.0) return res;   // no signal

    // amp_ref: reference amplitude for the erf fit range (may be updated by pol2)
    double amp_ref = peak_amp;

    // ── DUT0-2: pol2 fit around peak ─────────────────────────────────────────
    if (sensor_id != 3) {
        const double kHalfWin = 0.4;
        double t_lo_p2 = peak_time - kHalfWin;
        double t_hi_p2 = peak_time + kHalfWin;

        std::vector<double> px, py;
        px.reserve(16); py.reserve(16);
        for (int i = 0; i < n; ++i) {
            if (gx[i] >= t_lo_p2 && gx[i] <= t_hi_p2) {
                px.push_back(gx[i]);
                py.push_back(gy[i]);
            }
        }

        if (static_cast<int>(px.size()) >= 3) {
            double a0, a1, a2;
            if (FitPol2Analytical(px.data(), py.data(), static_cast<int>(px.size()), a0, a1, a2)) {
                if (a2 < 0.) {  // concave down = local maximum
                    double t_pk = -a1 / (2. * a2);
                    double v_pk =  a0 - a1 * a1 / (4. * a2);
                    if (v_pk > 0.) {
                        res.peakTime_Fit = static_cast<float>(t_pk);
                        res.ampMax_Fit   = static_cast<float>(v_pk);
                        amp_ref          = v_pk;
                    }
                }
            }
        }
    } else {
        // DUT3: peak from data
        res.peakTime_Fit = static_cast<float>(peak_time);
    }

    // ── erf fit on the rising edge (polarity-corrected space) ────────────────
    // V(t) = A · Freq((t − t₀) / σ),  A > 0,  σ > 0
    //
    // Find t_5%: last time before peak where corrected signal crosses 5% of amp_ref
    const double kFracLo    = 0.05;
    double       thresh_5   = kFracLo * amp_ref;
    double       t_lo_erf   = gx[0];

    for (int i = peak_idx; i >= 1; --i) {
        if (gy[i] >= thresh_5 && gy[i - 1] < thresh_5) {
            double dx = gx[i] - gx[i - 1];
            double dy = gy[i] - gy[i - 1];
            t_lo_erf = (std::abs(dy) > 1e-9)
                       ? gx[i - 1] + dx / dy * (thresh_5 - gy[i - 1])
                       : gx[i - 1];
            break;
        }
    }

    // Upper bound of erf window
    double t_hi_erf = (sensor_id != 3) ? peak_time : t_lo_erf + 5.0;
    if (t_hi_erf > gx[n - 1]) t_hi_erf = gx[n - 1];

    double dt_erf = t_hi_erf - t_lo_erf;
    if (dt_erf <= 0.) return res;

    // Collect samples within the erf window
    std::vector<double> ex, ey;
    ex.reserve(64); ey.reserve(64);
    for (int i = 0; i < n; ++i) {
        if (gx[i] >= t_lo_erf && gx[i] <= t_hi_erf) {
            ex.push_back(gx[i]);
            ey.push_back(gy[i]);
        }
    }
    if (static_cast<int>(ex.size()) < 4) return res;

    // Initial parameter estimates
    //   σ ≈ dt / 5.1  (5% ≈ t₀−1.64σ, 100% ≈ t₀+3.5σ → Δ ≈ 5.1σ)
    //   t₀ ≈ t_lo + 1.64·σ
    double sigma_init = dt_erf / 5.1;
    double t0_init    = t_lo_erf + 1.64 * sigma_init;

    static TF1* s_ferf = nullptr;
    if (!s_ferf)
        s_ferf = new TF1("_s_ferf_reuse", "[0]*TMath::Freq((x-[1])/[2])", 0., 1.);
    s_ferf->SetRange(t_lo_erf, t_hi_erf);
    s_ferf->SetParameter(0, amp_ref);
    s_ferf->SetParameter(1, t0_init);
    s_ferf->SetParameter(2, sigma_init);
    s_ferf->SetParLimits(0, 0.05 * amp_ref, 20. * amp_ref);
    s_ferf->SetParLimits(1, t_lo_erf - dt_erf, t_hi_erf);
    s_ferf->SetParLimits(2, 0.02, dt_erf * 3.);
    TGraph gerf(static_cast<int>(ex.size()), ex.data(), ey.data());
    gerf.Fit(s_ferf, "RQN");

    double A_fit     = s_ferf->GetParameter(0);
    double t0_fit    = s_ferf->GetParameter(1);
    double sigma_fit = std::abs(s_ferf->GetParameter(2));

    if (A_fit <= 0. || sigma_fit <= 0.) return res;

    // DUT3: ampMax_Fit = erf flat level |A|
    if (sensor_id == 3)
        res.ampMax_Fit = static_cast<float>(A_fit);

    // RT using configured thresholds (rt_low to rt_high)
    res.riseTime_Fit = static_cast<float>(
        TMath::Sqrt2() * sigma_fit *
        (TMath::ErfInverse(2.*rt_high - 1.) - TMath::ErfInverse(2.*rt_low - 1.)));

    // Leading edge: time at rt_low fraction of amplitude
    res.leadingEdge_Fit = static_cast<float>(
        t0_fit + TMath::Sqrt2() * sigma_fit * TMath::ErfInverse(2.*rt_low - 1.));

    // CFD at each configured threshold
    for (size_t k = 0; k < cfd_thresholds.size(); ++k) {
        double f   = cfd_thresholds[k] / 100.0;
        double arg = 2. * f - 1.;
        if (arg <= -1. || arg >= 1.) continue;
        res.timeCFD_Fit[k] = static_cast<float>(
            t0_fit + TMath::Sqrt2() * sigma_fit * TMath::ErfInverse(arg));
    }

    return res;
}

bool RunAnalysis(const AnalysisConfig &cfg, Long64_t eventStart = -1, Long64_t eventEnd = -1) {
  // Maximum file size for waveform plots: 4 GB
  const Long64_t MAX_PLOTS_FILE_SIZE = 4LL * 1024 * 1024 * 1024;  // 4 GB in bytes
  
  // Create waveform plots ROOT file if enabled
  TFile *waveformPlotsFile = nullptr;
  int waveformPlotsFileCounter = 0;

  string outname_base = cfg.output_dir()+'/';
  outname_base += to6digits(cfg.runnumber())+'/';
  outname_base += cfg.daq_name()+"/output/";
  
  auto openWaveformPlotsFile = [&](TFile *&outFile, int fileNum) {
    if (!cfg.waveform_plots_enabled) {
      return;
    }

    std::string baseFileName = cfg.waveform_plots_dir;
    std::string waveformPlotsFileName;
    
    if (fileNum == 0) {
      waveformPlotsFileName = BuildOutputPath(outname_base, "waveform_plots",
                                               baseFileName + ".root");
    } else {
      char suffix[32];
      std::snprintf(suffix, sizeof(suffix), "_%03d.root", fileNum);
      waveformPlotsFileName = BuildOutputPath(outname_base, "waveform_plots",
                                               baseFileName + suffix);
    }

    if (!EnsureParentDirectory(waveformPlotsFileName)) {
      std::cerr << "WARNING: Failed to create waveform plots output directory for "
                << waveformPlotsFileName << std::endl;
      return;
    }
    outFile = TFile::Open(waveformPlotsFileName.c_str(), "RECREATE");
    if (!outFile || outFile->IsZombie()) {
      std::cerr << "WARNING: Failed to create waveform plots output file "
                << waveformPlotsFileName << std::endl;
      std::cerr << "         Continuing without waveform plots output..." << std::endl;
      outFile = nullptr;
      return;
    }
    std::cout << "Waveform plots output enabled. Saving to: " << waveformPlotsFileName << std::endl;
    if (cfg.waveform_plots_only_signal && fileNum == 0) {
      std::cout << "  Only saving waveforms with detected signals (SNR > "
                << cfg.snr_threshold << ")" << std::endl;
    }
  };
  
  auto checkAndRotateWaveformPlotsFile = [&](TFile *&outFile) {
    if (!outFile || outFile->IsZombie()) {
      return;
    }

    Long64_t currentSize = outFile->GetSize();
    if (currentSize >= MAX_PLOTS_FILE_SIZE) {
      std::cout << "Waveform plots file size reached " << (currentSize / (1024.0 * 1024.0 * 1024.0))
                << " GB. Rotating to new file..." << std::endl;

      // Close current file
      std::string currentFileName = outFile->GetName();
      outFile->Close();
      delete outFile;
      outFile = nullptr;

      std::cout << "Saved waveform plots to: " << currentFileName << std::endl;

      // Open new file with incremented counter
      waveformPlotsFileCounter++;
      openWaveformPlotsFile(outFile, waveformPlotsFileCounter);
    }
  };

  openWaveformPlotsFile(waveformPlotsFile, waveformPlotsFileCounter);

  // Create quality check ROOT file
  TFile *qualityCheckFile = nullptr;
  if (cfg.waveform_plots_enabled) {
    // Use same naming scheme as waveform_plots_dir for quality check
    // If waveform_plots_dir is "waveform_plots", use "quality_check"
    // If waveform_plots_dir is "waveform_plots_chunk_0", use "quality_check_chunk_0"
    std::string qualityCheckBaseName = cfg.waveform_plots_dir;

    // Replace "waveform_plots" with "quality_check" in the base name
    size_t pos = qualityCheckBaseName.find("waveform_plots");
    if (pos != std::string::npos) {
      qualityCheckBaseName.replace(pos, std::string("waveform_plots").length(), "quality_check");
    } else {
      // Fallback: just use "quality_check" prefix
      qualityCheckBaseName = "quality_check_" + qualityCheckBaseName;
    }
    
    std::string qualityCheckFileName = BuildOutputPath(outname_base, "quality_check",
                                                       qualityCheckBaseName + ".root");
    if (!EnsureParentDirectory(qualityCheckFileName)) {
      std::cerr << "WARNING: Failed to create quality_check output directory for "
                << qualityCheckFileName << std::endl;
    } else {
      qualityCheckFile = TFile::Open(qualityCheckFileName.c_str(), "RECREATE");
      if (!qualityCheckFile || qualityCheckFile->IsZombie()) {
        std::cerr << "WARNING: Failed to create quality_check output file "
                  << qualityCheckFileName << std::endl;
        std::cerr << "         Continuing without quality_check output..." << std::endl;
        qualityCheckFile = nullptr;
      } else {
        std::cout << "Quality check output enabled. Saving to: " << qualityCheckFileName << std::endl;
      }
    }
  }

  // Build input path: output_dir/root/input_root
  std::string inputPath = BuildOutputPath(outname_base, "root", cfg.input_root());

  // Open input ROOT file
  TFile *inputFile = TFile::Open(inputPath.c_str(), "READ");
  if (!inputFile || inputFile->IsZombie()) {
    std::cerr << "ERROR: cannot open input ROOT file " << inputPath << std::endl;
    return false;
  }
  std::cout << "Reading input file: " << inputPath << std::endl;

  TTree *inputTree = dynamic_cast<TTree *>(inputFile->Get(cfg.input_tree().c_str()));
  if (!inputTree) {
    std::cerr << "ERROR: cannot find tree " << cfg.input_tree() << std::endl;
    inputFile->Close();
    return false;
  }

  // Pre-warm I/O: reduces the slow-start caused by cold disk reads
  inputTree->SetCacheSize(256LL * 1024 * 1024);  // 256 MB read-ahead cache
  inputTree->AddBranchToCache("*", true);
  inputTree->StopCacheLearningPhase();

  // Get number of entries
  Long64_t totalEntries = inputTree->GetEntries();
  if (totalEntries == 0) {
    std::cerr << "ERROR: input tree has no entries" << std::endl;
    inputFile->Close();
    return false;
  }

  // Determine event range to process
  Long64_t startEntry = (eventStart >= 0) ? eventStart : 0;
  Long64_t endEntry = (eventEnd >= 0) ? eventEnd : totalEntries;

  // Validate range
  if (startEntry < 0) startEntry = 0;
  if (endEntry > totalEntries) endEntry = totalEntries;
  if (startEntry >= endEntry) {
    std::cerr << "ERROR: invalid event range [" << startEntry << ", " << endEntry << ")" << std::endl;
    inputFile->Close();
    return false;
  }

  Long64_t nEntries = endEntry - startEntry;
  std::cout << "Processing event range [" << startEntry << ", " << endEntry << ") - "
            << nEntries << " events" << std::endl;

  const NsamplesPolicy policy = ResolveNsamplesPolicy(cfg.common.nsamples_policy);

  // Set up input branches
  int eventIdx = 0;
  int nChannels = 0;
  int nsamples = 0;
  std::vector<float> *timeAxis = nullptr;
  std::vector<std::vector<float> *> chPed(cfg.n_channels(), nullptr);
  std::vector<int> *nsamplesPerChannel = nullptr;

  inputTree->SetBranchAddress("event", &eventIdx);
  inputTree->SetBranchAddress("n_channels", &nChannels);
  inputTree->SetBranchAddress("nsamples", &nsamples);
  inputTree->SetBranchAddress("time_ns", &timeAxis);
  if (inputTree->GetBranch("nsamples_per_channel")) {
    inputTree->SetBranchAddress("nsamples_per_channel", &nsamplesPerChannel);
  }

  for (int ch = 0; ch < cfg.n_channels(); ++ch) {
    char bname[32];
    std::snprintf(bname, sizeof(bname), "ch%02d_ped", ch);
    if (inputTree->GetBranch(bname)) {
      inputTree->SetBranchAddress(bname, &chPed[ch]);
    }
  }

  // Build output path: output_dir/root/output_root
  std::string outputPath = BuildOutputPath(outname_base, "root", cfg.output_root());
  
  // Create directory if needed
  size_t lastSlash = outputPath.find_last_of('/');
  if (lastSlash != std::string::npos) {
    std::string dirPath = outputPath.substr(0, lastSlash);
    if (!CreateDirectoryIfNeeded(dirPath)) {
      std::cerr << "ERROR: failed to create output directory: " << dirPath << std::endl;
      inputFile->Close();
      return false;
    }
  }

  // Create output ROOT file
  TFile *outputFile = TFile::Open(outputPath.c_str(), "RECREATE");
  if (!outputFile || outputFile->IsZombie()) {
    std::cerr << "ERROR: cannot create output ROOT file " << outputPath << std::endl;
    inputFile->Close();
    return false;
  }
  std::cout << "Creating output file: " << outputPath << std::endl;

  TTree *outputTree = new TTree(cfg.output_tree().c_str(), "Analyzed waveform features");

  // Create output branches - per channel vectors
  int event = 0;
  std::vector<int> sensorID(cfg.n_channels());
  std::vector<int> sensorCol(cfg.n_channels());
  std::vector<int> sensorRow(cfg.n_channels());
  std::vector<bool> isHorizontal(cfg.n_channels());
  std::vector<bool> hasSignal(cfg.n_channels());
  std::vector<float> baseline(cfg.n_channels());
  std::vector<float> rmsNoise(cfg.n_channels());
  std::vector<float> rmsNoise_mV(cfg.n_channels(), 0.f);
  std::vector<float> noise1Point(cfg.n_channels());
  std::vector<float> ampMinBefore(cfg.n_channels());
  std::vector<float> ampMaxBefore(cfg.n_channels());
  std::vector<float> ampMax(cfg.n_channels());
  std::vector<float> ampMax_mV(cfg.n_channels(), 0.f);
  std::vector<float> charge(cfg.n_channels());
  std::vector<float> charge_mV(cfg.n_channels(), 0.f);
  std::vector<float> signalOverNoise(cfg.n_channels());
  std::vector<float> peakTime(cfg.n_channels());
  std::vector<float> riseTime(cfg.n_channels());
  std::vector<float> slewRate(cfg.n_channels());
  std::vector<float> slewRate_mV(cfg.n_channels(), 0.f);
  std::vector<float> jitterRMS(cfg.n_channels());

  // Multi-threshold count (needed by Fit variables below)
  const size_t nCFD    = cfg.cfd_thresholds.size();
  const size_t nLE     = cfg.le_thresholds.size();
  const size_t nCharge = cfg.charge_thresholds.size();

  // Fit-based output variables (_Fit suffix)
  std::vector<float> ampMax_Fit(cfg.n_channels(),      FitFeatures::kBad);
  std::vector<float> ampMax_Fit_mV(cfg.n_channels(),  0.f);
  std::vector<float> peakTime_Fit(cfg.n_channels(),    FitFeatures::kBad);
  std::vector<float> riseTime_Fit(cfg.n_channels(),    FitFeatures::kBad);
  std::vector<float> slewRate_Fit(cfg.n_channels(),    FitFeatures::kBad);
  std::vector<float> slewRate_Fit_mV(cfg.n_channels(), 0.f);
  std::vector<float> jitterRMS_Fit(cfg.n_channels(),   FitFeatures::kBad);
  std::vector<float> leadingEdge_Fit(cfg.n_channels(), FitFeatures::kBad);
  std::vector<std::vector<float>> timeCFD_Fit(
      cfg.n_channels(), std::vector<float>(nCFD, FitFeatures::kBad));

  // Initialize sensor and strip IDs from config
  for (int ch = 0; ch < cfg.n_channels(); ++ch) {
    sensorID[ch]  = cfg.sensor_ids[ch];
    sensorCol[ch] = cfg.sensor_cols[ch]; // Now holds Strip IDs
    sensorRow[ch] = cfg.sensor_rows[ch]; // Now holds Column IDs
    // Check sensor orientation
    isHorizontal[ch] = IsSensorHorizontal(sensorID[ch], cfg);
  }
  
  auto defineScalarBranches = [&]() {
    outputTree->Branch("nChannels", &nChannels);
    outputTree->Branch("event", &event);
    outputTree->Branch("sensorID", &sensorID);
    outputTree->Branch("sensorStrip", &sensorCol); // Rename to sensorStrip in Output? Or keep sensorCol?
                                                   // Keeping sensorCol for consistency with struct, but maybe better to rename to stripID in output?
                                                   // User said "unify to sensor_cols". So keep as sensorCol?
                                                   // Wait, previous code had "stripID" branch.
                                                   // If I remove stripID branch, compatibility might break?
                                                   // I will export sensorCol (which is strips) as "stripID" if I want to maintain branch structure,
                                                   // OR follow user request to "unify to sensor_cols".
                                                   // "Analysis process strip_id is 0-4 but... values are 28"
                                                   // If I change branch name to "sensorCol", user sees "sensorCol" with 0-4.
                                                   // If ID keep "stripID" branch name but fill with sensorCol.
                                                   // Let's use `sensorCol` branch name if user really wants to UNIFY.
                                                   // But maybe "stripID" is clearer for downstream.
                                                   // I will stick to "sensorCol" branch name to match the C++ structure update.
    outputTree->Branch("sensorCol", &sensorCol);
    outputTree->Branch("sensorRow", &sensorRow);
    // outputTree->Branch("stripID", &stripID); // Removed
    outputTree->Branch("isHorizontal", &isHorizontal);
    outputTree->Branch("hasSignal", &hasSignal);
    outputTree->Branch("baseline", &baseline);
    outputTree->Branch("rmsNoise",    &rmsNoise);
    outputTree->Branch("rmsNoise_mV", &rmsNoise_mV);
    outputTree->Branch("noise1Point", &noise1Point);
    outputTree->Branch("ampMinBefore", &ampMinBefore);
    outputTree->Branch("ampMaxBefore", &ampMaxBefore);
    outputTree->Branch("ampMax", &ampMax);
    outputTree->Branch("ampMax_mV", &ampMax_mV);
    outputTree->Branch("charge", &charge);
    outputTree->Branch("charge_mV", &charge_mV);
    outputTree->Branch("signalOverNoise", &signalOverNoise);
    outputTree->Branch("peakTime", &peakTime);
    outputTree->Branch("riseTime", &riseTime);
    outputTree->Branch("slewRate",    &slewRate);
    outputTree->Branch("slewRate_mV", &slewRate_mV);
    outputTree->Branch("jitterRMS",   &jitterRMS);
  };

  // Multi-threshold timing branches (per channel, per threshold)
  std::vector<std::vector<float>> timeCFD(cfg.n_channels(), std::vector<float>(nCFD));
  std::vector<std::vector<float>> jitterCFD(cfg.n_channels(), std::vector<float>(nCFD));
  std::vector<std::vector<float>> timeLE(cfg.n_channels(), std::vector<float>(nLE));
  std::vector<std::vector<float>> jitterLE(cfg.n_channels(), std::vector<float>(nLE));
  std::vector<std::vector<float>> totLE(cfg.n_channels(), std::vector<float>(nLE));
  std::vector<std::vector<float>> timeCharge(cfg.n_channels(), std::vector<float>(nCharge));
  auto defineTimingBranches = [&]() {
    for (int ch = 0; ch < cfg.n_channels(); ++ch) {
      for (size_t i = 0; i < nCFD; ++i) {
        outputTree->Branch(Form("ch%02d_timeCFD_%dpc", ch, cfg.cfd_thresholds[i]),
                          &timeCFD[ch][i]);
        outputTree->Branch(Form("ch%02d_jitterCFD_%dpc", ch, cfg.cfd_thresholds[i]),
                          &jitterCFD[ch][i]);
      }
      for (size_t i = 0; i < nLE; ++i) {
        outputTree->Branch(Form("ch%02d_timeLE_%.1fmV", ch, cfg.le_thresholds[i]),
                          &timeLE[ch][i]);
        outputTree->Branch(Form("ch%02d_jitterLE_%.1fmV", ch, cfg.le_thresholds[i]),
                          &jitterLE[ch][i]);
        outputTree->Branch(Form("ch%02d_totLE_%.1fmV", ch, cfg.le_thresholds[i]),
                          &totLE[ch][i]);
      }
      for (size_t i = 0; i < nCharge; ++i) {
        outputTree->Branch(Form("ch%02d_timeCharge_%dpc", ch, cfg.charge_thresholds[i]),
                          &timeCharge[ch][i]);
      }
    }
  };

  auto defineFitBranches = [&]() {
    outputTree->Branch("ampMax_Fit",      &ampMax_Fit);
    outputTree->Branch("ampMax_Fit_mV",   &ampMax_Fit_mV);
    outputTree->Branch("peakTime_Fit",    &peakTime_Fit);
    outputTree->Branch("riseTime_Fit",    &riseTime_Fit);
    outputTree->Branch("slewRate_Fit",    &slewRate_Fit);
    outputTree->Branch("slewRate_Fit_mV", &slewRate_Fit_mV);
    outputTree->Branch("jitterRMS_Fit",   &jitterRMS_Fit);
    outputTree->Branch("leadingEdge_Fit", &leadingEdge_Fit);
    for (int ch = 0; ch < cfg.n_channels(); ++ch) {
      for (size_t i = 0; i < nCFD; ++i) {
        outputTree->Branch(
            Form("ch%02d_timeCFD_Fit_%dpc", ch, cfg.cfd_thresholds[i]),
            &timeCFD_Fit[ch][i]);
      }
    }
  };

  defineScalarBranches();
  defineTimingBranches();
  defineFitBranches();

  // Determine DAQ index from daq_name (e.g. "daq00" -> 0, "daq01" -> 1)
  int daqIndex = -1;
  {
    const std::string &name = cfg.daq_name();
    size_t pos = name.find_first_of("0123456789");
    if (pos != std::string::npos) {
      try { daqIndex = std::stoi(name.substr(pos)); } catch (...) {}
    }
  }

  // Process all events in the specified range
  std::cout << "Analyzing " << nEntries << " events..." << std::endl;

  // Calculate progress reporting interval (report at least 10 times)
  Long64_t reportInterval = (nEntries < 10) ? 1 : nEntries / 10;
  bool nsamplesError = false;
  bool loggedNsamplesTrim = false;
  std::vector<float> trimmedAmpBuf;
  std::vector<float> trimmedTimeBuf;

  for (Long64_t i = 0; i < nEntries; ++i) {
    Long64_t entry = startEntry + i;

    if (i % reportInterval == 0 || i == nEntries - 1) {
      std::cout << "Processing entry " << entry << " (" << i << " / " << nEntries
                << " = " << (100 * i / nEntries) << "%)" << std::endl;
    }

    inputTree->GetEntry(entry);
    event = eventIdx;

    if (!timeAxis || timeAxis->empty()) {
      std::cerr << "WARNING: empty time axis at entry " << entry << std::endl;
      continue;
    }

    std::vector<int> effectiveSamples(cfg.n_channels(), 0);
    int minSamples = std::numeric_limits<int>::max();
    int maxSamples = 0;
    bool haveSamples = false;

    for (int ch = 0; ch < cfg.n_channels(); ++ch) {
      if (!chPed[ch] || chPed[ch]->empty()) {
        continue;
      }

      int chSamples = nsamples;
      if (nsamplesPerChannel &&
          ch < static_cast<int>(nsamplesPerChannel->size())) {
        chSamples = nsamplesPerChannel->at(ch);
      }

      chSamples =
          std::min(chSamples, static_cast<int>(chPed[ch]->size()));
      chSamples =
          std::min(chSamples, static_cast<int>(timeAxis->size()));

      if (chSamples <= 0) {
        continue;
      }

      haveSamples = true;
      effectiveSamples[ch] = chSamples;
      minSamples = std::min(minSamples, chSamples);
      maxSamples = std::max(maxSamples, chSamples);
    }

    if (!haveSamples) {
      continue;
    }

    const bool hasMismatch = maxSamples != minSamples;
    if (hasMismatch && policy == NsamplesPolicy::kStrict) {
      std::cerr << "ERROR: nsamples mismatch at entry " << entry
                << " (min " << minSamples << ", max " << maxSamples << ")"
                << std::endl;
      nsamplesError = true;
      break;
    }

    if (hasMismatch && !loggedNsamplesTrim &&
        policy == NsamplesPolicy::kPad) {
      std::cout << "INFO: nsamples mismatch detected at entry " << entry
                << ", trimming analysis to per-channel sample counts" << std::endl;
      loggedNsamplesTrim = true;
    }

    // Analyze each channel
    for (int ch = 0; ch < cfg.n_channels(); ++ch) {
      if (!chPed[ch] || chPed[ch]->empty()) {
        continue;
      }

      const int samplesToUse = effectiveSamples[ch];
      if (samplesToUse <= 0) {
        continue;
      }

      const bool needsTrim =
          samplesToUse != static_cast<int>(chPed[ch]->size()) ||
          static_cast<size_t>(samplesToUse) != timeAxis->size();

      const std::vector<float> *ampPtr = chPed[ch];
      const std::vector<float> *timePtr = timeAxis;

      if (needsTrim) {
        trimmedAmpBuf.assign(chPed[ch]->begin(),
                             chPed[ch]->begin() + samplesToUse);
        trimmedTimeBuf.assign(timeAxis->begin(),
                              timeAxis->begin() + samplesToUse);
        ampPtr = &trimmedAmpBuf;
        timePtr = &trimmedTimeBuf;
      }

      WaveformFeatures features = AnalyzeWaveform(*ampPtr, *timePtr, cfg, ch);

      // Fit-based features (DUT0-2: pol2 peak + erf edge; DUT3: erf only)
      ampMax_Fit[ch]      = FitFeatures::kBad;
      peakTime_Fit[ch]    = FitFeatures::kBad;
      riseTime_Fit[ch]    = FitFeatures::kBad;
      leadingEdge_Fit[ch] = FitFeatures::kBad;
      std::fill(timeCFD_Fit[ch].begin(), timeCFD_Fit[ch].end(), FitFeatures::kBad);

      if (features.hasSignal) {
        int sid = (ch < static_cast<int>(cfg.sensor_ids.size()))
                  ? cfg.sensor_ids[ch] : -1;
        int pol = (ch < static_cast<int>(cfg.signal_polarity.size()))
                  ? cfg.signal_polarity[ch] : 1;
        FitFeatures ff = ComputeFitFeatures(
            *ampPtr, *timePtr, features.baseline, pol, sid, cfg.cfd_thresholds,
            cfg.rise_time_low, cfg.rise_time_high);
        ampMax_Fit[ch]      = ff.ampMax_Fit;
        if (daqIndex >= 0 && daqIndex < 2 && ch < 16 && g_calib_pol1[daqIndex][ch].valid
            && ff.ampMax_Fit != FitFeatures::kBad)
          ampMax_Fit_mV[ch] = g_calib_pol1[daqIndex][ch].slope * ff.ampMax_Fit;
        else
          ampMax_Fit_mV[ch] = 0.f;
        peakTime_Fit[ch]    = ff.peakTime_Fit;
        riseTime_Fit[ch]    = ff.riseTime_Fit;
        if (ff.ampMax_Fit != FitFeatures::kBad && ff.riseTime_Fit > 0.f) {
          slewRate_Fit[ch]  = ff.ampMax_Fit
                              * static_cast<float>(cfg.rise_time_high - cfg.rise_time_low)
                              / ff.riseTime_Fit;
          if (daqIndex >= 0 && daqIndex < 2 && ch < 16 && g_calib_pol1[daqIndex][ch].valid)
            slewRate_Fit_mV[ch] = g_calib_pol1[daqIndex][ch].slope * slewRate_Fit[ch];
          else
            slewRate_Fit_mV[ch] = slewRate_Fit[ch];
          jitterRMS_Fit[ch] = (slewRate_Fit[ch] > 0.f)
                              ? features.rmsNoise / slewRate_Fit[ch]
                              : FitFeatures::kBad;
        } else {
          slewRate_Fit[ch]    = FitFeatures::kBad;
          slewRate_Fit_mV[ch] = 0.f;
          jitterRMS_Fit[ch]   = FitFeatures::kBad;
        }
        leadingEdge_Fit[ch] = ff.leadingEdge_Fit;
        for (size_t k = 0; k < nCFD && k < ff.timeCFD_Fit.size(); ++k)
          timeCFD_Fit[ch][k] = ff.timeCFD_Fit[k];
      }

      hasSignal[ch] = features.hasSignal;
      baseline[ch] = features.baseline;
      rmsNoise[ch] = features.rmsNoise;
      if (daqIndex >= 0 && daqIndex < 2 && ch < 16 && g_calib_pol1[daqIndex][ch].valid)
        rmsNoise_mV[ch] = g_calib_pol1[daqIndex][ch].slope * features.rmsNoise;
      else
        rmsNoise_mV[ch] = features.rmsNoise;  // no calibration: ADC = mV (slope=1)
      noise1Point[ch] = features.noise1Point;
      ampMinBefore[ch] = features.ampMinBefore;
      ampMaxBefore[ch] = features.ampMaxBefore;
      ampMax[ch] = features.ampMax;
      if (daqIndex >= 0 && daqIndex < 2 && ch < 16 &&
          g_calib_pol1[daqIndex][ch].valid && features.ampMax > 0.f) {
        ampMax_mV[ch] = g_calib_pol1[daqIndex][ch].Eval(features.ampMax);
      } else {
        ampMax_mV[ch] = 0.f;
      }
      charge[ch] = features.charge;
      if (daqIndex >= 0 && daqIndex < 2 && ch < 16 &&
          g_calib_pol1[daqIndex][ch].valid) {
        charge_mV[ch] = g_calib_pol1[daqIndex][ch].slope * features.charge;
      } else {
        charge_mV[ch] = 0.f;
      }
      signalOverNoise[ch] = features.signalOverNoise;
      peakTime[ch] = features.peakTime;
      riseTime[ch] = features.riseTime;
      slewRate[ch] = features.slewRate;
      if (daqIndex >= 0 && daqIndex < 2 && ch < 16 && g_calib_pol1[daqIndex][ch].valid)
        slewRate_mV[ch] = g_calib_pol1[daqIndex][ch].slope * features.slewRate;
      else
        slewRate_mV[ch] = features.slewRate;
      jitterRMS[ch] = features.jitterRMS;
      
      for (size_t i = 0; i < nCFD && i < features.timeCFD.size(); ++i) {
        timeCFD[ch][i] = features.timeCFD[i];
        jitterCFD[ch][i] = features.jitterCFD[i];
      }
      for (size_t i = 0; i < nLE && i < features.timeLE.size(); ++i) {
        timeLE[ch][i] = features.timeLE[i];	
	jitterLE[ch][i] = features.jitterLE[i];
        totLE[ch][i] = features.totLE[i];
      }
      for (size_t i = 0; i < nCharge && i < features.timeCharge.size(); ++i) {
        timeCharge[ch][i] = features.timeCharge[i];
      }

      // Save waveform plots if enabled
      if (waveformPlotsFile) {
        // Check if we should save this waveform
        bool shouldSave = !cfg.waveform_plots_only_signal || features.hasSignal;
        if (shouldSave) {
          SaveWaveformPlots(waveformPlotsFile, eventIdx, ch, *chPed[ch], *timeAxis, features, cfg);
        }
      }
    }

    // Create 2D histograms for each sensor showing signal amplitudes
    if (waveformPlotsFile || qualityCheckFile) {
      // Determine which sensors are present and how many strips each has
      std::map<int, std::vector<int>> sensorStrips;   // sensor ID -> list of strip IDs
      std::map<int, std::vector<int>> sensorChannels; // sensor ID -> list of channel indices
      std::map<int, std::vector<int>> sensorCols;     // sensor ID -> list of colum indices
      std::map<int, std::vector<int>> sensorRows;     // sensor ID -> list of row indices

      for (int ch = 0; ch < cfg.n_channels(); ++ch) {
        int sensorID = cfg.sensor_ids[ch];
        // In the unified naming:
        // cfg.sensor_cols holds Strip ID
        // cfg.sensor_rows holds Column ID
        int stripID = cfg.sensor_cols[ch]; 
        int colID   = cfg.sensor_rows[ch];
        
        sensorStrips[sensorID].push_back(stripID);
        sensorChannels[sensorID].push_back(ch);
        sensorCols[sensorID].push_back(colID);
        sensorRows[sensorID].push_back(colID);
      }

      // Create event directory if not exists (for waveformPlotsFile)
      char eventDirName[64];
      std::snprintf(eventDirName, sizeof(eventDirName), "event_%06d", eventIdx);
      TDirectory *eventDir = nullptr;
      if (waveformPlotsFile) {
        eventDir = waveformPlotsFile->GetDirectory(eventDirName);
        if (!eventDir) {
          eventDir = waveformPlotsFile->mkdir(eventDirName);
        }
      }

      // Store histograms for quality check canvas
      std::map<int, TH2F*> sensorHistograms;

      // Create histogram for each sensor
      for (const auto &sensorPair : sensorChannels) {
        int sensorID = sensorPair.first;
        const std::vector<int> &strips = sensorPair.second;
        const std::vector<int> &channels = sensorChannels[sensorID];
	
        // Find strip range
        int minStrip = *std::min_element(strips.begin(), strips.end());
        int maxStrip = *std::max_element(strips.begin(), strips.end());
        //int nStrips = maxStrip - minStrip + 1;

	// Check sensor orientation
        bool isHorizontal = IsSensorHorizontal(sensorID, cfg);
	
        TH2F *hist = nullptr;

        if (!isHorizontal) {
          // Vertical: X=Column (sensor_row), Y=Strip (sensor_col)
          hist = new TH2F(Form("sensor%02d_amplitude_map", sensorID),
                          Form("Event %d - Sensor %02d Amplitude Map;Column;Strip;Amplitude (V)",
                               eventIdx, sensorID),
                          2, 0, 2,  // X axis: column (sensor_rows)
                          5, 0, 5);  // Y axis: strip (sensor_cols)
	} else {
	  // Horizontal: X=Strip (sensor_col), Y=Column (sensor_row)
          hist = new TH2F(Form("sensor%02d_amplitude_map", sensorID),
                          Form("Event %d - Sensor %02d Amplitude Map;Strip;Column;Amplitude (V)",
                               eventIdx, sensorID),
                          5, 0, 5,  // X axis: strip (sensor_cols)
                          2, 0, 2);  // Y axis: column (sensor_rows)
        }
	
	std::vector<float> weight_col_pos;
	
        // Fill histogram with amplitude values
        for (size_t i = 0; i < channels.size(); ++i) {
          int ch = channels[i];
          int strip = sensorCol[ch]; // sensor_cols now holds strip ID
          int col = sensorRow[ch];   // sensor_rows now holds column ID
	  float amplitude = ampMax[ch];
	  
	  if (!isHorizontal) // Vertical: X=Column, Y=Strip
	    hist->Fill(col, strip, amplitude);
	  else // Horizontal: X=Strip, Y=Column
	    hist->Fill(strip, col, amplitude);
        }
	
        // Save to sensor directory within event (for waveformPlotsFile)
        if (waveformPlotsFile && eventDir) {
          TDirectory *sensorDir = eventDir->GetDirectory(Form("sensor%02d", sensorID));
          if (!sensorDir) {
            sensorDir = eventDir->mkdir(Form("sensor%02d", sensorID));
          }
          sensorDir->cd();	  
	  hist->Write(hist->GetName(), TObject::kOverwrite);
	}
        // Store histogram for quality check canvas	
	sensorHistograms[sensorID] = hist;
      }
      // Create quality check canvas with all sensors
      if (qualityCheckFile) {
        qualityCheckFile->cd();

        // Determine max sensor ID to size the canvas appropriately
        // This ensures sensors from different DAQs align correctly when merged with hadd
        int maxSensorID = 0;
        for (const auto &histPair : sensorHistograms) {
          if (histPair.first > maxSensorID) {
            maxSensorID = histPair.first;
          }
        }

        // Canvas needs maxSensorID+1 pads to accommodate sensor IDs 0 through maxSensorID
        // Use a minimum of 4 to support typical dual-DAQ setups (sensors 0,1,2,3)
        int numPads = std::max(4, maxSensorID + 1);

        // Create canvas with enough columns for all potential sensors
        char canvasName[64];
        std::snprintf(canvasName, sizeof(canvasName), "event_%06d_quality_check", eventIdx);
	TCanvas *canvas = new TCanvas(canvasName,
                                      Form("Event %d - All Sensors Quality Check", eventIdx),
                                      600 * numPads, 800);  // Width scales with pad count
        canvas->Divide(numPads, 1);  // Enough pads for all sensors

        // Draw each sensor in the pad corresponding to its sensor ID
        // This ensures sensor 0 goes to pad 1, sensor 1 to pad 2, etc.
        for (const auto &histPair : sensorHistograms) {
          int sensorID = histPair.first;
          canvas->cd(sensorID + 1);  // Sensor ID 0 -> pad 1, sensor ID 1 -> pad 2, etc.

          // Clone the histogram to avoid deletion issues
          TH2F *histClone = (TH2F*)histPair.second->Clone(Form("sensor%02d_qc_clone", sensorID));
          histClone->SetStats(0);  // Hide statistics box
          histClone->SetMaximum(5000);  // Fix z-axis maximum to 5000 for consistent comparison
          histClone->Draw("COLZ TEXT");  // Draw with color scale and text values
        }

        // Save canvas to quality check file
        canvas->Write(canvasName, TObject::kOverwrite);
        delete canvas;  // Canvas deletion will handle cloned histograms
      }

      // Clean up sensor histograms
      for (auto &pair : sensorHistograms) {
        delete pair.second;
      }

      if (waveformPlotsFile) {
        waveformPlotsFile->cd();
      }
    }

    // Check if waveform plots file needs rotation (after saving all plots for this event)
    if (waveformPlotsFile) {
      checkAndRotateWaveformPlotsFile(waveformPlotsFile);
    }

    outputTree->Fill();
  }

  if (nsamplesError) {
    if (waveformPlotsFile) {
      waveformPlotsFile->cd();
      waveformPlotsFile->Close();
      delete waveformPlotsFile;
      waveformPlotsFile = nullptr;
    }
    if (qualityCheckFile) {
      qualityCheckFile->cd();
      qualityCheckFile->Close();
      delete qualityCheckFile;
      qualityCheckFile = nullptr;
    }
    outputFile->Close();
    inputFile->Close();
    return false;
  }

  outputFile->cd();
  outputTree->Write();
  outputFile->Close();
  inputFile->Close();

  // Close waveform plots file if it was created
  if (waveformPlotsFile) {
    std::string finalFileName = waveformPlotsFile->GetName();
    waveformPlotsFile->cd();
    waveformPlotsFile->Close();
    delete waveformPlotsFile;
    std::cout << "Waveform plots output saved to " << finalFileName << std::endl;
    if (waveformPlotsFileCounter > 0) {
      std::cout << "  Total files created: " << (waveformPlotsFileCounter + 1)
                << " (split due to 4GB size limit)" << std::endl;
    }
  }

  // Close quality check file if it was created
  if (qualityCheckFile) {
    std::string finalFileName = qualityCheckFile->GetName();
    qualityCheckFile->cd();
    qualityCheckFile->Close();
    delete qualityCheckFile;
    std::cout << "Quality check output saved to " << finalFileName << std::endl;
  }

  std::string outputFullPath = BuildOutputPath(outname_base, "root", cfg.output_root());
  std::cout << "Analysis complete. Output written to " << outputFullPath << std::endl;
  return true;
}

void PrintUsage(const char *prog) {
  std::cout << "Analyze waveforms: Extract timing and amplitude features from ROOT file\n"
            << "Usage: " << prog << " [options]\n"
            << "Options:\n"
            << "  --config PATH          Load analysis settings from JSON file\n"
            << "  --input FILE           Override input ROOT file\n"
            << "  --output FILE          Override output ROOT file\n"
            << "  --event-range START:END  Process only events in range [START, END)\n"
            << "  --waveform-plots       Enable waveform plots output (saves detailed waveform plots)\n"
            << "  --waveform-plots-file NAME  Set waveform plots output ROOT file name (default: waveform_plots.root)\n"
            << "  --waveform-plots-all   Save all waveforms (default: only with signal)\n"
            << "  -h, --help             Show this help message\n";
}

} // namespace

int main(int argc, char **argv) {
  AnalysisConfig cfg;

  // Try to load default config
  std::string defaultPath = "converter_config.json";
  std::string err;
  if (LoadAnalysisConfigFromJson(defaultPath, cfg, &err)) {
    std::cout << "Loaded configuration from " << defaultPath << std::endl;
  }
  
  // Event range parameters
  Long64_t eventStart = -1;
  Long64_t eventEnd = -1;

  // Parse command line
  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    if (arg == "--help" || arg == "-h") {
      PrintUsage(argv[0]);
      return 0;
    } else if (arg == "--config") {
      if (i + 1 >= argc) {
        std::cerr << "ERROR: --config requires a value" << std::endl;
        return 1;
      }
      if (!LoadAnalysisConfigFromJson(argv[++i], cfg, &err)) {
        std::cerr << "ERROR: " << err << std::endl;
        return 1;
      }
      std::cout << "Loaded configuration from " << argv[i] << std::endl;
    } else if (arg == "--input") {
      if (i + 1 >= argc) {
        std::cerr << "ERROR: --input requires a value" << std::endl;
        return 1;
      }
      cfg.set_input_root(argv[++i]);
    } else if (arg == "--output") {
      if (i + 1 >= argc) {
        std::cerr << "ERROR: --output requires a value" << std::endl;
        return 1;
      }
      cfg.set_output_root(argv[++i]);
    } else if (arg == "--event-range") {
      if (i + 1 >= argc) {
        std::cerr << "ERROR: --event-range requires a value (START:END)" << std::endl;
        return 1;
      }
      std::string range = argv[++i];
      size_t colonPos = range.find(':');
      if (colonPos == std::string::npos) {
        std::cerr << "ERROR: --event-range format must be START:END" << std::endl;
        return 1;
      }
      try {
        eventStart = std::stoll(range.substr(0, colonPos));
        eventEnd = std::stoll(range.substr(colonPos + 1));
        std::cout << "Event range: [" << eventStart << ", " << eventEnd << ")" << std::endl;
      } catch (...) {
        std::cerr << "ERROR: invalid event range format" << std::endl;
        return 1;
      }
    } else if (arg == "--waveform-plots") {
      cfg.waveform_plots_enabled = true;
      std::cout << "Waveform plots output enabled" << std::endl;
    } else if (arg == "--waveform-plots-file") {
      if (i + 1 >= argc) {
        std::cerr << "ERROR: --waveform-plots-file requires a value" << std::endl;
        return 1;
      }
      cfg.waveform_plots_dir = argv[++i];
      // Remove .root extension if provided
      if (cfg.waveform_plots_dir.size() > 5 &&
          cfg.waveform_plots_dir.substr(cfg.waveform_plots_dir.size() - 5) == ".root") {
        cfg.waveform_plots_dir = cfg.waveform_plots_dir.substr(0, cfg.waveform_plots_dir.size() - 5);
      }
    } else if (arg == "--waveform-plots-all") {
      cfg.waveform_plots_only_signal = false;
      std::cout << "Will save all waveforms (not just signals)" << std::endl;
    } else {
      std::cerr << "ERROR: unknown option " << arg << std::endl;
      PrintUsage(argv[0]);
      return 1;
    }
  }

  try {
    if (!RunAnalysis(cfg, eventStart, eventEnd)) {
      return 2;
    }
  } catch (const std::exception &ex) {
    std::cerr << "Unhandled exception: " << ex.what() << std::endl;
    return 3;
  }

  return 0;
}
