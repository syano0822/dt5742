#include "analysis/multi_sensor_analysis.h"

#include <iostream>
#include <map>
#include <algorithm>
#include <cmath>

// ROOT includes
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"
#include "TDirectory.h"

SensorHitInfo SelectMaxAmplitudeStrip(const std::vector<SensorHitInfo>& hits) {
    SensorHitInfo max_hit;
    float max_amp = 0.0f;

    for (const auto& hit : hits) {
        if (hit.hasSignal && hit.ampMax > max_amp) {
            max_amp = hit.ampMax;
            max_hit = hit;
        }
    }

    return max_hit;
}

std::vector<CombinedEvent> ReadAndMatchEvents(const std::string& daq01_root,
                                               const std::string& daq02_root) {
    std::cout << "Opening ROOT files..." << std::endl;
    std::cout << "  DAQ01: " << daq01_root << std::endl;
    std::cout << "  DAQ02: " << daq02_root << std::endl;

    // Open ROOT files
    TFile* file1 = TFile::Open(daq01_root.c_str(), "READ");
    TFile* file2 = TFile::Open(daq02_root.c_str(), "READ");

    if (!file1 || file1->IsZombie()) {
        std::cerr << "ERROR: Cannot open DAQ01 file: " << daq01_root << std::endl;
        return {};
    }
    if (!file2 || file2->IsZombie()) {
        std::cerr << "ERROR: Cannot open DAQ02 file: " << daq02_root << std::endl;
        return {};
    }

    // Get Analysis trees
    TTree* tree1 = (TTree*)file1->Get("Analysis");
    TTree* tree2 = (TTree*)file2->Get("Analysis");

    if (!tree1) {
        std::cerr << "ERROR: Cannot find Analysis tree in DAQ01 file" << std::endl;
        return {};
    }
    if (!tree2) {
        std::cerr << "ERROR: Cannot find Analysis tree in DAQ02 file" << std::endl;
        return {};
    }

    Long64_t nEntries1 = tree1->GetEntries();
    Long64_t nEntries2 = tree2->GetEntries();
    Long64_t nEntries = std::min(nEntries1, nEntries2);

    std::cout << "DAQ01 entries: " << nEntries1 << std::endl;
    std::cout << "DAQ02 entries: " << nEntries2 << std::endl;
    std::cout << "Processing: " << nEntries << " events" << std::endl;

    // Variables for DAQ01
    int event1;
    std::vector<float>* baseline1 = nullptr;
    std::vector<float>* rmsNoise1 = nullptr;
    std::vector<float>* ampMax1 = nullptr;
    std::vector<float>* charge1 = nullptr;
    std::vector<float>* peakTime1 = nullptr;
    std::vector<bool>* hasSignal1 = nullptr;

    // Variables for DAQ02
    int event2;
    std::vector<float>* baseline2 = nullptr;
    std::vector<float>* rmsNoise2 = nullptr;
    std::vector<float>* ampMax2 = nullptr;
    std::vector<float>* charge2 = nullptr;
    std::vector<float>* peakTime2 = nullptr;
    std::vector<bool>* hasSignal2 = nullptr;

    // CFD timing branches for DAQ01 (16 channels)
    float cfd10_1[16], cfd20_1[16], cfd30_1[16];
    // CFD timing branches for DAQ02 (16 channels)
    float cfd10_2[16], cfd20_2[16], cfd30_2[16];

    // Set branch addresses for DAQ01
    tree1->SetBranchAddress("event", &event1);
    tree1->SetBranchAddress("baseline", &baseline1);
    tree1->SetBranchAddress("rmsNoise", &rmsNoise1);
    tree1->SetBranchAddress("ampMax", &ampMax1);
    tree1->SetBranchAddress("charge", &charge1);
    tree1->SetBranchAddress("peakTime", &peakTime1);
    tree1->SetBranchAddress("hasSignal", &hasSignal1);

    // Set CFD timing branches for DAQ01
    for (int ch = 0; ch < 16; ++ch) {
        tree1->SetBranchAddress(Form("ch%02d_timeCFD_10pc", ch), &cfd10_1[ch]);
        tree1->SetBranchAddress(Form("ch%02d_timeCFD_20pc", ch), &cfd20_1[ch]);
        tree1->SetBranchAddress(Form("ch%02d_timeCFD_30pc", ch), &cfd30_1[ch]);
    }

    // Set branch addresses for DAQ02
    tree2->SetBranchAddress("event", &event2);
    tree2->SetBranchAddress("baseline", &baseline2);
    tree2->SetBranchAddress("rmsNoise", &rmsNoise2);
    tree2->SetBranchAddress("ampMax", &ampMax2);
    tree2->SetBranchAddress("charge", &charge2);
    tree2->SetBranchAddress("peakTime", &peakTime2);
    tree2->SetBranchAddress("hasSignal", &hasSignal2);

    // Set CFD timing branches for DAQ02
    for (int ch = 0; ch < 16; ++ch) {
        tree2->SetBranchAddress(Form("ch%02d_timeCFD_10pc", ch), &cfd10_2[ch]);
        tree2->SetBranchAddress(Form("ch%02d_timeCFD_20pc", ch), &cfd20_2[ch]);
        tree2->SetBranchAddress(Form("ch%02d_timeCFD_30pc", ch), &cfd30_2[ch]);
    }

    // Read events and match
    std::vector<CombinedEvent> events;
    events.reserve(nEntries);

    int n_mismatches = 0;

    for (Long64_t entry = 0; entry < nEntries; ++entry) {
        tree1->GetEntry(entry);
        tree2->GetEntry(entry);

        // Check event number match
        if (event1 != event2) {
            if (n_mismatches < 10) {
                std::cerr << "WARNING: Event number mismatch at entry " << entry
                         << " (DAQ01=" << event1 << ", DAQ02=" << event2 << ")" << std::endl;
            }
            n_mismatches++;
            continue;
        }

        CombinedEvent evt;
        evt.event_number = event1;

        // Extract sensor 1, 2 from DAQ01 (channels 0-15)
        // Sensor 1: channels 0-7, strips 0-7
        // Sensor 2: channels 8-15, strips 0-7
        for (int ch = 0; ch < 16; ++ch) {
            int sensorID = (ch < 8) ? 1 : 2;
            int stripID = ch % 8;

            SensorHitInfo hit;
            hit.sensor_id = sensorID;
            hit.strip_id = stripID;
            hit.channel_global = ch;
            hit.hasSignal = (*hasSignal1)[ch];
            hit.baseline = (*baseline1)[ch];
            hit.rmsNoise = (*rmsNoise1)[ch];
            hit.ampMax = (*ampMax1)[ch];
            hit.charge = (*charge1)[ch];
            hit.peakTime = (*peakTime1)[ch];
            hit.timeCFD10 = cfd10_1[ch];
            hit.timeCFD20 = cfd20_1[ch];
            hit.timeCFD30 = cfd30_1[ch];

            evt.sensor_hits[sensorID - 1].push_back(hit);
            evt.has_sensor[sensorID - 1] = true;
        }

        // Extract sensor 3, 4 from DAQ02 (channels 0-15)
        // Sensor 3: channels 0-7, strips 0-7
        // Sensor 4: channels 8-15, strips 0-7
        for (int ch = 0; ch < 16; ++ch) {
            int sensorID = (ch < 8) ? 3 : 4;
            int stripID = ch % 8;

            SensorHitInfo hit;
            hit.sensor_id = sensorID;
            hit.strip_id = stripID;
            hit.channel_global = ch;
            hit.hasSignal = (*hasSignal2)[ch];
            hit.baseline = (*baseline2)[ch];
            hit.rmsNoise = (*rmsNoise2)[ch];
            hit.ampMax = (*ampMax2)[ch];
            hit.charge = (*charge2)[ch];
            hit.peakTime = (*peakTime2)[ch];
            hit.timeCFD10 = cfd10_2[ch];
            hit.timeCFD20 = cfd20_2[ch];
            hit.timeCFD30 = cfd30_2[ch];

            evt.sensor_hits[sensorID - 1].push_back(hit);
            evt.has_sensor[sensorID - 1] = true;
        }

        events.push_back(evt);

        if ((entry + 1) % 1000 == 0) {
            std::cout << "  Processed " << (entry + 1) << " / " << nEntries << " events\r" << std::flush;
        }
    }

    std::cout << std::endl;

    if (n_mismatches > 0) {
        std::cerr << "WARNING: Total event mismatches: " << n_mismatches << std::endl;
    }

    std::cout << "Successfully matched " << events.size() << " events" << std::endl;

    file1->Close();
    file2->Close();
    delete file1;
    delete file2;

    return events;
}

void GenerateAmplitudeMaps(const std::vector<CombinedEvent>& events,
                          TDirectory* outputDir,
                          const std::string& plotDir) {
    std::cout << "Generating amplitude maps..." << std::endl;

    // Create 2D histograms for each sensor (4 total)
    std::map<int, TH2F*> amplitude_maps;

    for (int sensorID = 1; sensorID <= 4; ++sensorID) {
        // Vertical orientation: X=1 bin, Y=strips (0-7), Z=amplitude
        TH2F* hist = new TH2F(Form("sensor%02d_amplitude_map", sensorID),
                             Form("Sensor %02d Amplitude Map;X;Strip;Amplitude (V)", sensorID),
                             1, 0, 1,        // X: single bin
                             8, -0.5, 7.5);  // Y: 8 strips

        amplitude_maps[sensorID] = hist;
    }

    // Fill histograms from all events
    for (const auto& evt : events) {
        for (int sensorID = 1; sensorID <= 4; ++sensorID) {
            const auto& hits = evt.sensor_hits[sensorID - 1];
            for (const auto& hit : hits) {
                amplitude_maps[sensorID]->Fill(0.5, hit.strip_id, hit.ampMax);
            }
        }
    }

    // Save to ROOT file and create PNG images
    outputDir->cd();
    for (auto& pair : amplitude_maps) {
        int sensorID = pair.first;
        TH2F* hist = pair.second;

        hist->Write();

        // Save as PNG
        TCanvas* c = new TCanvas(Form("c_sensor%02d_ampmap", sensorID), "", 800, 600);
        hist->Draw("COLZ");
        c->SaveAs(Form("%s/amplitude_map_sensor%02d.png", plotDir.c_str(), sensorID));
        delete c;
    }

    std::cout << "  Generated amplitude maps for 4 sensors" << std::endl;
}

void AnalyzeSensorBaseline(const std::vector<SensorHitInfo>& hits,
                          int sensorID,
                          TDirectory* outputDir,
                          const std::string& plotDir,
                          std::ofstream& summaryFile) {
    std::cout << "  Analyzing baseline for sensor " << sensorID << "..." << std::endl;

    // Create histogram
    TH1F* h_baseline = new TH1F(Form("sensor%02d_baseline", sensorID),
                                Form("Sensor %02d Baseline Distribution;Baseline (V);Counts", sensorID),
                                200, -0.1, 0.1);

    // Fill with all baseline values
    for (const auto& hit : hits) {
        h_baseline->Fill(hit.baseline);
    }

    // Calculate statistics
    double baseline_mean = h_baseline->GetMean();
    double baseline_rms = h_baseline->GetRMS();
    double signal_threshold = baseline_mean + 3.0 * baseline_rms;

    // Create plot with threshold line
    TCanvas* c = new TCanvas(Form("c_sensor%02d_baseline", sensorID), "", 800, 600);
    h_baseline->Draw();

    TLine* line = new TLine(signal_threshold, 0, signal_threshold, h_baseline->GetMaximum());
    line->SetLineColor(kRed);
    line->SetLineWidth(2);
    line->Draw();

    TLatex* label = new TLatex();
    label->SetTextColor(kRed);
    label->DrawLatexNDC(0.6, 0.8, Form("Threshold = %.4f V", signal_threshold));

    c->SaveAs(Form("%s/baseline_sensor%02d.png", plotDir.c_str(), sensorID));

    // Save to ROOT file
    outputDir->cd();
    h_baseline->Write();

    // Write to summary
    summaryFile << "=============================================================================\n";
    summaryFile << "Sensor " << sensorID << " Baseline Analysis:\n";
    summaryFile << "=============================================================================\n";
    summaryFile << "  Mean: " << baseline_mean << " V\n";
    summaryFile << "  RMS: " << baseline_rms << " V\n";
    summaryFile << "  Signal Threshold (mean + 3*RMS): " << signal_threshold << " V\n";
    summaryFile << std::endl;

    delete c;
}

void AnalyzeSensorAmplitude(const std::vector<SensorHitInfo>& hits,
                           int sensorID,
                           TDirectory* outputDir,
                           const std::string& plotDir,
                           std::ofstream& summaryFile) {
    std::cout << "  Analyzing amplitude for sensor " << sensorID << "..." << std::endl;

    // Create histogram (only hits with signal)
    TH1F* h_amplitude = new TH1F(Form("sensor%02d_amplitude", sensorID),
                                 Form("Sensor %02d Amplitude Distribution;Amplitude (V);Counts", sensorID),
                                 200, 0, 1.0);

    for (const auto& hit : hits) {
        if (hit.hasSignal && hit.ampMax > 0) {
            h_amplitude->Fill(hit.ampMax);
        }
    }

    // Fit with Landau
    TF1* landau_fit = new TF1(Form("landau_sensor%02d", sensorID), "landau", 0, 1.0);
    landau_fit->SetParameters(h_amplitude->GetMaximum(),
                             h_amplitude->GetMean(),
                             h_amplitude->GetRMS());

    h_amplitude->Fit(landau_fit, "SQR");

    // Extract fit parameters
    double mpv = landau_fit->GetParameter(1);
    double width = landau_fit->GetParameter(2);
    double chi2 = landau_fit->GetChisquare();
    int ndf = landau_fit->GetNDF();

    // Create plot
    TCanvas* c = new TCanvas(Form("c_sensor%02d_amplitude", sensorID), "", 800, 600);
    h_amplitude->Draw();
    landau_fit->SetLineColor(kRed);
    landau_fit->SetLineWidth(2);
    landau_fit->Draw("same");

    TLegend* legend = new TLegend(0.6, 0.7, 0.9, 0.9);
    legend->AddEntry(h_amplitude, "Data", "l");
    legend->AddEntry(landau_fit, "Landau Fit", "l");
    legend->AddEntry((TObject*)0, Form("MPV = %.4f V", mpv), "");
    legend->AddEntry((TObject*)0, Form("Width = %.4f V", width), "");
    legend->AddEntry((TObject*)0, Form("#chi^{2}/ndf = %.2f/%d", chi2, ndf), "");
    legend->Draw();

    c->SaveAs(Form("%s/amplitude_sensor%02d.png", plotDir.c_str(), sensorID));

    // Save to ROOT
    outputDir->cd();
    h_amplitude->Write();
    landau_fit->Write();

    // Write to summary
    summaryFile << "=============================================================================\n";
    summaryFile << "Sensor " << sensorID << " Amplitude Analysis:\n";
    summaryFile << "=============================================================================\n";
    summaryFile << "  Landau MPV: " << mpv << " V\n";
    summaryFile << "  Landau Width: " << width << " V\n";
    summaryFile << "  Chi-square/ndf: " << chi2 << "/" << ndf;
    if (ndf > 0) {
        summaryFile << " = " << (chi2 / ndf);
    }
    summaryFile << "\n" << std::endl;

    delete c;
    delete legend;
}

std::vector<TimingPairInfo> GenerateTimingPairs(const std::vector<CombinedEvent>& events) {
    std::cout << "Generating timing pairs..." << std::endl;

    std::vector<TimingPairInfo> pairs;

    // All sensor pair combinations
    const int sensor_pairs[][2] = {{1,2}, {1,3}, {1,4}, {2,3}, {2,4}, {3,4}};

    for (const auto& evt : events) {
        // Count sensors with hits
        int n_sensors_with_hits = 0;
        for (int i = 0; i < 4; ++i) {
            bool has_hit = false;
            for (const auto& hit : evt.sensor_hits[i]) {
                if (hit.hasSignal) {
                    has_hit = true;
                    break;
                }
            }
            if (has_hit) {
                n_sensors_with_hits++;
            }
        }

        // Only analyze events with 2+ sensors
        if (n_sensors_with_hits < 2) continue;

        // For each sensor pair
        for (const auto& pair_ids : sensor_pairs) {
            int s1 = pair_ids[0];
            int s2 = pair_ids[1];

            if (!evt.has_sensor[s1-1] || !evt.has_sensor[s2-1]) continue;

            // Select strip with maximum amplitude for each sensor
            auto hit1 = SelectMaxAmplitudeStrip(evt.sensor_hits[s1-1]);
            auto hit2 = SelectMaxAmplitudeStrip(evt.sensor_hits[s2-1]);

            if (!hit1.hasSignal || !hit2.hasSignal) continue;

            // Create timing pair
            TimingPairInfo pair_info;
            pair_info.sensor1_id = s1;
            pair_info.sensor2_id = s2;
            pair_info.event_number = evt.event_number;
            pair_info.sensor1_max_strip = hit1.strip_id;
            pair_info.sensor2_max_strip = hit2.strip_id;
            pair_info.sensor1_ampMax = hit1.ampMax;
            pair_info.sensor2_ampMax = hit2.ampMax;
            pair_info.delta_time_CFD10 = hit1.timeCFD10 - hit2.timeCFD10;
            pair_info.delta_time_CFD20 = hit1.timeCFD20 - hit2.timeCFD20;
            pair_info.delta_time_CFD30 = hit1.timeCFD30 - hit2.timeCFD30;

            pairs.push_back(pair_info);
        }
    }

    std::cout << "  Generated " << pairs.size() << " timing pairs" << std::endl;

    return pairs;
}

void AnalyzeTimingCorrelation(const std::vector<TimingPairInfo>& pairs,
                             int sensor1_id,
                             int sensor2_id,
                             int cfd_threshold_percent,
                             TDirectory* outputDir,
                             const std::string& plotDir,
                             std::ofstream& summaryFile) {
    std::string cfd_name = Form("CFD%d", cfd_threshold_percent);

    // Filter pairs for this sensor combination
    std::vector<TimingPairInfo> filtered_pairs;
    for (const auto& pair : pairs) {
        if (pair.sensor1_id == sensor1_id && pair.sensor2_id == sensor2_id) {
            filtered_pairs.push_back(pair);
        }
    }

    if (filtered_pairs.empty()) {
        std::cout << "  No timing pairs for sensors " << sensor1_id << " vs " << sensor2_id
                 << " (" << cfd_name << ")" << std::endl;
        return;
    }

    std::cout << "  Analyzing timing: Sensor " << sensor1_id << " vs " << sensor2_id
             << " (" << cfd_name << ") - " << filtered_pairs.size() << " pairs" << std::endl;

    // Create histograms
    TH2F* h2d_amp1 = new TH2F(
        Form("sensor%d_vs_%d_amp1_%s", sensor1_id, sensor2_id, cfd_name.c_str()),
        Form("Sensor %d vs %d (%s);Sensor %d Amplitude (V);#Deltat (ns)",
             sensor1_id, sensor2_id, cfd_name.c_str(), sensor1_id),
        100, 0, 1.0,   // X: amplitude
        100, -50, 50); // Y: time difference

    TH2F* h2d_amp2 = new TH2F(
        Form("sensor%d_vs_%d_amp2_%s", sensor1_id, sensor2_id, cfd_name.c_str()),
        Form("Sensor %d vs %d (%s);Sensor %d Amplitude (V);#Deltat (ns)",
             sensor1_id, sensor2_id, cfd_name.c_str(), sensor2_id),
        100, 0, 1.0, 100, -50, 50);

    TH1F* h1d_delta = new TH1F(
        Form("sensor%d_vs_%d_delta_%s", sensor1_id, sensor2_id, cfd_name.c_str()),
        Form("Sensor %d vs %d (%s);#Deltat (ns);Counts",
             sensor1_id, sensor2_id, cfd_name.c_str()),
        200, -50, 50);

    // Fill histograms
    for (const auto& pair : filtered_pairs) {
        float delta_t;
        if (cfd_threshold_percent == 10) {
            delta_t = pair.delta_time_CFD10;
        } else if (cfd_threshold_percent == 20) {
            delta_t = pair.delta_time_CFD20;
        } else {
            delta_t = pair.delta_time_CFD30;
        }

        h2d_amp1->Fill(pair.sensor1_ampMax, delta_t);
        h2d_amp2->Fill(pair.sensor2_ampMax, delta_t);
        h1d_delta->Fill(delta_t);
    }

    // Fit 1D with Gaussian
    TF1* gauss_fit = new TF1(Form("gauss_%d_%d_%s", sensor1_id, sensor2_id, cfd_name.c_str()),
                            "gaus", -50, 50);
    gauss_fit->SetParameters(h1d_delta->GetMaximum(),
                            h1d_delta->GetMean(),
                            h1d_delta->GetRMS());
    h1d_delta->Fit(gauss_fit, "SQR");

    double mean = gauss_fit->GetParameter(1);
    double sigma = gauss_fit->GetParameter(2);
    double chi2 = gauss_fit->GetChisquare();
    int ndf = gauss_fit->GetNDF();

    // Save 2D histogram 1
    TCanvas* c1 = new TCanvas(Form("c_%d_%d_amp1_%s", sensor1_id, sensor2_id, cfd_name.c_str()), "", 800, 600);
    h2d_amp1->Draw("COLZ");
    c1->SaveAs(Form("%s/timing_2d_amp1_sensor%d_vs_%d_%s.png",
                   plotDir.c_str(), sensor1_id, sensor2_id, cfd_name.c_str()));

    // Save 2D histogram 2
    TCanvas* c2 = new TCanvas(Form("c_%d_%d_amp2_%s", sensor1_id, sensor2_id, cfd_name.c_str()), "", 800, 600);
    h2d_amp2->Draw("COLZ");
    c2->SaveAs(Form("%s/timing_2d_amp2_sensor%d_vs_%d_%s.png",
                   plotDir.c_str(), sensor1_id, sensor2_id, cfd_name.c_str()));

    // Save 1D histogram with fit
    TCanvas* c3 = new TCanvas(Form("c_%d_%d_delta_%s", sensor1_id, sensor2_id, cfd_name.c_str()), "", 800, 600);
    h1d_delta->Draw();
    gauss_fit->SetLineColor(kRed);
    gauss_fit->SetLineWidth(2);
    gauss_fit->Draw("same");

    TLegend* legend = new TLegend(0.6, 0.7, 0.9, 0.9);
    legend->AddEntry(h1d_delta, "Data", "l");
    legend->AddEntry(gauss_fit, "Gaussian Fit", "l");
    legend->AddEntry((TObject*)0, Form("Mean = %.3f ns", mean), "");
    legend->AddEntry((TObject*)0, Form("#sigma = %.3f ns", sigma), "");
    legend->AddEntry((TObject*)0, Form("#chi^{2}/ndf = %.2f/%d", chi2, ndf), "");
    legend->Draw();

    c3->SaveAs(Form("%s/timing_1d_delta_sensor%d_vs_%d_%s.png",
                   plotDir.c_str(), sensor1_id, sensor2_id, cfd_name.c_str()));

    // Save to ROOT
    outputDir->cd();
    h2d_amp1->Write();
    h2d_amp2->Write();
    h1d_delta->Write();
    gauss_fit->Write();

    // Write to summary
    summaryFile << "Timing Correlation: Sensor " << sensor1_id << " vs " << sensor2_id
                << " (" << cfd_name << "):\n";
    summaryFile << "  Gaussian Mean: " << mean << " ns\n";
    summaryFile << "  Gaussian Sigma: " << sigma << " ns\n";
    summaryFile << "  Chi-square/ndf: " << chi2 << "/" << ndf;
    if (ndf > 0) {
        summaryFile << " = " << (chi2 / ndf);
    }
    summaryFile << "\n" << std::endl;

    delete c1;
    delete c2;
    delete c3;
    delete legend;
}
