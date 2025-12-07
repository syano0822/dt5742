#include <iostream>
#include <string>
#include <cstring>
#include <fstream>
#include <sys/stat.h>
#include <sys/types.h>
#include <chrono>
#include <ctime>

// ROOT includes
#include "TFile.h"
#include "TDirectory.h"
#include "TROOT.h"

// Analysis includes
#include "analysis/multi_sensor_analysis.h"

// Helper function to create directory
bool CreateDirectory(const std::string& path) {
#ifdef _WIN32
    return _mkdir(path.c_str()) == 0 || errno == EEXIST;
#else
    return mkdir(path.c_str(), 0755) == 0 || errno == EEXIST;
#endif
}

void PrintUsage(const char* prog) {
    std::cout << "Combined Multi-Sensor Analysis (Stage 4)\n";
    std::cout << "Usage: " << prog << " [options]\n";
    std::cout << "\n";
    std::cout << "Options:\n";
    std::cout << "  --daq01-root FILE    DAQ01 analyzed ROOT file (required)\n";
    std::cout << "  --daq02-root FILE    DAQ02 analyzed ROOT file (required)\n";
    std::cout << "  --output-dir DIR     Output directory (default: ./combined)\n";
    std::cout << "  --verbose            Enable verbose output\n";
    std::cout << "  -h, --help           Show this help message\n";
    std::cout << "\n";
    std::cout << "Example:\n";
    std::cout << "  " << prog << " \\\n";
    std::cout << "    --daq01-root /data/test07/daq01/output/root/waveforms_analyzed.root \\\n";
    std::cout << "    --daq02-root /data/test07/daq02/output/root/waveforms_analyzed.root \\\n";
    std::cout << "    --output-dir /data/test07/combined\n";
}

int main(int argc, char* argv[]) {
    // Parse command-line arguments
    std::string daq01_root;
    std::string daq02_root;
    std::string output_dir = "./combined";
    bool verbose = false;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-h" || arg == "--help") {
            PrintUsage(argv[0]);
            return 0;
        } else if (arg == "--daq01-root" && i + 1 < argc) {
            daq01_root = argv[++i];
        } else if (arg == "--daq02-root" && i + 1 < argc) {
            daq02_root = argv[++i];
        } else if (arg == "--output-dir" && i + 1 < argc) {
            output_dir = argv[++i];
        } else if (arg == "--verbose") {
            verbose = true;
        } else {
            std::cerr << "ERROR: Unknown option: " << arg << std::endl;
            PrintUsage(argv[0]);
            return 1;
        }
    }

    // Check required arguments
    if (daq01_root.empty() || daq02_root.empty()) {
        std::cerr << "ERROR: Both --daq01-root and --daq02-root are required" << std::endl;
        PrintUsage(argv[0]);
        return 1;
    }

    // Print configuration
    std::cout << "=============================================" << std::endl;
    std::cout << "Combined Multi-Sensor Analysis (Stage 4)" << std::endl;
    std::cout << "=============================================" << std::endl;
    std::cout << "DAQ01 ROOT file: " << daq01_root << std::endl;
    std::cout << "DAQ02 ROOT file: " << daq02_root << std::endl;
    std::cout << "Output directory: " << output_dir << std::endl;
    std::cout << "=============================================" << std::endl;
    std::cout << std::endl;

    // Create output directory
    if (!CreateDirectory(output_dir)) {
        std::cerr << "WARNING: Could not create output directory (may already exist)" << std::endl;
    }

    // Create plots subdirectory
    std::string plot_dir = output_dir + "/plots";
    if (!CreateDirectory(plot_dir)) {
        std::cerr << "WARNING: Could not create plots directory (may already exist)" << std::endl;
    }

    // Disable ROOT graphics (batch mode)
    gROOT->SetBatch(true);

    // Read and match events
    std::vector<CombinedEvent> events = ReadAndMatchEvents(daq01_root, daq02_root);

    if (events.empty()) {
        std::cerr << "ERROR: No events to analyze" << std::endl;
        return 1;
    }

    // Create output ROOT file
    std::string output_root_file = output_dir + "/combined_analysis.root";
    TFile* output_file = new TFile(output_root_file.c_str(), "RECREATE");
    if (!output_file || output_file->IsZombie()) {
        std::cerr << "ERROR: Cannot create output ROOT file: " << output_root_file << std::endl;
        return 1;
    }

    std::cout << "\nCreating ROOT file: " << output_root_file << std::endl;

    // Create directories in ROOT file
    TDirectory* dir_amplitude_maps = output_file->mkdir("AmplitudeMaps");
    TDirectory* dir_baseline = output_file->mkdir("BaselineAnalysis");
    TDirectory* dir_amplitude = output_file->mkdir("AmplitudeAnalysis");
    TDirectory* dir_timing = output_file->mkdir("TimingCorrelation");
    TDirectory* dir_timing_cfd10 = dir_timing->mkdir("CFD10");
    TDirectory* dir_timing_cfd20 = dir_timing->mkdir("CFD20");
    TDirectory* dir_timing_cfd30 = dir_timing->mkdir("CFD30");

    // Open summary text file
    std::string summary_file_path = output_dir + "/analysis_summary.txt";
    std::ofstream summary_file(summary_file_path);

    // Get current time for summary header
    auto now = std::chrono::system_clock::now();
    std::time_t now_time = std::chrono::system_clock::to_time_t(now);

    summary_file << "===============================================================================\n";
    summary_file << "Combined Multi-Sensor Analysis Summary\n";
    summary_file << "Generated: " << std::ctime(&now_time);
    summary_file << "===============================================================================\n";
    summary_file << "\n";
    summary_file << "Data Sources:\n";
    summary_file << "  DAQ01: " << daq01_root << "\n";
    summary_file << "  DAQ02: " << daq02_root << "\n";
    summary_file << "  Total Events: " << events.size() << "\n";
    summary_file << "\n";

    // Generate amplitude maps
    std::cout << "\n=============================================" << std::endl;
    std::cout << "Generating Amplitude Maps" << std::endl;
    std::cout << "=============================================" << std::endl;
    dir_amplitude_maps->cd();
    GenerateAmplitudeMaps(events, dir_amplitude_maps, plot_dir);

    // Baseline analysis for each sensor
    std::cout << "\n=============================================" << std::endl;
    std::cout << "Baseline Analysis" << std::endl;
    std::cout << "=============================================" << std::endl;
    summary_file << "===============================================================================\n";
    summary_file << "BASELINE ANALYSIS\n";
    summary_file << "===============================================================================\n";
    summary_file << "\n";

    dir_baseline->cd();
    for (int sensorID = 1; sensorID <= 4; ++sensorID) {
        // Collect all hits for this sensor
        std::vector<SensorHitInfo> sensor_hits;
        for (const auto& evt : events) {
            for (const auto& hit : evt.sensor_hits[sensorID - 1]) {
                sensor_hits.push_back(hit);
            }
        }
        AnalyzeSensorBaseline(sensor_hits, sensorID, dir_baseline, plot_dir, summary_file);
    }

    // Amplitude analysis for each sensor
    std::cout << "\n=============================================" << std::endl;
    std::cout << "Amplitude Analysis" << std::endl;
    std::cout << "=============================================" << std::endl;
    summary_file << "===============================================================================\n";
    summary_file << "AMPLITUDE ANALYSIS\n";
    summary_file << "===============================================================================\n";
    summary_file << "\n";

    dir_amplitude->cd();
    for (int sensorID = 1; sensorID <= 4; ++sensorID) {
        // Collect all hits for this sensor
        std::vector<SensorHitInfo> sensor_hits;
        for (const auto& evt : events) {
            for (const auto& hit : evt.sensor_hits[sensorID - 1]) {
                sensor_hits.push_back(hit);
            }
        }
        AnalyzeSensorAmplitude(sensor_hits, sensorID, dir_amplitude, plot_dir, summary_file);
    }

    // Generate timing pairs
    std::cout << "\n=============================================" << std::endl;
    std::cout << "Timing Correlation Analysis" << std::endl;
    std::cout << "=============================================" << std::endl;
    std::vector<TimingPairInfo> timing_pairs = GenerateTimingPairs(events);

    summary_file << "===============================================================================\n";
    summary_file << "TIMING CORRELATIONS\n";
    summary_file << "===============================================================================\n";
    summary_file << "\n";

    // Sensor pair combinations
    const int sensor_pairs[][2] = {{1,2}, {1,3}, {1,4}, {2,3}, {2,4}, {3,4}};
    const int cfd_thresholds[] = {10, 20, 30};

    // Analyze timing correlations for each CFD threshold
    for (int cfd_idx = 0; cfd_idx < 3; ++cfd_idx) {
        int cfd_threshold = cfd_thresholds[cfd_idx];
        TDirectory* dir_cfd = nullptr;

        if (cfd_threshold == 10) {
            dir_cfd = dir_timing_cfd10;
        } else if (cfd_threshold == 20) {
            dir_cfd = dir_timing_cfd20;
        } else {
            dir_cfd = dir_timing_cfd30;
        }

        dir_cfd->cd();

        std::cout << "\nCFD" << cfd_threshold << " Timing Analysis:" << std::endl;
        summary_file << "-----------------------------------------------------------------------------\n";
        summary_file << "CFD" << cfd_threshold << " Timing Analysis:\n";
        summary_file << "-----------------------------------------------------------------------------\n";

        for (const auto& pair_ids : sensor_pairs) {
            int s1 = pair_ids[0];
            int s2 = pair_ids[1];

            AnalyzeTimingCorrelation(timing_pairs, s1, s2, cfd_threshold,
                                    dir_cfd, plot_dir, summary_file);
        }

        summary_file << "\n";
    }

    // Close files
    summary_file << "===============================================================================\n";
    summary_file << "Analysis Complete\n";
    summary_file << "===============================================================================\n";
    summary_file.close();

    output_file->Write();
    output_file->Close();
    delete output_file;

    std::cout << "\n=============================================" << std::endl;
    std::cout << "Analysis Complete!" << std::endl;
    std::cout << "=============================================" << std::endl;
    std::cout << "Output ROOT file: " << output_root_file << std::endl;
    std::cout << "Summary text file: " << summary_file_path << std::endl;
    std::cout << "Plots directory: " << plot_dir << std::endl;
    std::cout << "=============================================" << std::endl;

    return 0;
}
