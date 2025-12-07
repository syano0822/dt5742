#pragma once

#include <vector>
#include <string>
#include <fstream>

// Forward declarations for ROOT classes
class TDirectory;
class TH1F;
class TH2F;
class TF1;

// Per-strip hit information
struct SensorHitInfo {
    int sensor_id;        // 1, 2, 3, or 4
    int strip_id;         // 0-7
    int channel_global;   // Original channel index in ROOT tree
    bool hasSignal;
    float baseline;
    float rmsNoise;
    float ampMax;
    float charge;
    float peakTime;
    float timeCFD10;
    float timeCFD20;
    float timeCFD30;

    // Constructor with default values
    SensorHitInfo() : sensor_id(0), strip_id(0), channel_global(0), hasSignal(false),
                      baseline(0.0f), rmsNoise(0.0f), ampMax(0.0f), charge(0.0f),
                      peakTime(0.0f), timeCFD10(0.0f), timeCFD20(0.0f), timeCFD30(0.0f) {}
};

// Event structure combining both DAQs
struct CombinedEvent {
    int event_number;
    std::vector<SensorHitInfo> sensor_hits[4];  // 4 sensors, each with up to 8 strips
    bool has_sensor[4];  // Track which sensors have data

    // Constructor
    CombinedEvent() : event_number(0) {
        for (int i = 0; i < 4; ++i) {
            has_sensor[i] = false;
        }
    }
};

// Timing pair for correlation analysis
struct TimingPairInfo {
    int sensor1_id;
    int sensor2_id;
    int event_number;
    int sensor1_max_strip;
    int sensor2_max_strip;
    float sensor1_ampMax;
    float sensor2_ampMax;
    float delta_time_CFD10;  // sensor1 - sensor2
    float delta_time_CFD20;
    float delta_time_CFD30;

    // Constructor
    TimingPairInfo() : sensor1_id(0), sensor2_id(0), event_number(0),
                       sensor1_max_strip(0), sensor2_max_strip(0),
                       sensor1_ampMax(0.0f), sensor2_ampMax(0.0f),
                       delta_time_CFD10(0.0f), delta_time_CFD20(0.0f), delta_time_CFD30(0.0f) {}
};

// Function declarations

/**
 * Read and match events from two DAQ ROOT files
 * @param daq01_root Path to DAQ01 analyzed ROOT file
 * @param daq02_root Path to DAQ02 analyzed ROOT file
 * @return Vector of combined events
 */
std::vector<CombinedEvent> ReadAndMatchEvents(const std::string& daq01_root,
                                               const std::string& daq02_root);

/**
 * Generate amplitude maps for all sensors
 * @param events Vector of combined events
 * @param outputDir ROOT directory to save histograms
 * @param plotDir Directory to save PNG images
 */
void GenerateAmplitudeMaps(const std::vector<CombinedEvent>& events,
                          TDirectory* outputDir,
                          const std::string& plotDir);

/**
 * Analyze baseline distribution for a sensor
 * @param hits Vector of hits for this sensor
 * @param sensorID Sensor ID (1-4)
 * @param outputDir ROOT directory to save histograms
 * @param plotDir Directory to save PNG images
 * @param summaryFile Output stream for text summary
 */
void AnalyzeSensorBaseline(const std::vector<SensorHitInfo>& hits,
                          int sensorID,
                          TDirectory* outputDir,
                          const std::string& plotDir,
                          std::ofstream& summaryFile);

/**
 * Analyze amplitude distribution for a sensor with Landau fit
 * @param hits Vector of hits for this sensor
 * @param sensorID Sensor ID (1-4)
 * @param outputDir ROOT directory to save histograms
 * @param plotDir Directory to save PNG images
 * @param summaryFile Output stream for text summary
 */
void AnalyzeSensorAmplitude(const std::vector<SensorHitInfo>& hits,
                           int sensorID,
                           TDirectory* outputDir,
                           const std::string& plotDir,
                           std::ofstream& summaryFile);

/**
 * Generate timing pairs from events (only events with 2+ sensors with hits)
 * @param events Vector of combined events
 * @return Vector of timing pairs
 */
std::vector<TimingPairInfo> GenerateTimingPairs(const std::vector<CombinedEvent>& events);

/**
 * Analyze timing correlation between two sensors for a specific CFD threshold
 * @param pairs Vector of timing pairs (filtered for this sensor pair)
 * @param sensor1_id First sensor ID
 * @param sensor2_id Second sensor ID
 * @param cfd_threshold_percent CFD threshold (10, 20, or 30)
 * @param outputDir ROOT directory to save histograms
 * @param plotDir Directory to save PNG images
 * @param summaryFile Output stream for text summary
 */
void AnalyzeTimingCorrelation(const std::vector<TimingPairInfo>& pairs,
                             int sensor1_id,
                             int sensor2_id,
                             int cfd_threshold_percent,
                             TDirectory* outputDir,
                             const std::string& plotDir,
                             std::ofstream& summaryFile);

/**
 * Select the strip with maximum amplitude from a vector of hits
 * @param hits Vector of hits for a sensor
 * @return Hit with maximum amplitude
 */
SensorHitInfo SelectMaxAmplitudeStrip(const std::vector<SensorHitInfo>& hits);
