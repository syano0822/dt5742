#pragma once

#include <chrono>
#include <cstdint>
#include <deque>
#include <fstream>
#include <string>
#include <vector>

#include "config/monitor_config.h"

// Event statistics
struct EventStats {
  uint32_t latest_event_number = 0;
  uint32_t total_events_read = 0;
  uint32_t event_gaps_detected = 0;
  uint32_t corrupted_events = 0;
  std::chrono::steady_clock::time_point start_time;
  std::chrono::steady_clock::time_point last_update_time;

  void UpdateEventNumber(uint32_t new_event_number);
};

// QA status
enum class QAStatus { OK, WARNING, ERROR };

// Waveform QA results
struct WaveformQA {
  float baseline_mean = 0.0f;
  float baseline_rms = 0.0f;
  float waveform_min = 0.0f;
  float waveform_max = 0.0f;
  float noise_estimate = 0.0f;

  QAStatus baseline_status = QAStatus::OK;
  QAStatus range_status = QAStatus::OK;
  QAStatus noise_status = QAStatus::OK;

  bool HasIssues() const;
  std::string GetStatusString() const;
};

// QA summary statistics
struct QASummary {
  uint32_t total_checked = 0;
  uint32_t ok_count = 0;
  uint32_t warning_count = 0;
  uint32_t error_count = 0;

  float avg_baseline = 0.0f;
  float avg_noise = 0.0f;

  void Update(const WaveformQA& qa);
};

// Rate calculator with moving window average
class RateCalculator {
public:
  explicit RateCalculator(int window_seconds = 10);

  void RecordEvent(uint32_t event_number, std::chrono::steady_clock::time_point time);
  double GetRate() const;

private:
  std::deque<std::pair<std::chrono::steady_clock::time_point, uint32_t>> history_;
  std::chrono::seconds window_size_;
};

// File monitor for incremental reading
class FileMonitor {
public:
  explicit FileMonitor(const std::string& file_path);

  bool Open();
  bool IsOpen() const;
  bool CheckNewData();
  std::streampos GetPosition() const;
  std::ifstream& GetStream();

private:
  std::string file_path_;
  std::ifstream file_;
  std::streampos last_position_;
};

// QA checker
class QAChecker {
public:
  explicit QAChecker(const MonitorConfig& config);

  WaveformQA PerformChecks(const std::vector<float>& waveform);

private:
  const MonitorConfig& config_;

  QAStatus CheckBaseline(float baseline_mean);
  QAStatus CheckRange(float waveform_min, float waveform_max);
  QAStatus CheckNoise(float noise_rms);
};

// Display manager
class DisplayManager {
public:
  DisplayManager();

  void PrintStatus(const EventStats& stats, const RateCalculator& rate_calc,
                   const QASummary& qa_summary);
  void PrintFinalSummary(const EventStats& stats, const QASummary& qa_summary);
  void PrintWarning(uint32_t event_number, const WaveformQA& qa);

  bool ShouldUpdate(const EventStats& stats, int events_since_update);

private:
  std::chrono::steady_clock::time_point last_display_update_;
  int display_update_interval_ms_ = 1000;

  std::string FormatRate(double rate) const;
  std::string FormatDuration(std::chrono::seconds duration) const;
  std::string GetCurrentTime() const;
};

// Real-time monitor main class
class RealtimeMonitor {
public:
  explicit RealtimeMonitor(const MonitorConfig& config);

  bool Initialize();
  void Run();
  void Stop();

private:
  MonitorConfig config_;
  FileMonitor file_monitor_;
  EventStats stats_;
  QASummary qa_summary_;
  RateCalculator rate_calc_;
  DisplayManager display_;
  QAChecker qa_checker_;
  std::ofstream log_file_;
  bool running_ = true;
  int events_since_last_update_ = 0;

  void ProcessNewEvents();
  void PerformQACheck(const std::vector<float>& waveform, uint32_t event_num);
  void LogWarning(uint32_t event_number, const WaveformQA& qa);
};

// Utility functions
bool FileExists(const std::string& path);
