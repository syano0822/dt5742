#include "monitor/realtime_monitor.h"
#include "utils/file_io.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <sys/stat.h>
#include <thread>

// Utility function to check if file exists
bool FileExists(const std::string& path) {
  struct stat buffer;
  return (stat(path.c_str(), &buffer) == 0);
}

// EventStats implementation
void EventStats::UpdateEventNumber(uint32_t new_event_number) {
  if (latest_event_number > 0 && new_event_number != latest_event_number + 1) {
    event_gaps_detected++;
  }
  latest_event_number = new_event_number;
  total_events_read++;
  last_update_time = std::chrono::steady_clock::now();
}

// WaveformQA implementation
bool WaveformQA::HasIssues() const {
  return baseline_status != QAStatus::OK || range_status != QAStatus::OK ||
         noise_status != QAStatus::OK;
}

std::string WaveformQA::GetStatusString() const {
  std::stringstream ss;
  if (baseline_status != QAStatus::OK) {
    ss << "Baseline deviation = " << std::abs(baseline_mean - 3500.0f) << "V";
  }
  if (range_status != QAStatus::OK) {
    ss << (ss.str().empty() ? "" : ", ");
    ss << "Signal out of range (min=" << waveform_min << "V, max=" << waveform_max << "V)";
  }
  if (noise_status != QAStatus::OK) {
    ss << (ss.str().empty() ? "" : ", ");
    ss << "Excessive noise (RMS=" << noise_estimate << "V)";
  }
  return ss.str();
}

// QASummary implementation
void QASummary::Update(const WaveformQA& qa) {
  total_checked++;

  if (!qa.HasIssues()) {
    ok_count++;
  } else if (qa.baseline_status == QAStatus::ERROR || qa.range_status == QAStatus::ERROR ||
             qa.noise_status == QAStatus::ERROR) {
    error_count++;
  } else {
    warning_count++;
  }

  // Update running averages
  float alpha = 1.0f / total_checked;
  avg_baseline = avg_baseline * (1.0f - alpha) + qa.baseline_mean * alpha;
  avg_noise = avg_noise * (1.0f - alpha) + qa.noise_estimate * alpha;
}

// RateCalculator implementation
RateCalculator::RateCalculator(int window_seconds)
    : window_size_(window_seconds) {}

void RateCalculator::RecordEvent(uint32_t event_number,
                                 std::chrono::steady_clock::time_point time) {
  history_.push_back({time, event_number});

  // Remove events outside the window
  auto cutoff = time - window_size_;
  while (!history_.empty() && history_.front().first < cutoff) {
    history_.pop_front();
  }
}

double RateCalculator::GetRate() const {
  if (history_.size() < 2) {
    return 0.0;
  }

  auto time_span = std::chrono::duration<double>(history_.back().first - history_.front().first).count();
  if (time_span < 0.001) {
    return 0.0;
  }

  uint32_t event_span = history_.back().second - history_.front().second;
  return event_span / time_span;
}

// FileMonitor implementation
FileMonitor::FileMonitor(const std::string& file_path)
    : file_path_(file_path), last_position_(0) {}

bool FileMonitor::Open() {
  file_.open(file_path_, std::ios::binary);
  if (!file_.is_open()) {
    return false;
  }
  last_position_ = 0;
  return true;
}

bool FileMonitor::IsOpen() const {
  return file_.is_open();
}

bool FileMonitor::CheckNewData() {
  if (!file_.is_open()) {
    return false;
  }

  file_.seekg(0, std::ios::end);
  std::streampos current_size = file_.tellg();

  if (current_size > last_position_) {
    file_.seekg(last_position_);
    return true;
  }

  return false;
}

std::streampos FileMonitor::GetPosition() const {
  return last_position_;
}

std::ifstream& FileMonitor::GetStream() {
  return file_;
}

// QAChecker implementation
QAChecker::QAChecker(const MonitorConfig& config) : config_(config) {}

WaveformQA QAChecker::PerformChecks(const std::vector<float>& waveform) {
  WaveformQA qa;

  if (waveform.empty()) {
    return qa;
  }

  // Calculate baseline (mean of first N samples)
  int n_pedestal = std::min(config_.qa_pedestal_samples, static_cast<int>(waveform.size()));
  float sum = 0.0f;
  for (int i = 0; i < n_pedestal; i++) {
    sum += waveform[i];
  }
  qa.baseline_mean = sum / n_pedestal;

  // Calculate baseline RMS (noise estimate)
  float variance = 0.0f;
  for (int i = 0; i < n_pedestal; i++) {
    float diff = waveform[i] - qa.baseline_mean;
    variance += diff * diff;
  }
  qa.baseline_rms = std::sqrt(variance / n_pedestal);
  qa.noise_estimate = qa.baseline_rms;

  // Find min/max of waveform
  auto [min_it, max_it] = std::minmax_element(waveform.begin(), waveform.end());
  qa.waveform_min = *min_it;
  qa.waveform_max = *max_it;

  // Perform QA checks
  qa.baseline_status = CheckBaseline(qa.baseline_mean);
  qa.range_status = CheckRange(qa.waveform_min, qa.waveform_max);
  qa.noise_status = CheckNoise(qa.baseline_rms);

  return qa;
}

QAStatus QAChecker::CheckBaseline(float baseline_mean) {
  float deviation = std::abs(baseline_mean - config_.qa_baseline_target);

  if (deviation > config_.qa_baseline_tolerance * 2.0f) {
    return QAStatus::ERROR;
  } else if (deviation > config_.qa_baseline_tolerance) {
    return QAStatus::WARNING;
  }
  return QAStatus::OK;
}

QAStatus QAChecker::CheckRange(float waveform_min, float waveform_max) {
  if (waveform_min < config_.qa_signal_min || waveform_max > config_.qa_signal_max) {
    return QAStatus::ERROR;
  }

  // Check for suspiciously flat signals
  if (std::abs(waveform_max - waveform_min) < 1.0f) {
    return QAStatus::WARNING;
  }

  return QAStatus::OK;
}

QAStatus QAChecker::CheckNoise(float noise_rms) {
  if (noise_rms > config_.qa_noise_threshold * 2.0f) {
    return QAStatus::ERROR;
  } else if (noise_rms > config_.qa_noise_threshold) {
    return QAStatus::WARNING;
  }
  return QAStatus::OK;
}

// DisplayManager implementation
DisplayManager::DisplayManager()
    : last_display_update_(std::chrono::steady_clock::now()) {}

std::string DisplayManager::GetCurrentTime() const {
  auto now = std::chrono::system_clock::now();
  auto time = std::chrono::system_clock::to_time_t(now);
  std::stringstream ss;
  ss << std::put_time(std::localtime(&time), "%H:%M:%S");
  return ss.str();
}

std::string DisplayManager::FormatRate(double rate) const {
  std::stringstream ss;
  if (rate < 1.0) {
    ss << std::fixed << std::setprecision(1) << (rate * 60.0) << " evt/min";
  } else {
    ss << std::fixed << std::setprecision(1) << rate << " evt/s";
  }
  return ss.str();
}

std::string DisplayManager::FormatDuration(std::chrono::seconds duration) const {
  auto hours = std::chrono::duration_cast<std::chrono::hours>(duration);
  duration -= hours;
  auto minutes = std::chrono::duration_cast<std::chrono::minutes>(duration);
  duration -= minutes;
  auto seconds = duration;

  std::stringstream ss;
  ss << std::setfill('0') << std::setw(2) << hours.count() << ":"
     << std::setw(2) << minutes.count() << ":"
     << std::setw(2) << seconds.count();
  return ss.str();
}

bool DisplayManager::ShouldUpdate(const EventStats& stats, int events_since_update) {
  auto now = std::chrono::steady_clock::now();
  auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(
      now - last_display_update_).count();

  return (elapsed >= display_update_interval_ms_) || (events_since_update >= 10);
}

void DisplayManager::PrintStatus(const EventStats& stats, const RateCalculator& rate_calc,
                                  const QASummary& qa_summary) {
  auto now = std::chrono::steady_clock::now();
  auto runtime = std::chrono::duration_cast<std::chrono::seconds>(now - stats.start_time);

  double rate = rate_calc.GetRate();
  std::string rate_str = (rate > 0.0) ? FormatRate(rate) : "calculating...";

  std::cout << "\r[" << GetCurrentTime() << "] "
            << "Event: " << stats.latest_event_number << " | "
            << "Rate: " << rate_str << " | "
            << "Total: " << stats.total_events_read << " | "
            << "QA: OK=" << qa_summary.ok_count
            << " WARN=" << qa_summary.warning_count
            << " ERR=" << qa_summary.error_count << " | "
            << "Runtime: " << FormatDuration(runtime)
            << std::flush;

  last_display_update_ = now;
}

void DisplayManager::PrintWarning(uint32_t event_number, const WaveformQA& qa) {
  std::string severity = "WARNING";
  if (qa.baseline_status == QAStatus::ERROR || qa.range_status == QAStatus::ERROR ||
      qa.noise_status == QAStatus::ERROR) {
    severity = "ERROR";
  }

  std::cout << "\n[" << severity << "] Event " << event_number << ": "
            << qa.GetStatusString() << std::endl;
}

void DisplayManager::PrintFinalSummary(const EventStats& stats, const QASummary& qa_summary) {
  auto runtime = std::chrono::duration_cast<std::chrono::seconds>(
      stats.last_update_time - stats.start_time);
  double avg_rate = (runtime.count() > 0) ? (stats.total_events_read / static_cast<double>(runtime.count())) : 0.0;

  std::cout << "\n\n";
  std::cout << "═════════════════════════════════════════════════════\n";
  std::cout << "         Monitoring Session Summary\n";
  std::cout << "═════════════════════════════════════════════════════\n";
  std::cout << "  Total Events:       " << stats.total_events_read << "\n";
  std::cout << "  Latest Event:       " << stats.latest_event_number << "\n";
  std::cout << "  Event Gaps:         " << stats.event_gaps_detected << "\n";
  std::cout << "  Runtime:            " << FormatDuration(runtime) << "\n";
  std::cout << "  Average Rate:       " << std::fixed << std::setprecision(1) << avg_rate << " evt/s\n";
  std::cout << "\n";

  if (qa_summary.total_checked > 0) {
    std::cout << "  QA Checks:          " << qa_summary.total_checked << "\n";
    std::cout << "    OK:               " << qa_summary.ok_count
              << " (" << std::fixed << std::setprecision(1)
              << (100.0 * qa_summary.ok_count / qa_summary.total_checked) << "%)\n";
    std::cout << "    Warnings:         " << qa_summary.warning_count << "\n";
    std::cout << "    Errors:           " << qa_summary.error_count << "\n";
    std::cout << "    Avg Baseline:     " << std::fixed << std::setprecision(1)
              << qa_summary.avg_baseline << " V\n";
    std::cout << "    Avg Noise:        " << std::fixed << std::setprecision(1)
              << qa_summary.avg_noise << " V RMS\n";
  }

  std::cout << "═════════════════════════════════════════════════════\n\n";
}

// RealtimeMonitor implementation
RealtimeMonitor::RealtimeMonitor(const MonitorConfig& config)
    : config_(config),
      file_monitor_(config.input_file),
      rate_calc_(config.rate_window_seconds),
      qa_checker_(config) {
  display_.display_update_interval_ms_ = config.display_update_interval_ms;
}

bool RealtimeMonitor::Initialize() {
  // Wait for file to exist
  while (!FileExists(config_.input_file) && running_) {
    std::cout << "Waiting for DAQ to start (" << config_.input_file << ")...\r" << std::flush;
    std::this_thread::sleep_for(std::chrono::milliseconds(config_.polling_interval_ms));
  }

  if (!running_) {
    return false;
  }

  // Open file
  if (!file_monitor_.Open()) {
    std::cerr << "\nError: Cannot open file " << config_.input_file << std::endl;
    return false;
  }

  // Open log file if enabled
  if (config_.log_warnings) {
    log_file_.open(config_.log_file, std::ios::app);
    if (log_file_.is_open()) {
      auto now = std::chrono::system_clock::now();
      auto time = std::chrono::system_clock::to_time_t(now);
      log_file_ << "\n[" << std::put_time(std::localtime(&time), "%Y-%m-%d %H:%M:%S")
                << "] Monitor started, file: " << config_.input_file << std::endl;
    }
  }

  stats_.start_time = std::chrono::steady_clock::now();
  stats_.last_update_time = stats_.start_time;

  std::cout << "\nMonitoring started. Press Ctrl+C to stop.\n\n";
  return true;
}

void RealtimeMonitor::Run() {
  while (running_) {
    // Check for new data
    if (file_monitor_.CheckNewData()) {
      ProcessNewEvents();
    }

    // Sleep until next poll
    std::this_thread::sleep_for(std::chrono::milliseconds(config_.polling_interval_ms));
  }
}

void RealtimeMonitor::Stop() {
  running_ = false;
}

void RealtimeMonitor::ProcessNewEvents() {
  auto& file = file_monitor_.GetStream();

  while (file.good() && running_) {
    std::streampos current_pos = file.tellg();

    // Read header
    ChannelHeader header;
    if (!ReadHeader(file, header)) {
      // Incomplete header, rewind and wait for next poll
      file.clear();
      file.seekg(current_pos);
      break;
    }

    // Validate event size
    if (header.eventSize <= HEADER_BYTES) {
      std::cerr << "\nError: Invalid event size at event " << stats_.total_events_read << std::endl;
      stats_.corrupted_events++;
      break;
    }

    // Read waveform
    uint32_t payload_bytes = header.eventSize - HEADER_BYTES;
    int nsamples = payload_bytes / sizeof(float);
    std::vector<float> waveform(nsamples);

    file.read(reinterpret_cast<char*>(waveform.data()), payload_bytes);
    if (!file.good() && !file.eof()) {
      // Incomplete payload, rewind and wait for next poll
      file.clear();
      file.seekg(current_pos);
      break;
    }

    // Update statistics
    stats_.UpdateEventNumber(header.eventCounter);
    rate_calc_.RecordEvent(header.eventCounter, std::chrono::steady_clock::now());
    events_since_last_update_++;

    // Perform QA check (sampled)
    if (config_.qa_enabled && stats_.total_events_read % config_.qa_sampling_interval == 0) {
      PerformQACheck(waveform, header.eventCounter);
    }

    // Update display
    if (display_.ShouldUpdate(stats_, events_since_last_update_)) {
      display_.PrintStatus(stats_, rate_calc_, qa_summary_);
      events_since_last_update_ = 0;
    }
  }

  // Update file position for next poll
  file_monitor_.last_position_ = file.tellg();
}

void RealtimeMonitor::PerformQACheck(const std::vector<float>& waveform, uint32_t event_num) {
  WaveformQA qa = qa_checker_.PerformChecks(waveform);
  qa_summary_.Update(qa);

  if (qa.HasIssues()) {
    display_.PrintWarning(event_num, qa);
    LogWarning(event_num, qa);
  }
}

void RealtimeMonitor::LogWarning(uint32_t event_number, const WaveformQA& qa) {
  if (!log_file_.is_open()) {
    return;
  }

  std::string severity = "WARNING";
  if (qa.baseline_status == QAStatus::ERROR || qa.range_status == QAStatus::ERROR ||
      qa.noise_status == QAStatus::ERROR) {
    severity = "ERROR";
  }

  auto now = std::chrono::system_clock::now();
  auto time = std::chrono::system_clock::to_time_t(now);

  log_file_ << "[" << std::put_time(std::localtime(&time), "%Y-%m-%d %H:%M:%S")
            << "] Event " << event_number << ": " << severity << " - "
            << qa.GetStatusString() << std::endl;
}
