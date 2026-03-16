#ifndef LOGGER_H
#define LOGGER_H

#include <iostream>
#include <sstream>
#include <memory>
#include <chrono>
#include <iomanip>
#include <ctime>

enum class LogLevel { trace, debug, info, warning, error, fatal };

class LogStream : public std::ostream {
  class LogBuf : public std::streambuf {
    LogLevel level_;
    std::string buffer_;
  protected:
    int overflow(int c) override {
      if (c != EOF) {
        buffer_ += static_cast<char>(c);
        if (c == '\n') {
          flush_line();
        }
      }
      return c;
    }
    int sync() override {
      if (!buffer_.empty()) flush_line();
      return 0;
    }
    void flush_line() {
      auto now = std::chrono::system_clock::now();
      auto t = std::chrono::system_clock::to_time_t(now);
      std::tm tm{};
      localtime_r(&t, &tm);
      static const char* names[] = {"trace","debug","info","warning","error","fatal"};
      std::cerr << std::put_time(&tm, "%Y-%m-%d %H:%M:%S")
                << ": <" << names[static_cast<int>(level_)] << "> "
                << buffer_;
      if (!buffer_.empty() && buffer_.back() != '\n') std::cerr << '\n';
      buffer_.clear();
    }
  public:
    LogBuf(LogLevel level) : level_(level) {}
    ~LogBuf() { sync(); }
  };
  LogBuf buf_;
public:
  LogStream(LogLevel level) : std::ostream(&buf_), buf_(level) {
    setf(std::ios::fixed, std::ios::floatfield);
    precision(6);
  }
};

class Logger {
public:
  std::unique_ptr<std::ostream> trace() { return std::make_unique<LogStream>(LogLevel::trace); }
  std::unique_ptr<std::ostream> debug() { return std::make_unique<LogStream>(LogLevel::debug); }
  std::unique_ptr<std::ostream> info()  { return std::make_unique<LogStream>(LogLevel::info); }
  std::unique_ptr<std::ostream> warn()  { return std::make_unique<LogStream>(LogLevel::warning); }
  std::unique_ptr<std::ostream> error() { return std::make_unique<LogStream>(LogLevel::error); }
  std::unique_ptr<std::ostream> fatal() { return std::make_unique<LogStream>(LogLevel::fatal); }
};

extern Logger logger;

#endif
