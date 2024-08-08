/**
 * Boost logger meant to look like an output stream
 */
#include "logger.h"

int LogBuf::overflow(int c) {
  if (c != EOF) {
    m_buffer += static_cast<char>(c);
    if (c == '\n') {
      BOOST_LOG_SEV(m_logger, m_level) << m_buffer;
      m_buffer.clear();
    }
  }
  return c;
}

int LogBuf::sync() {
  if (!m_buffer.empty()) {
    BOOST_LOG_SEV(m_logger, m_level) << m_buffer;
    m_buffer.clear();
  }
  return 0;
}

LogBuf::LogBuf(boost::log::sources::severity_logger<boost::log::trivial::severity_level>& logger,
  boost::log::trivial::severity_level level)
  : m_logger(logger), m_level(level) {}

LogBuf::~LogBuf() {
  sync();
}


LogStream::LogStream(boost::log::sources::severity_logger<boost::log::trivial::severity_level>& logger,
  boost::log::trivial::severity_level level)
  : std::ostream(&m_buf), m_buf(logger, level) {}

Logger::Logger() {}

std::unique_ptr<std::ostream> Logger::trace() { return std::make_unique<LogStream>(m_logger, boost::log::trivial::trace); }
std::unique_ptr<std::ostream> Logger::debug() { return std::make_unique<LogStream>(m_logger, boost::log::trivial::debug); }
std::unique_ptr<std::ostream> Logger::info() { return std::make_unique<LogStream>(m_logger, boost::log::trivial::info); }
std::unique_ptr<std::ostream> Logger::warn() { return std::make_unique<LogStream>(m_logger, boost::log::trivial::warning); }
std::unique_ptr<std::ostream> Logger::error() { return std::make_unique<LogStream>(m_logger, boost::log::trivial::error); }
std::unique_ptr<std::ostream> Logger::fatal() { return std::make_unique<LogStream>(m_logger, boost::log::trivial::fatal); }
