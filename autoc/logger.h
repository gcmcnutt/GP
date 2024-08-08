#ifndef LOGGER_H
#define LOGGER_H

/**
 * Boost logger meant to look like an output stream
 */
#include <boost/log/trivial.hpp>
#include <boost/log/sources/severity_logger.hpp>
#include <sstream>

class LogBuf : public std::streambuf {
private:
  boost::log::sources::severity_logger<boost::log::trivial::severity_level>& m_logger;
  boost::log::trivial::severity_level m_level;
  std::string m_buffer;

protected:
  virtual int overflow(int c);

  virtual int sync();

public:
  LogBuf(boost::log::sources::severity_logger<boost::log::trivial::severity_level>& logger,
    boost::log::trivial::severity_level level);

  ~LogBuf();
};


class LogStream : public std::ostream {
private:
  LogBuf m_buf;

public:
  LogStream(boost::log::sources::severity_logger<boost::log::trivial::severity_level>& logger,
    boost::log::trivial::severity_level level);
};


class Logger {
private:
  boost::log::sources::severity_logger<boost::log::trivial::severity_level> m_logger;

public:
  Logger();

  std::unique_ptr<std::ostream> trace();
  std::unique_ptr<std::ostream> debug();
  std::unique_ptr<std::ostream> info();
  std::unique_ptr<std::ostream> warn();
  std::unique_ptr<std::ostream> error();
  std::unique_ptr<std::ostream> fatal();
};

extern Logger logger;

#endif