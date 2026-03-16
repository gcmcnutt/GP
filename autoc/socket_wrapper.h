#pragma once

#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <netdb.h>
#include <unistd.h>

#include <cerrno>
#include <cstring>
#include <memory>
#include <stdexcept>
#include <string>

// A connected TCP socket (client or accepted connection).
// Move-only, RAII: closes the fd in destructor.
class TcpSocket {
public:
    TcpSocket() : fd_(-1) {}

    explicit TcpSocket(int fd) : fd_(fd) {}

    ~TcpSocket() { close_fd(); }

    // Move
    TcpSocket(TcpSocket&& other) noexcept : fd_(other.fd_) { other.fd_ = -1; }
    TcpSocket& operator=(TcpSocket&& other) noexcept {
        if (this != &other) {
            close_fd();
            fd_ = other.fd_;
            other.fd_ = -1;
        }
        return *this;
    }

    // No copy
    TcpSocket(const TcpSocket&) = delete;
    TcpSocket& operator=(const TcpSocket&) = delete;

    // Connect to a remote host:port.
    void connect(const std::string& host, unsigned short port) {
        struct addrinfo hints{};
        hints.ai_family   = AF_INET;
        hints.ai_socktype = SOCK_STREAM;

        struct addrinfo* res = nullptr;
        int rc = ::getaddrinfo(host.c_str(), std::to_string(port).c_str(), &hints, &res);
        if (rc != 0) {
            throw std::runtime_error(std::string("TcpSocket::connect getaddrinfo: ") + gai_strerror(rc));
        }

        int fd = ::socket(res->ai_family, res->ai_socktype, res->ai_protocol);
        if (fd < 0) {
            ::freeaddrinfo(res);
            throw std::runtime_error(std::string("TcpSocket::connect socket: ") + std::strerror(errno));
        }

        if (::connect(fd, res->ai_addr, res->ai_addrlen) < 0) {
            ::freeaddrinfo(res);
            ::close(fd);
            throw std::runtime_error(std::string("TcpSocket::connect connect: ") + std::strerror(errno));
        }

        ::freeaddrinfo(res);
        close_fd();
        fd_ = fd;
    }

    // Read exactly n bytes, looping on partial reads.
    void read(void* buf, std::size_t n) {
        auto* ptr = static_cast<char*>(buf);
        std::size_t remaining = n;
        while (remaining > 0) {
            ssize_t bytes = ::recv(fd_, ptr, remaining, MSG_WAITALL);
            if (bytes < 0) {
                throw std::runtime_error(std::string("TcpSocket::read recv: ") + std::strerror(errno));
            }
            if (bytes == 0) {
                throw std::runtime_error("TcpSocket::read: connection closed by peer");
            }
            ptr       += bytes;
            remaining -= static_cast<std::size_t>(bytes);
        }
    }

    // Write exactly n bytes, looping on partial writes.
    void write(const void* buf, std::size_t n) {
        const auto* ptr = static_cast<const char*>(buf);
        std::size_t remaining = n;
        while (remaining > 0) {
            ssize_t bytes = ::send(fd_, ptr, remaining, MSG_NOSIGNAL);
            if (bytes < 0) {
                throw std::runtime_error(std::string("TcpSocket::write send: ") + std::strerror(errno));
            }
            if (bytes == 0) {
                throw std::runtime_error("TcpSocket::write: connection closed by peer");
            }
            ptr       += bytes;
            remaining -= static_cast<std::size_t>(bytes);
        }
    }

    int fd() const { return fd_; }

private:
    int fd_;

    void close_fd() {
        if (fd_ >= 0) {
            ::close(fd_);
            fd_ = -1;
        }
    }
};

// Server-side TCP acceptor. Binds immediately on construction.
// Pass port=0 to let the OS choose an ephemeral port.
// Move-only, RAII: closes the listening fd in destructor.
class TcpAcceptor {
public:
    explicit TcpAcceptor(unsigned short port = 0) : fd_(-1), port_(0) {
        int fd = ::socket(AF_INET, SOCK_STREAM, 0);
        if (fd < 0) {
            throw std::runtime_error(std::string("TcpAcceptor: socket: ") + std::strerror(errno));
        }

        int opt = 1;
        if (::setsockopt(fd, SOL_SOCKET, SO_REUSEADDR, &opt, sizeof(opt)) < 0) {
            ::close(fd);
            throw std::runtime_error(std::string("TcpAcceptor: setsockopt: ") + std::strerror(errno));
        }

        struct sockaddr_in addr{};
        addr.sin_family      = AF_INET;
        addr.sin_addr.s_addr = INADDR_ANY;
        addr.sin_port        = htons(port);

        if (::bind(fd, reinterpret_cast<struct sockaddr*>(&addr), sizeof(addr)) < 0) {
            ::close(fd);
            throw std::runtime_error(std::string("TcpAcceptor: bind: ") + std::strerror(errno));
        }

        if (::listen(fd, SOMAXCONN) < 0) {
            ::close(fd);
            throw std::runtime_error(std::string("TcpAcceptor: listen: ") + std::strerror(errno));
        }

        // Discover the actual port (important when port==0).
        struct sockaddr_in bound{};
        socklen_t len = sizeof(bound);
        if (::getsockname(fd, reinterpret_cast<struct sockaddr*>(&bound), &len) < 0) {
            ::close(fd);
            throw std::runtime_error(std::string("TcpAcceptor: getsockname: ") + std::strerror(errno));
        }

        fd_    = fd;
        port_  = ntohs(bound.sin_port);
    }

    ~TcpAcceptor() { close_fd(); }

    // Move
    TcpAcceptor(TcpAcceptor&& other) noexcept : fd_(other.fd_), port_(other.port_) {
        other.fd_   = -1;
        other.port_ = 0;
    }
    TcpAcceptor& operator=(TcpAcceptor&& other) noexcept {
        if (this != &other) {
            close_fd();
            fd_         = other.fd_;
            port_       = other.port_;
            other.fd_   = -1;
            other.port_ = 0;
        }
        return *this;
    }

    // No copy
    TcpAcceptor(const TcpAcceptor&) = delete;
    TcpAcceptor& operator=(const TcpAcceptor&) = delete;

    // The port this acceptor is bound to (useful when constructed with port=0).
    unsigned short port() const { return port_; }

    // Block until a client connects; returns the connected socket.
    std::unique_ptr<TcpSocket> accept() {
        struct sockaddr_in client{};
        socklen_t len = sizeof(client);
        int client_fd = ::accept(fd_, reinterpret_cast<struct sockaddr*>(&client), &len);
        if (client_fd < 0) {
            throw std::runtime_error(std::string("TcpAcceptor::accept: ") + std::strerror(errno));
        }
        return std::make_unique<TcpSocket>(client_fd);
    }

private:
    int           fd_;
    unsigned short port_;

    void close_fd() {
        if (fd_ >= 0) {
            ::close(fd_);
            fd_ = -1;
        }
    }
};
