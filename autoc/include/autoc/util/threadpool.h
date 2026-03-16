/**
 * Thread pool that spawns minisim child processes and communicates via TCP sockets.
 */
#pragma once

#include <vector>
#include <queue>
#include <functional>
#include <memory>
#include <atomic>
#include <mutex>
#include <condition_variable>
#include <thread>
#include <chrono>
#include <cstdlib>
#include <cstring>
#include <unistd.h>
#include <sys/wait.h>

#include "autoc/autoc.h"
#include "autoc/util/config.h"
#include "autoc/util/logger.h"
#include "autoc/util/socket_wrapper.h"

using namespace std;

class ThreadPool {
private:
  std::vector<std::thread> threads;
  std::vector<std::unique_ptr<WorkerContext>> worker_contexts;
  std::queue<std::function<void(WorkerContext&)>> tasks;
  std::mutex queue_mutex;
  std::condition_variable condition;
  bool stop;
  std::atomic<int> active_tasks{ 0 };
  std::atomic<int> tasks_in_queue{ 0 };

  void worker(int id, AutocConfig& extraCfg) {
    auto& context = *worker_contexts[id];
    context.workerId = id;

    // Bind a listening socket (port=0 for ephemeral)
    TcpAcceptor acceptor(extraCfg.minisimPortOverride);
    unsigned short port_ = acceptor.port();

    if (extraCfg.minisimPortOverride > 0) {
      *logger.info() << "Now manually launch sim on port " << port_ << endl;
    }
    else {
      std::string subprocess_path = extraCfg.minisimProgram;
      *logger.info() << "Launching: [" << id << "] " << subprocess_path << " " << port_ << endl;

      // Deterministic stagger to avoid simultaneous launches
      unsigned int delayMs = static_cast<unsigned int>((id * 431) % 2000);
      if (delayMs > 0) {
        std::this_thread::sleep_for(std::chrono::milliseconds(delayMs));
      }

      // Fork + exec the minisim process
      pid_t pid = fork();
      if (pid < 0) {
        *logger.error() << "fork() failed for worker " << id << ": " << strerror(errno) << endl;
        return;
      }
      if (pid == 0) {
        // Child process
        setenv("AUTOC_DETERMINISTIC", "1", 1);
        std::string pidStr = std::to_string(getppid());
        std::string idStr = std::to_string(id);
        std::string portStr = std::to_string(port_);
        execl(subprocess_path.c_str(), subprocess_path.c_str(),
              pidStr.c_str(), idStr.c_str(), portStr.c_str(), nullptr);
        // execl only returns on error
        std::cerr << "execl failed for worker " << id << ": " << strerror(errno) << std::endl;
        _exit(1);
      }
      // Parent: child is running, we don't wait on it (it'll exit when socket closes)
      context.childPid = pid;
    }

    // Accept connection from minisim
    context.socket = acceptor.accept();

    while (true) {
      std::function<void(WorkerContext&)> task;
      {
        std::unique_lock<std::mutex> lock(queue_mutex);
        condition.wait(lock, [this] { return stop || !tasks.empty(); });
        if (stop && tasks.empty()) return;
        task = std::move(tasks.front());
        tasks.pop();
        tasks_in_queue--;
        active_tasks++;
      }
      task(context);
    }
  }

public:
  ThreadPool(AutocConfig& extraCfg) : stop(false) {
    worker_contexts.reserve(extraCfg.evalThreads);
    for (int i = 0; i < extraCfg.evalThreads; ++i) {
      worker_contexts.push_back(std::make_unique<WorkerContext>());
      threads.emplace_back(&ThreadPool::worker, this, i, std::ref(extraCfg));
    }
  }

  ~ThreadPool() {
    {
      std::unique_lock<std::mutex> lock(queue_mutex);
      stop = true;
    }
    condition.notify_all();
    for (auto& thread : threads) {
      thread.join();
    }
  }

  template<class F>
  void enqueue(F&& f) {
    {
      std::unique_lock<std::mutex> lock(queue_mutex);
      tasks.emplace([this, f = std::forward<F>(f)](WorkerContext& context) {
        f(context);
        active_tasks--;
        condition.notify_all();
        });
      tasks_in_queue++;
    }
    condition.notify_one();
  }

  void wait_for_tasks() {
    std::unique_lock<std::mutex> lock(queue_mutex);
    condition.wait(lock, [this] {
      return (tasks_in_queue == 0) && (active_tasks == 0);
      });
  }
};
