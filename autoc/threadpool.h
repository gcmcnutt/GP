/**
 * This threadpool spawns child processes corresponding to a task thread, and uses sockets for communication
 */
#pragma once

#include <boost/asio.hpp>
#include <boost/thread.hpp>
#include <vector>
#include <queue>
#include <functional>
#include <memory>
#include <atomic>

#include "autoc.h"
#include "logger.h"

using namespace std;

class ThreadPool {
private:
  std::vector<boost::thread> threads;
  std::vector<std::unique_ptr<WorkerContext>> worker_contexts;
  std::queue<std::function<void(WorkerContext&)>> tasks;
  boost::mutex queue_mutex;
  boost::condition_variable condition;
  bool stop;
  boost::asio::io_service io_service;
  std::atomic<int> active_tasks{ 0 };
  std::atomic<int> tasks_in_queue{ 0 };

  void worker(int id, ExtraConfig& extraCfg) {
    auto& context = *worker_contexts[id];

    // Launch subprocess with port as argument
    tcp::acceptor acceptor_(io_service, tcp::endpoint(tcp::v4(), extraCfg.minisimPortOverride));
    unsigned short port_ = acceptor_.local_endpoint().port();
    if (extraCfg.minisimPortOverride > 0) {
      *logger.info() << "Now manually launch sim on port " << port_ << endl;
    }
    else {
      std::string subprocess_path = extraCfg.minisimProgram;
      std::vector<std::string> args = { std::to_string(id), std::to_string(port_) };
      *logger.info() << "Launching: [" << id << "] " << subprocess_path << " " << port_ << endl;

      context.child_process = boost::process::child(
        subprocess_path,
        args
        // bp::std_out > bp::null,
        // bp::std_err > bp::null,
      );
      context.child_process.detach();
    }

    // Set up socket connection
    context.socket = std::make_unique<boost::asio::ip::tcp::socket>(io_service);
    acceptor_.accept(*context.socket);

    while (true) {
      std::function<void(WorkerContext&)> task;
      {
        boost::unique_lock<boost::mutex> lock(queue_mutex);
        condition.wait(lock, [this] { return stop || !tasks.empty(); });
        if (stop && tasks.empty()) return;
        task = std::move(tasks.front());
        tasks.pop();
        tasks_in_queue--;
        active_tasks++;
      }
      task(context);  // Pass the context to the task
    }
  }

public:
  ThreadPool(ExtraConfig& extraCfg) : stop(false) {
    worker_contexts.reserve(extraCfg.evalThreads);
    for (int i = 0; i < extraCfg.evalThreads; ++i) {
      worker_contexts.push_back(std::make_unique<WorkerContext>());
      threads.emplace_back(&ThreadPool::worker, this, i, extraCfg);
    }
  }

  ~ThreadPool() {
    {
      boost::unique_lock<boost::mutex> lock(queue_mutex);
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
      boost::unique_lock<boost::mutex> lock(queue_mutex);
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
    boost::unique_lock<boost::mutex> lock(queue_mutex);
    condition.wait(lock, [this] {
      return (tasks_in_queue == 0) && (active_tasks == 0);
      });
  }
};