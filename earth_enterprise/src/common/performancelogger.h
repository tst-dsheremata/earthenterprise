// Copyright 2018 The Open GEE Contributors
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef PERFORMANCELOGGER_H
#define PERFORMANCELOGGER_H

#include <iostream>

#include <assert.h>
#include <pthread.h>
#include <string>
#include <time.h>

#include "common/timeutils.h"
#include "common/khThread.h"  // defines khMutex


#ifdef LOG_PERFORMANCE


namespace performance_logger {

// We need a mutex for performance logging, but we've intstrumented the
// khMutex class hierarchy, and we don't want to see log entries for logging
// operations, so we have a non-instrumented implementation of the mutex
// classes here. :P
class khMutexBase {
 public:
  pthread_mutex_t mutex;

  void lock(void) { Lock(); }
  void Lock(void);
  void unlock(void) { Unlock(); }
  void Unlock(void);
  bool trylock(void) { return TryLock(); }
  bool TryLock(void);
};

// Simple non-recursive, non-checking mutex
class khMutex : public khMutexBase
{
 public:
  khMutex(void);
  ~khMutex(void);
};

class khLockGuard {
  khMutexBase &mutex;
 public:
  khLockGuard(khMutexBase &mutex_) : mutex(mutex_) {
    mutex.Lock();
  }
  ~khLockGuard(void) {
    mutex.Unlock();
  }
};


/*
 * Singleton class for logging event performance. This class is intended for
 * performance debugging.
 */
class PerformanceLogger {
  public:
    static PerformanceLogger& instance() {
      static PerformanceLogger _instance;

      return _instance;
    }

    // Initialize member variables in constructor:
    PerformanceLogger()
    : write_mutex()
    {}

    void logTiming(
        const std::string & operation, // The operation being timed
        const std::string & object,    // The object that the operation is performed on
        const timespec startTime,      // The start time of the operation
        const timespec endTime,        // The end time of the operation
        const size_t size = 0);        // The size of the object, if applicable
    // void logThread(/* Add params as needed*/) {
    //   // Format thread data
    //   doNotify("", /* the formatted data */, "threadfile.csv");
    // }
  private:
    khMutex write_mutex;
    PerformanceLogger() {}
    void doNotify(std::string data, std::string fileName);
};

/*
 * A convenience class for timing a block of code. Timing begins when an
 * instance of this class is created and ends when the user calls "end" or the
 * instance goes out of scope.
 */
template <class PerfLoggerCls>
class BlockPerformanceLogger {
  public:
    BlockPerformanceLogger(
        const std::string & operation,
        const std::string & object,
        const size_t size = 0) :
      operation(operation),
      object(object),
      size(size),
      startTime(getime::getMonotonicTime()),
      ended(false) {}
    void end() {
      if (!ended) {
        ended = true;
        const timespec endTime = getime::getMonotonicTime();
        PerfLoggerCls::instance().logTiming(operation, object, startTime, endTime, size);
      }
    }
    ~BlockPerformanceLogger() {
      end();
    }
  private:
    const std::string operation;
    const std::string object;
    const size_t size;
    const timespec startTime;
    bool ended;
};

}

// Programmers should use the macros below to time code instead of using the
// classes above directly. This makes it easy to use compile time flags to
// exclude the time code.
//
// For the macros below, the variadic arguments include the object name
// (required) and the size (optional).
//
// If you use multiple performance logging statements in the same scope you
// must give each one a unique name.
#define BEGIN_PERF_LOGGING(name, op, ...) \
  performance_logger::BlockPerformanceLogger<performance_logger::PerformanceLogger> name(op, __VA_ARGS__)
#define END_PERF_LOGGING(name) name.end()

#else

#define BEGIN_PERF_LOGGING( ... )
#define END_PERF_LOGGING( ... )

#endif  // LOG_PERFORMANCE

#endif // PERFORMANCELOGGER_H
