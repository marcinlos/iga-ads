#include "ads/executor/galois.hpp"
#include <galois/runtime/Statistics.h>

namespace ads {

    galois_executor::galois_executor(int threads) {
        galois::runtime::setStatFile("/dev/null");
        thread_count(threads);
    }

    void galois_executor::thread_count(int threads) {
        galois::setActiveThreads(threads);
    }
}
