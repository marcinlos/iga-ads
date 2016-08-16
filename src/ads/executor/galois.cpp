#include "ads/executor/galois.hpp"

namespace ads {

    galois_executor::galois_executor(int threads) {
        thread_count(threads);
    }

    void galois_executor::thread_count(int threads) {
        Galois::setActiveThreads(threads);
    }
}
