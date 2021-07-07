// SPDX-FileCopyrightText: 2015 - 2021 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

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

}  // namespace ads
