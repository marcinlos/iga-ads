// SPDX-FileCopyrightText: 2015 - 2023 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#include "ads/executor/galois.hpp"

#include "ads/config.hpp"

#ifdef ADS_USE_GALOIS

#    include <galois/runtime/Statistics.h>

namespace ads {

galois_executor::galois_executor(int threads) {
    galois::runtime::setStatFile("/dev/null");
    thread_count(threads);
}

void galois_executor::thread_count(int threads) {
    galois::setActiveThreads(threads);
}

}  // namespace ads

#endif  // defined(ADS_USE_GALOIS)
