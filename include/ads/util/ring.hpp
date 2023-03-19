// SPDX-FileCopyrightText: 2015 - 2023 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef ADS_UTIL_RING_HPP
#define ADS_UTIL_RING_HPP

#include <vector>

namespace ads {
namespace util {
template <typename T>

class ring {
private:
    std::vector<T> buffer_;
    int size_;
    int first_ = 0;

public:
    template <typename... Args>
    explicit ring(int size, Args&&... args)
    : size_(size) {
        for (int i = 0; i < size; ++i) {
            buffer_.emplace_back(std::forward<Args>(args)...);
        }
    }

    int size() const { return size_; }

    T& operator[](int k) {
        int idx = (first_ + k) % size_;
        return buffer_[idx];
    }

    const T& operator[](int k) const {
        int idx = (first_ + k) % size_;
        return buffer_[idx];
    }

    void rotate(int n = 1) { first_ = (first_ + size_ - n) % size_; }
};

}  // namespace util
}  // namespace ads

#endif  // ADS_UTIL_RING_HPP
