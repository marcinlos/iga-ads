#ifndef ADS_CONFIG_HPP
#define ADS_CONFIG_HPP

#include <string_view>


namespace ads {

struct version_info {
    std::string_view full;
    int major;
    int minor;
    int patch;
};

auto version() -> version_info;

}

#endif // ADS_CONFIG_HPP
