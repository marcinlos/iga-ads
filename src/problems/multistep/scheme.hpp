#ifndef ADS_PROBLEMS_MULTISTEP_SCHEME_HPP_
#define ADS_PROBLEMS_MULTISTEP_SCHEME_HPP_

#include <vector>
#include <iterator>
#include <iostream>
#include <sstream>


namespace ads {

struct scheme {
    int s;
    std::vector<double> as;
    std::vector<double> bs;
};

std::ostream& operator << (std::ostream& out, const scheme& s);

ads::scheme get_scheme(const std::string& name);

scheme parse_scheme(const std::string& text);

}

#endif /* ADS_PROBLEMS_MULTISTEP_SCHEME_HPP_ */
