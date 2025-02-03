/*
 * Written by Scott Norton
 * 20 Dec 2022
 */

// Copyright © 2025 Merck & Co., Inc., Rahway, NJ, USA and its affiliates. All rights reserved.
// This file is part of DRAGoN.
//
// This source code is licensed under the MIT License found in the
// LICENSE file in the root directory of this source tree.

#ifndef DRUGSEQ_UTIL_H
#define DRUGSEQ_UTIL_H

#include <string>
#include <vector>
#include <chrono>
#include <iomanip>
#include <ranges>

template <>
struct std::hash<std::pair<std::string, std::string>> {
    unsigned long long operator()(const std::pair<std::string, std::string>& x) const {
        const std::hash<std::string> hasher {};
        return (hasher(x.first) << 1) ^ hasher(x.second);
    }
};

namespace drugseq {
    // Join a vector of strings with a delimiter
    template<typename InputIt>
    std::string strjoin(InputIt begin, InputIt end, std::string const& j) {
        std::string result {};
        for (auto it = begin;; result += j) {
            result += *it;
            if (++it == end) {
                break;
            }
        }
        return result;
    }

    static inline std::string strjoin(std::vector<std::string> const& v, std::string const& j) {
        return strjoin(v.cbegin(), v.cend(), j);
    }

    template <std::ranges::forward_range RangeT, typename CharT = char>
        requires std::is_same_v<std::ranges::range_value_t<RangeT>, std::basic_string<CharT>>
    std::basic_string<CharT> strjoin(RangeT&& rng, std::basic_string<CharT> const& delim) {
        std::basic_string<CharT> result {};
        for (std::basic_string<CharT> sep {}; auto word : rng) {
            result += sep + word;
            sep = delim;
        }
        return result;
    }

    template <std::ranges::forward_range RangeT, typename CharT = char>
        requires std::is_same_v<std::ranges::range_value_t<RangeT>, std::basic_string<CharT>>
    std::basic_string<CharT> strjoin(RangeT&& rng, const CharT *delim) {
        std::basic_string<CharT> result {};
        for (const CharT *sep {""}; auto word : rng) {
            result += sep + word;
            sep = delim;
        }
        return result;
    }

    // Join a vector of strings with a space delimeter. Escapes double-quotes and backslashes. If a space exists in a substring, it will be double-quoted (unescaped).
    template<typename InputIt>
    std::string shlexjoin(InputIt begin, InputIt end) {
        std::stringstream result {};
        for (auto it = begin;; result << ' ') {
            const std::string& s = *it++;
            if (s.find(' ') != std::string::npos) {
                result << std::quoted(s);
            } else {
                result << s;
            }
            if (it == end) {
                break;
            }
        }
        return result.str();
    }

    static inline std::string shlexjoin(std::vector<std::string> const& v) {
        return shlexjoin(v.cbegin(), v.cend());
    }

    // Split a string into a vector of strings
    template<typename OutputIt>
    OutputIt strsplit(OutputIt out, std::string const& s, std::string const& d) {
        std::size_t pos = 0, nextpos;
        nextpos = s.find(d);
        while (pos != std::string::npos) {
            if (nextpos == std::string::npos) {
                *out++ = s.substr(pos);
                break;
            }
            *out++ = s.substr(pos, nextpos - pos);
            // consume adjacent delimiters
            do {
                pos = nextpos + d.length();
                nextpos = s.find(d, pos);
            } while (nextpos != std::string::npos && nextpos == pos);
        }
        return out;
    }

    static inline std::vector<std::string> strsplit(std::string const& s, std::string const& d) {
        std::vector<std::string> ret = {};
        strsplit(std::back_inserter(ret), s, d);
        return ret;
    }

    // IOMANIP object to print the current local system time
    template<typename _CharT = char>
    std::_Put_time<_CharT> put_time() {
        const std::time_t t_c = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
        return std::put_time(std::localtime(&t_c), "%F %T");
    }
}

#endif //DRUGSEQ_UTIL_H
