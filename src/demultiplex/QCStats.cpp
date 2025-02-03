/*
 * Written by Scott Norton
 * 12-Dec-2022
 */

// Copyright © 2025 Merck & Co., Inc., Rahway, NJ, USA and its affiliates. All rights reserved.
// This file is part of DRAGoN.
//
// This source code is licensed under the MIT License found in the
// LICENSE file in the root directory of this source tree.

#include <string>
#include <vector>
#include <unordered_map>
#include <iostream>
#include <boost/log/trivial.hpp>
#include "ReadPair.hpp"
#include "BarcodeMatcher.hpp"
#include "QCStats.hpp"
#include <nlohmann/json.hpp>

Statistics::Statistics(BarcodeMatcher const& matcher, const int max_file, const int &max_mismatch) :
    total(0),
    written_counts(max_file + 1),
    well_counts(matcher.get_numbc() + 1),
    mismatch_counts(matcher.get_numbc() * (max_mismatch + 1)),
    fail_counts((matcher.get_numbc() + 1) * QCFailReason::MAX),
    unmatched_counts(),
    matcher(matcher),
    max_file(max_file),
    max_mismatch(max_mismatch)
{}

void Statistics::record(const ReadPair &rp, const QCFailReason &fail, const int &idx, const int &mismatch, const int &fileidx) {
    std::size_t _idx = idx;
    std::size_t _fileidx = fileidx;
    bool isNoBC = (idx < 0);
    if (isNoBC) {
        _idx = matcher.get_numbc();
        _fileidx = max_file;
    }
    ++total;
    ++fail_counts.at(_idx * QCFailReason::MAX + fail);
    ++written_counts.at(_fileidx);
    ++well_counts.at(_idx);
    if (!isNoBC) {
        ++mismatch_counts.at(_idx * (max_mismatch + 1) + mismatch);
    } else {
        ++unmatched_counts[rp.bcseq];
    }
}

nlohmann::ordered_json Statistics::to_json() const {
    BOOST_LOG_TRIVIAL(debug) << "to_json() invoked";
    nlohmann::ordered_json data = {
        {"read_count", total},
        {"well_counts", nlohmann::json::object()},
        {"written_counts", written_counts},
        {"mismatch_counts", nlohmann::json::array()},
        {"fail_counts", nlohmann::json::object()},
        {"unmatched_counts", unmatched_counts}
    };
    std::vector<std::string> const& bcnames = matcher.get_bcnames();
    int const& numbc = matcher.get_numbc();

    BOOST_LOG_TRIVIAL(debug) << "get well_counts";
    nlohmann::ordered_json& well_counts_j = data["well_counts"];
    for (int i = 0; i < numbc; ++i) {
        well_counts_j[bcnames.at(i)] = well_counts.at(i);
    }
    well_counts_j["noBC"] = well_counts.at(numbc);

    BOOST_LOG_TRIVIAL(debug) << "get mismatch_counts, max_mismatch = " << max_mismatch;
    nlohmann::ordered_json& mismatch_counts_j = data["mismatch_counts"];
    for (int j = 0; j <= max_mismatch; ++j) {
        nlohmann::ordered_json& entry = mismatch_counts_j.emplace_back(nlohmann::ordered_json::object());
        for (int i = 0; i < numbc; ++i) {
            entry[bcnames.at(i)] = mismatch_counts.at(i * (max_mismatch + 1) + j);
        }
    }

    BOOST_LOG_TRIVIAL(debug) << "get fail_counts";
    nlohmann::ordered_json& fail_counts_j = data["fail_counts"];
    for (int j = QCFailReason::PASS; j != QCFailReason::MAX; ++j) {
        nlohmann::ordered_json& entry = fail_counts_j[QCFailNames[j]] = nlohmann::ordered_json::object();
        for (int i = 0; i < numbc; ++i) {
            entry[bcnames.at(i)] = fail_counts.at(i * QCFailReason::MAX + j);
        }
        entry["noBC"] = fail_counts.at(numbc * QCFailReason::MAX + j);
    }

    BOOST_LOG_TRIVIAL(debug) << "get unmatched_counts";
    nlohmann::ordered_json& unmatched_counts_j = data["unmatched_counts"];

    BOOST_LOG_TRIVIAL(debug) << "to_json() done";
    return data;
}

std::ostream &operator<<(std::ostream &strm, Statistics const &self) {
    // Stream as JSON
    strm << std::setw(4) << self.to_json() << std::flush;
    BOOST_LOG_TRIVIAL(debug) << "flushed";
    return strm;
}
