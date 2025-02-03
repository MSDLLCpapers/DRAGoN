
// Copyright © 2025 Merck & Co., Inc., Rahway, NJ, USA and its affiliates. All rights reserved.
// This file is part of DRAGoN.
//
// This source code is licensed under the MIT License found in the
// LICENSE file in the root directory of this source tree.

#include <sstream>
#include <iomanip>
#include <limits>
#include <array>
#include <cmath>
#include <drugseq/Gtf.hpp>

namespace drugseq {

std::istream& operator>>(std::istream& strm, GtfRecord& me) {
    std::string line {};
    std::array<std::string, 9> fields {};
    do {
        if (!std::getline(strm, line)) {
            return strm;
        }
    } while (line.empty() || line[0] == '#');
    line = line.substr(0, line.find('#'));
    std::stringstream line_ss (line);
    for (int i = 0; i < 9; ++i) {
        if (!std::getline(line_ss, fields[i], '\t') && i != 8) {
            throw std::out_of_range("invalid gtf entry");
        }
    }
    me.seqname = fields[0];
    me.source = fields[1];
    me.feature = fields[2];
    me.start = std::stoll(fields[3]) - 1;  // represent as 0-based internally
    me.end = std::stoll(fields[4]) - 1;  // represent as 0-based internally
    me.score = fields[5] == "." ? std::numeric_limits<double>::quiet_NaN() : std::stod(fields[5]);
    me.strand = fields[6] == "." ? "" : fields[6];
    me.frame = fields[7] == "." ? 0xFF : static_cast<unsigned char>(std::stoul(fields[7]));
    me.attr.clear();
    std::string attr_s;
    std::stringstream attr_ss(fields[8]);
    while (std::getline(attr_ss, attr_s, ';')) {
        std::string key;
        std::stringstream key_ss(attr_s);
        key_ss >> key;
        char sp;
        key_ss.read(&sp, 1);
        if (key_ss.peek() == '"') {
            std::string value;
            key_ss >> std::quoted(value);
            me.attr[key] = value;
        } else {
            double value;
            key_ss >> value;
            me.attr[key] = value;
        }
    }
    return strm;
}

std::ostream& operator<<(std::ostream& strm, GtfRecord const& me) {
    strm << me.seqname << "\t" << me.source << "\t" << me.feature << "\t" << me.start << "\t" << me.end << "\t";
    if (std::isnan(me.score)) {
        strm << ".";
    } else {
        strm << me.score;
    }
    strm << "\t";
    if (me.strand.empty()) {
        strm << ".";
    } else {
        strm << me.strand;
    }
    strm << "\t";
    if (me.frame == 0xFF) {
        strm << ".";
    } else {
        strm << me.frame;
    }
    strm << "\t";
    for (const auto [key, value] : me.attr) {
        try {
            strm << key << " " << std::quoted(std::get<std::string>(value)) << "; ";
        } catch (std::bad_variant_access) {
            strm << key << " " << std::get<double>(value) << "; ";
        }
    }
    return strm;
}
}
