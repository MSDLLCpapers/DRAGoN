#ifndef LIBDRUGSEQ_GTF_HPP
#define LIBDRUGSEQ_GTF_HPP

// Written by Scott Norton
// 12 Sep 2023

// Copyright © 2025 Merck & Co., Inc., Rahway, NJ, USA and its affiliates. All rights reserved.
// This file is part of DRAGoN.
//
// This source code is licensed under the MIT License found in the
// LICENSE file in the root directory of this source tree.

#include <string>
#include <iostream>
#include <fstream>
#include <map>
#include <variant>
#include <vector>
#include <ranges>

namespace drugseq {
    // Represents a single record in a GTF file
    struct GtfRecord {
        std::string seqname;
        std::string source;
        std::string feature;
        long long start;
        long long end;
        double score;
        std::string strand;
        unsigned char frame;
        std::map<std::string, std::variant<std::string, double>> attr;

        GtfRecord() = default;
        friend std::istream& operator>>(std::istream& strm, GtfRecord& me);
        friend std::ostream& operator<<(std::ostream& strm, GtfRecord const& me);
    };

    // Specialized GTF reader
    class Gtf {
        std::vector<GtfRecord> records;
        std::ifstream ifstream;
    public:
        explicit Gtf(const std::string& filename) : ifstream(filename) {}
        Gtf(const Gtf& other) = delete;
        Gtf(Gtf&& other) : ifstream(std::move(other.ifstream)) {}

        // A non-storing version of load (defined below)
        auto lazy_load(const std::string& feature = "gene") {
            auto is_view = std::ranges::istream_view<GtfRecord>(ifstream);
            auto filter_view = std::ranges::filter_view(is_view, [=](const GtfRecord& rec) -> bool {
                return rec.feature == feature;
            });
            return filter_view;
        }

        auto lazy_load(const std::string& feature, const std::string& metafeature) {
            auto is_view = std::ranges::istream_view<GtfRecord>(ifstream);
            auto filter_view = std::ranges::filter_view(is_view, [=](const GtfRecord& rec) -> bool {
                return rec.feature == feature || rec.feature == metafeature;
            });
            return filter_view;
        }

        Gtf& load(const std::string& feature = "gene") {
            // would like to use std::ranges::to but that's not available
            // in g++ yet...
            records.clear();
            for (auto x : lazy_load(feature)) {
                records.push_back(x);
            }
            return *this;
        }

        // Extracts a range with gene_ids
        std::vector<std::string> get_gene_names(std::string const& attr_name = "gene_id") const {
            auto view = std::ranges::transform_view(records, [&](const GtfRecord& rec) -> std::string {
                try {
                    return std::get<std::string>(rec.attr.at(attr_name));
                } catch (std::out_of_range e) {
                    std::cerr << "terminate called after throwing an instance of 'std::out_of_range'\n"
                                "what():  map::at\n"
                                "offending record: " << rec << "\n";
                    std::exit(1);
                }
            });
            return {view.begin(), view.end()};
        }

        const std::vector<GtfRecord>& get_records() const {return records;}
        operator bool() const {
            return !records.empty();
        }

        bool operator!() const {
            return records.empty();
        }
    };
}

#endif //LIBDRUGSEQ_GTF_HPP
