/*
 * File:   lrucache.hpp
 * Author: Alexander Ponomarev
 *
 * Created on June 20, 2013, 5:09 PM
 *
 * https://raw.githubusercontent.com/lamerman/cpp-lru-cache/master/include/lrucache.hpp
 * Issues: https://github.com/lamerman/cpp-lru-cache/issues
 *
 * Modified by Scott Norton on January 5, 2023
 *
 * Used in this project under a BSD3 license (see LICENSE_THIRD_PARTY)
 *
 */

#ifndef _LRUCACHE_HPP_INCLUDED_
#define    _LRUCACHE_HPP_INCLUDED_

#include <unordered_map>
#include <list>
#include <cstddef>
#include <stdexcept>

namespace cache {

template<typename key_t, typename value_t, size_t _max_size = 0>
class lru_cache {
public:
    typedef typename std::pair<key_t, value_t> key_value_pair_t;
    typedef typename std::list<key_value_pair_t>::iterator list_iterator_t;

    // If max_size == 0, cache size is unlimited
    lru_cache() = default;

    void put(const key_t& key, const value_t& value) {
        auto it = _cache_items_map.find(key);
        // This is the most-recently-used key, so put it in front
        _cache_items_list.push_front(key_value_pair_t(key, value));

        // If this is an update operation, invalidate the existing entry
        if (it != _cache_items_map.end()) {
            _cache_items_list.erase(it->second);
            _cache_items_map.erase(it);
        }
        // Map the key to the front of the list
        _cache_items_map[key] = _cache_items_list.begin();

        // If this makes the cache too big, evict the least-recently-used entry
        if (_max_size != 0 && _cache_items_map.size() > _max_size) {
            _cache_items_map.erase(_cache_items_list.back().first);
            _cache_items_list.pop_back();
        }
    }

    const value_t& get(const key_t& key) {
        auto it = _cache_items_map.find(key);
        if (it == _cache_items_map.end()) {
            throw std::range_error("There is no such key in cache");
        } else {
            // Accessed item is now the most-recently-used, so bump it to the front
            // This does not invalidate it->second
            _cache_items_list.splice(_cache_items_list.begin(), _cache_items_list, it->second);
            return it->second->second;
        }
    }

    bool exists(const key_t& key) const {
        return _cache_items_map.find(key) != _cache_items_map.end();
    }

    size_t size() const {
        return _cache_items_map.size();
    }

private:
    std::list<key_value_pair_t> _cache_items_list;
    std::unordered_map<key_t, list_iterator_t> _cache_items_map;
};

} // namespace cache

#endif    /* _LRUCACHE_HPP_INCLUDED_ */
