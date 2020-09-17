#pragma once
#include <vector>
#include <stdexcept>

#define lchild(x) ((x << 1) + 1)
#define rchild(x) ((x << 1) + 2)
#define parent(x) (x ? (x-1) >> 1 : 0)

// If desired we could implement a custom compare function 
// for comparing non-numerical keys

// For creaign a max heap, use
// max_heap_comp<K> as the Compare template arg
template<typename K>
struct max_heap_comp {
    bool operator()(const K &a, const K &b) const {
        return a > b;
    }
};

// For min heap
template<typename K>
struct min_heap_comp {
    bool operator()(const K &a, const K &b) const {
        return a < b;
    }
};

template<typename K, typename V, typename Compare>
class BinaryHeap
{
private:
    std::vector<std::pair<K, V>> m_heap;
    Compare m_comparator;

    void BubbleUp(size_t idx) {
        size_t p = parent(idx);
        while (m_comparator(m_heap[idx].first, m_heap[p].first)) {
            std::swap(m_heap[idx], m_heap[p]);
            idx = p;
            p = parent(p);
        }
    }

    void BubbleDown(size_t idx) {
        size_t n = m_heap.size();
        size_t l = lchild(idx);
        size_t r = rchild(idx);
        size_t extreme = idx;
        while (l < n) {
            if (m_comparator(m_heap[l].first, m_heap[idx].first)) {
                extreme = l;
            }

            if (r < n && m_comparator(m_heap[r].first, m_heap[l].first)) {
                extreme = r;
            }
            
            if (extreme == idx) {
                break;
            }
            else {
                std::swap(m_heap[idx], m_heap[extreme]);
                idx = extreme;
                l = lchild(idx);
                r = rchild(idx);
            }
        }
    }

public:
    BinaryHeap() { }
    BinaryHeap(const std::vector<K> &keys, const std::vector<V> &vals) {
        if (keys.size() != vals.size()) {
            throw "Must have same number of keys and values.";
        }
        else {
            for (size_t i = 0; i < keys.size(); i++) {
                Insert(keys[i], vals[i]);
            }
        }
    }

    void Insert(const K &key, const V &val) {
        size_t n = m_heap.size();
        m_heap.push_back(std::make_pair(key, val));
        BubbleUp(n);
    }

    const V &GetTop() const {
        if (m_heap.size() > 0) {
            return m_heap[0].second;
        }
        else {
            throw std::length_error("Heap empty");
        }
    }

    V Pop() {
        size_t n = m_heap.size();
        if (n > 0) {
            V result = m_heap[0].second;
            std::swap(m_heap[0], m_heap[n - 1]);
            m_heap.pop_back();
            BubbleDown(0);
            return result;
        }
        else {
            throw std::length_error("Heap empty");
        }
    }

    void ChangeKey(size_t idx, const K &new_key) {
        if (m_comparator(new_key, m_heap[idx].first)) {
            m_heap[idx].first = new_key;
            BubbleUp(idx);
        }
        else {
            m_heap[idx].first = new_key;
            BubbleDown(idx);
        }
    }
};