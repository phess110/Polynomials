#pragma once
#include <vector>
#include <stdexcept>

#define lchild(x) ((x << 1) + 1)
#define rchild(x) ((x << 1) + 2)
#define parent(x) (x ? (x-1) >> 1 : 0)

// If desired we could implement a custom compare function 
// for comparing non-numerical keys

/*! For creating a max heap, use
     max_heap_comp<K> as the Compare template arg
*/
template<typename K>
struct max_heap_comp {
    bool operator()(const K &a, const K &b) const {
        return a > b;
    }
};

//! For min heap
template<typename K>
struct min_heap_comp {
    bool operator()(const K &a, const K &b) const {
        return a < b;
    }
};

/*! 
    \brief Templated class for a binary heap data structure
    \tparam K the key data type
    \tparam V the value data type
    \tparam Compare A struct containing a callable comparision operator

    \details The heap is implemented using a vector of key-value pairs
*/
template<typename K, typename V, typename Compare>
class BinaryHeap
{
private:
    std::vector<std::pair<K, V>> m_heap;
    Compare m_comparator;


    /*!
        Maintain heap property by bubbling the node at index i upward
    */
    void BubbleUp(size_t idx) {
        size_t p = parent(idx);
        while (m_comparator(m_heap[idx].first, m_heap[p].first)) {
            std::swap(m_heap[idx], m_heap[p]);
            idx = p;
            p = parent(p);
        }
    }

    /*!
        Maintain heap property by bubbling the node at index i downward
    */
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

    /*!
        \brief Constructs a binary heap given two lists of keys and corresponding values
        \throw user-error If the two lists do not have the same number of values
    */
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

    /*!
        \brief Inserts the given key-value pair into the heap
        \param [in] key 
        \param [in] val
    */
    void Insert(const K &key, const V &val) {
        size_t n = m_heap.size();
        m_heap.push_back(std::make_pair(key, val));
        BubbleUp(n);
    }

    /*!
        \return The element at the top of the heap
        \throw std::length_error If the heap is empty
    */
    const V &GetTop() const {
        if (m_heap.size() > 0) {
            return m_heap[0].second;
        }
        else {
            throw std::length_error("Heap empty");
        }
    }

    /*!
        \brief Removes the top element from the heap and returns its value
        \return The element at the top of the heap
        \throw std::length_error If the heap is empty
    */
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

    /*!
        \param [in] idx Location of the value to be changed
        \param [in] new_key The value the key is to be updated to
    */
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