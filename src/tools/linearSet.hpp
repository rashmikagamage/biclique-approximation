
#ifndef LINEARSET
#define LINEARSET
#include <tuple>
#include <utility>
#include <cassert>
#include <cstring>

class LinearSet {
private:
    uint32_t * vSet;
    uint32_t * fIndex;
    uint32_t sz;

public:
    LinearSet() { sz = 0; }
    LinearSet(uint32_t sz_) {
        resize(sz_);
    }
    LinearSet(const LinearSet& other) : sz(other.sz), vSet(new uint32_t[other.sz]), fIndex(new uint32_t[other.sz]) {
        std::copy(other.vSet, other.vSet + sz, vSet);
        std::copy(other.fIndex, other.fIndex + sz, fIndex);
    }

    // Copy assignment operator
    LinearSet& operator=(const LinearSet& other) {
        if (this != &other) {
            LinearSet temp(other);  // Use copy constructor
            std::swap(vSet, temp.vSet);
            std::swap(fIndex, temp.fIndex);
            std::swap(sz, temp.sz);
        }
        return *this;
    }
    void resize(uint32_t sz_) {
        sz = sz_;
        vSet = new uint32_t[sz];
        fIndex = new uint32_t[sz];
        for(uint32_t i = 0; i < sz; i++) {
            vSet[i] = fIndex[i] = i;
        }
    }

    ~LinearSet() { delete [] fIndex; delete [] vSet; }
    uint32_t * begin() {
        return vSet;
    }
    uint32_t* end() {
        return vSet + sz;  // Points to one past the last element
    }
    uint32_t operator [] (uint32_t i) {
        // if(i >= g->maxSize()) {
        //     printf("error index\n"); return -1;
        // }
        return vSet[i];
    }
    // changing a vertex to a new position in vSet
    // u is the vertex to be changed, p is the position to be changed
    //if you want get vertex lcoation with id 5 then use fIndex[5] to get the location in vSet
    void changeTo(uint32_t u, uint32_t p) {
        uint32_t pU = fIndex[u];
        std::swap(fIndex[u], fIndex[vSet[p]]);
        std::swap(vSet[pU], vSet[p]);
    }
    //get the u's location in vSet
    uint32_t idx(uint32_t u) {
        return fIndex[u];
    }
//pU and p are indexes 
    void changeToByPos(uint32_t pU, uint32_t p) {
        std::swap(fIndex[vSet[pU]], fIndex[vSet[p]]);
        std::swap(vSet[pU], vSet[p]);
    }

    void copy(uint32_t * p, uint32_t r, uint32_t st = 0) {
        assert(st <= r);
        assert(r <= sz);

        if(r - st >= 4) {
            memcpy(p, vSet + st, sizeof(uint32_t) * (r - st));
        }
        else {
            for(uint32_t i = st; i < r; i++) {
                p[i - st] = vSet[i];
            }
        }
    }
    void clear() {
        delete[] vSet;   // Deallocate memory for vSet
        delete[] fIndex; // Deallocate memory for fIndex
        vSet = nullptr;  // Set to nullptr to avoid dangling pointers
        fIndex = nullptr; // Set to nullptr to avoid dangling pointers
        sz = 0;          // Reset size to zero
    }
    uint32_t size() { return sz; }
};

#endif