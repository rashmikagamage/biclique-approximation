#ifndef HOPSTOTCH
#define HOPSTOTCH

#ifndef _WIN32
#include <x86intrin.h>
#else
#include <intrin.h>
#endif

#include <immintrin.h>
#include <stdint.h>

#include <cstdint>
#include <cstdio>
#include <cstring>
#include <utility>
#ifdef NDEBUG
#undef NDEBUG
#endif
#include <cassert>

#define EMPTY 0xffffffff

constexpr uint32_t L2_CACHE_LINE_SIZE = 64;
constexpr uint32_t H = L2_CACHE_LINE_SIZE / sizeof(uint32_t);

#ifdef __AVX512__
#define setFunction _mm512_set1_epi32
#define cmpFunction _mm512_cmpeq_epi32_mask
#define loadFunction _mm512_loadu_si512
#endif

#ifdef __AVX512__
const __m256i eightEMPTY = _mm256_set1_epi32(EMPTY);
const __m512i sixteenEMPTY = _mm512_set1_epi32(EMPTY);
#endif

class hopstotchHash {
   public:
    uint32_t *v;
    uint32_t n;  // Number of elements
    uint32_t roundN, preZero;
    uint32_t hashKey = 1e9 + 7;

   public:
#ifndef __SSE__
#define __SSE__
#endif

#ifdef __SSE__
#ifndef __AVX512__
    bool containSIMD(uint32_t u) {
        uint32_t hashV = hash(u);
        uint32_t p = hashV;

        __m256i eightU = _mm256_set1_epi32(u);

        auto address = reinterpret_cast<const __m256i *>(v + p);
        __m256i eightNextV = _mm256_loadu_si256(address);
        __m256i cmpRes = _mm256_cmpeq_epi32(eightU, eightNextV);
        auto msk = _mm256_movemask_epi8(cmpRes);
#ifdef _WIN32
        if (__popcnt(msk) > 0) return true;
#else
        if (_popcnt32(msk) > 0) return true;
#endif
        eightNextV = _mm256_loadu_si256(address + 1);
        cmpRes = _mm256_cmpeq_epi32(eightU, eightNextV);
        msk = _mm256_movemask_epi8(cmpRes);
#ifdef _WIN32
        if (__popcnt(msk) > 0) return true;
#else
        if (_popcnt32(msk) > 0) return true;
#endif

        return false;
    }
#else
    bool containSIMD(uint32_t u) {
        uint32_t hashV = hash(u);
        uint32_t p = hashV;

        __m256i eightU = _mm256_set1_epi32(u);

        __m256i eightNextV = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(v + p));
        __mmask8 cmpRes = _mm256_cmpeq_epi32_mask(eightU, eightNextV);
        if (cmpRes) return true;

        eightNextV = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(v + p + 8));
        cmpRes = _mm256_cmpeq_epi32_mask(eightU, eightNextV);
        if (cmpRes) return true;

        return false;
    }
#endif
#endif

    hopstotchHash(uint32_t initialCapacity = 16) : n(0) {
        initialize(initialCapacity);
    }

    ~hopstotchHash() {
        delete[] v;
    }

    void initialize(uint32_t initialCapacity) {
        preZero = 0;
        while ((1u << preZero) <= 2 * initialCapacity - 1) preZero++;
        preZero = 32 - preZero;

        roundN = 1u << (32 - preZero);
        if (roundN < H) roundN = H;
        v = new uint32_t[roundN + H];

        memset(v, EMPTY, sizeof(uint32_t) * (roundN + H));
    }

    void build(uint32_t *adjList, uint32_t n_) {
        n = n_;
        if (n > 0) {
            preZero = 0;
            while ((1u << preZero) <= 2 * n - 1) preZero++;
            preZero = 32 - preZero;

            roundN = 1u << (32 - preZero);
            if (roundN < H) roundN = H;
            delete[] v;  // Delete the old array if one exists
            v = new uint32_t[roundN + H];

            memset(v, EMPTY, sizeof(uint32_t) * (roundN + H));
            for (uint32_t i = 0; i < n; i++) {
                uint32_t u = adjList[i];
                if (!insert(u)) {
                    bool f = reBuild();
                    if (!f) {
                        printf("error build hash table\n");
                        fflush(stdout);
                        exit(-1);
                    } else
                        insert(u);
                }
            }
            for (uint32_t i = 0; i < H - 1; i++) {
                v[i + roundN] = v[i];
            }

            for (uint32_t i = 0; i < n; i++) {
                uint32_t u = adjList[i];
                if (!contain(u)) {
                    printf("error build hash table 2\n");
                    fflush(stdout);
                    exit(-1);
                }
            }
        }
    }

    bool insert(uint32_t u) {
        if (n >= roundN / 2) {  // Load factor threshold to trigger rehash
            rehash(roundN * 2);
        }

        uint32_t hashV = hash(u);
        uint32_t p = hashV;
        uint32_t i = 0;
        for (; i < roundN; i++) {
            if (v[p] == EMPTY) break;
            p = (p + 1) % roundN;
        }

        while ((p - hashV + roundN) % roundN >= H) {
            uint32_t t = (p - H + 1 + roundN) % roundN;
            bool f = false;

            for (uint32_t i = t; i != p; i = (i + 1) % roundN) {
                if ((p - hash(v[i]) + roundN) % roundN < H) {
                    v[p] = v[i];
                    v[i] = EMPTY;
                    p = i;
                    f = true;
                    break;
                }
            }

            if (!f) return false;
        }

        v[p] = u;
        n++;
        return true;
    }

    bool containNormal(uint32_t u) {
        uint32_t hashV = hash(u);

        // Iterate over the neighborhood range
        for (uint32_t i = 0; i < H; i++) {
            uint32_t j = (hashV + i) % roundN;  // Ensure wrapping around the table
            if (v[j] == u) {                    // Check if the element is found
                return true;
            }
            if (v[j] == EMPTY) {  // Stop early if an empty slot is found
                return false;
            }
        }

        return false;  // Return false if the element is not found within the neighborhood
    }

    bool contain(uint32_t u) {
        return containNormal(u);  // Use this for debugging and simplifying
    }

    bool reBuild() {
        uint32_t *tmpV = v;
        v = new uint32_t[roundN];
        memset(v, EMPTY, sizeof(uint32_t) * roundN);

        uint32_t rebuildTimes = 0;
        while (rebuildTimes < 100) {
            bool f = true;
            for (uint32_t i = rebuildTimes; i < roundN + rebuildTimes; i++) {
                uint32_t j = i % roundN;
                if (tmpV[j] != EMPTY) {
                    if (!insert(tmpV[j])) {
                        f = false;
                        break;
                    }
                }
            }

            if (f)
                break;
            else {
                rebuildTimes++;
                memset(v, EMPTY, sizeof(uint32_t) * roundN);
            }
        }

        delete[] tmpV;

        if (rebuildTimes == 100) return false;

        return true;
    }

    void rehash(uint32_t newCapacity) {
        uint32_t *oldV = v;
        uint32_t oldRoundN = roundN;

        initialize(newCapacity);

        for (uint32_t i = 0; i < oldRoundN; ++i) {
            if (oldV[i] != EMPTY) {
                insert(oldV[i]);
            }
        }

        delete[] oldV;
    }

    // Hash function
    uint32_t hash(uint32_t t) {
        uint32_t tmp = hashKey * t;
        uint32_t ans = tmp >> preZero;
        return ans;
    }
};

#endif
