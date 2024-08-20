#include <cstdlib>
#include <iostream>
#include <vector>

class CuckooHash {
   private:
    int size;
    int threshold;
    std::vector<int> table1, table2;

    inline int hash1(int key) const {
        return key & (size - 1);  // Using bitwise AND instead of modulo
    }

    inline int hash2(int key) const {
        return (key >> 5) & (size - 1);  // Using bitwise operations
    }

    void rehash() {
        std::vector<int> oldTable1 = table1;
        std::vector<int> oldTable2 = table2;

        size *= 2;
        table1.assign(size, -1);
        table2.assign(size, -1);

        for (int key : oldTable1) {
            if (key != -1) {
                insert(key);
            }
        }
        for (int key : oldTable2) {
            if (key != -1) {
                insert(key);
            }
        }
    }

   public:
    CuckooHash(int size) : size(size), threshold(size), table1(size, -1), table2(size, -1) {}

    bool insert(int key) {
        int pos1 = hash1(key);
        if (table1[pos1] == -1) {
            table1[pos1] = key;
            return true;
        }

        int pos2 = hash2(key);
        if (table2[pos2] == -1) {
            table2[pos2] = key;
            return true;
        }

        int currentKey = key;
        for (int i = 0; i < threshold; ++i) {  // Reduced threshold for fewer iterations
            std::swap(currentKey, table1[pos1]);
            pos2 = hash2(currentKey);
            if (table2[pos2] == -1) {
                table2[pos2] = currentKey;
                return true;
            }
            std::swap(currentKey, table2[pos2]);
            pos1 = hash1(currentKey);
        }

        // Rehash if needed
        rehash();
        return insert(key);
    }

    bool contain(int key) const {
        int pos1 = hash1(key);
        if (table1[pos1] == key) {
            return true;
        }

        int pos2 = hash2(key);
        if (table2[pos2] == key) {
            return true;
        }

        return false;
    }

    bool remove(int key) {
        int pos1 = hash1(key);
        if (table1[pos1] == key) {
            table1[pos1] = -1;
            return true;
        }

        int pos2 = hash2(key);
        if (table2[pos2] == key) {
            table2[pos2] = -1;
            return true;
        }

        return false;
    }
};
