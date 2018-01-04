#ifndef SOA_H
#define SOA_H

#include "../common/smart_assert.h"

#include <string>
#include <vector>

struct Chunk
{
    Chunk(std::size_t blockSize);
    ~Chunk();
    void* Allocate(std::size_t numBlocks);
    void Clear();

    char* pData_;
    unsigned int nextAvailableData_;
    unsigned int blockSize_;
    static const unsigned int kChunkSize = 8 * 1024 * 1024;

private:
    Chunk(const Chunk&);
    Chunk& operator=(const Chunk&);
};

class FixedSizeObjectAllocator
{
public:
    FixedSizeObjectAllocator(std::size_t blockSize);
    ~FixedSizeObjectAllocator();
    
    void* Allocate(std::size_t numBlocks);

    void Clear();

private:
    std::size_t RoundUpTo16(std::size_t size) {
		if ((size % 16) == 0) return size;
        std::size_t r = 16 - (size % 16);
        return size + r;
    }

private:
    std::vector<Chunk*> chunkList_;
    std::size_t firstAvailableChunkIdx_;
    std::size_t blockSize_;
};

#endif // SOA_H
