#include "soa.h"

Chunk::Chunk(std::size_t blockSize)
{
    pData_ = new char[kChunkSize];
    nextAvailableData_ = 0;
    blockSize_ = blockSize;
}

Chunk::~Chunk()
{
    delete[] pData_;
}

void
Chunk::Clear()
{
    nextAvailableData_ = 0;
}

void*
Chunk::Allocate(std::size_t numBlocks)
{
    if (nextAvailableData_ + blockSize_ * numBlocks > kChunkSize) return 0;
    char* p = pData_ + nextAvailableData_;
    nextAvailableData_ += blockSize_ * numBlocks;
    return static_cast<void*>(p);
}

FixedSizeObjectAllocator::FixedSizeObjectAllocator(std::size_t blockSize)
{
    blockSize_ = RoundUpTo16(blockSize);
    firstAvailableChunkIdx_ = 0;
    chunkList_.push_back(new Chunk(blockSize_));
}

FixedSizeObjectAllocator::~FixedSizeObjectAllocator()
{
    std::vector<Chunk*>::iterator iter;
    for (iter = chunkList_.begin(); iter != chunkList_.end(); ++iter) {
        delete (*iter);
    }
}

void
FixedSizeObjectAllocator::Clear()
{
    std::vector<Chunk*>::iterator iter;
    for (iter = chunkList_.begin(); iter != chunkList_.end(); ++iter) {
        (**iter).Clear();
    }
    firstAvailableChunkIdx_ = 0;
}

void*
FixedSizeObjectAllocator::Allocate(std::size_t numBlocks)
{
    void* p = chunkList_[firstAvailableChunkIdx_]->Allocate(numBlocks);
    if (!p) {
        if (firstAvailableChunkIdx_ == chunkList_.size() - 1) {
            chunkList_.push_back(new Chunk(blockSize_));
        }
        ++firstAvailableChunkIdx_;
        p = chunkList_[firstAvailableChunkIdx_]->Allocate(numBlocks);
    }
    r_assert(p);
    return p;
}
