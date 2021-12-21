#pragma once

#include <iostream>

struct AllocationData
{
    int totalAllocated = 0;
    int totalFreed = 0;

    int GetTotalAllocated(){
        return this->totalAllocated;
    }
    int GetTotalFreed(){
        return this->totalFreed;
    }
    int GetCurrentUsage(){
        return this->totalAllocated - this->totalFreed;
    }
};

static AllocationData allocationData;

void* operator new(size_t size){
    allocationData.totalAllocated += size;
    return malloc(size);
}

void operator delete(void* memory, size_t size) noexcept{
    allocationData.totalFreed += size;
    return free(memory);
}