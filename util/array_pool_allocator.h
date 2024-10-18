#include <parlay/alloc.h>


class abstract_type_allocator {
    virtual void* alloc() = 0;
    virtual void free(void* ptr) = 0;
};

template<typename T>
class type_allocator_wrapper : public abstract_type_allocator {
    public: parlay::type_allocator<T> allocator;
    void* alloc() {
        return (void*) allocator.alloc();
    }
    void free(void* ptr) {
        allocator.free((T*) ptr);
    }
};

template<typename T>
class ArrayPoolAllocator {
private:
    std::vector<abstract_type_allocator*> allocators;
public:
    T* alloc(uint8_t len) {
        T* array_ptr = (T*) allocators[len-1]->alloc();
        return array_ptr;
    }
    void free(T* ptr, uint8_t len) {
        allocators[len-1]->free((void*)ptr);
    }
    ArrayPoolAllocator() {
        allocators.push_back(new type_allocator_wrapper<T[1]>);
        allocators.push_back(new type_allocator_wrapper<T[2]>);
        allocators.push_back(new type_allocator_wrapper<T[3]>);
        allocators.push_back(new type_allocator_wrapper<T[4]>);
        allocators.push_back(new type_allocator_wrapper<T[5]>);
        allocators.push_back(new type_allocator_wrapper<T[6]>);
        allocators.push_back(new type_allocator_wrapper<T[7]>);
        allocators.push_back(new type_allocator_wrapper<T[8]>);
        allocators.push_back(new type_allocator_wrapper<T[9]>);
        allocators.push_back(new type_allocator_wrapper<T[10]>);
        allocators.push_back(new type_allocator_wrapper<T[11]>);
        allocators.push_back(new type_allocator_wrapper<T[12]>);
        allocators.push_back(new type_allocator_wrapper<T[13]>);
        allocators.push_back(new type_allocator_wrapper<T[14]>);
        allocators.push_back(new type_allocator_wrapper<T[15]>);
        allocators.push_back(new type_allocator_wrapper<T[16]>);
        allocators.push_back(new type_allocator_wrapper<T[17]>);
        allocators.push_back(new type_allocator_wrapper<T[18]>);
        allocators.push_back(new type_allocator_wrapper<T[19]>);
        allocators.push_back(new type_allocator_wrapper<T[20]>);
        allocators.push_back(new type_allocator_wrapper<T[21]>);
        allocators.push_back(new type_allocator_wrapper<T[22]>);
        allocators.push_back(new type_allocator_wrapper<T[23]>);
        allocators.push_back(new type_allocator_wrapper<T[24]>);
        allocators.push_back(new type_allocator_wrapper<T[25]>);
        allocators.push_back(new type_allocator_wrapper<T[26]>);
        allocators.push_back(new type_allocator_wrapper<T[27]>);
        allocators.push_back(new type_allocator_wrapper<T[28]>);
        allocators.push_back(new type_allocator_wrapper<T[29]>);
        allocators.push_back(new type_allocator_wrapper<T[30]>);
        allocators.push_back(new type_allocator_wrapper<T[31]>);
        allocators.push_back(new type_allocator_wrapper<T[32]>);
        allocators.push_back(new type_allocator_wrapper<T[33]>);
        allocators.push_back(new type_allocator_wrapper<T[34]>);
        allocators.push_back(new type_allocator_wrapper<T[35]>);
        allocators.push_back(new type_allocator_wrapper<T[36]>);
        allocators.push_back(new type_allocator_wrapper<T[37]>);
        allocators.push_back(new type_allocator_wrapper<T[38]>);
        allocators.push_back(new type_allocator_wrapper<T[39]>);
        allocators.push_back(new type_allocator_wrapper<T[40]>);
        allocators.push_back(new type_allocator_wrapper<T[41]>);
        allocators.push_back(new type_allocator_wrapper<T[42]>);
        allocators.push_back(new type_allocator_wrapper<T[43]>);
        allocators.push_back(new type_allocator_wrapper<T[44]>);
        allocators.push_back(new type_allocator_wrapper<T[45]>);
        allocators.push_back(new type_allocator_wrapper<T[46]>);
        allocators.push_back(new type_allocator_wrapper<T[47]>);
        allocators.push_back(new type_allocator_wrapper<T[48]>);
        allocators.push_back(new type_allocator_wrapper<T[49]>);
        allocators.push_back(new type_allocator_wrapper<T[50]>);
    }
};
