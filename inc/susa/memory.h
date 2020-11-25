/*
 * This file is part of Susa.
 *
 * Susa is free software: you can redistribute it and/or modify
 * it under the terms of the Lesser GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * at your option) any later version.
 *
 * Susa is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * Lesser GNU General Public License for more details.
 *
 * You should have received a copy of the Lesser GNU General Public License
 * along with Susa. If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * @file memory.h
 * @brief The <i>memory</i> management and logging (declaration and definition).
 * @author Behrooz Kamary
 * 
 * @defgroup MEMO Memory Management 
 */

#ifndef SUSA_MEMORY_H
#define SUSA_MEMORY_H

#include <atomic>
#include <cstdlib>
#include <susa/debug.h>

namespace susa
{
/**
* @brief The <i>memory_tracker</i> is a singleton class to track allocated and deallocated memory.
*
* @ingroup MEMO
*
*/
class memory_tacker 
{
    private:
        memory_tacker() = default;
        ~memory_tacker() = default;
        std::atomic_ullong ull_used_memory;

    public:
        memory_tacker(const memory_tacker&) = delete;
        memory_tacker(memory_tacker&&) = delete;
        memory_tacker& operator=(const memory_tacker&) = delete;
        memory_tacker& operator=(memory_tacker&&) = delete;

        //! Returns the instance of <i>memory_tacker</i>
        static memory_tacker& instance()
        {
            static memory_tacker INSTANCE;
            return INSTANCE;
        }

        //! Add a number of bytes
        void add(size_t nbytes)
        {
            ull_used_memory.fetch_add(nbytes);
        }

        //! Subtract a number of bytes
        void sub(size_t nbytes)
        {
            ull_used_memory.fetch_sub(nbytes);
        }

        //! Read the allocated number of bytes
        size_t read()
        {
            return ull_used_memory.load();
        }
};

/**
* @brief The allocator_log class inherits std::allocator to manage and log the dynamic memory allocation/deallocation.
*
* @ingroup MEMO
*
*/
template <typename T> class allocator_log : public std::allocator<T>
{
    public:
        using base = std::allocator<T>;
        using pointer = typename std::allocator_traits<base>::pointer;
        using size_type = typename std::allocator_traits<base>::size_type;
        using value_type = T;


        template <typename U> struct rebind
        {
            using type = allocator_log<U>;
        };

        /**
         * Allocates a memory space.
         *
         * @param size the number of elements
         * @param hint the hint
         */
        pointer allocate(size_type size, const void *hint = 0)
        {
            memory_tacker::instance().add(size * sizeof(T));
            pointer ptr = base::allocate(size, hint);
            SUSA_LOG_INF("allocate " << size * sizeof(T) << " bytes. pointer (" << ptr << ")");
            return ptr;
        }

        /**
         * Releases an allocated memory space.
         *
         * @param ptr the allocated pointer
         * @param size the size of allocated memory
         */
        void deallocate(pointer ptr, size_type size)
        {
            SUSA_LOG_INF("deallocate " << size * sizeof(T) << " bytes. pointer (" << ptr << ")");
            memory_tacker::instance().sub(size * sizeof(T));
            return base::deallocate(ptr, size);
        }

        allocator_log() noexcept : base() { SUSA_LOG_INF("constructor"); }

        allocator_log(const allocator_log& arg) noexcept : base(arg) {}

        template <class U> allocator_log(const allocator_log<U>& arg) noexcept : base(arg) {}

        ~allocator_log() noexcept {}
};

/**
* @brief The <i>memory</i> class uses std::malloc() and std::free() to manage the dynamic memory.
*
* <i>memory</i> is a base memory manager class for Susa types.
*
* @ingroup MEMO
*
*/
template <class T> class memory
{
    protected:
        T*     _matrix;
        size_t sizet_bytes;
        size_t sizet_objects;

        /**
         * Allocates a memory space without initialization.
         *
         * @param sizet_size the number of T elements.
         */
        void allocate(size_t sizet_size);

        //! Release the allocated memory.
        void deallocate();

    public:
        //! Constructor
        memory();

        //! Destructor
        virtual ~memory() noexcept;

        //! Copy assignment constructor
        memory& operator=(const memory <T>& mat_arg);

        //! Copy constructor
        memory(const memory <T>& mat_arg);

        //! Copy Constructor for rvalues
        memory(memory <T> &&mat_arg) noexcept;


        //! Returns the number of allocated objects.
        size_t size() const;
};

template <class T> memory<T>::memory()
{
    _matrix        = nullptr;
    sizet_bytes    = 0;
    sizet_objects  = 0;
}

template <class T> memory<T>::~memory() noexcept
{
    deallocate();
}

template <class T> memory<T>::memory(const memory <T>& mat_arg)
{
    _matrix           = nullptr;
    sizet_objects     = 0;
    sizet_bytes       = 0;

    size_t sizet_size = mat_arg.sizet_objects;

    if (sizet_size != 0)
    {
        this->allocate(sizet_size);

        for (size_t sizet_index = 0; sizet_index < sizet_size; sizet_index++)
        {
            this->_matrix[sizet_index] = mat_arg._matrix[sizet_index];
        }
    }
    else
    {
        this->deallocate();
    }

}


template <class T> memory<T>& memory<T>::operator=(const memory <T>& mat_arg)
{

    size_t sizet_size = mat_arg.sizet_objects;

    if (this != &mat_arg)
    {
        if (sizet_size != 0)
        {

            if (this->sizet_objects != sizet_size)
            {
                this->allocate(sizet_size);
            }

            for (size_t sizet_index = 0; sizet_index < sizet_size; sizet_index++)
            {
                this->_matrix[sizet_index] = mat_arg._matrix[sizet_index];
            }
        }
        else
        {
            this->deallocate();
        }
    }

    return *this;
}

template <class T> void memory<T>::allocate(size_t sizet_size)
{

    if (sizet_objects == sizet_size) return;

    sizet_objects = sizet_size;
    sizet_bytes   = sizet_size * sizeof(T);

    if (_matrix == nullptr)
    {
        void* block = std::malloc(sizet_bytes);
        SUSA_ASSERT_MESSAGE(block != nullptr, "memory allocation failed.");
        if (block == nullptr) std::exit(EXIT_FAILURE);
        _matrix = static_cast<T*>(block);
    }
    else if (sizet_size == 0)
    {
        deallocate();
    }
    else
    {
        void* block = std::realloc((void*)_matrix, sizet_bytes);

        SUSA_ASSERT_MESSAGE(block != nullptr, "memory allocation failed.");

        if (block == nullptr)
        {
            deallocate();
            std::exit(EXIT_FAILURE);
        }

        _matrix = static_cast<T*>(block);
    }
}

template <class T> void memory<T>::deallocate()
{
    if (_matrix != nullptr)
    {
        std::free(_matrix);
        _matrix       = nullptr;
        sizet_objects = 0;
        sizet_bytes   = 0;
    }
}

template <class T> size_t memory <T>::size() const
{
    return sizet_objects;
}


template <class T> memory <T>::memory(memory&& mat_arg) noexcept
{

  this->sizet_objects = mat_arg.sizet_objects;
  this->sizet_bytes   = mat_arg.sizet_bytes;

  this->_matrix         = mat_arg._matrix;
  mat_arg._matrix       = nullptr;
  mat_arg.sizet_objects = 0;
  mat_arg.sizet_bytes   = 0;

}

}
#endif
