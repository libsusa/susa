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
 * @brief The <i>memory</i> template definition and declaration.
 * @author Behrooz Kamary Aliabadi
 * @version 1.0.0
 */

#ifndef SUSA_MEMORY_H
#define SUSA_MEMORY_H

#include <cstdlib>
#include <susa/debug.h>

namespace susa
{
  /**
   * @brief The <i>memory</i> class.
   *
   * This class may not be used as is. It is a memory manager class for the types.
   *
   * @ingroup TYPES
   *
   */
    template <class T> class memory
    {
        protected:
            T* _matrix;
            size_t sizet_bytes;
            size_t sizet_objects;

        public:
            //! Constructor
            memory();

            //! Destructor
            ~memory();

            /**
             * Allocates a memory space without initialization.
             *
             * @param sizet_size the number of T elements.
             */
            void allocate(size_t sizet_size);

            //! Release the allocated memory.
            void deallocate();

            //! Returns the number of allocated objects.
            size_t size() const;
    };

    template <class T> memory<T>::memory()
    {
        _matrix        = nullptr;
        sizet_bytes    = 0;
        sizet_objects  = 0;
    }

    template <class T> memory<T>::~memory()
    {
        deallocate();
    }

    template <class T> void memory<T>::allocate(size_t sizet_size)
    {
        sizet_objects = sizet_size;
        sizet_bytes   = sizet_size * sizeof(T);

        if (_matrix == nullptr)
        {
            void* block = std::malloc(sizet_size * sizeof(T));
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
            void* block = std::realloc(_matrix, sizet_size * sizeof(T));

            if (block == nullptr)
            {
                deallocate();
                SUSA_ASSERT_MESSAGE(block != nullptr, "memory allocation failed.");
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
            _matrix      = nullptr;
            sizet_objects = 0;
            sizet_bytes   = 0;
        }
    }

    template <class T> size_t memory <T>::size() const
    {
        return sizet_objects;
    }
}
#endif
