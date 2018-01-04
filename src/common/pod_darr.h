#ifndef POD_DARR_H
#define POD_DARR_H

#include "aux_tools.h"

template <typename T>
class PODArray
{
public:
    PODArray(const idx size = 1024) 
        : m_used(0),
          m_alloc(0),
          m_data(NULL) {
              check_and_realloc(size);
          }
    
    void destroy() {
        if (m_alloc) {
            sfree(m_data);
            m_data = NULL;
        }
        m_used = m_alloc = 0;
    }

    ~PODArray() {
        destroy();
    }

public:
    const T* begin() const {
        return m_data;
    }

    T* begin() {
        return const_cast<T*>(
                static_cast<const PODArray<T>&>(*this).begin()
                );
    }

    const T* end() const {
        return m_data + m_used;
    }

    T* end() {
        return const_cast<T*>(
                static_cast<const PODArray<T>&>(*this).end()
                );
    }

    const T& front() const {
        return m_data[0];
    }
    
    T& front() {
        return const_cast<T&>(
                static_cast<const PODArray<T>&>(*this).front()
                );
    }

    const T& back() const {
        return m_data[m_used - 1];
    }

    T& back() {
        return const_cast<T&>(
                static_cast<const PODArray<T>&>(*this).back()
                );
    }

    const T& operator[](const idx i) const {
        return m_data[i];
    }

    T& operator[](const idx i) {
        return const_cast<T&>(
                static_cast<const PODArray<T>&>(*this)[i]
                );
    }

    idx size() const {
        return m_used;
    }

    idx size() {
        return static_cast<const PODArray<T>&>(*this).size();
    }

    const T* data() const {
        return m_data;
    }

    T* data() {
        return const_cast<T*>(
                static_cast<const PODArray<T>&>(*this).data()
                );
    }

    void clear() {
        m_used = 0;
    }

    bool empty() const {
        return !m_used;
    }

    void push_back(const T* p, const idx size) {
        check_and_realloc(size);
        memcpy(m_data + m_used, p, sizeof(T) * size);
        m_used += size;
    }

    void push_back(const T& t) {
        push_back(&t, 1);
    }

    void pop_back() {
        --m_used;
    }

    void resize(const idx new_size) {
		realloc_data(new_size);
		m_used = new_size;
    }
	
	void reserve(const idx reserve_size) {
		realloc_data(reserve_size);
	}

private:
    void check_and_realloc(const idx added_size) {
		if (m_used + added_size > m_alloc) {
			idx new_alloc = m_alloc ? m_alloc : 1024;
			while (m_used + added_size > new_alloc) new_alloc *= 2;
			realloc_data(new_alloc);
		}
	}
		
	void realloc_data(const idx new_alloc) {
		if (m_alloc < new_alloc) {
			dsnew(new_data, T, new_alloc);
			memcpy(new_data, m_data, sizeof(T) * m_used);
			sfree(m_data);
			m_data = new_data;
			m_alloc = new_alloc;
		}
	}

private:
    idx     m_used;
    idx     m_alloc;
    T*      m_data;
};

#endif // POD_DARR_H
