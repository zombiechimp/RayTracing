//
//  Vector.cpp
//  Algoritam_pracenja_zrake
//
//  Created by Tea Jakić on 01/06/2018.
//  Copyright © 2018 Tea Jakić. All rights reserved.
//

#include <array>
#include <memory>
#include <algorithm>
#include <numeric>
#include <sstream>
#include <math.h>

template <typename T, size_t N>
class Vector
{
public:
    std::array<T, N> m_elements;
    
    //Default constructor
    Vector() :
    m_elements({})
    {
    }
    
    //Overloaded
    Vector(std::array<T, N> elements) :
    m_elements(elements)
    {
        
    }
    
    //Copy-constructor
    Vector(const Vector& other) :
    m_elements(other.m_elements)
    {
    }
    
    //Move constructor
    Vector(Vector&& other) noexcept:
    m_elements(other.m_elements)
    {
    }
    
    Vector& operator=(Vector other)
    {
        
        using std::swap;
        this->swap(other);
        
        return *this;
    }
    
    Vector& operator+=(const Vector&rhs)
    {
        for (int i = 0; i < N; i++)
        {
            m_elements[i] += rhs.m_elements[i];
        }
        
        return *this;
    }
    
    Vector& operator-=(const Vector&rhs)
    {
        for (int i = 0; i < N; i++)
        {
            m_elements[i] -= rhs.m_elements[i];
        }
        
        return *this;
    }
    
    Vector& operator*=(T const& scalar)
    {
        std::transform(m_elements.begin(), m_elements.end(), m_elements.begin(), [&](T el) -> T {return el * scalar; });
        
        return *this;
    }
    
    Vector& operator*=(Vector& rhs)
    {
        if (N != 3) {
            throw 11;
        }
        
        std::array<T, N> result;
        
        result[0] = m_elements[1] * rhs[2] - m_elements[2] * rhs[1];
        result[1] = -(m_elements[0] * rhs[2] - m_elements[2] * rhs[0]);
        result[2] = m_elements[0] * rhs[1] - m_elements[1] * rhs[0];
        
        this->m_elements = result;
        
        return *this;
    }
    
    static size_t get_num_dimesions()
    {
        return N;
    }
    
    T& operator[] (const int index)
    {
        return m_elements[index];
    }
    
    double scalar_product(Vector<T, N> other)
    {
        double res = 0;
        
        for (auto i = 0; i < N; ++i)
            res += other[i] * (*this)[i];
        
        return res;
    }
    
    //TODO: Properties for stream
    void print(std::ostream& out) const
    {
        std::copy(m_elements.cbegin(), m_elements.cend(), std::ostream_iterator<T>(out, " "));
    }
    
    Vector<T, N - 1> fromHomogeneous();
    
    double cosine(const Vector<T, N>& other)
    {
        return ((*this)*other) / (magnitude()*other.magnitude());
    }
    
    Vector normalize()
    {
        double norm = magnitude();
        if(!magnitude()) return *this;
        std::transform(m_elements.begin(), m_elements.end(), m_elements.begin(), [&](T el) -> T {return el / norm; });
        return *this;
    }
    
    Vector getNewNormalized()
    {
        double norm = magnitude();
        std::array<T, N> result;
        std::transform(m_elements.begin(), m_elements.end(), result.begin(), [&](T el) -> T {return el / norm; });
        
        return Vector<T, N>(result);
    }
    
    double magnitude()
    {
        double sumSquares = std::inner_product(m_elements.begin(), m_elements.end(), m_elements.begin(), 0.0);
        return sqrt(sumSquares);
    }
    
    std::array<T, N> toArray()
    {
        return m_elements;
    }
    
private:
    
    void swap(Vector& other)
    {
        std::swap(m_elements, other.m_elements);
    }
};

template <typename T, size_t N>
Vector<T, N - 1> Vector<T, N>::fromHomogeneous()
{
    std::array<T, N - 1> newElems;
    for (auto i = 0; i < N - 1; i++)
    {
        newElems[i] = m_elements[i] / ((double)m_elements[N - 1]);
    }
    
    return Vector<T, N - 1>(newElems);
}

template<class T, size_t N>
inline Vector<T, N> operator+(Vector<T, N> lhs, const Vector<T, N>& rhs)
{
    lhs += rhs;
    return lhs;
}


template<class T, size_t N>
inline Vector<T, N> operator*(T const& scalar, Vector<T, N> v)
{
    return v *= scalar;
}

template<class T, size_t N>
inline Vector<T, N> operator*(Vector<T, N> v, T const & scalar)
{
    return v *= scalar;
}

template<class T, size_t N>
inline Vector<T, N> operator*(Vector<T, N> v, Vector<T, N> v2)
{
    return v *= v2;
}


template<class T, size_t N>
inline Vector<T, N> operator-(Vector<T, N> lhs, const Vector<T, N>& rhs)
{
    lhs -= rhs;
    return lhs;
}


template<class T, size_t N>
std::ostream& operator<<(std::ostream& o, const Vector<T, N>& vec)
{
    vec.print(o);
    return o;
}

template<class T, size_t N>
std::istream& operator>>(std::istream& is, Vector<T, N>& obj)
{
    std::string line;
    std::getline(is, line);
    
    std::istringstream ss(line);
    
    T element;
    size_t el_count = 0;
    
    
    while (ss >> element) {
        if (el_count >= N) {
            throw std::invalid_argument("invalid element count");
        }
        
        obj[el_count++] = element;
    }
    
    return is;
}


typedef Vector<double, 3> Vector3d;
typedef Vector<double, 2> Vector2d;
