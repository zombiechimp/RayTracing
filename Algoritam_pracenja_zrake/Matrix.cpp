//
//  Matrix.cpp
//  Algoritam_pracenja_zrake
//
//  Created by Tea Jakić on 01/06/2018.
//  Copyright © 2018 Tea Jakić. All rights reserved.
//

#include <array>
#include <iostream>
#include <memory>
#include <algorithm>
#include <numeric>
#include <vector>
#include <exception>
#include <string>
#include <sstream>

/*
 N: rows count
 M: cols count
 */


template <typename T>
class Matrix
{
public:
    std::vector<std::vector<T>> mElements;
    
    int rowCount;
    int columnCount;
    
    //Default constructor
    Matrix(int n, int m) :
    rowCount(n),
    columnCount(m),
    mElements(n, std::vector<T>(m)) {}
    
    //Overloaded
    Matrix(std::vector<std::vector<T>> elements) :
    mElements(elements),
    rowCount(elements.size()),
    columnCount(elements[0].size())
    {
    }
    
    //Copy-constructor
    Matrix(const Matrix& other) :
    mElements(other.mElements),
    rowCount(other.rowCount),
    columnCount(other.columnCount)
    {
    }
    
    //Move constructor
    Matrix(Matrix&& other) :
    mElements(other.mElements),
    rowCount(other.rowCount),
    columnCount(other.columnCount)
    {
    }
    
    Matrix& operator=(Matrix other)
    {
        using std::swap;
        this->swap(other);
        
        return *this;
    }
    
    Matrix& operator+=(const Matrix&rhs)
    {
        for (int i = 0; i < rowCount; i++)
        {
            for (int j = 0; j < columnCount; j++)
                mElements[i][j] += rhs.mElements[i][j];
        }
        
        return *this;
    }
    
    Matrix& operator-=(const Matrix&rhs)
    {
        for (int i = 0; i < rowCount; i++)
        {
            for (int j = 0; j < columnCount; j++)
                mElements[i][j] -= rhs.mElements[i][j];
        }
        
        return *this;
    }
    
    Matrix& operator*=(T const& scalar)
    {
        for (int i = 0; i < rowCount; i++)
        {
            for (int j = 0; j < columnCount; j++)
                mElements[i][j] *= scalar;
        }
        
        return *this;
    }
    
    
    Matrix& operator*=(Matrix& rhs)
    {
        if (mElements[0].size() != rhs.rowCount)
        {
            throw 11;
        }
        
        std::vector<std::vector<T>> result(mElements.size(), std::vector<T>(rhs.columnCount));
        
        for (int i = 0; i < rowCount; i++)
        {
            for (int j = 0; j < rhs.columnCount; j++)
            {
                for (int k = 0; k < columnCount; k++)
                {
                    result[i][j] += mElements[i][k] * rhs[k][j];
                }
            }
        }
        
        this->mElements = result;
        
        return *this;
    }
    
    std::vector<double>& operator[] (const int index)
    {
        return mElements[index];
    }
    
    Matrix<T> transposed()
    {
        Matrix<T> transpose(columnCount, rowCount);
        
        for (int i = 0; i < rowCount; i++)
            for (int j = 0; j < columnCount; j++)
                transpose[i][j] = mElements[j][i];
        
        return transpose;
    }
    
    Matrix<T> eliminatedRow(int ind)
    {
        Matrix<T> eliminated(rowCount - 1, columnCount);
        for (int i = 0; i < rowCount; ++i)
        {
            for (int j = 0; j < columnCount; ++j)
            {
                if (i < ind)
                {
                    eliminated[i][j] = mElements[i][j];
                }
                else if (i > ind)
                {
                    eliminated[i - 1][j] = mElements[i][j];
                }
            }
        }
        
        return eliminated;
    }
    
    Matrix<T> eliminatedColumn(int ind)
    {
        Matrix<T> eliminated(rowCount, columnCount - 1);
        
        for (int i = 0; i < rowCount; i++)
        {
            for (int j = 0; j < columnCount; ++j)
            {
                if (j < ind) {
                    eliminated[i][j] = mElements[i][j];
                }
                else if (j > ind)
                {
                    eliminated[i][j - 1] = mElements[i][j];
                }
            }
        }
        
        return eliminated;
    }
    
    Matrix<T> cofactorMatrix()
    {
        Matrix<T> cofactors(rowCount, columnCount);
        
        for (int i = 0; i < rowCount; ++i)
        {
            for (int j = 0; j < columnCount; ++j)
            {
                // i % 2 will be 0 for rows that have the initial sign as +
                // j % 2 will be 0 for cols that have initial plus, if i % 2 is 0
                int sign = ((j % 2) ? -1 : 1) * ((i % 2) ? -1 : 1);
                
                double cofactor = sign * eliminatedRow(i).eliminatedColumn(j).determinant();
                cofactors[i][j] = cofactor;
            }
        }
        
        
        return cofactors;
    }
    
    Matrix<T> inverse()
    {
        if (rowCount != columnCount)
            throw 11;
        
        Matrix<T> adjoint = cofactorMatrix().transposed();
        double det = determinant();
        
        if (det == 0) {
            throw std::domain_error("inverse not defined for given matrix");
        }
        
        return adjoint * (1 / det);
    }
    
    std::vector<std::vector<T>> toVector()
    {
        return mElements;
    }
    
    double determinant() {
        
        if (rowCount != columnCount) {
            throw 11;
        }
        
        
        size_t dim = rowCount;
        if (dim == 1)
        {
            return mElements[0][0];
        }
        
        if (dim == 2) {
            return mElements[0][0] * mElements[1][1] - mElements[0][1] * mElements[1][0];
        }
        if(dim == 3){
            return mElements[0][0] * (mElements[1][1]*mElements[2][2] - mElements[1][2] * mElements[2][1]) - mElements[0][1] * (mElements[1][0]*mElements[2][2] - mElements[2][1] * mElements[0][2]) + mElements[0][2] * (mElements[1][0]*mElements[1][2] - mElements[1][1] * mElements[0][2]);
        }
        
        double determinant = 0;
        int elimination_sign = 1;
        for (int elimination_index = 0; elimination_index < rowCount; elimination_index++)
        {
            Matrix<T> elRow = eliminatedRow(0);
            Matrix<T> elRowCol = elRow.eliminatedColumn(elimination_index);
            double subDet = elRowCol.determinant();
            determinant += subDet * mElements[0][elimination_index] * elimination_sign;
            
            elimination_sign *= -1;
        }
        
        return determinant;
    }
    
    //TODO: Properties for stream
    void print(std::ostream& out) const
    {
        for (int i = 0; i < rowCount; i++) {
            std::copy(mElements[i].cbegin(), mElements[i].cend(), std::ostream_iterator<T>(out, " "));
            out << std::endl;
        }
    }
    
    
private:
    
    void swap(Matrix& other)
    {
        std::swap(mElements, other.mElements);
    }
};

template<class T>
inline Matrix<T> operator+(Matrix<T> lhs, const Matrix<T>& rhs)
{
    lhs += rhs;
    return lhs;
}


template<class T>
inline Matrix<T> operator*(T const& scalar, Matrix<T> mtx)
{
    return mtx *= scalar;
}

template<class T>
inline Matrix<T> operator*(Matrix<T> mtx, T const & scalar)
{
    return mtx *= scalar;
}

template<class T>
inline Matrix<T> operator*(Matrix<T> lhs, Matrix<T> rhs)
{
    return lhs *= rhs;
}


template<class T>
inline Matrix<T> operator-(Matrix<T> lhs, const Matrix<T>& rhs)
{
    lhs -= rhs;
    return lhs;
}


template<class T>
std::ostream& operator<<(std::ostream& o, const Matrix<T>& mtx)
{
    mtx.print(o);
    return o;
}


template<class T>
std::istream& operator>>(std::istream& is, Matrix<T>& obj)
{
    std::vector<std::string> delimitedTokens(obj.rowCount, "");
    
    std::string s;
    std::getline(is, s);
    std::string delimiter = "|";
    
    size_t pos = 0;
    size_t tokenCount = 0;
    std::string token;
    while ((pos = s.find(delimiter)) != std::string::npos)
    {
        token = s.substr(0, pos);
        s.erase(0, pos + delimiter.length());
        
        delimitedTokens[tokenCount] = token;
        tokenCount++;
    }
    
    delimitedTokens[tokenCount] = s;
    tokenCount++;
    
    if (tokenCount != obj.rowCount)
    {
        throw std::invalid_argument("invalid number of rows");
    }
    
    for (int i = 0; i < obj.rowCount; ++i)
    {
        std::istringstream iss(delimitedTokens[i]);
        
        T element;
        size_t elCount = 0;
        while (iss >> element)
        {
            if (elCount < (size_t)obj.columnCount) {
                obj[i][elCount] = element;
                elCount++;
            }
            else
                throw std::invalid_argument("invalid number of columns");
        }
    }
    
    return is;
}

