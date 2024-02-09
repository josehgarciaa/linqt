#ifndef CHEB_MOMENTS1D_HPP
#define CHEB_MOMENTS1D_HPP

#include "chebyshev_moments.hpp"

class Moments1D : public Moments {
public:
    explicit Moments1D(size_t m0 = 0) 
        : _numMoms(m0) {
        this->MomentVector(vector_t(_numMoms, 0.0));
    }

    //Read from an external file
    explicit Moments1D(const std::string& momfilename);

    // Getters
    size_t MomentNumber() const { return _numMoms; }
    size_t HighestMomentNumber() const { return _numMoms; }

    // Operator overloading for direct access
    Moments::value_t& operator()(size_t m0) { return this->MomentVector(m0); }

    // Setters
    void MomentNumber(size_t numMoms);

    // File IO
    void saveIn(const std::string& filename);

    // Transformation
    void ApplyJacksonKernel(double broad);

    // Input/Output
    void Print() const;

private:
    size_t _numMoms;
};

#endif