#ifndef CHEB_MOMENTSTD_HPP
#define CHEB_MOMENTSTD_HPP

#include "chebyshev_moments.hpp"



class MomentsTD : public Moments {
public:
    static constexpr double HBAR = 0.6582119624; // Planck constant in eV.fs

    // Using delegating constructors to reduce code repetition
    MomentsTD() : MomentsTD(1, 1) {}

    explicit MomentsTD(size_t m, size_t n) 
        : _numMoms(m), _maxTimeStep(n), _timeStep(0), _dt(0) {
        this->MomentVector(Moments::vector_t(m * n, 0.0));
    }

    explicit MomentsTD(const std::string& momfilename); // Implementation needed

    // Getters
    size_t MomentNumber() const { return _numMoms; }
    size_t HighestMomentNumber() const { return _numMoms; }
    size_t CurrentTimeStep() const { return _timeStep; }
    size_t MaxTimeStep() const { return _maxTimeStep; }
    double TimeDiff() const { return _dt; }
    double ChebyshevFreq() const { return HalfWidth() / CUTOFF / HBAR; }
    double ChebyshevFreq_0() const { return BandCenter() / HBAR; }

    // Setters
    void MomentNumber(size_t mom) { _numMoms = mom; }
    void MaxTimeStep(size_t maxTimeStep) { _maxTimeStep = maxTimeStep; }
    void IncreaseTimeStep() { ++_timeStep; }
    void ResetTime() { _timeStep = 0; }
    void TimeDiff(double dt) { _dt = dt; }

    // Operations
    int Evolve(const vector_t& Phi); // Implementation needed
    void ApplyJacksonKernel(double broad); // Implementation needed

    // File IO
    void saveIn(const std::string& filename); // Implementation needed
    void Print() const; // Implementation needed

    // Operators
    Moments::value_t& operator()(size_t m, size_t n) {
        return this->MomentVector()[m * MaxTimeStep() + n];
    }

private:
    size_t _numMoms, _maxTimeStep, _timeStep;
    double _dt;
};
#endif
