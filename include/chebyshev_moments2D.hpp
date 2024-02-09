#ifndef CHEB_MOMENTS2D_HPP
#define CHEB_MOMENTS2D_HPP

#include "chebyshev_moments.hpp"


class Moments2D : public Moments {
public:
    std::array<size_t, 2> numMoms = {0, 0};

    Moments2D() = default;

    explicit Moments2D(size_t m0, size_t m1) : numMoms{m0, m1} {
        this->MomentVector(vector_t(numMoms[0] * numMoms[1], 0.0));
    }

    explicit Moments2D(const std::string& momfilename);

    Moments2D(const Moments2D& mom, size_t m0, size_t m1) {
        this->getMomentsParams(mom);
        numMoms = {m0, m1};
        this->MomentVector(vector_t(numMoms[0] * numMoms[1], 0.0));
    }

    // Getters
    std::array<int, 2> MomentNumber() const { return numMoms; }
    int HighestMomentNumber(int i) const { return numMoms[i]; }
    int HighestMomentNumber() const { return std::max(numMoms[0], numMoms[1]); }

    // Setters
    void MomentNumber(int mom0, int mom1) { numMoms = {mom0, mom1}; }

    // Operators
    Moments::value_t& operator()(int m0, int m1) {
        return this->MomentVector()[m0 * numMoms[1] + m1];
    }

    // Transformations
    void ApplyJacksonKernel(double b0, double b1);

    // Input/Output
    void saveIn(const std::string& filename);
    void Print() const;

	void AddSubMatrix( Moments2D& sub , const int mL, const int mR)
	{
		for(int m0=0; m0<sub.HighestMomentNumber(0); m0++)
		for(int m1=0; m1<sub.HighestMomentNumber(1); m1++)
			this->operator()(mL+m0,mR+m1) += sub(m0,m1);
	} 

	void InsertSubMatrix( Moments2D& sub , const int mL, const int mR)
	{
		for(int m0=0; m0<sub.HighestMomentNumber(0); m0++)
		for(int m1=0; m1<sub.HighestMomentNumber(1); m1++)
			this->operator()(mL+m0,mR+m1) = sub(m0,m1);
	} 
private:
	array<int, 2> numMoms;
};

  class MomentsTD : public Moments
  {
	  public:
	  const double HBAR = 0.6582119624 ;//planck constant in eV.fs
	  
	  MomentsTD():
	  _numMoms(1), _maxTimeStep(1),
	  _timeStep(0),
	  _dt(0)
	  {};
	  
	  MomentsTD( const size_t m, const size_t n ): 
	  _numMoms(m), _maxTimeStep(n),
	  _timeStep(0),
	  _dt(0)
	  { this->MomentVector( Moments::vector_t(m*n, 0.0) );    };
	  
	  MomentsTD( std::string momfilename );
	  
	  //GETTERS
	  inline
	  size_t MomentNumber() const { return _numMoms;};
	  
	  inline
	  size_t HighestMomentNumber() const { return _numMoms;};
	  
	  inline
	  size_t CurrentTimeStep() const { return _timeStep;};

	  inline 
	  size_t MaxTimeStep() const { return _maxTimeStep; };
	  
	  inline
	  double TimeDiff() const   { return _dt; };
	  	  
	  inline
	  double ChebyshevFreq() const   { return HalfWidth()/chebyshev::CUTOFF/HBAR; };

	  inline
	  double ChebyshevFreq_0() const   { return BandCenter()/HBAR; };
	  
	  //SETTERS
	  
	  void MomentNumber(const size_t mom);
	  
	  void MaxTimeStep(const  size_t maxTimeStep )  {  _maxTimeStep = maxTimeStep; };

	  inline
	  void IncreaseTimeStep(){ _timeStep++; };

	  inline
	  void ResetTime(){ _timeStep=0; };
	  
	  inline
	  void TimeDiff(const double dt ) { _dt = dt; };
	  
	  
	  int Evolve(  vector_t& Phi);
	  
	  //OPERATORS
	  inline
	  Moments::value_t& operator()(const size_t m, const size_t n)
	  {
		  return this->MomentVector( m*MaxTimeStep() + n );
	  };
	  
	  //Transformation
	  void ApplyJacksonKernel( const double broad );
	  
	  //COSTFUL FUNCTIONS
	  void saveIn(std::string filename);
	  
	  void Print();

  private:
    size_t _numMoms, _maxTimeStep, _timeStep;
    double _dt;
  };

#endif
