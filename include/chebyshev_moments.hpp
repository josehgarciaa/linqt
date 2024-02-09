// With contributions made by Angel D. Prieto S.
#ifndef CHEBYSHEV_MOMENTS
#define CHEBYSHEV_MOMENTS


#include <complex>
#include <vector>
#include <string>
#include <array>

#include "sparse_matrix.hpp" //contain SparseMatrixType
#include <cassert>			 //needed for assert
#include <fstream>   		 //For ifstream and ofstream
#include <limits>    		 //Needed for dbl::digits10
#include "linear_algebra.hpp"
#include "vector_list.hpp"
#include "special_functions.hpp"
#include "chebyshev_coefficients.hpp"

namespace chebyshev 
{
	const double CUTOFF = 0.99;

class Moments
{
	public:
	typedef std::complex<double>  value_t;
	typedef std::vector< value_t > vector_t;

	//default constructor
	Moments():
	_pNHAM(0),system_label(""),system_size(0),
	band_width(0),band_center(0){};

	//GETTERS
	void getMomentsParams( Moments& mom)
	{
		this->SetHamiltonian( mom.Hamiltonian() ) ; 
		this->SystemLabel( mom.SystemLabel());
		this->BandWidth( mom.BandWidth() );
		this->BandCenter( mom.BandCenter() );
	};	

	//In C++ all members declared within the class definition are 
	//automatically considered inline	

	//Getters
	size_t SystemSize() const { return system_size; };
	string SystemLabel() const { return system_label; };
	double BandWidth() const { return band_width; };
	double HalfWidth() const { return BandWidth()/2.0; };
	double BandCenter() const { return band_center; };
	double ScaleFactor() const { return  2.0*chebyshev::CUTOFF/BandWidth(); };
	double ShiftFactor() const { return -2.0*chebyshev::CUTOFF*BandCenter()/BandWidth(); };
	vector_t& MomentVector() { return mu ;}
	value_t& MomentVector(const int i){return  mu[i]; };
	Moments::vector_t& Chebyshev0(){ return ChebV0; } 
	Moments::vector_t& Chebyshev1(){ return ChebV1; } 
	SparseMatrixType& Hamiltonian()
	{ 
		return *_pNHAM; 
	};

	//SETTERS
	void SystemSize(const int dim)  { system_size = dim; };
	void SystemLabel(string label)  { system_label = label; };
	void BandWidth( const double x)  { band_width = x; };
	void BandCenter(const double x) { band_center = x; };
	void MomentVector(const vector_t _mu ) { mu= _mu;}
	void SetInitVectors( const vector_t& T0 );
	void SetInitVectors( SparseMatrixType &OP ,const vector_t& T0 );
	void SetHamiltonian( SparseMatrixType& NHAM )
	{ 
		if ( this->SystemSize() == 0 ) //Use the rank of the hamiltonian as system size
			this->SystemSize( NHAM.rank() );	
		assert( NHAM.rank() == this->SystemSize()  );
		_pNHAM = &NHAM; 
	};

	//Heavy functions
	int  Rescale2ChebyshevDomain()
	{
		this->Hamiltonian().Rescale(this->ScaleFactor(),this->ShiftFactor());
		return 0;
	};

	void SetAndRescaleHamiltonian(SparseMatrixType& NHAM)
	{ 
		this->SetHamiltonian(NHAM);
		this->Rescale2ChebyshevDomain();
	};

	double Rescale2ChebyshevDomain(const double energ)
	{ 
		return this->ScaleFactor()*energ +this->ShiftFactor(); 
	};

    int JacksonKernelMomCutOff( const double broad );
    double JacksonKernel(const double m,  const double Mom );


	private:
	SparseMatrixType* _pNHAM;

	Moments::vector_t ChebV0,ChebV1,OPV;
	std::string system_label;
	size_t system_size;
	double band_width,band_center;
	vector_t mu;	
};


/*
class Vectors : public Moments {
public:
    using vectorList_t = VectorList<Moments::value_t>;

    explicit Vectors(size_t nMoms = 0, size_t dim = 0) : Chebmu(nMoms, dim), numVecs(nMoms) {}

    template<typename MomentType>
    Vectors(MomentType& mom) : Chebmu(mom.HighestMomentNumber(), mom.SystemSize()), numVecs(mom.HighestMomentNumber()) {
        this->getMomentsParams(mom);
        if (!mom.Chebyshev0().empty()) {
            this->SetInitVectors(mom.Chebyshev0());
        }
    }

    // Getters
    size_t NumberOfVectors() const { return numVecs; }
    size_t Size() const { return static_cast<size_t>(this->SystemSize()) * numVecs; }
    double SizeInGB() const {
        return sizeof(value_t) * Size() * std::pow(2.0, -30.0);
    }
    size_t HighestMomentNumber() const { return Chebmu.ListSize(); }
    vectorList_t& List() { return Chebmu; }
    Moments::vector_t& Vector(size_t m0) { return Chebmu.ListElem(m0); }
    Moments::value_t& operator()(size_t m0) { return Chebmu(m0, 0); }

    // Setters
    void SetNumberOfVectors(size_t x) { numVecs = x; }

    // Operations
    int IterateAll();
    int EvolveAll(double DeltaT, double Omega0);
    int Multiply(SparseMatrixType& OP);

    double MemoryConsumptionInGB();

    Vectors( MomentsTD& mom ): Chebmu(mom.HighestMomentNumber(), mom.SystemSize() )
	{ 
		this->getMomentsParams(mom);
	};

	void getMomentsParams( Moments& mom)
	{
		this->SetHamiltonian( mom.Hamiltonian() ) ; 
		this->SystemLabel( mom.SystemLabel());
		this->BandWidth( mom.BandWidth() );
		this->BandCenter( mom.BandCenter() );
		if ( mom.Chebyshev0().size() == this->SystemSize() )
			this->SetInitVectors( mom.Chebyshev0() );
	}

	inline
	size_t Size() const
	{
		return  (long unsigned int)this->SystemSize()*
				(long unsigned int)this->NumberOfVectors();
	}


private:
    vectorList_t Chebmu;
    size_t numVecs = 0;
	Moments::vector_t OPV;

};
*/
#endif 


