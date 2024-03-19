#include "chebyshev_moments.hpp"


int chebyshev::Vectors_sliced::IterateAllSliced(int s )
{
  size_t segment_size = ( s == num_sections_ ? last_section_size_ : section_size_ ),
    segment_start = s * section_size_,
    DIM = this->SystemSize();

	//The vectorss Chebyshev0() and Chebyshev1() are assumed to have
	// been initialized
       linalg::extract_segment( Chebyshev0(), DIM, segment_start, Vector(0), segment_size);//Not parallelized. Easy with omp
	for(int m=1; m < this->NumberOfVectors(); m++ )
	{
	  linalg::extract_segment( Chebyshev1(), DIM, segment_start, Vector(m),  segment_size );//Not parallelized. Easy with omp
	  this->Hamiltonian().Multiply(2.0,Chebyshev1(),-1.0,Chebyshev0());
	  Chebyshev0().swap(Chebyshev1());
	}
	return 0;
};


int chebyshev::Vectors_sliced::MultiplySliced( SparseMatrixType &OP, int s)
{
  size_t segment_size = ( s == num_sections_ ? last_section_size_ : section_size_ ),
    segment_start = s * section_size_,
    DIM = this->SystemSize();
  
	assert( OP.rank() == this->SystemSize() );
	if( OPV_.size()!= OP.rank() )
	       OPV_ = Moments::vector_t ( OP.rank() );


	for(int i=0; i<OP.rank(); i++)//not parallelized; with omp/ eigen this is straightforward;
	  OPV_[i] = 0.0;


	
	for(int m=0; m < this->NumberOfVectors(); m++ )
	{
	  linalg::introduce_segment(Chebmu_.ListElem(m), segment_size, OPV_, DIM, segment_start);//Suboptimal. Ideally, there would be a blockOP x cheb_segment_vector
		OP.Multiply(  OPV_, Chebmu_.ListElem(m) );
	}

	return 0;
};



/*
int chebyshev::Vectors::IterateAll( )
{	
	//The vectorss Chebyshev0() and Chebyshev1() are assumed to have
	// been initialized
	linalg::copy( this->Chebyshev0() ,this->Vector(0) );
	for(int m=1; m < this->NumberOfVectors(); m++ )
	{
		linalg::copy( Chebyshev1() , this->Vector(m) );
		this->Hamiltonian().Multiply(2.0,Chebyshev1(),-1.0,Chebyshev0());
		Chebyshev0().swap(Chebyshev1());
	}
	return 0;
};



int chebyshev::Vectors::Multiply( SparseMatrixType &OP )
{
	assert( OP.rank() == this->SystemSize() );
	if( this->OPV.size()!= OP.rank() )
		this->OPV = Moments::vector_t ( OP.rank() );
	
	for(size_t m=0; m < this->NumberOfVectors(); m++ )
	{
		linalg::copy( this->Chebmu.ListElem(m), this->OPV ); 
		OP.Multiply(  this->OPV, this->Chebmu.ListElem(m) );
	}

	return 0;
};
*/

int chebyshev::Vectors::EvolveAll(const double DeltaT, const double Omega0)
{
	const auto dim = this->SystemSize();
	const auto numVecs = this->NumberOfVectors();

	if( this->Chebyshev0().size()!= dim )
		this->Chebyshev0() = Moments::vector_t(dim,Moments::value_t(0)); 

	if( this->Chebyshev1().size()!= dim )
		this->Chebyshev1() = Moments::vector_t(dim,Moments::value_t(0)); 
	//From now on this-> will be discarded in Chebyshev0() and Chebyshev1()

	const auto I = Moments::value_t(0, 1);
	const double x = Omega0*DeltaT;
	for(size_t m=0; m < this->NumberOfVectors(); m++ )
	{
		auto& myVec = this->Vector(m);
		
		int n = 0;
		double Jn = besselJ(n,x);
		linalg::copy(myVec , Chebyshev0());
		linalg::scal(0, myVec); //Set to zero
		linalg::axpy( Jn , Chebyshev0(), myVec);

		double Jn1 = besselJ(n+1,x);	
		this->Hamiltonian().Multiply(Chebyshev0(), Chebyshev1());
		linalg::axpy(-value_t(2) * I * Jn1, Chebyshev1(), myVec);
		
		auto nIp =-I;
		while( 0.5*(std::abs(Jn)+std::abs(Jn1) ) > 1e-15)
		{
			nIp*=-I ;
			Jn  = Jn1;
			Jn1 = besselJ(n, x);
			this->Hamiltonian().Multiply(2.0, Chebyshev1(), -1.0, Chebyshev0());
			linalg::axpy(value_t(2) * nIp * value_t(Jn1), Chebyshev0(), myVec);
			Chebyshev0().swap(Chebyshev1());
			n++;
		}
	}
  return 0;
};


double chebyshev::Vectors::MemoryConsumptionInGB()
{
	return SizeInGB()+2.0*( (double)this->SystemSize() )*pow(2.0,-30.0) ;
}
	
