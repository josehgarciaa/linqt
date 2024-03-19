
#include "linear_algebra.hpp"

void linalg::scal(const complex<double>& a, vector< complex<double> >& x)
{
	cblas_zscal(x.size(), &a, &x[0], 1);
	return ;
}


void linalg::axpy(const complex<double>& a, const vector< complex<double> >& x, vector< complex<double> >& y)
{
	assert( x.size() == y.size() );
	cblas_zaxpy(x.size(), &a, &x[0], 1, &y[0], 1);
	return ;
}

void linalg::copy(const vector< complex<double> >&x,vector< complex<double> >& y)
{
	assert( x.size() == y.size() );
	cblas_zcopy(x.size(), &x[0], 1, &y[0], 1);
}

complex<double> linalg::vdot(const vector< complex<double> >& x,const vector< complex<double> >& y)
{
	assert( x.size() == y.size() );
	complex<double> dotc;
	cblas_zdotc_sub(x.size(), &x[0], 1, &y[0], 1, &dotc);
	return dotc;
}

double linalg::nrm2(const vector< complex<double> >& x)
{
	return  cblas_dznrm2 (x.size(),&x[0],1);
};


void linalg::scal(const int dim, complex<double> a, complex<double> *x)
{
	cblas_zscal(dim, &a, x, 1);
	return ;
}

void linalg::axpy(const int dim, complex<double> a, const complex<double> *x, complex<double> *y)
{
	cblas_zaxpy(dim, &a, x, 1, y, 1);
	return ;
}

void linalg::copy(const int dim, const complex<double> *x, complex<double> *y)
{
	cblas_zcopy(dim, x, 1, y, 1);
}

complex<double> linalg::vdot(const int dim, const complex<double> *x, complex<double> *y)
{
	complex<double> dotc;
	cblas_zdotc_sub(dim, x, 1, y, 1, &dotc);
	return dotc;
}


double linalg::nrm2(const int dim, const complex<double> *x)
{
	return  cblas_dznrm2 (dim,x,1);
};

void linalg::batch_vdot(const int dim,const int batchSize,const complex<double>* leftb,const complex<double>* rightb,complex<double>* output)
{
	complex<double> alpha=1.0, beta=1.0; //This gives the quantity <L|R>* because there is no pure conjugation in MKL
	cblas_zgemm (CblasRowMajor, CblasNoTrans, CblasConjTrans,batchSize, batchSize , dim, &alpha, leftb,dim ,rightb,dim,&beta,output,batchSize);

	const long tot_dim = batchSize*batchSize;
	#pragma omp parallel for 
	for(long int j =0 ; j < tot_dim ; j++)
		output[j]= std::conj(output[j]);
	return ;
}



void linalg::extract_segment(vector< complex<double> >&x, size_t size_x, size_t start_x,  vector< complex<double> >& y, size_t size_y ){//size_x >> size_y
  for(size_t i = 0; i < size_y; i++)//Terrible: not parallelized. Trivial in Eigen/OMP, can't find a block_copy subroutine in cblas.
    x[start_x + i] = y[i];
};

void linalg::introduce_segment(vector< complex<double> >&x, size_t size_x, vector< complex<double> >& y, size_t size_y, size_t start_y ){//size_y >> size_y

  for(size_t i = 0; i< size_x; i++)//Terrible: not parallelized. Trivial in Eigen/OMP, can't find a block_copy subroutine in cblas.
    y[start_y + i] = x[i];

};  
