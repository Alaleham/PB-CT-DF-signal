// Implementation of the fftwf3drc class

#include "fftwf3drc.h"

void Fftwf3drc::Cleanup()
{
	if (aplan != 0) { fftwf_destroy_plan(aplan); aplan = 0; }
	if (bplan != 0) { fftwf_destroy_plan(bplan); bplan = 0; }
	if (pin != 0) { fftwf_free(pin); pin = 0; }
	if (pout != 0) { fftwf_free(pout); pout = 0; }
	nz = ny = nx = 0;
}


void Fftwf3drc::GetRealXArray3D(xar::XArray3D<float>& aaa)
{
	if ((int)aaa.GetDim1() != nz || (int)aaa.GetDim2() != ny || (int)aaa.GetDim3() != nx)
		throw std::runtime_error("dimensions of the output real array are different from the internal one in Fftwf3drc");
	int m = 0;
	float fnorm = 1.0f / float(GetNr());
	for (index_t k1 = 0; k1 < nz; k1++)
		for (index_t j = 0; j < ny; j++)
			for (index_t i = 0; i < nx; i++)
				aaa[k1][j][i] = pin[m++] * fnorm;
}


void Fftwf3drc::SetRealXArray3D(xar::XArray3D<float> aaa)
{
	if ((int)aaa.GetDim1() != nz || (int)aaa.GetDim2() != ny || (int)aaa.GetDim3() != nx)
		throw std::runtime_error("dimensions of the input real array are different from the internal one in Fftwf3drc");
	int m = 0;
	for (index_t k = 0; k < nz; k++)
		for (index_t j = 0; j < ny; j++)
			for (index_t i = 0; i < nx; i++)
				pin[m++] = aaa[k][j][i];
}


void Fftwf3drc::GetComplexXArray3D(xar::XArray3D<xar::fcomplex>& aaa)
{
	int nx2 = GetNx2();
	if ((int)aaa.GetDim1() != nz || (int)aaa.GetDim2() != ny || (int)aaa.GetDim3() != nx2)
		throw std::runtime_error("dimensions of the output complex array are different from the internal one in Fftwf3drc");
	int m = 0;
	for (index_t k = 0; k < nz; k++)
		for (index_t j = 0; j < ny; j++)
			for (index_t i = 0; i < nx2; i++)
			{ 
				aaa[k][j][i] = xar::fcomplex(pout[m][0], pout[m][1]); 
				m++; 
			}
}


void Fftwf3drc::SetComplexXArray3D(xar::XArray3D<xar::fcomplex> aaa)
{
	int nx2 = GetNx2();
	if ((int)aaa.GetDim1() != nz || (int)aaa.GetDim2() != ny || (int)aaa.GetDim3() != nx2)
		throw std::runtime_error("dimensions of the input complex array are different from the internal one in Fftwf3drc");
	int m = 0;
	for (index_t k = 0; k < nz; k++)
		for (index_t j = 0; j < ny; j++)
			for (index_t i = 0; i < nx2; i++)
			{
				pout[m][0] = aaa[k][j][i].real();
				pout[m][1] = aaa[k][j][i].imag();
				m++;
			}
}


void Fftwf3drc::PrintRealArray(const char* message)
{
	float* pin = GetReal();
	printf(message);
	int m = 0;
	for (index_t k = 0; k < nz; k++)
		for (index_t j = 0; j < ny; j++)
			for (index_t i = 0; i < nx; i++)
				printf("\npin[%zd,%zd,%zd] = %g", k, j, i, pin[m++]);
}


void Fftwf3drc::PrintComplexArray(const char* message)
{
	printf(message);
	int m = 0;
	int nx2 = GetNx2();
	for (index_t k = 0; k < nz; k++)
		for (index_t j = 0; j < ny; j++)
			for (index_t i = 0; i < nx2; i++)
			{
				printf("\npout[%zd,%zd,%zd] = (%g, %g)", k, j, i, pout[m][0], pout[m][1]); 
				m++; 
			}
}

