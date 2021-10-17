// PbiDF.cpp : This file contains code for extraction of dark-field signal from in-line phase-contrast images
//

#include <chrono>
#include <omp.h>
#include "XA_ini.h"
#include "IXAHWave.h"
#include "XArray2D.h"
#include "XA_data.h"
#include "XA_file.h"
#include "XA_fft2.h"
#include "XA_tie.h"
#include "XA_born.h"

using namespace xar;

#define PBIDF
#ifdef PBIDF

int main()
{
	// start the execution timer
	std::chrono::system_clock::time_point start_time = std::chrono::system_clock::now();

	try
	{
		printf("\nStarting PbiDF program ...");

		char cline[1024], ctitle[1024], cparam[1024], cparam1[1024], cparam2[1024], cparam3[1024], cparam4[1024], cparam5[1024];
		printf("\nReading PbiDF.txt input parameter file ...");
		FILE* ff0 = fopen("PbiDF.txt", "rt");
		if (!ff0) throw std::exception("Error: cannot open parameter file PbiDF.txt.");

		fgets(cline, 1024, ff0); // 1st line - comment

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // file name base for input files
		if (sscanf(cline, "%s %[^\n]s", ctitle, cparam) != 2) throw std::exception("Error reading input file name base from input parameter file.");
		string strfilepathin = cparam;
		printf("\nInput file name base = %s", strfilepathin.c_str());

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // file name base for output Born-diff files
		if (sscanf(cline, "%s %[^\n]s", ctitle, cparam) != 2) throw std::exception("Error reading output file name base from input parameter file.");
		string strfilepathout = cparam;
		printf("\nOutput file name base = %s", strfilepathout.c_str());

		bool bTIEout(true);
		fgets(cline, 1024, ff0); strtok(cline, "\n"); // optional file name base for output TIE-hom files
		if (sscanf(cline, "%s %[^\n]s", ctitle, cparam) < 2) throw std::exception("Error reading optional TIE-hom output file name base from input parameter file.");
		string strfilepathout1 = cparam;
		if (strcmp(strfilepathout1.c_str(), "N\0"))
			printf("\nOutput TIE-hom file name base = %s", strfilepathout1.c_str());
		else
		{
			bTIEout = false;
			printf("\nThere will be no output of TIE-hom files");
		}

		bool bI1TIEout(true);
		fgets(cline, 1024, ff0); strtok(cline, "\n"); // optional file name base for output reprojected TIE-hom files
		if (sscanf(cline, "%s %[^\n]s", ctitle, cparam) < 2) throw std::exception("Error reading optional reprojected TIE-hom output file name base from input parameter file.");
		string strfilepathout2 = cparam;
		if (strcmp(strfilepathout2.c_str(), "N\0"))
			printf("\nOutput reprojected TIE-hom file name base = %s", strfilepathout2.c_str());
		else
		{
			bI1TIEout = false;
			printf("\nThere will be no output of TIE-hom files");
		}

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // number of CT projection angles, stride
		if (sscanf(cline, "%s %s %s", ctitle, cparam, cparam1) != 3) throw std::exception("Error reading the number of CT projection angles & stride parameters from input parameter file.");
		int nangles = atoi(cparam); // needs to be 'int' for OpenMP 
		int istride = atoi(cparam1);
		printf("\nNumber of CT projection angles = %d", nangles);

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // CT angle range in degrees
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading CT angle range parameter from input parameter file.");
		double angle_range = atof(cparam);
		printf("\nCT angle range = %g", angle_range);
		double angle_step = angle_range / nangles, angle;

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // parameters for trimming or padding the input image (to bring the dims to powers of 2)
		if (sscanf(cline, "%s %s %s %s %s %s", ctitle, cparam1, cparam2, cparam3, cparam4, cparam5) != 6) throw std::exception("Error reading image trim/pad parameters from input parameter file.");
		int iXLeft = atoi(cparam1);
		int iXRight = atoi(cparam2);
		int iYTop = atoi(cparam3);
		int iYBottom = atoi(cparam4);
		double fPadVal = atof(cparam5);
		printf("\nImage trim(-)/pad(+) parameters: iXLeft = %d, iXRight = %d, iYTop = %d, iYBottom = %d, PadValue = %g", iXLeft, iXRight, iYTop, iYBottom, fPadVal);
		if (iXLeft > 0 && (iXRight < 0 || iYTop < 0 || iYBottom < 0)) throw std::exception("All trim/pad parameters must be either non-negative or non-positive.");
		else if (iXLeft < 0 && (iXRight > 0 || iYTop > 0 || iYBottom > 0)) throw std::exception("All trim/pad parameters must be either non-negative or non-positive.");
		else if (iXRight > 0 && (iYTop < 0 || iYBottom < 0)) throw std::exception("All trim/pad parameters must be either non-negative or non-positive.");
		else if (iXRight < 0 && (iYTop > 0 || iYBottom > 0)) throw std::exception("All trim/pad parameters must be either non-negative or non-positive.");
		else if (iYTop > 0 && iYBottom < 0) throw std::exception("All trim/pad parameters must be either non-negative or non-positive.");
		else if (iYTop < 0 && iYBottom > 0) throw std::exception("All trim/pad parameters must be either non-negative or non-positive.");

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // wavelength in microns (if 0, will be read from file)
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading wavelength parameter from input parameter file.");
		double wl = atof(cparam);
		printf("\nWavelength (microns) = %g", wl);

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // defocus distance in microns
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading defocus distance parameter from input parameter file.");
		double defocus = atof(cparam);
		printf("\nDefocus distance (microns) = %g", defocus);

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // delta/beta for TIE-hom and for Born-hom
		if (sscanf(cline, "%s %s %s", ctitle, cparam, cparam1) != 3) throw std::exception("Error reading delta/beta parameters from input parameter file.");
		double delta2beta = atof(cparam);
		double delta2beta1 = atof(cparam1);
		printf("\nDelta/beta for TIE-hom = %g, delta/beta for Born-hom = %g", delta2beta, delta2beta1);

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // X and Y widths of Gauss filter in detector plane (microns)
		if (sscanf(cline, "%s %s %s", ctitle, cparam, cparam1) != 3) throw std::exception("Error reading Gaussian filter width parameters from input parameter file.");
		double GaussX = atof(cparam);
		double GaussY = atof(cparam1);
		printf("\nWidth of Gaussian filter in detector plane (pixels): X-width = %g, Y-width = %g", GaussX, GaussY);

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // regularization parameter for 1st Born
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading regularization parameter for 1st Born from input parameter file.");
		double alpha = atof(cparam);
		printf("\nRegularization parameter for 1st Born = %g", alpha);

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // calculate Born-hom(0) or TIE-hom+Born-hom(1)
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading Born-hom / Tie+Born-hom / Fresnel-hom parameter from input parameter file.");
		int iMode = atoi(cparam);
		switch (iMode)
		{
		case 0:
			printf("\nThis program will calculate Born-hom retrieval increment to TIE-hom");
			break;
		case 1:
			printf("\nThis program will calculate TIE+Born-hom retrieval");
			break;
		case 2:
			printf("\nThis program will calculate Fresnel-hom retrieval");
			break;
		default:
			throw std::exception("Born-hom / TIE+Born-hom / Fresnel-hom parameter values can only be 0, 1 or 2.");
		}
		fgets(cline, 1024, ff0); strtok(cline, "\n"); // number of worker threads
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading number of worker threads from input parameter file.");
		int nThreads = atoi(cparam);
		printf("\nNumber of worker threads = %d", nThreads);

		fclose(ff0);

		string infile_i; // input file with an in-line projection image
		string outfile_j; // output file with SAXS image
		string outfile1_j; // output file with TIE image
		string outfile2_j; // output file with reprojected TIE image

		XArray2D<double> xaobjtie; // TIE-hom retrieved intensity array
		XArray2D<double> xaint, xaint0, imageI; // auxillary real array
		XArray2D<dcomplex> xacamp; // auxillary complex array

		index_t kk = 2;

		// create formatting string to add properly formatted indexes at the end of the output file names
		size_t i_dot = strfilepathin.rfind('.'), o_dot = strfilepathout.rfind('.'), o_dot1 = strfilepathout1.rfind('.'), o_dot2 = strfilepathout2.rfind('.');
		size_t nfieldB_length;
		char ndig[8];
		string myformat(""), myformat1("");
		if (nangles > 1)
		{
			nfieldB_length = 1 + size_t(log10(double(nangles - 1))); //maximum number of digits corresponding to angles in the input file name
			sprintf(ndig, "%zd", nfieldB_length); //convert the calculated maximum number of digits corresponding to angles into a string, e.g. 3 into "3"
			myformat += "%0" + string(ndig) + "d"; //construct format string for inserting 0-padded angle indexes into file names - see usage below
			nfieldB_length = 1 + size_t(log10(double(nangles / istride - 1))); //maximum number of digits corresponding to angles in the output file name
			sprintf(ndig, "%zd", nfieldB_length); //convert the calculated maximum number of digits corresponding to angles into a string, e.g. 3 into "3"
			myformat1 += "%0" + string(ndig) + "d"; //construct format string for inserting 0-padded angle indexes into file names - see usage below
		}

		omp_set_num_threads(nThreads);
#pragma omp parallel default(none) private(xaobjtie, xaint, xaint0, xacamp, infile_i, outfile_j, outfile1_j)
#pragma omp for schedule(dynamic) nowait
		for (int i = 0; i < nangles; i += istride)
		{
			try // this is needed in OMP mode in order to catch exceptions within the worker threads
			{
				int j = i / istride;
				angle = angle_step * double(i);
				printf("\nAngle = %f", angle);

				// generate input and output file names
				char buffer[128], buffer1[128];
				infile_i = strfilepathin;
				outfile_j = strfilepathout;
				sprintf(buffer, myformat.data(), i);
				sprintf(buffer1, myformat1.data(), j);
				infile_i.insert(i_dot, buffer);
				outfile_j.insert(o_dot, buffer1);
				if (bTIEout)
				{
					outfile1_j = strfilepathout1;
					outfile1_j.insert(o_dot1, buffer1);
				}
				if (bI1TIEout)
				{
					outfile2_j = strfilepathout2;
					outfile2_j.insert(o_dot2, buffer1);
				}

				// read the in-line projection image from input file
				printf("\nReading input file %s ...", infile_i.c_str());
				XArData::ReadFileGRD(xaint, infile_i.c_str(), wl);
				printf("\nWavelength = %g", ((IXAHWave2D*)xaint.GetHeadPtr())->GetWl());
				XArray2DMove<double> xamoveint(xaint);
				if (iXRight > 0 || iXLeft > 0 || iYTop > 0 || iYBottom > 0) xamoveint.Pad(index_t(iYTop), index_t(iYBottom), index_t(iXLeft), index_t(iXRight), fPadVal);
				else if (iXRight < 0 || iXLeft < 0 || iYTop < 0 || iYBottom < 0) xamoveint.Trim(index_t(-iYTop), index_t(-iYBottom), index_t(-iXLeft), index_t(-iXRight));
				double l2dim1 = log2(xaint.GetDim1()), l2dim2 = log2(xaint.GetDim2());
				if (int(l2dim1) != l2dim1 || int(l2dim2) != l2dim2) throw std::exception("Dimensions of the input image after pad/trim are not integer powers of 2.");
				// remove any negative values
				if (xaint.Norm(xar::eNormMin) <= 0)  
				{
					for (int i1 = 0; i1 < xaint.GetDim1(); i1++)
						for (int j1 = 0; j1 < xaint.GetDim2(); j1++) {
							if (xaint[i1][j1] <= 0) {
								xaint[i1][j1] = 0.00001f;
							}
						}
				}

				// If input images have small negative values, this can be fixed differently in the future, but currently we take a safe bet
				if (xaint.Norm(xar::eNormMin) <= 0)
					throw std::exception("Input image intensity distribution contains some negative or zero values.");

				if (iMode == 0 || iMode == 1)
				{
					xaint0 = xaint; // save the trimmed original image for later use
					// do TIE-hom phase retrieval
					XA_2DTIE<double> xatie;
					xatie.DP(xaint, delta2beta, defocus);
					xaobjtie = xaint; // save the TIE-hom retrieved intensity for later use
					if (bTIEout)
					{
						printf("\nWriting TIE output file %s ...", outfile1_j.c_str());
						XArData::WriteFileGRD(xaint, outfile1_j.c_str(), eGRDBIN);
					}

					// homogenise the object-plane complex amplitude and do forward propagation
					xacamp.Resize(xaint.GetDim1(), xaint.GetDim2());
					dcomplex* arrC = &xacamp.front();
					double* arrI = &xaint.front();
					double amp;
					for (index_t k = 0; k < xacamp.size(); k++)
					{
						amp = pow(arrI[k], 0.5f);
						arrC[k] = std::polar<double>(amp, delta2beta * log(amp));
					}
					xacamp.SetHeadPtr(xaint.GetHeadPtr() ? xaint.GetHeadPtr()->Clone() : 0);
					XArray2DFFT<double> xafft2(xacamp);
					xafft2.Fresnel(defocus); // propagate forward to the image plane
					Abs2(xacamp, xaint);

					// Optionally Gauss-filter the re-propagated TIE-hom image in the detector plane
					//XArray2DFFTRe<double> xafilt(xaint);
					//xafilt.FilterGauss(GaussY / 2.0, GaussX / 2.0);
					if (bI1TIEout)
					{
						printf("\nWriting reprojected TIE output file %s ...", outfile2_j.c_str());
						XArData::WriteFileGRD(xaint, outfile2_j.c_str(), eGRDBIN);
					}

					// Calculate the difference between the original image and DP-repropagated image
					//void FilterGauss(double dblSigmaY, double dblSigmaX);
					double* arrI0 = &xaint0.front();
					arrI = &xaint.front();
					for (index_t k = 0; k < xacamp.size(); k++) arrI[k] = arrI0[k] - arrI[k] + 1.0f;
					//for (index_t k = 0; k < xacamp.size(); k++) arrI[k] = arrI0[k] - std::norm(arrC[k]) + 1.0f;
					imageI = xaint;
					// Do 1st Born on the difference between the original image and DP-repropagated image
					double fIin = xaint.Norm(xar::eNormAver);
					XA_2DBorn<double> xaborn;
					XA_2DBorn<double> xaborn2;
					xaborn.BornSC(imageI, defocus, delta2beta1, alpha, false);
					xaborn2.BornSCIter(xaint, imageI, defocus, delta2beta1, alpha, 1.e-14, kk);
					arrI = &xaint.front();
					double* arrItie = &xaobjtie.front();
					
					if (iMode == 1)	// we add TIE-hom to mu_Born here
						for (index_t k = 0; k < xaint.size(); k++) arrI[k] = arrItie[k] + arrI[k] - 1.0f;
				}

				if (iMode == 2)
				{
					XA_2DBorn<double> xaborn;
					xaborn.BornSC(xaint, defocus, delta2beta, alpha, true);
				}

				// trim back and write the result to output file
				if (iXRight > 0 || iXLeft > 0 || iYTop > 0 || iYBottom > 0) xamoveint.Trim(index_t(iYTop), index_t(iYBottom), index_t(iXLeft), index_t(iXRight));
				printf("\nWriting SAXS output file %s ...", outfile_j.c_str());
				XArData::WriteFileGRD(xaint, outfile_j.c_str(), eGRDBIN);
			}
			catch (std::exception& E)
			{
				printf("\n\n!!!Exception: %s\n", E.what());
				exit(1);
			}
		}
	}
	catch (std::exception& E)
	{
		printf("\n\n!!!Exception: %s\n", E.what());
	}

	std::chrono::system_clock::time_point end_time = std::chrono::system_clock::now();
	printf("\n\nMain program finished. Execution time = %I64d s.\n", std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time).count());

	//printf("\nPress any key to exit..."); getchar();

	return 0;
}


#endif // PBIDF