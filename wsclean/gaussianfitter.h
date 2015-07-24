#ifndef GAUSSIAN_FITTER_H
#define GAUSSIAN_FITTER_H

#include "matrix2x2.h"

#include <cmath>
#include <cstring>
#include <iostream>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multifit_nlin.h>

#include "uvector.h"

class GaussianFitter
{
public:
	GaussianFitter() : _posConstrained(0.0) { }
	
	void SetPosConstrained(double positionOffsetConstrained)
	{
		_posConstrained = positionOffsetConstrained;
	}
	
	void Fit2DGaussianCentred(const double* image, size_t width, size_t height, double beamEst, double& beamMaj, double& beamMin, double& beamPA)
	{
		size_t prefSize = std::max<size_t>(10, std::ceil(beamEst*10.0));
		if(prefSize%2 != 0) ++prefSize;
		if(prefSize < width || prefSize < height)
		{
			size_t boxWidth  = std::min(prefSize, width);
			size_t boxHeight = std::min(prefSize, height);
			size_t nIter = 0;
			bool boxWasLargeEnough;
			do {
				fit2DGaussianCentredInBox(image, width, height, beamEst, beamMaj, beamMin, beamPA, boxWidth, boxHeight);
				
				boxWasLargeEnough =
					(beamMaj*4.0 < boxWidth || width>=boxWidth) &&
					(beamMaj*4.0 < boxHeight || height>=boxHeight);
				if(!boxWasLargeEnough)
				{
					prefSize = std::max<size_t>(10, std::ceil(beamMaj*10.0));
					if(prefSize%2 != 0) ++prefSize;
				}
				++nIter;
			} while(!boxWasLargeEnough && nIter < 5);
		}
		else {
			fit2DGaussianCentred(image, width, height, beamEst, beamMaj, beamMin, beamPA);
		}
	}
	
	void Fit2DGaussianFull(const double* image, size_t width, size_t height, double& val, double& posX, double& posY, double& beamMaj, double& beamMin, double& beamPA, double* floorLevel = 0)
	{
		size_t prefSize = std::max<size_t>(10, std::ceil(beamMaj*10.0));
		if(prefSize%2 != 0) ++prefSize;
		if(prefSize < width || prefSize < height)
		{
			size_t xStart  = std::max<int>(0, int(round(posX)) - int(prefSize)/2);
			size_t xEnd    = std::min(width, size_t(round(posX)) + prefSize/2);
			size_t yStart  = std::max<int>(0, int(round(posY)) - int(prefSize)/2);
			size_t yEnd    = std::min(height, size_t(round(posY)) + prefSize/2);
			size_t nIter = 0;
			bool boxWasLargeEnough;
			do {
				fit2DGaussianFullInBox(image, width, height, val, posX, posY, beamMaj, beamMin, beamPA, floorLevel, xStart, xEnd, yStart, yEnd);
				
				size_t boxWidth = xEnd - xStart;
				size_t boxHeight = yEnd - yStart;
				boxWasLargeEnough =
					(beamMaj*4.0 < boxWidth || width>=boxWidth) &&
					(beamMaj*4.0 < boxHeight || height>=boxHeight);
				if(!boxWasLargeEnough)
				{
					prefSize = std::max<size_t>(10, std::ceil(beamMaj*10.0));
					if(prefSize%2 != 0) ++prefSize;
				}
				++nIter;
			} while(!boxWasLargeEnough && nIter < 5);
		}
		else {
			fit2DGaussianFull(image, width, height, val, posX, posY, beamMaj, beamMin, beamPA, floorLevel);
		}
	}
	
private:
	const double* _image;
	size_t _width, _height, _scaleFactor;

	void fit2DGaussianCentredInBox(const double* image, size_t width, size_t height, double beamEst, double& beamMaj, double& beamMin, double& beamPA, size_t boxWidth, size_t boxHeight)
	{
		size_t startX = (width-boxWidth)/2;
		size_t startY = (height-boxHeight)/2;
		ao::uvector<double> smallImage(boxWidth*boxHeight);
		for(size_t y=startY; y!=(height+boxHeight)/2; ++y)
		{
			memcpy(&smallImage[(y-startY)*boxWidth], &image[y*width + startX], sizeof(double)*boxWidth);
		}
		
		fit2DGaussianCentred(&smallImage[0], boxWidth, boxHeight, beamEst, beamMaj, beamMin, beamPA);
	}
	
	void fit2DGaussianCentred(const double* image, size_t width, size_t height, double beamEst, double& beamMaj, double& beamMin, double& beamPA)
	{
		_width = width;
		_height = height;
		_image = image;
		_scaleFactor = (width + height)/2;
		
		const gsl_multifit_fdfsolver_type *T = gsl_multifit_fdfsolver_lmsder;
		gsl_multifit_fdfsolver *solver = gsl_multifit_fdfsolver_alloc (T, _width*_height, 3);
		
		gsl_multifit_function_fdf fdf;
		fdf.f = &fitting_func_centered;
		fdf.df = &fitting_deriv_centered;
		fdf.fdf = &fitting_both_centered;
		fdf.n = _width*_height;
		fdf.p = 3;
		fdf.params = this;
		
		// Using the FWHM formula for a Gaussian:
		const long double sigmaToBeam = 2.0L * sqrtl(2.0L * logl(2.0L));
		double initialValsArray[3] = {
			beamEst/(_scaleFactor*double(sigmaToBeam)),
			beamEst/(_scaleFactor*double(sigmaToBeam)),
			0.0
		};
		gsl_vector_view initialVals = gsl_vector_view_array (initialValsArray, 3);
		gsl_multifit_fdfsolver_set (solver, &fdf, &initialVals.vector);

		int status;
		size_t iter = 0;
		do {
			iter++;
			status = gsl_multifit_fdfsolver_iterate (solver);
			
			if(status)
				break;
			
			status = gsl_multifit_test_delta(solver->dx, solver->x, 1e-7, 1e-7);
			
		} while (status == GSL_CONTINUE && iter < 500);
		
		double
			sx = gsl_vector_get (solver->x, 0),
			sy = gsl_vector_get (solver->x, 1),
			beta = gsl_vector_get (solver->x, 2);
			
		gsl_multifit_fdfsolver_free(solver);
		
		convertShapeParameters(sx, sy, beta, beamMaj, beamMin, beamPA);
	}
	
	static int fitting_func_centered(const gsl_vector *xvec, void *data, gsl_vector *f)
	{
		GaussianFitter& fitter=*static_cast<GaussianFitter*>(data);
		double sx = gsl_vector_get(xvec, 0);
		double sy = gsl_vector_get(xvec, 1);
		double beta = gsl_vector_get(xvec, 2);
		const size_t width = fitter._width, height = fitter._height;
		int xMid = width/2, yMid = height/2;
		double scale = 1.0/fitter._scaleFactor;
		
		size_t dataIndex = 0;
		double errSum = 0.0;
		for(size_t yi=0; yi!=height; ++yi)
		{
			double y = (yi - yMid) * scale;
			for(size_t xi=0; xi!=width; ++xi)
			{
				double x = (xi - xMid) * scale;
				double e = err_centered(fitter._image[dataIndex], x, y, sx, sy, beta);
				errSum += e*e;
				gsl_vector_set(f, dataIndex, e);
				++dataIndex;
			}
		}
		//std::cout << "sx=" << sx << ", sy=" << sy << ", beta=" << beta << ", err=" << errSum << '\n';
		return GSL_SUCCESS;
	}
	
	static int fitting_deriv_centered(const gsl_vector *xvec, void *data, gsl_matrix *J)
	{
		GaussianFitter& fitter=*static_cast<GaussianFitter*>(data);
		double sx = gsl_vector_get(xvec, 0);
		double sy = gsl_vector_get(xvec, 1);
		double beta = gsl_vector_get(xvec, 2);
		const size_t width = fitter._width, height = fitter._height;
		int xMid = width/2, yMid = height/2;
		double scale = 1.0 / fitter._scaleFactor;
		
		size_t dataIndex = 0;
		for(size_t yi=0; yi!=height; ++yi)
		{
			double y = (yi - yMid)*scale;
			for(size_t xi=0; xi!=width; ++xi)
			{
				double x = (xi - xMid)*scale;
				double expTerm = exp(-x*x/(2.0*sx*sx) - beta*x*y/(sx*sy) - y*y/(2.0*sy*sy));
				double dsx = (beta*x*y/(sx*sx*sy)+x*x/(sx*sx*sx)) * expTerm;
				double dsy = (beta*x*y/(sy*sy*sx)+y*y/(sy*sy*sy)) * expTerm;
				double dbeta = -x*y/(sx*sy) * expTerm;
				gsl_matrix_set(J, dataIndex, 0, dsx);
				gsl_matrix_set(J, dataIndex, 1, dsy);
				gsl_matrix_set(J, dataIndex, 2, dbeta);
				++dataIndex;
			}
		}
		return GSL_SUCCESS;
	}
	
	static int fitting_both_centered(const gsl_vector *x, void *data, gsl_vector *f, gsl_matrix *J)
	{
		fitting_func_centered(x, data, f);
		fitting_deriv_centered(x, data, J);
		return GSL_SUCCESS;
	}
	
	static double err_centered(double val, double x, double y, double sx, double sy, double beta)
	{
		return exp(-x*x/(2.0*sx*sx) - beta*x*y/(sx*sy) - y*y/(2.0*sy*sy)) - val;
	}
	
	static double err_full(double val, double v, double x, double y, double sx, double sy, double beta)
	{
		return exp(-x*x/(2.0*sx*sx) - beta*x*y/(sx*sy) - y*y/(2.0*sy*sy))*v - val;
	}
	
	void fit2DGaussianFullInBox(const double* image, size_t width, size_t height, double& val, double& posX, double& posY, double& beamMaj, double& beamMin, double& beamPA, double* floorLevel, size_t xStart, size_t xEnd, size_t yStart, size_t yEnd)
	{
		size_t boxWidth = xEnd - xStart;
		size_t boxHeight = yEnd - yStart;
		ao::uvector<double> smallImage(boxWidth*boxHeight);
		for(size_t y=yStart; y!=yEnd; ++y)
		{
			memcpy(&smallImage[(y-yStart)*boxWidth], &image[y*width + xStart], sizeof(double)*boxWidth);
		}
		
		posX -= xStart;
		posY -= yStart;
		fit2DGaussianFull(&smallImage[0], boxWidth, boxHeight, val, posX, posY, beamMaj, beamMin, beamPA, floorLevel);
		posX += xStart;
		posY += yStart;
	}
	
	void fit2DGaussianFull(const double* image, size_t width, size_t height, double& val, double& posX, double& posY, double& beamMaj, double& beamMin, double& beamPA, double* floorLevel)
	{
		_width = width;
		_height = height;
		_image = image;
		_scaleFactor = (width + height)/2;
		
		if(floorLevel == 0)
			fit2DGaussianFull(image, width, height, val, posX, posY, beamMaj, beamMin, beamPA);
		else
			fit2DGaussianFullWithFloor(image, width, height, val, posX, posY, beamMaj, beamMin, beamPA, *floorLevel);
	}
	
	void fit2DGaussianFull(const double* image, size_t width, size_t height, double& val, double& posX, double& posY, double& beamMaj, double& beamMin, double& beamPA)
	{
		const gsl_multifit_fdfsolver_type *T = gsl_multifit_fdfsolver_lmsder;
		gsl_multifit_fdfsolver *solver = gsl_multifit_fdfsolver_alloc (T, _width*_height, 6);
		
		gsl_multifit_function_fdf fdf;
		fdf.f = &fitting_func_full;
		fdf.df = &fitting_deriv_full;
		fdf.fdf = &fitting_both_full;
		fdf.n = _width*_height;
		fdf.p = 6;
		fdf.params = this;
		
		// Using the FWHM formula for a Gaussian:
		const long double sigmaToBeam = 2.0L * sqrtl(2.0L * logl(2.0L));
		_xInit = -(posX-width/2)/_scaleFactor;
		_yInit = -(posY-height/2)/_scaleFactor;
		double initialValsArray[6] = {
			val,
			_xInit,
			_yInit,
			beamMaj/(_scaleFactor*double(sigmaToBeam)),
			beamMaj/(_scaleFactor*double(sigmaToBeam)),
			0.0
		};
		gsl_vector_view initialVals = gsl_vector_view_array (initialValsArray, 6);
		gsl_multifit_fdfsolver_set (solver, &fdf, &initialVals.vector);

		int status;
		size_t iter = 0;
		do {
			iter++;
			status = gsl_multifit_fdfsolver_iterate (solver);
			
			if(status)
				break;
			
			status = gsl_multifit_test_delta(solver->dx, solver->x, 1e-7, 1e-7);
			
		} while (status == GSL_CONTINUE && iter < 500);
		
		val = gsl_vector_get (solver->x, 0);
		posX = -1.0*gsl_vector_get (solver->x, 1)*_scaleFactor + width/2;
		posY = -1.0*gsl_vector_get (solver->x, 2)*_scaleFactor + height/2;
		double
			sx = gsl_vector_get (solver->x, 3),
			sy = gsl_vector_get (solver->x, 4),
			beta = gsl_vector_get (solver->x, 5);
			
		gsl_multifit_fdfsolver_free(solver);
		
		convertShapeParameters(sx, sy, beta, beamMaj, beamMin, beamPA);
	}
	
	void fit2DGaussianFullWithFloor(const double* image, size_t width, size_t height, double& val, double& posX, double& posY, double& beamMaj, double& beamMin, double& beamPA, double& floorLevel)
	{
		const gsl_multifit_fdfsolver_type *T = gsl_multifit_fdfsolver_lmsder;
		gsl_multifit_fdfsolver *solver = gsl_multifit_fdfsolver_alloc (T, _width*_height, 7);
		
		gsl_multifit_function_fdf fdf;
		fdf.f = &fitting_func_floor;
		fdf.df = &fitting_deriv_floor;
		fdf.fdf = &fitting_both_floor;
		fdf.n = _width*_height;
		fdf.p = 7;
		fdf.params = this;
		
		// Using the FWHM formula for a Gaussian:
		const long double sigmaToBeam = 2.0L * sqrtl(2.0L * logl(2.0L));
		_xInit = -(posX-width/2)/_scaleFactor;
		_yInit = -(posY-height/2)/_scaleFactor;
		double initialValsArray[7] = {
			val,
			_xInit,
			_yInit,
			beamMaj/(_scaleFactor*double(sigmaToBeam)),
			beamMaj/(_scaleFactor*double(sigmaToBeam)),
			0.0,
			0.0
		};
		gsl_vector_view initialVals = gsl_vector_view_array (initialValsArray, 7);
		gsl_multifit_fdfsolver_set (solver, &fdf, &initialVals.vector);

		int status;
		size_t iter = 0;
		do {
			iter++;
			status = gsl_multifit_fdfsolver_iterate (solver);
			
			if(status)
				break;
			
			status = gsl_multifit_test_delta(solver->dx, solver->x, 1e-7, 1e-7);
			
		} while (status == GSL_CONTINUE && iter < 500);
		
		val = gsl_vector_get (solver->x, 0);
		posX = -1.0*gsl_vector_get (solver->x, 1)*_scaleFactor + width/2;
		posY = -1.0*gsl_vector_get (solver->x, 2)*_scaleFactor + height/2;
		double
			sx = gsl_vector_get (solver->x, 3),
			sy = gsl_vector_get (solver->x, 4),
			beta = gsl_vector_get (solver->x, 5);
		floorLevel = gsl_vector_get(solver->x, 6);
			
		gsl_multifit_fdfsolver_free(solver);
		
		convertShapeParameters(sx, sy, beta, beamMaj, beamMin, beamPA);
	}
	
	static int fitting_func_full(const gsl_vector *xvec, void *data, gsl_vector *f)
	{
		GaussianFitter& fitter=*static_cast<GaussianFitter*>(data);
		double
			v = gsl_vector_get (xvec, 0),
			xc = gsl_vector_get (xvec, 1),
			yc = gsl_vector_get (xvec, 2),
			sx = gsl_vector_get (xvec, 3),
			sy = gsl_vector_get (xvec, 4),
			beta = gsl_vector_get (xvec, 5);
		const size_t width = fitter._width, height = fitter._height;
		int xMid = width/2, yMid = height/2;
		double scale = 1.0/fitter._scaleFactor;
		
		size_t dataIndex = 0;
		double errSum = 0.0;
		for(int yi=0; yi!=int(height); ++yi)
		{
			double yS = yc + (yi - yMid) * scale;
			for(int xi=0; xi!=int(width); ++xi)
			{
				double xS = xc + (xi - xMid) * scale;
				double e = err_full(fitter._image[dataIndex], v, xS, yS, sx, sy, beta);
				errSum += e*e;
				gsl_vector_set(f, dataIndex, e);
				++dataIndex;
			}
		}
		//std::cout << "v=" << v << ", x=" << xc << ", y=" << yc << ", sx=" << sx << ", sy=" << sy << ", beta=" << beta << ", err=" << errSum << '\n';
		return GSL_SUCCESS;
	}
	
	static int fitting_deriv_full(const gsl_vector *xvec, void *data, gsl_matrix *J)
	{
		GaussianFitter& fitter=*static_cast<GaussianFitter*>(data);
		const double scale = 1.0 / fitter._scaleFactor;
		double
			v = gsl_vector_get (xvec, 0),
			xc = gsl_vector_get (xvec, 1),
			yc = gsl_vector_get (xvec, 2),
			sx = gsl_vector_get (xvec, 3),
			sy = gsl_vector_get (xvec, 4),
			beta = gsl_vector_get (xvec, 5);
		if(fitter._posConstrained!=0.0 && (std::fabs(xc-fitter._xInit)>fitter._posConstrained*scale || std::fabs(yc-fitter._yInit)>fitter._posConstrained*scale))
		{
			std::cout << "GSL_EDOM\n";
			return GSL_EDOM;
		}
		const size_t width = fitter._width, height = fitter._height;
		int xMid = width/2, yMid = height/2;
		
		size_t dataIndex = 0;
		for(int yi=0; yi!=int(height); ++yi)
		{
			double y = yc + (yi - yMid)*scale;
			for(int xi=0; xi!=int(width); ++xi)
			{
				double x = xc + (xi - xMid)*scale;
				double expTerm = exp(-x*x/(2.0*sx*sx) - beta*x*y/(sx*sy) - y*y/(2.0*sy*sy));
				double dv = expTerm;
				expTerm *= v;
				double dx = (-beta*y/(sx*sy) - x/(sx*sx)) * expTerm;
				double dy = (-beta*x/(sy*sx) - y/(sy*sy)) * expTerm;
				double dsx = (beta*x*y/(sx*sx*sy)+x*x/(sx*sx*sx)) * expTerm;
				double dsy = (beta*x*y/(sy*sy*sx)+y*y/(sy*sy*sy)) * expTerm;
				double dbeta = -x*y/(sx*sy) * expTerm;
				gsl_matrix_set(J, dataIndex, 0, dv);
				gsl_matrix_set(J, dataIndex, 1, dx);
				gsl_matrix_set(J, dataIndex, 2, dy);
				gsl_matrix_set(J, dataIndex, 3, dsx);
				gsl_matrix_set(J, dataIndex, 4, dsy);
				gsl_matrix_set(J, dataIndex, 5, dbeta);
				++dataIndex;
			}
		}
		//std::cout << "diff: v=" << v << ", x=" << xc << ", y=" << yc << ", sx=" << sx << ", sy=" << sy << ", beta=" << beta << '\n';
		return GSL_SUCCESS;
	}
	
	static int fitting_both_full(const gsl_vector *x, void *data, gsl_vector *f, gsl_matrix *J)
	{
		fitting_func_full(x, data, f);
		fitting_deriv_full(x, data, J);
		return GSL_SUCCESS;
	}
	
	
	
	static int fitting_func_floor(const gsl_vector *xvec, void *data, gsl_vector *f)
	{
		GaussianFitter& fitter=*static_cast<GaussianFitter*>(data);
		const double scale = 1.0/fitter._scaleFactor;
		double
			v = gsl_vector_get (xvec, 0),
			xc = gsl_vector_get (xvec, 1),
			yc = gsl_vector_get (xvec, 2),
			sx = gsl_vector_get (xvec, 3),
			sy = gsl_vector_get (xvec, 4),
			beta = gsl_vector_get (xvec, 5),
			fl = gsl_vector_get (xvec, 6);
		if(fitter._posConstrained!=0.0 && (std::fabs(xc-fitter._xInit)>fitter._posConstrained*scale || std::fabs(yc-fitter._yInit)>fitter._posConstrained*scale))
			return GSL_EDOM;
		const size_t width = fitter._width, height = fitter._height;
		int xMid = width/2, yMid = height/2;
		
		size_t dataIndex = 0;
		double errSum = 0.0;
		for(int yi=0; yi!=int(height); ++yi)
		{
			double yS = yc + (yi - yMid) * scale;
			for(int xi=0; xi!=int(width); ++xi)
			{
				double xS = xc + (xi - xMid) * scale;
				double e = err_full(fitter._image[dataIndex], v, xS, yS, sx, sy, beta) + fl;
				errSum += e*e;
				gsl_vector_set(f, dataIndex, e);
				++dataIndex;
			}
		}
		//std::cout << "v=" << v << ", x=" << xc << ", y=" << yc << ", sx=" << sx << ", sy=" << sy << ", beta=" << beta << ", err=" << errSum << '\n';
		return GSL_SUCCESS;
	}
	
	static int fitting_deriv_floor(const gsl_vector *xvec, void *data, gsl_matrix *J)
	{
		GaussianFitter& fitter=*static_cast<GaussianFitter*>(data);
		double
			v = gsl_vector_get (xvec, 0),
			xc = gsl_vector_get (xvec, 1),
			yc = gsl_vector_get (xvec, 2),
			sx = gsl_vector_get (xvec, 3),
			sy = gsl_vector_get (xvec, 4),
			beta = gsl_vector_get (xvec, 5);
		const size_t width = fitter._width, height = fitter._height;
		int xMid = width/2, yMid = height/2;
		double scale = 1.0 / fitter._scaleFactor;
		
		size_t dataIndex = 0;
		for(int yi=0; yi!=int(height); ++yi)
		{
			double y = yc + (yi - yMid)*scale;
			for(int xi=0; xi!=int(width); ++xi)
			{
				double x = xc + (xi - xMid)*scale;
				double expTerm = exp(-x*x/(2.0*sx*sx) - beta*x*y/(sx*sy) - y*y/(2.0*sy*sy));
				double dv = expTerm;
				expTerm *= v;
				double dx = (-beta*y/(sx*sy) - x/(sx*sx)) * expTerm;
				double dy = (-beta*x/(sy*sx) - y/(sy*sy)) * expTerm;
				double dsx = (beta*x*y/(sx*sx*sy)+x*x/(sx*sx*sx)) * expTerm;
				double dsy = (beta*x*y/(sy*sy*sx)+y*y/(sy*sy*sy)) * expTerm;
				double dbeta = -x*y/(sx*sy) * expTerm;
				double dfl = 1.0;
				gsl_matrix_set(J, dataIndex, 0, dv);
				gsl_matrix_set(J, dataIndex, 1, dx);
				gsl_matrix_set(J, dataIndex, 2, dy);
				gsl_matrix_set(J, dataIndex, 3, dsx);
				gsl_matrix_set(J, dataIndex, 4, dsy);
				gsl_matrix_set(J, dataIndex, 5, dbeta);
				gsl_matrix_set(J, dataIndex, 6, dfl);
				++dataIndex;
			}
		}
		//std::cout << "diff: v=" << v << ", x=" << xc << ", y=" << yc << ", sx=" << sx << ", sy=" << sy << ", beta=" << beta << '\n';
		return GSL_SUCCESS;
	}
	
	static int fitting_both_floor(const gsl_vector *x, void *data, gsl_vector *f, gsl_matrix *J)
	{
		fitting_func_floor(x, data, f);
		fitting_deriv_floor(x, data, J);
		return GSL_SUCCESS;
	}
	
	static double dfdsx(double x, double y, double sx, double sy, double beta)
	{
		// f[x,y,sx,sy,beta]:=exp(-x*x/(2.0*sx*sx) - beta*x*y/(sx*sy) - y*y/(2.0*sy*sy));
		// diff(f[x,y,sx,sy,beta],sx);
		return (beta*x*y/(sx*sx*sy)+x*x/(sx*sx*sx)) * exp(-x*x/(2.0*sx*sx) - beta*x*y/(sx*sy) - y*y/(2.0*sy*sy));
	}
	
	static double dfdbeta(double x, double y, double sx, double sy, double beta)
	{
		return x*y/(sx*sy) * exp(-x*x/(2.0*sx*sx) - beta*x*y/(sx*sy) - y*y/(2.0*sy*sy));
	}
	
	static double dfdx(double x, double y, double sx, double sy, double beta)
	{
		// g[x,y,v,dx,dy,sx,sy,beta]:=v*exp(-(x+dx)*(x+dx)/(2.0*sx*sx) - beta*(x+dx)*(y+dy)/(sx*sy) - (y+dy)*(y+dy)/(2.0*sy*sy));
		// diff(g[x,y,v,dx,dy,sx,sy,beta],dx);
		return (-beta*y/(sx*sy) - x/(sx*sx)) *
			exp(-x*x/(2.0*sx*sx) - beta*x*y/(sx*sy) - y*y/(2.0*sy*sy));
	}
	
	static double dfdv(double x, double y, double sx, double sy, double beta)
	{
		// diff(g[x,y,v,dx,dy,sx,sy,beta],v);
		return exp(-x*x/(2.0*sx*sx) - beta*x*y/(sx*sy) - y*y/(2.0*sy*sy));
	}
	
	void convertShapeParameters(double sx, double sy, double beta, double& ellipseMaj, double& ellipseMin, double& ellipsePA)
	{
		const long double sigmaToBeam = 2.0L * sqrtl(2.0L * logl(2.0L));
		const double betaFact = 1.0 - beta*beta;
		double cov[4];
		cov[0] = sx*sx / betaFact;
		cov[1] = beta * sx*sy / betaFact;
		cov[2] = cov[1];
		cov[3] = sy*sy / betaFact;
		
		double e1, e2, vec1[2], vec2[2];
		Matrix2x2::EigenValuesAndVectors(cov, e1, e2, vec1, vec2);
		ellipseMaj = sqrt(std::fabs(e1)) * sigmaToBeam * _scaleFactor;
		ellipseMin = sqrt(std::fabs(e2)) * sigmaToBeam * _scaleFactor;
		ellipsePA = atan2(vec1[0], vec1[1]);
	}
	
	double _xInit, _yInit, _posConstrained;
};

#endif
