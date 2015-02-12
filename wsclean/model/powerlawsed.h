#ifndef POWER_LAW_SED_H
#define POWER_LAW_SED_H

#include "spectralenergydistribution.h"

class PowerLawSED : public SpectralEnergyDistribution
{
public:
	PowerLawSED() : _referenceFrequency(0.0)
	{
		for(size_t p=0; p!=4; ++p) _factors[p] = 0.0;
	}
	
	virtual PowerLawSED* Clone() const
	{
		return new PowerLawSED(*this);
	}
	
	virtual std::string ToString() const
	{
		std::ostringstream str;
		str << "PowerLawSED, reference frequency=" << _referenceFrequency*1e-6 << " MHz, brightness=["
			<< _factors[0] << ", " << _factors[1] << ", " << _factors[2] << ", " << _factors[3]
			<< "], SI terms: {";
		if(!_terms.empty())
		{
			str << _terms[0];
			for(size_t i=1; i!=_terms.size(); ++i)
				str << ", " << _terms[i];
		}
		str << "}";
		return str.str();
	}
	
	virtual long double FluxAtFrequencyFromIndex(long double frequencyHz, size_t pIndex) const
	{
		return NonLinearPowerLawFitter::Evaluate(frequencyHz/_referenceFrequency, _terms) * _factors[pIndex];
	}
	
	virtual long double IntegratedFlux(long double startFrequency, long double endFrequency, PolarizationEnum polarization) const
	{
		size_t pIndex = Polarization::StokesToIndex(polarization);
		long double sum = 0.0;
		for(size_t i=0; i!=101; ++i)
		{
			long double frequency = startFrequency + (endFrequency - startFrequency)*double(i)/100.0;
			sum += FluxAtFrequencyFromIndex(frequency, pIndex);
		}
		return sum / 101.0;
	}
	
	virtual long double AverageFlux(long double startFrequency, long double endFrequency, PolarizationEnum polarization) const
	{
		return IntegratedFlux(startFrequency, endFrequency, polarization);
	}
	
	virtual bool operator<(const SpectralEnergyDistribution &other) const
	{
		return other.FluxAtFrequencyFromIndex(_referenceFrequency, 0) < FluxAtFrequencyFromIndex(_referenceFrequency, 0);
	}
	
	virtual void operator*=(double factor)
	{
		for(size_t p=0; p!=4; ++p)
			_factors[p] *= factor;
	}
	
	virtual void operator+=(const SpectralEnergyDistribution &other)
	{
		throw std::runtime_error("operator+= not yet implemented for power law sed");
	}
	
	void SetData(double referenceFrequency, const double* brightnessVector, const std::vector<double>& siTerms)
	{
		_referenceFrequency = referenceFrequency;
		double refBrightness = brightnessVector[0];
		if(refBrightness <= 0.0)
			refBrightness = 1.0;
		_terms.resize(siTerms.size()+1);
		_terms[0] = NonLinearPowerLawFitter::FactorToTerm0(refBrightness, siTerms[0]);
		for(size_t i=0; i!=siTerms.size(); ++i)
			_terms[i+1] = siTerms[i];
		for(size_t p=0; p!=4; ++p)
			_factors[p] = brightnessVector[p] / refBrightness;
		std::cout << ToString() << '\n';
	}
private:
	double _referenceFrequency;
	double _factors[4];
	std::vector<double> _terms;
};

#endif
