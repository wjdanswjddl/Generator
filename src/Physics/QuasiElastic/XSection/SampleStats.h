#pragma once
// Computes running values of sample statistics.
// Based on http://tinyurl.com/mean-var-onl-alg
class SampleStats {

  public:

    SampleStats() : fN(0), fMean( 0. ), fM2( 0. ) {}

    void AddValue(double value) {
      ++fN;
      double delta = value - fMean;
      fMean += delta / fN;
      double delta2 = value - fMean;
      fM2 += delta * delta2;
    }

    inline double Mean() const { return fMean; }
    //inline double FinitePopulationVariance() const { return fM2 / fN; }
    inline double SampleVariance() const { return fM2 / (fN - 1); }
    inline double SampleStdDeviation() const
      { return std::sqrt( this->SampleVariance() ); }
    inline int SampleSize() const { return fN; }
    inline double StdErrorOnMean() const
      { return this->SampleStdDeviation() / std::sqrt(fN); }

  protected:

    int fN; // sample size
    double fMean;
    double fM2;
};
