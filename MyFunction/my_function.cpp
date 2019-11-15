#include "my_function.h"

#define _USE_MATH_DEFINES
#include <math.h>

// ------------------------------------------------------------------------------------------------
TMyProblem::TMyProblem()
{
  mIsInitialized = false;
  mDimension = 1;
}

// ------------------------------------------------------------------------------------------------
int TMyProblem::SetConfigPath(const std::string& configPath)
{
  return IProblem::OK;
}

// ------------------------------------------------------------------------------------------------
int TMyProblem::SetDimension(int dimension)
{
  if (dimension > 0 && dimension <= mMaxDimension)
  {
    mDimension = dimension;
    return IProblem::OK;
  }
  else
    return IProblem::ERROR;
}

// ------------------------------------------------------------------------------------------------
int TMyProblem::GetDimension() const
{
  return mDimension;
}

// ------------------------------------------------------------------------------------------------
int TMyProblem::Initialize()
{
  if (mDimension > 0)
  {
    mIsInitialized = true;
    return IProblem::OK;
  }
  else
    return IProblem::ERROR;
}

// ------------------------------------------------------------------------------------------------
void TMyProblem::GetBounds(double* lower, double *upper)
{
  if (mIsInitialized)
    for (int i = 0; i < mDimension; i++)
    {
      lower[i] = -1.7;
      upper[i] = 7.1;
    }
}

// ------------------------------------------------------------------------------------------------
int TMyProblem::GetOptimumValue(double& value) const
{
  if (!mIsInitialized)
    return IProblem::UNDEFINED;

  value = 0.0;
  return IProblem::OK;
}

// ------------------------------------------------------------------------------------------------
int TMyProblem::GetOptimumPoint(double* point) const
{
  if (!mIsInitialized)
    return IProblem::UNDEFINED;

  for (int i = 0; i < mDimension; i++)
    point[i] = 0.0;
  return IProblem::OK;
}

// ------------------------------------------------------------------------------------------------
int TMyProblem::GetNumberOfFunctions() const
{
  return 1;
}

// ------------------------------------------------------------------------------------------------
int TMyProblem::GetNumberOfConstraints() const
{
  return 0;
}

// ------------------------------------------------------------------------------------------------
int TMyProblem::GetNumberOfCriterions() const
{
  return 1;
}

// ------------------------------------------------------------------------------------------------
double TMyProblem::CalculateFunctionals(const double* x, int fNumber)
{
  double sum = 0.;
  for (int j = 0; j < mDimension; j++)
    sum += x[j] * log10(abs(x[j]) * 2) - exp(abs(x[j]) * 1.1 - x[j])*sin(M_PI * x[j] / 3);
  return sum;
}

// ------------------------------------------------------------------------------------------------
TMyProblem::~TMyProblem()
{

}

// ------------------------------------------------------------------------------------------------
LIB_EXPORT_API IProblem* create()
{
  return new TMyProblem();
}

// ------------------------------------------------------------------------------------------------
LIB_EXPORT_API void destroy(IProblem* ptr)
{
  delete ptr;
}
// - end of file ----------------------------------------------------------------------------------