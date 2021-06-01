#include "Tasks.h"
#include <iostream>
#include "../sample_src/grishagin_function.hpp"

// ------------------------------------------------------------------------------------------------
MyProblem::MyProblem()
{
}

// ------------------------------------------------------------------------------------------------
int MyProblem::SetConfigPath(const std::string& configPath)
{
  //mConfigPath = std::string(configPath);
  return IProblem::OK;
}

// ------------------------------------------------------------------------------------------------
int MyProblem::SetDimension(int dimension)
{
  /*if (dimension > 0 && dimension <= mMaxDimension)
  {
    mDimension = dimension;
    return IProblem::OK;
  }
  else
    return IProblem::ERROR;*/
  dim = dimension;
  return IProblem::OK;
}

// ------------------------------------------------------------------------------------------------
int MyProblem::GetDimension() const
{
  return BigProblem->GetDimension();
}

// ------------------------------------------------------------------------------------------------
int MyProblem::Initialize()
{
  BigProblem = new TGrishaginProblem(10);
  std::cout << "Grishagin Problem 1" << std::endl;
  std::cout << "GrishaginProblem ( 0.5, 0.5 ) = " << BigProblem->ComputeFunction({ 0.5, 0.5 }) << std::endl;
  std::cout << "GrishaginProblem Derivatives ( 0.5, 0.5 ) = {" <<
    BigProblem->ComputeFunctionDerivatives({ 0.5, 0.5 })[0] << ", " <<
    BigProblem->ComputeFunctionDerivatives({ 0.5, 0.5 })[1] << "}" << std::endl;
  std::cout << std::endl;
  return IProblem::OK;
}

// ------------------------------------------------------------------------------------------------
void MyProblem::GetBounds(double* lower, double *upper)
{
  std::vector<double> l(dim);
  std::vector<double> u(dim);
  BigProblem->GetBounds(l, u);
  for (int i = 0; i < dim; i++) {
    lower[i] = l[i];
    upper[i] = u[i];
  }
}

// ------------------------------------------------------------------------------------------------
int MyProblem::GetOptimumValue(double& value) const
{
  value = BigProblem->GetOptimumValue();
  return IProblem::OK;
}

// ------------------------------------------------------------------------------------------------
int MyProblem::GetOptimumPoint(double* point) const
{
  std::vector<double> x = BigProblem->GetOptimumPoint();
  for (int i = 0; i < dim; i++)
    point[i] = x[i];

  return IProblem::OK;
}

// ------------------------------------------------------------------------------------------------
int MyProblem::GetNumberOfFunctions() const
{
  return 1;
}

// ------------------------------------------------------------------------------------------------
int MyProblem::GetNumberOfConstraints() const
{
  return 0;
}

// ------------------------------------------------------------------------------------------------
int MyProblem::GetNumberOfCriterions() const
{
  return 1;
}

// ------------------------------------------------------------------------------------------------
double MyProblem::CalculateFunctionals(const double* x, int fNumber)
{
  std::vector<double> f;
  for (int i = 0; i < dim; i++)
    f.push_back(x[i]);

  return BigProblem->ComputeFunction(f);
}

// ------------------------------------------------------------------------------------------------
MyProblem::~MyProblem()
{

}

// ------------------------------------------------------------------------------------------------
LIB_EXPORT_API IProblem* create()
{
  return new MyProblem();
}

// ------------------------------------------------------------------------------------------------
LIB_EXPORT_API void destroy(IProblem* ptr)
{
  delete ptr;
}
// - end of file ----------------------------------------------------------------------------------