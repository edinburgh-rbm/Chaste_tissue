#include "PopulationConstants.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "DefaultCellProliferativeType.hpp"

boost::shared_ptr<AbstractCellProperty> p_wild_type(new DefaultCellProliferativeType());
boost::shared_ptr<AbstractCellProperty> p_diff_type(new DifferentiatedCellProliferativeType());
double wild_type_lambda = 0.12;
double diff_type_lambda = 0.12;
double mixed_type_lambda = 0.2;

