
#include "MatteoCellCycleModel.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "Debug.hpp"

MatteoCellCycleModel::MatteoCellCycleModel()
    : AbstractSimpleCellCycleModel(),
      mMinCellCycleDuration(12.0), // Hours
      mMaxCellCycleDuration(14.0)  // Hours
{
}

/*called every time step and tells if a cell is aready to divide (using the cell cycle model)*/

bool MatteoCellCycleModel::ReadyToDivide()
{
   bool ready = mpCell->GetCellData()->GetItem("divide");
    mpCell->GetCellData()->SetItem("divide", 0);
    return ready;
    

}

MatteoCellCycleModel::MatteoCellCycleModel(const MatteoCellCycleModel& rModel)
   : AbstractSimpleCellCycleModel(rModel),
     mMinCellCycleDuration(rModel.mMinCellCycleDuration),
     mMaxCellCycleDuration(rModel.mMaxCellCycleDuration)
{
}

AbstractCellCycleModel* MatteoCellCycleModel::CreateCellCycleModel()
{
    return new MatteoCellCycleModel(*this);
}

void MatteoCellCycleModel::SetCellCycleDuration()
{
}

double MatteoCellCycleModel::GetMinCellCycleDuration()
{
    return mMinCellCycleDuration;
}

void MatteoCellCycleModel::SetMinCellCycleDuration(double minCellCycleDuration)
{
    mMinCellCycleDuration = minCellCycleDuration;
}

double MatteoCellCycleModel::GetMaxCellCycleDuration()
{
    return mMaxCellCycleDuration;
}

void MatteoCellCycleModel::SetMaxCellCycleDuration(double maxCellCycleDuration)
{
    mMaxCellCycleDuration = maxCellCycleDuration;
}

double MatteoCellCycleModel::GetAverageTransitCellCycleTime()
{
    return 0.5*(mMinCellCycleDuration + mMaxCellCycleDuration);
}

double MatteoCellCycleModel::GetAverageStemCellCycleTime()
{
    return 0.5*(mMinCellCycleDuration + mMaxCellCycleDuration);
}

void MatteoCellCycleModel::OutputCellCycleModelParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<MinCellCycleDuration>" << mMinCellCycleDuration << "</MinCellCycleDuration>\n";
    *rParamsFile << "\t\t\t<MaxCellCycleDuration>" << mMaxCellCycleDuration << "</MaxCellCycleDuration>\n";

    // Call method on direct parent class
    AbstractSimpleCellCycleModel::OutputCellCycleModelParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(MatteoCellCycleModel)
