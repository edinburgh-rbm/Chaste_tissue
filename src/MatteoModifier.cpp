
#include "MatteoModifier.hpp"
#include "RandomNumberGenerator.hpp"
#include "CellLabel.hpp"
#include "Debug.hpp"

template<unsigned DIM>
MatteoModifier<DIM>::MatteoModifier()
    : AbstractCellBasedSimulationModifier<DIM>()
{
}

template<unsigned DIM>
MatteoModifier<DIM>::~MatteoModifier()
{
}

template<unsigned DIM>
void MatteoModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
   

    //do every 10 units of time
    if (fmod(SimulationTime::Instance()->GetTime(),10) == 0)
    {
        // Store each cell's current fitness
        UpdateCellData(rCellPopulation);

        // Randomly pick one cell to divide
        std::map<CellPtr, double> map;
        double prev = rCellPopulation.Begin()->GetCellData()->GetItem("fitness");
        for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
        {
            if (*cell_iter == *(rCellPopulation.Begin()))
            {
                map[*cell_iter] = prev;
            }
            else
            {
                map[*cell_iter] = prev + cell_iter->GetCellData()->GetItem("fitness");
                prev +=cell_iter->GetCellData()->GetItem("fitness");
            }
        }
        double total_fitness = prev;
        for (std::map<CellPtr, double>::iterator iter = map.begin(); iter!=map.end(); ++iter)
        {
            iter->second /= total_fitness;
        }

        double r = RandomNumberGenerator::Instance()->ranf();
        for (std::map<CellPtr, double>::iterator iter = map.begin(); iter!=map.end(); ++iter)
        {
	std::map<CellPtr, double>::iterator other_iter = iter;
	++other_iter;
            if ((iter->second < r) && (other_iter->second >= r))
            {
                TellCellToDivide(iter->first);
                break;
            }
        }
    }


}

template<unsigned DIM>
void MatteoModifier<DIM>::TellCellToDivide(CellPtr pCell)
{
pCell->GetCellData()->SetItem("divide", 1);
}

template<unsigned DIM>
void MatteoModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
}

template<unsigned DIM>
void MatteoModifier<DIM>::UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
double b = 10;
double c = 5;
double delta = 0.01;
    // Iterate over cell population
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        bool is_defector = cell_iter->template HasCellProperty<CellLabel>();
        double payoff = 0;
        std::set<unsigned> neighbours = rCellPopulation.GetNeighbouringLocationIndices(*cell_iter);
        for (std::set<unsigned>::iterator iter = neighbours.begin();
             iter != neighbours.end();
             ++iter)
        {
            CellPtr p_neighbour = rCellPopulation.GetCellUsingLocationIndex(*iter);
            bool is_neighbour_defector = p_neighbour->template HasCellProperty<CellLabel>();

            double contribution = b - c;
            if (is_defector && !is_neighbour_defector)
            {
                contribution = b;
            }
            else if (!is_defector && is_neighbour_defector)
            {
                contribution = -c;
            }
            else if (is_defector && is_neighbour_defector)
            {
                contribution = 0;
            }
            payoff += contribution;
        }
        // Get the fitness of this cell
        double cell_fitness = pow(1 + delta, payoff);
        
        // Store the cell's fitness in CellData
        cell_iter->GetCellData()->SetItem("fitness", cell_fitness);
    }
}

template<unsigned DIM>
void MatteoModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class MatteoModifier<1>;
template class MatteoModifier<2>;
template class MatteoModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(MatteoModifier)

