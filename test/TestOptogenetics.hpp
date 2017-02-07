

#ifndef TESTOPTOGENETICS_HPP_
#define TESTOPTOGENETICS_HPP_


#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "OffLatticeSimulation.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "MatteoForce.hpp"
#include "ConstantTargetAreaModifier.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "WildTypeCellMutationState.hpp"
#include "SmartPointers.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "NoCellCycleModel.hpp"

/*
 * The next header file defines a simple subcellular reaction network model that includes the functionality
 * for solving each cell's Delta/Notch signalling ODE system at each time step, using information about neighbouring
 * cells through the {{{CellData}}} class.
 */
#include "MatteoSrnModel.hpp"
/*
 * The next header defines the simulation class modifier corresponding to the Delta-Notch SRN model.
 * This modifier leads to the {{{CellData}}} cell property being updated at each timestep to deal with Delta-Notch signalling.
 */
#include "DeltaNotchTrackingModifier.hpp"

/* Having included all the necessary header files, we proceed by defining the test class.
 */
class TestOptogenetics : public AbstractCellBasedTestSuite
{
public:

    void TestVertexBasedMonolayerWithDeltaNotch() throw (Exception)
    {
        /* We include the next line because Vertex simulations cannot be run in parallel */
        EXIT_IF_PARALLEL;

        /* First we create a regular vertex mesh. */
        HoneycombVertexMeshGenerator generator(10, 10);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);

        for (unsigned elem_index=0; elem_index<p_mesh->GetNumElements(); elem_index++)
        {
            NoCellCycleModel* p_cc_model = new NoCellCycleModel();
            p_cc_model->SetDimension(2);

            /* We choose to initialise the concentrations to random levels in each cell. */
            std::vector<double> initial_conditions;
            initial_conditions.push_back(RandomNumberGenerator::Instance()->ranf());
            initial_conditions.push_back(RandomNumberGenerator::Instance()->ranf());
            MatteoSrnModel* p_srn_model = new MatteoSrnModel();
            p_srn_model->SetInitialConditions(initial_conditions);

            CellPtr p_cell(new Cell(p_state, p_cc_model, p_srn_model));
            p_cell->SetCellProliferativeType(p_diff_type);
            double birth_time = -RandomNumberGenerator::Instance()->ranf()*12.0;
            p_cell->SetBirthTime(birth_time);
            cells.push_back(p_cell);
        }

        /* Using the vertex mesh and cells, we create a cell-based population object, and specify which results to
         * output to file. */
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        /* We are now in a position to create and configure the cell-based simulation object, pass a force law to it,
         * and run the simulation. We can make the simulation run for longer to see more patterning by increasing the end time. */
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestOptogenetics");
        simulator.SetSamplingTimestepMultiple(10);
        simulator.SetEndTime(11.0);

        /* Then, we define the modifier class, which automatically updates the values of Delta and Notch within the cells in {{{CellData}}} and passes it to the simulation.*/
        MAKE_PTR(DeltaNotchTrackingModifier<2>, p_modifier);
        simulator.AddSimulationModifier(p_modifier);

        MAKE_PTR(MatteoForce<2>, p_force);
        simulator.AddForce(p_force);

        /* This modifier assigns target areas to each cell, which are required by the {{{MatteoHondaForce}}}.
         */
        MAKE_PTR(ConstantTargetAreaModifier<2>, p_growth_modifier);
        simulator.AddSimulationModifier(p_growth_modifier);
        simulator.Solve();
    }
};

#endif /*TESTOPTOGENETICS_HPP_*/
