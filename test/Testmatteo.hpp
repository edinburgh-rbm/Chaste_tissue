#ifndef TESTRUNNINGDIFFERENTIALADHESIONSIMULATIONSTUTORIAL_HPP_
#define TESTRUNNINGDIFFERENTIALADHESIONSIMULATIONSTUTORIAL_HPP_

#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "ToroidalHoneycombVertexMeshGenerator.hpp"
#include "Toroidal2dVertexMesh.hpp"
#include "CellsGenerator.hpp"
#include "MatteoCellCycleModel.hpp"
#include "MatteoForce.hpp"
#include "CellLabel.hpp"
#include "PopulationConstants.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "OffLatticeSimulation.hpp"
#include "SmartPointers.hpp"
#include "FakePetscSetup.hpp"
#include "FarhadifarForce.hpp"
#include "ConstantTargetAreaModifier.hpp"
#include "MatteoModifier.hpp"
#include "CellLabelWriter.hpp"

class Testmatteo : public AbstractCellBasedTestSuite
{
public:

    void TestVertexBasedDifferentialAdhesionSimulation() throw (Exception)
    {
        /* First we create a regular vertex mesh. Here we choose to set the value of the cell rearrangement threshold. */
        HoneycombVertexMeshGenerator generator(20, 20);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();
        //ToroidalHoneycombVertexMeshGenerator generator(5, 5);
        //Toroidal2dVertexMesh* p_mesh = generator.GetToroidalMesh();
        p_mesh->SetCellRearrangementThreshold(0.1);

        /* We then create some cells using the helper class {{{CellsGenerator}}}. Note that in this simulation
         * the cells are all differentiated, and thus no cell division occurs; if we wished, we could modify
         * the three lines below in a straightforward manner to incorporate cell proliferation and investigate
         * the effect of this on the cell sorting process. */
        std::vector<CellPtr> cells;
        CellsGenerator<MatteoCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumElements(), std::vector<unsigned>(), p_diff_type);

       for(unsigned i=0; i<cells.size();i++)
         {
           cells[i]->GetCellData()->SetItem("divide", 0);
           cells[i]->GetCellData()->SetItem("fitness", 1);
         }

        /* Using the vertex mesh and cells, we create a cell-based population object, and specify which results to
         * output to file. */
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.AddCellWriter<CellLabelWriter>();

        /* We randomly label some cells using the cell property {{{CellLabel}}}. We begin by creating a shared pointer to
         * this cell property using the helper singleton {{{CellPropertyRegistry}}}. We then loop over the cells and label
         * each cell independently with probability 0.5. Note that since the cells have been passed to the
         * {{{VertexBasedCellPopulation}}} object, the vector {{{cells}}} above is now empty, so we must use the
         * {{{Iterator}}} to loop over cells. */
         boost::shared_ptr<AbstractCellProperty> p_label(CellPropertyRegistry::Instance()->Get<CellLabel>());
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            if (RandomNumberGenerator::Instance()->ranf() < 0.1)
            {
                cell_iter->AddCellProperty(p_label);
            }
        }

        /* We are now in a position to create and configure the cell-based simulation object.
         * We can make the simulation run for longer to see more cell sorting by increasing the end time. */
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestMatteo");
        simulator.SetSamplingTimestepMultiple(10);
        simulator.SetEndTime(11.0);

        /* Next we create the differential adhesion force law. This builds upon the model of Nagai, Honda and co-workers
         * encounted in the TestRunningVertexBasedSimulationsTutorial by allowing different values of the adhesion
         * energy parameters depending on the types of two neighbouring cells. Here we interpret the 'type' of a cell
         * as whether or not it has the cell property {{{CellLabel}}}; it would be straightforward to create a similar
         * force law that took account of a cell's mutation state, for example. Having created the force law, we set the
         * values of the parameters. If the adhesion energy for two neighbouring homotypic cells is less than that of two
         * heterotypic cells, then we may expect cell sorting to occur, in which the cells of each type will tend to locally
         * aggregate over time. */
        MAKE_PTR(MatteoForce<2>, p_force);
        simulator.AddForce(p_force);

        /* A {{{NagaiHondaForceDifferentialAdhesionForce}}} assumes that each cell has been assigned a target area.
         * The {{{ConstantTargetAreaModifier}}} will assign and update the target areas of all cells.
         */
        MAKE_PTR(ConstantTargetAreaModifier<2>, p_growth_modifier);
        simulator.AddSimulationModifier(p_growth_modifier);

        MAKE_PTR(MatteoModifier<2>, p_modifier);
        simulator.AddSimulationModifier(p_modifier);

        /* Finally, we run the simulation. */
        simulator.Solve();
    }

    /*
     * EMPTYLINE
     *
     * To visualize the results, use Paraview. See the UserTutorials/VisualizingWithParaview tutorial for more information.
     *
     * Load the file {{{/tmp/$USER/testoutput/TestVertexBasedDifferentialAdhesionSimulation/results_from_time_0/results.pvd}}}.
     */
};

#endif /*TESTRUNNINGDIFFERENTIALADHESIONSIMULATIONSTUTORIAL_HPP_*/
