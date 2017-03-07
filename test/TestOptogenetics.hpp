

#ifndef TESTOPTOGENETICS_HPP_
#define TESTOPTOGENETICS_HPP_

#include <cxxtest/TestSuite.h>
#include <boost/program_options.hpp>
#include <boost/format.hpp>
#include "CheckpointArchiveTypes.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "OffLatticeSimulation.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "MatteoForce.hpp"
#include "RandomMotionForce.hpp"
#include "ConstantTargetAreaModifier.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "WildTypeCellMutationState.hpp"
#include "SmartPointers.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "NoCellCycleModel.hpp"
#include "CommandLineArguments.hpp"
#include "PopulationConstants.hpp"
#include "CellProliferativeTypesWriter.hpp"

namespace po = boost::program_options;

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
class TestOptogenetics : public AbstractCellBasedTestSuite {
protected:
  po::variables_map args;

  void ProcessCommandLineArguments() throw (Exception) {
    /* Process Command Line Arguments */
    po::options_description description("Chaste Tissue Optogenetics Test Usage");

    description.add_options()
	("help,h", "Display this help message")
	("lambda,l", po::value<double>()->default_value(0.12),"Rigidity of wild-type boundary")
	("diff,m", po::value<double>()->default_value(0.12), "Rigidity of differentiated boundary")
	("mixed,x", po::value<double>()->default_value(0.0), "Rigidity of mixed boundary")
	("proportion,p", po::value<double>()->default_value(0.1), "Proportion of population that is mutant")
	("number,n", po::value<unsigned>()->default_value(16), "sqrt(number of cells)")
        ("noise,z", po::value<double>()->default_value(0.05), "Noise parameter")
	("dt,d", po::value<double>()->default_value(1.0/200.0), "Simulation time step")
	("sample,s", po::value<unsigned>()->default_value(200), "Sampling time step multiple")
	("time,t", po::value<double>()->default_value(10.0), "Simulation end time");

    int argc = *(CommandLineArguments::Instance()->p_argc);
    TS_ASSERT_LESS_THAN(0, argc); // argc should always be 1 or greater
    char** argv = *(CommandLineArguments::Instance()->p_argv);
    assert(argv != NULL);

    po::store(po::command_line_parser(argc, argv).options(description).run(), args);
    po::notify(args);

    if(args.count("help")){
      std::cout << description;
      exit(0);
    }

    wild_type_lambda = args["lambda"].as<double>();
    diff_type_lambda = args["diff"].as<double>();
    mixed_type_lambda = args["mixed"].as<double>();
  }

public:

  void TestVertexBasedMonolayerWithDeltaNotch() throw (Exception) {
    /* We include the next line because Vertex simulations cannot be run in parallel */
    EXIT_IF_PARALLEL;

    ProcessCommandLineArguments();

    /* First we create a regular vertex mesh. */
    HoneycombVertexMeshGenerator generator(args["number"].as<unsigned>(), args["number"].as<unsigned>());
    MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();
    p_mesh->SetCellRearrangementThreshold(0.1);

    std::vector<CellPtr> cells;
    MAKE_PTR(WildTypeCellMutationState, p_state);

    for (unsigned elem_index=0; elem_index<p_mesh->GetNumElements(); elem_index++) {
      NoCellCycleModel* p_cc_model = new NoCellCycleModel();
      p_cc_model->SetDimension(2);

      /* We choose to initialise the concentrations to random levels in each cell. */
      std::vector<double> initial_conditions;
      initial_conditions.push_back(RandomNumberGenerator::Instance()->ranf());
      initial_conditions.push_back(RandomNumberGenerator::Instance()->ranf());
      MatteoSrnModel* p_srn_model = new MatteoSrnModel();
      p_srn_model->SetInitialConditions(initial_conditions);

      CellPtr p_cell(new Cell(p_state, p_cc_model, p_srn_model));
      double birth_time = -RandomNumberGenerator::Instance()->ranf()*12.0;
      p_cell->SetBirthTime(birth_time);

      // Set a target area rather than setting a growth modifier. (the modifiers don't work correctly as making very long G1 phases)
      p_cell->GetCellData()->SetItem("target area", 1.0);

      if (RandomNumberGenerator::Instance()->ranf() < args["proportion"].as<double>()) {
        p_cell->SetCellProliferativeType(p_diff_type);
      } else {
        p_cell->SetCellProliferativeType(p_wild_type);
      }

      cells.push_back(p_cell);
    }

    /* Using the vertex mesh and cells, we create a cell-based population object, and specify which results to
     * output to file. */
    VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);
    cell_population.AddCellWriter<CellProliferativeTypesWriter>();

    double noise = args["noise"].as<double>();

    /* We are now in a position to create and configure the cell-based simulation object, pass a force law to it,
     * and run the simulation. We can make the simulation run for longer to see more patterning by increasing the end time. */
    OffLatticeSimulation<2> simulator(cell_population);

    simulator.SetOutputDirectory(boost::str(boost::format("Optogenetics-l%1%-m%2%-x%3%-z%4%") % wild_type_lambda % diff_type_lambda % mixed_type_lambda % noise));

    /* set up the timing */
    simulator.SetDt(args["dt"].as<double>());
    simulator.SetSamplingTimestepMultiple(args["sample"].as<unsigned>());
    simulator.SetEndTime(args["time"].as<double>());

    /* Then, we define the modifier class, which automatically updates the values of Delta and Notch within the cells in {{{CellData}}} and passes it to the simulation.*/
    MAKE_PTR(DeltaNotchTrackingModifier<2>, p_modifier);
    simulator.AddSimulationModifier(p_modifier);

    MAKE_PTR(MatteoForce<2>, p_force);
    simulator.AddForce(p_force);

    // Add some noise to avoid local minimum
    MAKE_PTR(RandomMotionForce<2>, p_random_force);
    p_random_force->SetMovementParameter(noise);
    simulator.AddForce(p_random_force);

    /* This modifier assigns target areas to each cell, which are required by the {{{MatteoHondaForce}}}.
     */
    MAKE_PTR(ConstantTargetAreaModifier<2>, p_growth_modifier);
    simulator.AddSimulationModifier(p_growth_modifier);
    simulator.Solve();
  }
};

#endif /*TESTOPTOGENETICS_HPP_*/
