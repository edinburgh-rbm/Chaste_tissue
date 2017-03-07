
#include "MatteoForce.hpp"
#include "PopulationConstants.hpp"

template<unsigned DIM>
MatteoForce<DIM>::MatteoForce()
   : FarhadifarForce<DIM>()
     {
}

template<unsigned DIM>
MatteoForce<DIM>::~MatteoForce()
{
}


template<unsigned DIM>
double MatteoForce<DIM>::GetLineTensionParameter(unsigned elem_index, Node<DIM>* pNodeA, Node<DIM>* pNodeB, VertexBasedCellPopulation<DIM>& rVertexCellPopulation)
{
    // Find the indices of the elements owned by each node
    std::set<unsigned> elements_containing_nodeA = pNodeA->rGetContainingElementIndices();
    std::set<unsigned> elements_containing_nodeB = pNodeB->rGetContainingElementIndices();

    // Find common elements
    std::set<unsigned> shared_elements;
    std::set_intersection(elements_containing_nodeA.begin(),
                          elements_containing_nodeA.end(),
                          elements_containing_nodeB.begin(),
                          elements_containing_nodeB.end(),
                          std::inserter(shared_elements, shared_elements.begin()));

    // Check that the nodes have a common edge
    assert(!shared_elements.empty());

    double tension;

    if (shared_elements.size() == 1) {
	CellPtr c = rVertexCellPopulation.GetCellUsingLocationIndex(elem_index);
	if (c->GetCellProliferativeType() == p_wild_type) {
	    tension = wild_type_lambda;
	} else {
	    tension = diff_type_lambda;
	}
    } else {
	unsigned n_wild = 0, n_diff = 0;

	for (std::set<unsigned>::iterator iter = shared_elements.begin(); iter != shared_elements.end(); ++iter) {
	    unsigned i = *(iter);
	    CellPtr c = rVertexCellPopulation.GetCellUsingLocationIndex(i);
	    if (c->GetCellProliferativeType() == p_wild_type) {
		n_wild++;
	    } else {
		n_diff++;
	    }
	}
	
	assert(n_wild + n_diff == 2);
	if (n_wild == 2) {
	    tension = wild_type_lambda;
	} else if (n_diff == 2) {
	    tension = diff_type_lambda;
	} else {
	    tension = mixed_type_lambda;
	}

	// If the edge corresponds to a single element, then the cell is on the boundary
	// if not on the boundary it will be visited twice.
	tension /= 2;
    }

    return tension;
}


template<unsigned DIM>
void MatteoForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    // Call method on direct parent class
    FarhadifarForce<DIM>::OutputForceParameters(rParamsFile);
}

// Explicit instantiation
template class MatteoForce<1>;
template class MatteoForce<2>;
template class MatteoForce<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(MatteoForce)
