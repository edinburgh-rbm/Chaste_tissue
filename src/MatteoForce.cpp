
#include "MatteoForce.hpp"

int ntypes;
double **costs;
double **demographics;

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
    CellPtr c = rVertexCellPopulation.GetCellUsingLocationIndex(elem_index);
    int colour = c->GetCellData()->GetItem("cell type");

    if (shared_elements.size() == 1) {
      tension = costs[colour][colour];
    } else {
      int ocolour = colour;

      for (std::set<unsigned>::iterator iter = shared_elements.begin(); iter != shared_elements.end(); ++iter) {
	unsigned i = *(iter);
	CellPtr oc = rVertexCellPopulation.GetCellUsingLocationIndex(i);
	int o = oc->GetCellData()->GetItem("cell type");
	// in this loop, we will visit the current cell as well, so
	// only remember colour if it is different
	if (o != colour)
	  ocolour = o;

	tension = costs[colour][ocolour];
	
	// If the edge corresponds to a single element, then the cell is on the boundary
	// if not on the boundary it will be visited twice.
	tension /= 2;
      }

    }
    return tension;
}

template<unsigned DIM>
void MatteoForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
  *rParamsFile << "\t\t\t<Costs>" << ntypes << " " << ntypes << "\n";
  for (int i=0; i<ntypes; i++) {
    for (int j=0; j<ntypes; j++) {
      *rParamsFile << costs[i][j];
      if (j < ntypes-1) {
	*rParamsFile << " ";
      } else {
	*rParamsFile << "\n";
      }
    }
  }
  *rParamsFile << "</Costs>\n";

  *rParamsFile << "\t\t\t<Demographics>" << ntypes << " 1\n";		  
  for (int i=0; i<ntypes; i++) {
    *rParamsFile << demographics[i][0] << "\n";
  }
  *rParamsFile << "</Demographics>\n";

  FarhadifarForce<DIM>::OutputForceParameters(rParamsFile);
}

// Explicit instantiation
template class MatteoForce<1>;
template class MatteoForce<2>;
template class MatteoForce<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(MatteoForce)
