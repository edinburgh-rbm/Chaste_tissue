
#ifndef MatteoFORCE_HPP_
#define MatteoFORCE_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include "Exception.hpp"

#include "FarhadifarForce.hpp"
#include "VertexBasedCellPopulation.hpp"

#include <iostream>


extern double **costs;
extern double **demographics;
extern int ntypes;

/**
 * A force class for use in Vertex-based simulations. This force is based on the
 * Energy function proposed by Matteo et al in  Curr. Biol., 2007, 17, 2095-2104.
 */


template<unsigned DIM>
class MatteoForce : public FarhadifarForce<DIM>
{
friend class TestForces;

private:

    friend class boost::serialization::access;
    /**
     * Boost Serialization method for archiving/checkpointing.
     * Archives the object and its member variables.
     *
     * @param archive  The boost archive.
     * @param version  The current version of this class.
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<FarhadifarForce<DIM> >(*this);
    }
  
public:

    /**
     * Constructor.
     */
    MatteoForce();
  
    /**
     * Destructor.
     */
    virtual ~MatteoForce();

    /**
     * Get the line tension parameter for the edge between two given nodes.
     *
     * @param pNodeA one node
     * @param pNodeB the other node
     * @param rVertexCellPopulation reference to the cell population
     *
     * @return the line tension parameter for this edge.
     */
    virtual double GetLineTensionParameter(unsigned elem_index, Node<DIM>* pNodeA, Node<DIM>* pNodeB, VertexBasedCellPopulation<DIM>& rVertexCellPopulation);

    void OutputForceParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(MatteoForce)

#endif /*MatteoFORCE_HPP_*/
