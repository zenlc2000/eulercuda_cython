# eulertour.pyx
#
# import low-level C declerations

cimport pyeulertour


def findEulerDevice(EulerVertex[:] d_ev,
                unsigned int[:] d_l, 
                unsigned int[:] d_e, 
                unsigned int vcount,
                EulerEdge[:] d_ee,
                unsigned int ecount,
                CircuitEdge[:] d_cg_edge, 
                unsigned int[:] cg_edgeCount,
                unsigned int[:] cg_vertexCount,
                unsigned int kmerLength):
                
    return pyeulertour.findEulerDevice(
            d_ev,
            d_l,
            d_e,
            vcount,
            d_ee,
            ecount,
            d_cg_edge,
            cg_edgeCount,
            cg_vertexCount,
            kmerLength)
                
                
def executeSwipeDevice(EulerVertex[:] d_ev,
                unsigned int[:] d_e, 
                unsigned int vcount, 
                EulerEdge[:] d_ee, 
                unsigned int ecount, 
                CircuitEdge[:] d_cg_edge,
                unsigned int cg_edgeCount , 
                unsigned int[:] d_tree,
                unsigned int treeCount):
                
    return pyeulertour.executeSwipeDevice(
            d_ev,
            d_e,
            vcount,
            d_ee,
            ecount,
            d_cg_edge,
            cg_edgeCount ,
            d_tree,
            treeCount)
                

def markContigStart(EulerEdge[:] d_ee, 
                unsigned char[:] d_contigStart, 
                unsigned int ecount):
                
    return pyeulertour.markContigStart(
        d_ee,
        d_contigStart,
        ecount)
                 
                 
                 

