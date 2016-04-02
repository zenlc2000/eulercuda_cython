# eulertour.pyx
#
# import low-level C declerations

cimport eulertour


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
                
                return eulertour.findEulerDevice(
                <EulerVertex *> &d_ev[0],
                <unsigned int *> &d_l[0], 
                <unsigned int *> &d_e[0], 
                vcount,
                <EulerEdge *> &d_ee[0],
                ecount,
                <CircuitEdge **> &d_cg_edge[0], 
                <unsigned int *> &cg_edgeCount[0],
                <unsigned int *> &cg_vertexCount[0],
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
                
                return eulertour.executeSwipeDevice(
                <EulerVertex *> &d_ev[0],
                <unsigned int *> &d_e[0], 
                vcount, 
                <EulerEdge *> &d_ee[0], 
                ecount, 
                <CircuitEdge * > &d_cg_edge[0],
                cg_edgeCount , 
                <unsigned int *> &d_tree[0],
                treeCount)
                

def markContigStart(EulerEdge[:] d_ee, 
                unsigned char[:] d_contigStart, 
                unsigned int ecount):
                
                return eulertour.markContigStart(<EulerEdge *> &d_ee[0], <unsigned char *> &d_contigStart[0],ecount)
                 
                 
                 

