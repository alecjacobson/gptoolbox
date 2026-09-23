import numpy as np
import igl

S = np.array(S).reshape(-1)
GV = np.array(GV);
GI = np.array(GI);
mc = igl.marching_cubes(S,GV,GI,isovalue)

class MarchingCubesSparseResult:
    def __init__(self,V,F):
        self.V = V
        self.F = F
 
res = MarchingCubesSparseResult(mc[0], mc[1])

