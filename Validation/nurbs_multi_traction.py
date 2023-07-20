import sys, os
from pathlib import Path

current_dir = Path(__file__).parent
site_packages_dir = (current_dir / '..' / 'build' / 'lib' / 'python3' / 'site-packages').resolve()
sys.path.insert(0, str(site_packages_dir))
print(f'Adding {site_packages_dir} to sys.path')

import Sofa
import SofaCaribou
import numpy as np
#import Caribou.Geometry
import Caribou
from Caribou.Topology import SplinePatch


# TODO NURBS Geometry data

nodes = np.array([[0., 0.], [1., 0.], [2., 0.], [3., 0.],
                  [0., 1.], [1., 1.], [2., 1.], [3., 1.],
                  [0., 2.], [1., 2.], [2., 2.], [3., 2.],
                  [0., 3.], [1., 3.], [2., 3.], [3., 3.]])
nodes = np.array(nodes, dtype=np.float64)
weights = np.array([1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0])
weights =  np.array(weights, dtype=np.float64)

i_indices = np.array([[0, 1, 2, 4, 5, 6, 8, 9, 10],
                      [1, 2, 3, 5, 6, 7, 9, 10, 11],
                      [4, 5, 6, 8, 9, 10, 12, 13, 14],
                      [5, 6, 7, 9, 10, 11, 13,14, 15]])
i_indices = np.array(i_indices, dtype=np.uint64)

knot_ranges = np.array([[0.0, 0.0, 0.5, 0.5],
                        [0.5, 0.0, 1.0, 0.5],
                        [0.0, 0.5, 0.5, 1.0],
                        [0.5, 0.5, 1.0, 1.0]])
knot_ranges = np.array(knot_ranges, dtype=np.float64)

knot1 = np.array([0., 0., 0., 0.5, 1., 1., 1.])
knot1 = np.array(knot1, dtype=np.float64)

knot2 = np.array([0., 0., 0., 0.5, 1., 1., 1.])
knot2 = np.array(knot2, dtype=np.float64)
patch = SplinePatch(nodes, weights, i_indices, knot1, knot2, knot_ranges)
#patch = SplinePatch(nodes, weights, i_indices, knot1, knot2, knotranges)

#print("I am in Canti Python")
#print(patch.knot_1().tolist())

#print("positions")
#print(patch.all_positions())

#print("weights")
#print(patch.all_weights())

#print("indices")
#print(patch.indices())

element_type = "NurbsSurf_2D"

FORCES = [0, 1000, 0]

material = "PlaneStress"


class ControlFrame(Sofa.Core.Controller):

    def __init__(self, node):
        Sofa.Core.Controller.__init__(self)
        self.root = self.CreateGraph(node)

    def CreateGraph(self, root):

        root.addObject('DefaultVisualManagerLoop')
        root.addObject('DefaultAnimationLoop')
        root.addObject('VisualStyle', displayFlags="showForceFields showBehaviorModels")
        root.addObject('RequiredPlugin',
                       pluginName="SofaOpenglVisual SofaBaseMechanics SofaBaseTopology SofaSparseSolver SofaImplicitOdeSolver SofaTopologyMapping SofaBoundaryCondition SofaEngine")
        root.gravity = [0, 0, 0]
        sofa_node = root.addChild("sofa_node")
        sofa_node.addObject('StaticSolver', newton_iterations="25", relative_correction_tolerance_threshold="1e-15",
                            relative_residual_tolerance_threshold="1e-10", printLog="1")
        sofa_node.addObject('SparseLDLSolver', template="CompressedRowSparseMatrixMat3x3d")
        self.sofa_mo = sofa_node.addObject('MechanicalObject', name="mo", template = "Vec2d", position=patch.all_positions().tolist())
        sofa_node.addObject('CaribouSplineTopology', name='topo', template = element_type,
                             indices=patch.indices().tolist(), knot_1 = patch.knot_1().tolist(),
                             knot_2 = patch.knot_2().tolist(), weights = patch.all_weights().tolist(), knot_spans = knot_ranges.tolist())
#        sofa_node.addObject('CaribouMass', density = 0.1)
        sofa_node.addObject('FixedConstraint', indices=[0, 1, 2, 3])
        sofa_node.addObject('TractionSplineForcefield', traction=FORCES, boundary = 3)
        sofa_node.addObject(material + "Material", young_modulus="10000", poisson_ratio="0.3")
        sofa_node.addObject('ElasticSplineForcefield', template = element_type, printLog=True)

        self.sofa_rest_position = np.array(self.sofa_mo.position.value.copy().tolist())
        print("Nurbs Initial Positions ", self.sofa_rest_position)

        return root

    def onSimulationInitDoneEvent(self, event):
        self.sofa_rest_position = np.array(self.sofa_mo.position.value.copy().tolist())

    def onAnimateBeginEvent(self, event):

        pass

    def onAnimateEndEvent(self, event):
        sofa_current_positions = np.array(self.sofa_mo.position.value.copy().tolist())
        displacement = sofa_current_positions - self.sofa_rest_position
        print("Displacement is ")
        print(displacement)


def createScene(node):
    node.addObject(ControlFrame(node))


# Choose in your script to activate or not the GUI
USE_GUI = True


def main():
    import SofaRuntime
    import Sofa.Gui
    SofaRuntime.importPlugin("SofaOpenglVisual")
    SofaRuntime.importPlugin("SofaImplicitOdeSolver")
    SofaRuntime.importPlugin("SofaLoader")

    root = Sofa.Core.Node("root")
    createScene(root)
    Sofa.Simulation.init(root)

    if not USE_GUI:
        for iteration in range(10):
            Sofa.Simulation.animate(root, root.dt.value)
    else:
        Sofa.Gui.GUIManager.Init("myscene", "qglviewer")
        Sofa.Gui.GUIManager.createGUI(root, __file__)
        Sofa.Gui.GUIManager.SetDimension(1080, 1080)
        Sofa.Gui.GUIManager.MainLoop(root)
        Sofa.Gui.GUIManager.closeGUI()


# Function used only if this script is called from a python environment
if __name__ == '__main__':
    main()
