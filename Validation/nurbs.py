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


ELEMENT_TYPE = "Tetrahedron"
ELEMENT_APPROXIMATION_DEGREE = 1
MATERIAL_MODEL = "NeoHookean"
FORCES = [0, -4000, 0]
# TODO improve the manual permutation for matching the redefinition of the hexahedron

nodes = np.array([[0., 0.], [1., 0.], [2., 0.], [0., 1.], [1., 1.],
                  [2., 1.], [0., 2.], [1., 2.], [2., 2.]])
nodes = np.array(nodes, dtype=np.float64)
weights = np.array([1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0])
weights =  np.array(weights, dtype=np.float64)

i_indices = np.array([0,1,2,3,4,5,6,7,8])
i_indices = np.array(i_indices, dtype=np.uint64)

knot_ranges = np.array([0.0, 0.0, 1.0, 1.0])
knot_ranges = np.array(knot_ranges, dtype=np.float64)

knot1 = np.array([0., 0., 0., 1., 1., 1.])
knot1 = np.array(knot1, dtype=np.float64)

knot2 = np.array([0., 0., 0., 1., 1., 1.])
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

FORCES = [0, 5, 0]

material = "PlaneStress"
#if MATERIAL_MODEL == "SaintVenantKirchhoff" or MATERIAL_MODEL == "NeoHookean":
#    material = MATERIAL_MODEL
#else:
#    raise ValueError('The material model is not implemented yet.')


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

        sofa_node.addObject('FixedConstraint', indices=[0, 1, 2])
        sofa_node.addObject('ConstantForceField', totalForce=FORCES, indices=[6,7,8])
        sofa_node.addObject(material + "Material", young_modulus="3000", poisson_ratio="0.3")
        sofa_node.addObject('ElasticSplineForcefield', template = element_type, printLog=True)

        self.sofa_rest_position = np.array(self.sofa_mo.position.value.copy().tolist())
        print("Nurbs Initial Positions " sofa_rest_position)

#        fenics_node = root.addChild("fenics_node")
#        fenics_node.addObject('StaticSolver', newton_iterations="25", relative_correction_tolerance_threshold="1e-15",
#                              relative_residual_tolerance_threshold="1e-10", printLog="1")
#        fenics_node.addObject('SparseLDLSolver', template="CompressedRowSparseMatrixMat3x3d")
#        self.fenics_mo = fenics_node.addObject('MechanicalObject', name="mo", position=mesh.points.tolist())
#        fenics_node.addObject('CaribouTopology', name='topology', template=element_fenics,
#                              indices=indices_fenics.tolist())
#        fenics_node.addObject('BoxROI', name="fixed_roi", box="-7.5 -7.5 -0.9 7.5 7.5 0.1")
#        fenics_node.addObject('FixedConstraint', indices="@fixed_roi.indices")
#        fenics_node.addObject('BoxROI', name="top_roi", box="-7.5 -7.5 79.9 7.5 7.5 80.1")
#        fenics_node.addObject('ConstantForceField', totalForce=FORCES, indices="@top_roi.indices")
#        fenics_node.addObject('FEniCS_Material', template=element_fenics, young_modulus="3000",
#                              poisson_ratio="0.3", C01=0.7, C10=-0.55, k=0.001, material_name=material, path="/home/..")

#        fenics_node.addObject('HyperelasticForcefield_FEniCS', printLog=True)

        return root

    def onSimulationInitDoneEvent(self, event):
        self.sofa_rest_position = np.array(self.sofa_mo.position.value.copy().tolist())
#        self.fenics_rest_position = np.array(self.fenics_mo.position.value.copy().tolist())

    def onAnimateBeginEvent(self, event):

        pass

    def onAnimateEndEvent(self, event):
        sofa_current_positions = np.array(self.sofa_mo.position.value.copy().tolist())

#        errors = []
#        for sofa_current_point, sofa_initial_point in zip(sofa_current_positions, self.sofa_rest_position):
#            if np.linalg.norm(sofa_current_point - sofa_initial_point) != 0:
#                errors.append(np.linalg.norm(sofa_current_point - fenics_current_point) / np.linalg.norm(
#                    sofa_current_point - sofa_initial_point))

        mean_error = np.mean(np.array(errors))

        displacement = sofa_current_positions - self.sofa_rest_position
#        print(displacement[:, 1].min())
        print("Displacement is ")
        print(displacement)

        print(f"Relative Mean Error: {100 * mean_error} %")


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
