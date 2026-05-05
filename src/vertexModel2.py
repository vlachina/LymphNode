from tyssue import PlanarGeometry, Sheet, History
from tyssue import config
from tyssue.behaviors import EventManager
from tyssue.behaviors.sheet import basic_events
from tyssue.draw import sheet_view
from tyssue.dynamics import effectors, model_factory
from tyssue.solvers.viscous import EulerSolver
import copy
import tyssue.draw.plt_draw as draw

from IPython.display import Image
from tyssue.draw.plt_draw import plot_forces

import src.brownianMotion as brownianMotion
import src.inputMechanicalParametersModel2 as inputMechanicalParametersModel2

from tyssue.geometry.base_geometry import update_periodic_dcoords  # Use the full path

import logging
log = logging.getLogger(__name__)



def initialize(numCellRows = 40):
    ## Defining energy contributions
    # https://tyssue.readthedocs.io/en/latest/_modules/tyssue/dynamics/effectors.html
    energyContributions_model = model_factory([
        #brownianMotion.BrownianMotion,
        effectors.FaceAreaElasticity,
        effectors.LineTension,
        effectors.LengthElasticity,
        #effectors.PerimeterElasticity,
        #effectors.CellAreaElasticity,
        #effectors.FaceContractility,
        #effectors.BarrierElasticity
        #effectors.LineViscosity
        #effectors.BorderElasticity
        ])

    ## Size of the patch
    #numCellRows = 40
    noiseCellShape = 0.2

    # noise = 0 -> hexagonal pattern
    # noise = 1 -> random voronoi
    cellMap = Sheet.planar_sheet_2d('tissue',
        nx=numCellRows, # approximate number of cells on the x axis
        ny=numCellRows, # approximate number of cells along the y axis
        distx=1, # distance between 2 cells along x
        disty=1, # distance between 2 cells along y
        noise=noiseCellShape)

    cellMap.remove(cellMap.cut_out([[1, numCellRows], [1, numCellRows]]), trim_borders=True)
    cellMap.reset_index()
    cellMap.reset_topo()


    ## Definition of the geometry of the sheet
    # PlanarGeometry: Geometry methods for 2D planar cell arangements
    # SheetGeometry: Geometry definitions for 2D sheets in 3D
    # BulkGeometry: Geometry functions for 3D cell arangements
    geom  = PlanarGeometry

    # Update geometry with the patch
    geom.update_all(cellMap)

    # Visualize the sheet
    #fig, ax = sheet_view(cellMap, mode="quick", figsize=(10, 10))

    ## Connect cells with energy contributions
    cellMap.update_specs(energyContributions_model.specs)

    return [cellMap, geom, energyContributions_model]


def on_topo_change(sheet):
    print('Topology changed!')


def solveEuler(cellMap, geom, energyContributions_model, endTime):

    ## Init history object
    # The History object records all the time steps
    history_cellMap = History(cellMap)

    ## Manager Initialization
    manager = EventManager("face", )
    #manager.append(basic_events.auto_reconnect)

    ## Init solver
    solver1 = EulerSolver(cellMap,
        geom,
        energyContributions_model,
        manager=manager,
        bounds=(
            -cellMap.edge_df.length.median()/10,
            cellMap.edge_df.length.median()/10
            ),
        history=history_cellMap,
        auto_reconnect=False)

    manager.update()

    ## Run the solver
    res1 = solver1.solve(tf=endTime, dt=1, on_topo_change=on_topo_change, topo_change_args=(solver1.eptm,))

    ## Deep copy to return it and being able to modify maintaining the previous one
    cellMap_new = copy.deepcopy(cellMap)
    history_new = copy.deepcopy(history_cellMap)

    return [cellMap_new, geom, energyContributions_model, history_new, solver1]


def solveEuler_with_energy_tracking(cellMap, geom, energyContributions_model, endTime):
    history_cellMap = History(cellMap)
    manager = EventManager("face")

    ## Init solver

    solver1 = EulerSolver(cellMap,
        geom,
        energyContributions_model,
        manager=manager,
        bounds=(
            -cellMap.edge_df.length.median() / 10,
            cellMap.edge_df.length.median() / 10
            ),
        history=history_cellMap,
        auto_reconnect=False)

    manager.update()
    energy_over_time = []  # List to store energy values

    for t in range(endTime):
        res1 = solver1.solve(tf=endTime, dt=1, on_topo_change=on_topo_change, topo_change_args=(solver1.eptm,))
        current_energy = energyContributions_model.compute_energy(cellMap) # Compute the current total energy
        energy_over_time.append(current_energy)  # Store the energy value

        # Your existing code to update the cellMap and record history, if necessary

    cellMap_new = copy.deepcopy(cellMap)
    history_new = copy.deepcopy(history_cellMap)

    return [cellMap_new, geom, energyContributions_model, history_new, energy_over_time, solver1]  # Return the list of energy values

from scipy.optimize import minimize

def energy_plotting(cellmap, energyContributions_model, geom):
    gradients = energyContributions_model.compute_gradient(cellmap, components=True)
    gradients = {label: (srce, trgt) for label, (srce, trgt)
                 in zip(energyContributions_model.labels, gradients)}
    fig, ax = plot_forces(cellmap, geom, energyContributions_model, ['z', 'y'], scaling=0.1)
    fig.set_size_inches(10, 12)



