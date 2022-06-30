"""Landlab component that calculates temporal and spatial changes in river
bed elevation and grain size distribution. Also, this component predicts 
fractional or total sediment transport based on the specified bed load 
transport model. Hydraulic properties are obtained from an external flow 
component. We recommend coupling with the overland flow component from 
Adams et al, 2017.

.. codeauthor:: Angel Monsalve, Sam Anderson

Examples
--------

>>> import numpy as np
>>> import copy
>>> from landlab import RasterModelGrid
>>> from landlab.components import RiverBedDynamics
>>> from landlab import imshow_grid
>>> from landlab.grid.mappers import map_mean_of_link_nodes_to_link

Create a grid on which to calculate sediment transport

>>> grid = RasterModelGrid((5, 5))

The grid will need some data to run the bedload transport component. 
To check the names of the fields that provide input to the bedload transport 
component, use the *input_var_names* class property.

>>> RiverBedDynamics.input_var_names
('bed_surface__grainSizeDistribution_location',
 'surface_water__depth',
 'surface_water__velocity',
 'surface_water__velocity_previous_time_step',
 'topographic__elevation',
 'topographic__elevation_original')

Create fields of data for each of these input variables. When running a 
complete simulation some of these variables will be created by the flow model.
Notice that surface water depth and velocity are required at links. However, 
specifying these variables at nodes is easier and then we can map the fields 
onto links. By doing so, we don't have to deal with links numbering. When this 
component is coupled to overland flow there is no need to map fields because it 
is done automatically within the component.

We create the bed surface grain size (GSD) distribution location. We assume 
that there are two different GSD within the watershed (labeled as 0 and 1)
>>> grid.at_node['bed_surface__grainSizeDistribution_location'] = np.array([
... 0, 1., 1., 1., 0, 
... 0, 1., 1., 1., 0, 
... 0, 1., 1., 1., 0, 
... 0, 1., 1., 1., 0, 
... 0, 1., 1., 1., 0,])

Now we create the topography data

>>> grid.at_node['topographic__elevation'] = np.array([
... 1.09, 1.08, 1.08, 1.08, 1.09,
... 1.08, 1.07, 1.07, 1.07, 1.08,
... 1.07, 1.06, 1.06, 1.06, 1.07,
... 1.06, 1.05, 1.05, 1.05, 1.06,
... 1.05, 1.04, 1.03, 1.04, 1.05,])

We copy this elevation to define the original bed surface elevation
>>> grid.at_node['topographic__elevation_original'] = \
    copy.deepcopy(grid.at_node['topographic__elevation']) 

We set the boundary conditions
>>> grid.set_watershed_boundary_condition(grid.at_node['topographic__elevation'])

And check the node status
>>> grid.status_at_node
array([4, 4, 4, 4, 4, 4, 0, 0, 0, 4, 4, 0, 0, 0, 4, 4, 0, 0, 0, 4, 4, 4, 1,
       4, 4], dtype=uint8)
Which tell us that there is one outlet located on the 22nd node

the topography data can be display using
>>> imshow_grid(grid,'topographic__elevation')

Now we add some water into the watershed. In this case is specified in nodes

>>> grid.at_node['surface_water__depth'] = np.array([
... 0.102, 0.102, 0.102, 0.102, 0.102,
... 0.102, 0.102, 0.102, 0.102, 0.102,
... 0.102, 0.102, 0.102, 0.102, 0.102,
... 0.102, 0.102, 0.102, 0.102, 0.102,
... 0.102, 0.102, 0.102, 0.102, 0.102,])

Now we give the water a velocity. 
>>> grid.at_node['surface_water__velocity'] = np.array([
... 0.25, 0.25, 0.25, 0.25, 0.25,
... 0.25, 0.25, 0.25, 0.25, 0.25,
... 0.25, 0.25, 0.25, 0.25, 0.25,
... 0.25, 0.25, 0.25, 0.25, 0.25,
... 0.25, 0.25, 0.25, 0.25, 0.25,])

Notice that we are trying to specify a vector at a node using a single value. 
In this example this was made on purpose to highlight this fact. The component 
will assume that this is the magnitude and, in this case, it will have 
different components depending the location in the grid. When overland flow is
used, there is no need to specify a velocity because it is a result of the
component.

Now (for the purpose of the example) we assume that is identical to the 
previous time.
>>> grid.at_node['surface_water__velocity_previous_time_step'] = \
    grid.at_node['surface_water__velocity']

By default, when creating our grid we used a spacing of 1 m in the x and y 
directions. Therefore, the discharge is 2.5 m**3/s.

Now we map nodes into links when it is required
>>> grid['link']['surface_water__depth'] = \
    map_mean_of_link_nodes_to_link(grid,'surface_water__depth')
>>> grid['link']['surface_water__velocity'] = \
    map_mean_of_link_nodes_to_link(grid,'surface_water__velocity')
>>> grid['link']['surface_water__velocity_previous_time_step'] = \
    map_mean_of_link_nodes_to_link(grid,'surface_water__velocity_previous_time_step')

We assign a GSD to each location. The structure of this array is:
First column contains the grain sizes in milimiters
Second column is location 0 in 'bed_grainSizeDistribution__location'
Third column is location 1 in 'bed_grainSizeDistribution__location', and so on
>>> gsd = np.array([[32, 100, 100],[16, 25, 50],[8, 0, 0]])

A time step must be provided. Usually this will be an output of overland flow,
but we can override it by specifying a custom value.
>>> timeStep = 1 # time step in seconds

Instantiate the `RiverBedDynamics` component to work on the grid, and run it.

>>> RBD = RiverBedDynamics(grid , gsd = gsd, dt = timeStep,bedloadEq = 'Parker1990')
>>> RBD.run_one_step()

After executing the bedloadTransport componet, new fields have been added to 
the grid. Use the *output_var_names* property to see the names of the fields 
that have been changed.

>>> RBD.output_var_names
('bed_surface__geometricMeanSize',
 'bed_surface__grainSizeDistribution',
 'bed_surface__medianSize',
 'bed_surface__sandFraction',
 'bed_surface__standardDevSize',
 'sediment_transport__bedloadGSD',
 'sediment_transport__bedloadRate',
 'sediment_transport__netBedload',
 'surface_water__shearStress')

The `sediment_transport__bedloadRate` field is defined at nodes.
>>> RBD.var_loc('sediment_transport__netBedload')
'node'

>>> grid.at_node['sediment_transport__netBedload'] # doctest: +NORMALIZE_WHITESPACE
array([  0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
         0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
        -4.50423253e-08,   0.00000000e+00,  -4.50423253e-08,
         0.00000000e+00,   0.00000000e+00,  -4.50423253e-08,
         0.00000000e+00,  -4.50423253e-08,   0.00000000e+00,
         0.00000000e+00,  -4.50423253e-08,   1.14595731e-04,
        -4.50423253e-08,   0.00000000e+00,   0.00000000e+00,
         0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
         0.00000000e+00])

The 'surface_water__shearStress' field is defined at links.

>>> RBD.var_loc('surface_water__shearStress')
'link'

>>> grid.at_link['surface_water__shearStress'] # doctest: +NORMALIZE_WHITESPACE
array([ 10.002783,   0.      ,   0.      ,  10.002783,  10.002783,
        10.002783,  10.002783,  10.002783,  10.002783,  10.002783,
         0.      ,   0.      ,  10.002783,  10.002783,  10.002783,
        10.002783,  10.002783,  10.002783,  10.002783,   0.      ,
         0.      ,  10.002783,  10.002783,  10.002783,  10.002783,
        10.002783,  10.002783,  10.002783,   0.      ,   0.      ,
        10.002783,  10.002783,  10.002783,  20.005566,  10.002783,
        10.002783,  10.002783,  10.002783,  10.002783,  10.002783])

Considering the link at the watershed exit, link Id 33, we can obtain the bed
load transport rate

>>> grid.at_link['sediment_transport__bedloadRate'][33]
0.00011474671940125335

Therefore, the bed load transport rate according to Parker 1991 surface-based
equation is 1.147 * 10^-4 m2/s

The GSD at this place is
>>> grid.at_link['sediment_transport__bedloadGSD'][33]
array([ 0.43305913,  0.56694087])

Which in cummulative percentage is equivalent to
D mm    % Finer
32      100.00
16      56.70
8       0.00

We can also check the  bed load grain size distribution in all nodes

>>> grid.at_node['sediment_transport__bedloadGSD']
array([[ 0.41234914,  0.58765086],
       [ 0.35174623,  0.64825377],
       [ 0.15411137,  0.34588863],
       [ 0.35174623,  0.64825377],
       [ 0.41234914,  0.58765086],
       [ 0.46791178,  0.53208822],
       [ 0.33876249,  0.66123751],
       [ 0.14112763,  0.35887237],
       [ 0.33876249,  0.66123751],
       [ 0.46791178,  0.53208822],
       [ 0.46791178,  0.53208822],
       [ 0.33876249,  0.66123751],
       [ 0.14112763,  0.35887237],
       [ 0.33876249,  0.66123751],
       [ 0.46791178,  0.53208822],
       [ 0.46791178,  0.53208822],
       [ 0.33876249,  0.66123751],
       [ 0.21643048,  0.28356952],
       [ 0.33876249,  0.66123751],
       [ 0.46791178,  0.53208822],
       [ 0.41234914,  0.58765086],
       [ 0.30822274,  0.69177726],
       [ 0.35764978,  0.64235022],
       [ 0.30822274,  0.69177726],
       [ 0.41234914,  0.58765086]])


After the flow acted on the bed and sediment transport occured we can check 
the new topographic elevation field

>>> grid.at_node['topographic__elevation']  # doctest: +NORMALIZE_WHITESPACE
array([ 1.09      ,  1.08      ,  1.08      ,  1.08      ,  1.09      ,
        1.08      ,  1.07000007,  1.07      ,  1.07000007,  1.08      ,
        1.07      ,  1.06000007,  1.06      ,  1.06000007,  1.07      ,
        1.06      ,  1.05000007,  1.04983077,  1.05000007,  1.06      ,
        1.05      ,  1.04      ,  1.02983077,  1.04      ,  1.05      ])


"""
from landlab import Component, FieldError
import numpy as np
import pandas as pd
import scipy.constants
from scipy.interpolate import interp1d
import copy
import os
import shutil

class RiverBedDynamics(Component):

    """Landlab component that predicts the evolution of a river bed 
    considering changes in elevation and grain size distribution in response to
    bed load transport according to the Exner equation and the transfer 
    function of Toro-Ecobar et al., (1996).
    
    To estimate temporal and spatial changes in river bed properties, this 
    component predicts the bedload transport rate and fractional tranport at 
    each link using the unsteady shear stress. Time varying hydraulics 
    variables are obtained from a surface flow flow, for instance OverlandFlow.

    The primary method of this class is :func:`run_one_step`.

    References
    ----------
    **Required Software Citation(s) Specific to this Component**

    Not required but recommended
    
    Adams, J., Gasparini, N., Hobley, D., Tucker, G., Hutton, E., Nudurupati,
    S., Istanbulluoglu, E. (2017). The Landlab v1. 0 OverlandFlow component:
    a Python tool for computing shallow-water flow across watersheds.
    Geoscientific Model Development  10(4), 1645.
    https://dx.doi.org/10.5194/gmd-10-1645-2017

    **Additional References**

    G. Parker (1990) Surface-based bedload transport relation for gravel 
    rivers, Journal of Hydraulic Research, 28:4, 417-436, 
    DOI: 10.1080/00221689009499058
    
    Wilcock, P. R., & Crowe, J. C. (2003). Surface-based transport model for 
    mixed-size sediment. Journal of hydraulic engineering, 129(2), 120-128.
    DOI: 10.1061/(ASCE)0733-9429(2003)129:2(120)
    
    Meyer-Peter, E. and Müller, R., 1948, Formulas for Bed-Load Transport, 
    Proceedings, 2nd Congress, International Association of Hydraulic Research, 
    Stockholm: 39-64.
    
    Fernandez Luque, R. and R. van Beek, 1976, Erosion and transport of 
    bedload sediment, Journal of Hydraulic Research, 14(2): 127-144.
    https://doi.org/10.1080/00221687609499677
    
    Mueller, E. R., J. Pitlick, and J. M. Nelson (2005), Variation in the 
    reference Shields stress for bed load transport in gravelbed streams and 
    rivers, Water Resour. Res., 41, W04006, doi:10.1029/2004WR003692
    
    Carlos M. Toro-Escobar, Chris Paola & Gary Parker (1996) Transfer function 
    for the deposition of poorly sorted gravel in response to streambed 
    aggradation, Journal of Hydraulic Research, 34:1, 35-53, 
    DOI: 10.1080/00221689609498763

    """

    _name = "RiverBedDynamics"

    _unit_agnostic = True

    _info = {
        'bed_surface__fixedElevation': {
          	"dtype": int,
          	"intent": "in",
          	"optional": True,
          	"units": "-",
          	"mapping": 'node',
          	"doc": "Sets a node as a fixed elevation",
        },
        'bed_surface__grainSizeDistribution_location': {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "-",
            "mapping": 'node',
            "doc": "Accounts for the initial grain size distribution that can \
                vary in space. Each value corresponds to a given GSD location.",
        },
        'bed_surface__fixed_grainSizeDistribution': {
            "dtype": int,
            "intent": "in",
            "optional": True,
            "units": "-",
            "mapping": 'node',
            "doc": "Determine locations that do not change its GSD in time",
        },            
        'surface_water__depth': {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "m",
            "mapping": 'link',
            "doc": "Depth of water on the surface at nodes. Usually it is loaded \
                from the flow model",
        },
        "surface_water__velocity": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "m/s",
            "mapping": 'link',
            "doc": "Average velocity of the surface water. Usually it is loaded \
                from the flow model",
        },
        "surface_water__velocity_previous_time_step": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "m/s",
            "mapping": 'link',
            "doc": "Average velocity of the surface water in the previous time \
                step. Usually it is loaded from the flow model",
        },                     
        'topographic__elevation': {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "m",
            "mapping": 'node',
            "doc": "Land surface topographic elevation. Usually loaded from \
                flow model. ",
        },		
        "sediment_transport__bedloadGSD": {
            	"dtype": float,
            	"intent": "out",
            	"optional": False,
            	"units": "-",
            	"mapping": 'link',
            	"doc": "Bed load grain size distribution",
        },    
        'sediment_transport__bedloadRate': {
            	"dtype": float,
            	"intent": "out",
            	"optional": False,
            	"units": "m^2/s",
            	"mapping": 'link',
            	"doc": "Volumetric bed load transport rate per unit width",
        },        
        "bed_surface__geometricMeanSize": {
            	"dtype": float,
            	"intent": "out",
            	"optional": False,
            	"units": "mm",
            	"mapping": 'link',
            	"doc": "Bed surface geometric mean grain size",
        },
        "bed_surface__grainSizeDistribution": {
            	"dtype": float,
            	"intent": "out",
            	"optional": False,
            	"units": " ",
            	"mapping": 'link',
            	"doc": "Bed surface grain size fraction",
        },     
        "bed_surface__medianSize": {
            	"dtype": float,
            	"intent": "out",
            	"optional": False,
            	"units": "mm",
            	"mapping": 'node',
            	"doc": "Bed surface median grain size",
        },        
        "bed_surface__sandFraction": {
            	"dtype": float,
            	"intent": "out",
            	"optional": False,
            	"units": "-",
            	"mapping": 'node',
            	"doc": "Bed surface sand content",
        },         
        "bed_surface__standardDevSize": {
            	"dtype": float,
            	"intent": "out",
            	"optional": False,
            	"units": "mm",
            	"mapping": 'node',
            	"doc": "Bed surface standard deviation grain size",
        },
        'sediment_transport__imposed_bedloadGSD': {
          	"dtype": float,
          	"intent": "in",
          	"optional": True,
          	"units": " ",
          	"mapping": 'link',
          	"doc": "Sets the sediment transport GSD where \
                  sediment supply is imposed",
        }, 
        'sediment_transport__imposed_sediment_supply': {
          	"dtype": float,
          	"intent": "in",
          	"optional": True,
          	"units": "m^2/s",
          	"mapping": 'link',
          	"doc": "Sets the sediment transport rate per unit width where \
                  sediment supply is imposed",
        },             
        'sediment_transport__netBedload': {
          	"dtype": float,
          	"intent": "out",
          	"optional": False,
          	"units": "m^3/s",
          	"mapping": 'node',
          	"doc": "Net sediment transport rate at a node",
        },                           
        "surface_water__shearStress": {
            	"dtype": float,
            	"intent": "out",
            	"optional": False,
            	"units": " Pa ",
            	"mapping": 'node',
            	"doc": "Shear stress applied by the surface water on the bed surface",
        },
 
    }

    def __init__(
        self,
        grid,
        gsd,                        # Initialize the gsd array.
        rho = 1000,                 # Fluid density (kg/m**3)
        rho_s = 2650,               # Sediment density (kg/m**3)
        bedloadEq = 'MPM',          # Used to select the equation
        variableCriticalShearStress = False, # Critical shear stress varies 
                                    # with bed slope
        useHydraulicsRadiusInShearStress = False, # Used the hydraulic radius
                                    # instead of water depth to calculate
                                    # shear stresses
        lambda_p = 0.35,            # Sediment porosity
        outletBoundaryCondition = "zeroGradient",
        evolveBed = True,           # Bed surface elevation and GSD changes in 
                                    # response to bed load transport.
        updateBedSurfaceGSD = True, # Bed surface GSD is updated considering the 
                                    # the interaction with bed load rate and bed 
                                    # load GSD.
        trackStratigraphy = False,  # Stores in the disk the GSD of each sediment 
                                    # layer at each node
        nCyclesToProcessStratigraphy = 10,    # Enter to the read/write 
                                    # stratigraphy routine every this number 
                                    # of time steps.
        newSurfaceLayerThickness = 1, # When reaching this thickness the deposi-
                                    # ted surface layer becomes a new subsurface
                                    # layer 
        output_vector = False,      # Shear stress and velocity vector can be
                                    # exported.
        dt = 1,                     # Time step (s). When coupled to OverlandFlow
                                    # is dynamically adjusted
        alpha = 1.0,                # An upwinding coefficient for a central 
                                    # difference scheme used when updating the 
                                    # bed GSD
    ):
        
        """Calculate the evolution of a river bed based on bedload transport 
        and fractional rates on links. An external flow hydraulics solver is 
        required to predict flow variables in time and space, for example, 
        OverlandFlow. The shear stress used in sediment transport equations 
        accounts for variability in time and space of water depth and flow 
        velocity. 
        
        This component adjusts topographic elevation and grain size
        distribution at each node within a landlab grid.

        Parameters
        ----------
        grid :  RasterModelGrid
            A landlab grid. 
        gsd :   Array float, mandatory, 
            Grain size distribution. Must contain as many GSD as different 
            indexes are in GSDLocation 
        rho :   float, optional
            Density of fluid, default set to water density of 1000 kg/m^3
        rho_s : float, optional
            Density of sediment, default set to sediment density of 2650 kg/m^3         
        bedloadEq :    String, optional
            Select the bedload transport equation, use 
            MPM for Meyer Peter and Muller
            FLvB for Fernandez Luque & van Beek (1976)
            Parker1990 for Parker 1990 
            WilcockAndCrowe for Wilcock and Crowe 2003
            default is MPM 
        variableCriticalShearStress: Boolean, optional
            If True the critical shear stress in MPM or FLvB
            will be obtained using Mueller et al. (2005) equation.
        useHydraulicsRadiusInShearStress: Boolean, optional
            If True shear stress calculations will use the hydraulics ratio
            default is False, so it uses water depth
        lambda_p : float, optional
            Sediment porosity, default value is 0.35 [ ]
        outletBoundaryCondition : str
            Boundary condition at watershed outlet. 
            Default condition is zeroGradient (Maps outlet upstream node).
            Other options are:
            fixedValue: Does not change value at outlet through the run 
        evolveBed : Boolean, optional
            If True the bed evolves according to local bed load transport rate
            and GSD. If False, bed load transport rate and GSD are calculated
            and output.
        updateBedSurfaceGSD : Boolean, optional
            If True the bed GSD can evolve according to the interaction between
            the current surface GSD and the bed load GSD and rate.
        trackStratigraphy : Boolean, optional
            If True the component will store the GSDD of every layer at each node.
            This is computationally demand because needs to read/write data at
            every time step. Recommended only when cycles of erosion and deposition
            are frequent are very important.    
        nCyclesToProcessStratigraphy : int, optional
            If trackStratigraphy is True data will be read and stores every 
            "nCyclesToProcessStratigraphy" times steps. Must be larger or equal to 1
        dt: float, optional
            Time step in seconds. When this component is coupled to a flow model
            it is dynamically updated.
        alpha : float, optional
            It is an upwinding coefficient for a central difference scheme when 
            updating the bed GSD - default value is 1.0 - a value of 0.5 generates
            an central differences scheme.
        """
        super().__init__(grid)
        
        self._g = scipy.constants.g # Acceleration due to gravity (m/s**2).
        self._rho = rho
        self._rho_s = rho_s
        self._R = (rho_s-rho)/rho                                   
        self._tauStar_rsgo = 0.0386 # Reference dimensionless shear stress for 
                                    # the median size
        self._beta = 0.0951         # Coefficient for the hiding function
        self._gsd = gsd
        self._bedloadEq = bedloadEq
        self._variableCriticalShearStress = variableCriticalShearStress
        self._useHydraulicsRadiusInShearStress = useHydraulicsRadiusInShearStress
        self._lambda_p = lambda_p
        self._alpha = alpha
        self._outletBoundaryCondition = outletBoundaryCondition
        
        self._gsd_location = \
            self._grid.at_node['bed_surface__grainSizeDistribution_location']
            
        self._grid._dt = dt
        
        # True if it is the very first iteration
        self._firstIteration = True
            
        # Initialize the bed surface grain size properties. At the beginnig and
        # only for the first time step it uses the GSD from gsd.xlsx. Later it 
        # uses the calculated one
        self.defineInitialBedProperties()
        
        self.link_at_border_cells()
               
        # This flag changes to True if there are links with fixed bedload GSD
        self._fixedLinkCalculate = False
        
        # This flag is used to activate or deactivate the bed evolution part
        # of the component.
        self._evolveBed = evolveBed
        
         # This flag is used to activate or deactivate the bed GSD updating part
        # of the component.       
        self._updateBedSurfaceGSD = updateBedSurfaceGSD
        
        # This flag is used to activate or deactivate option to store the GSD
        # of individual layers in each node. 
        self._trackStratigraphy = trackStratigraphy
        
        # This value is used to merge all deposited layers in a new subsurface
        # layer
        self._newSurfaceLayerThickness = newSurfaceLayerThickness
        
        # When newSurfaceLayerThickness is reached this flag is used to force 
        # recording and reading data
        self._computeStratigraphy  = False
                
        # When newSurfaceLayerThickness is reached this flag is used to update  
        # subsurface data
        self._updateSubSurface = False
        self._updateErodedSubSurface = False
        self._subsurfaceChanged = False
        
        # If true exports shear stress and velocity as vectors on nodes at each 
        # time step
        self._output_vector = output_vector
                       
        self._nLinks = self._grid.number_of_links
        
        # Sets initial values to enter into the write and read stratigraphy routine
        self._stratigraphyCycle = 0
        self._nCyclesToProcessStratigraphy = nCyclesToProcessStratigraphy
                
        self._po = np.array([
            0.6684,0.7639,0.8601,0.9096,0.9615,1,1.055,1.108,1.197,1.302,1.407,
            1.529,1.641,1.702,1.832,1.937,2.044,2.261,2.499,2.732,2.993,3.477,
            4.075,4.469,5.016,6.158,7.821,10.06,14.38,19.97,25.79,38.57,68.74,
            91.95,231.2,2320
            ])
        self._oo = np.array([
            1.011,1.011,1.01,1.008,1.004,0.9997,0.9903,0.9789,0.9567,0.9273,
            0.8964,0.8604,0.8287,0.8123,0.7796,0.7554,0.7326,0.6928,0.6585,
            0.6345,0.615,0.5877,0.564,0.5523,0.5395,0.5209,0.5045,0.4917,
            0.479,0.4712,0.4668,0.462,0.4578,0.4564,0.4541,0.4527
            ])
        self._so = np.array([
            0.8157,0.8157,0.8182,0.8233,0.8333,0.8439,0.8621,0.8825,0.9214,
            0.9723,1.025,1.083,1.13,1.153,1.196,1.225,1.25,1.287,1.313,1.333,
            1.352,1.38,1.403,1.414,1.426,1.444,1.458,1.469,1.48,1.486,1.49,
            1.493,1.497,1.498,1.499,1.5
            ])
        
        # Creating grid fields at time zero
        try:
            self._grid['link']['sediment_transport__bedloadRate'] = grid.add_zeros(
                'sediment_transport__bedloadRate', at='link',
                units=self._info['sediment_transport__bedloadRate']["units"])
        except FieldError:
            self._grid['link']['sediment_transport__bedloadRate'] = \
                grid.at_link['sediment_transport__bedloadRate']
            self._grid['link']['sediment_transport__bedloadRate'].fill(0.0)  
        
        try:
            self._grid['node']['bed_surface__fixedElevation'] = grid.add_zeros(
                'bed_surface__fixedElevation', at='node',
                units=self._info['bed_surface__fixedElevation']["units"])
        except FieldError:
            self._grid['node']['bed_surface__fixedElevation'] = \
                grid.at_node['bed_surface__fixedElevation']
            self._grid['node']['bed_surface__fixedElevation'].fill(0.0)  
        
        try:
            self._grid['link']['sediment_transport__imposed_sediment_supply'] = grid.add_zeros(
                'sediment_transport__imposed_sediment_supply', at='link',
                units=self._info['sediment_transport__imposed_sediment_supply']["units"])
        except FieldError:
            self._grid['link']['sediment_transport__imposed_sediment_supply'] = \
                grid.at_link['sediment_transport__imposed_sediment_supply']
            self._grid['link']['sediment_transport__imposed_sediment_supply'].fill(0.0)
        
        try:
            self._grid['node']['bed_surface__fixed_grainSizeDistribution'] = grid.add_zeros(
                'bed_surface__fixed_grainSizeDistribution', at='node',
                units=self._info['bed_surface__fixed_grainSizeDistribution']["units"])
        except FieldError:
            self._grid['node']['bed_surface__fixed_grainSizeDistribution'] = \
                grid.at_node['bed_surface__fixed_grainSizeDistribution']
            self._grid['node']['bed_surface__fixed_grainSizeDistribution'].fill(0.0)        
        
        try:
            self._grid['link']['sediment_transport__imposed_bedloadGSD'] = grid.add_zeros(
                'sediment_transport__imposed_bedloadGSD', at='link',
                units=self._info['sediment_transport__imposed_bedloadGSD']["units"])
        except FieldError:
            self._grid['link']['sediment_transport__imposed_bedloadGSD'] = \
                grid.at_link['sediment_transport__imposed_bedloadGSD']
            self._grid['link']['sediment_transport__imposed_bedloadGSD'].fill(0.0)             
        
        # Define faces normal vector
        self._normal = -(self._grid.link_dirs_at_node)
        
        # Makes a copy of the original bed surface elevation and maps into links
        self._grid['node']['topographic__elevation_original'] = \
            copy.deepcopy(self._grid['node']['topographic__elevation'])
        self._grid['link']['topographic__elevation_original'] = \
            self._grid.map_mean_of_link_nodes_to_link(self._grid['node']['topographic__elevation_original'])
        self._grid['link']['topographic__elevation_subsurface'] = \
            self._grid.map_mean_of_link_nodes_to_link(self._grid['node']['topographic__elevation_original'])
        self._grid['link']['bed_surface__surface_thickness_new_layer'] = \
            self._grid['link']['topographic__elevation_original'] * 0
            
        # Defines some variables for storing stratigraphy
        # location x and y
        self._x = 0.5 * (self._grid.node_x[self._grid.node_at_link_head] + \
                          self._grid.node_x[self._grid.node_at_link_tail])
        self._y = 0.5 * (self._grid.node_y[self._grid.node_at_link_head] + \
                          self._grid.node_y[self._grid.node_at_link_tail])    
        # folders location
        self._cwd = os.getcwd() 
        self._StratigraphyTempFilesPath = 'stratigraphyTempFiles'
        self._StratigraphyRawDataPath = 'stratigraphyRawData'
        
    def defineInitialBedProperties(self):              
        """ Conducts the initial configuration of bed grain size distribution 
        properties. Reads input data and fills the required variables.
        
        This definition is conducted only during the first time step. Then, the
        calculated or updated bed information, which will be already in the 
        required format, will be used.
        
        """
        nNodes = self._grid.number_of_nodes
        
        D = self._gsd[:,0]      # Grain sizes
        f = self._gsd[:,1:]     # Grain sizes frequency cumulative
        
        # Number of locations with different GSD during the first time step
        nGSDLocations = self._gsd.shape[1]-1  
        
        # Calculates the sand fraction in each location
        if np.min(D)<2:     # Grain smaller than 2 mm are sand
            FS_0 = self.calculateSandFraction()
            # Flag is set as True to adjust the GSD and calculate a sand-free 
            # GSD when using Parker 1990 eq
            self._adjustGSD_flag = True               
        else: 
            FS_0 = np.zeros(self._gsd.shape[1]-1)
            # Flag is set as False, there is no need to adjust the GSD
            self._adjustGSD_flag = False              
   
        # If Parker Eq is used we remove sand content
        if self._adjustGSD_flag is True:    
            
            # We add 2mm into GSD and update the GSD
            id2mm = np.argmin(D>=2)     # Location where 2 mm will be placed
   
            if self._bedloadEq == 'Parker1990':     
                # Adds 2 mm and removes sand and smaller grains                     
                D = np.concatenate([D[0:id2mm],[2]])        
                f = np.concatenate([f[0:id2mm,:],[FS_0]]) 
                gravelFraction = (f[0,:] - f[-1,:])/100
                # Calculates cumulative grain sizes frequency
                f_0 = np.flip(-np.diff(f/100,axis=0)/gravelFraction,axis=0)
                f = np.abs(np.flip(np.cumsum( f_0 ,axis=0),axis = 0)*100) 
                f = np.vstack((f,np.zeros([1,f.shape[1]])))           
                FS_0 = FS_0 * 0
            else:
                # Adds 2 mm but does not remove sand and smaller grains
                D = np.concatenate([D[0:id2mm],[2],D[id2mm:]])          
                f = np.concatenate([f[0:id2mm,:],[FS_0],f[id2mm:,:]])
        
        # Copies the grain sizes after removing sand if Parker Eq. is used.
        # All grains will be copies if other equation is used.
        self._DOrig = D                                     
        
        # Grain sizes frequency - Now is not cumulative anymore and has the 
        # same dimensions than the equivalent D or DEq
        f = np.abs(-np.diff(f/100,axis=0))                          
        DEq = (D[0:-1] * D[1:])**0.5 # Equivalent grain sizes                              
        FS = np.zeros_like(self._gsd_location)  # Sand fraction at each node             
      
        # Bed grain sizes frequency in each node
        f_DEq = np.zeros([nNodes,D.shape[0]-1])
        for i in np.arange(0,nGSDLocations):
            (id_gsdLocation,) = np.where(self._gsd_location==i)
            FS[id_gsdLocation] = FS_0[i]/100
            f_DEq[id_gsdLocation,:] = f[:,i]
                
        # Calculating D50 at each node - Grain sizes in Psi scale
        Psi_D = np.flip(np.log2(D),axis=0)  
        # Frequency of each grain size in each node
        f_D = np.hstack((np.zeros([self._gsd_location.shape[0],1]),\
                         np.cumsum(np.flip(f_DEq,axis=1),axis=1))) 
        
        # Finds the index of the grain size smaller than 50% in each node
        i0 = np.argmin(f_D<=0.5,axis=1) - 1 
        # Finds the index of the grain size larger than 50% in each node
        i1 = np.argmax(f_D>0.5,axis=1)
        nodesList = np.arange(0,nNodes)
        
        Psi_D50 = Psi_D[i0] + \
            ( (Psi_D[i1]-Psi_D[i0]) / (f_D[nodesList,i1]-f_D[nodesList,i0]) ) \
                * (0.5 - f_D[nodesList,i0])
        D50 = 2**Psi_D50 # Median grain size in each node 
         
        # Calculating the geometric mean and standard deviation
        Psi_DEq = np.log2(DEq) # Equivalent grain sizes in Psi scale
        Psi_DEq_mean = np.sum(f_DEq*Psi_DEq,axis=1)
        geoMean = 2**Psi_DEq_mean
        # Equivalent grain sizes in Psi scale at each node
        Psi_DEq = np.tile(Psi_DEq,(nNodes,1)) 
        Psi_DEq_mean = np.reshape(Psi_DEq_mean,[nNodes,1])
        geoStd = 2**np.sqrt(np.sum( ((Psi_DEq - Psi_DEq_mean)**2) \
                                   * f_DEq, axis=1))
       
        # We save all data into the grid so it can be shared with different 
        # functions and components        
        # Bed grain sizes frequency in each node
        self._grid['node']['bed_surface__grainSizeDistribution'] = f_DEq
        self._grid['node']['bed_surface__grainSizeDistribution_Original'] = copy.deepcopy(f_DEq)  
        self._grid['node']['bed_subsurface__grainSizeDistribution'] = copy.deepcopy(f_DEq) 
        self._grid['node']['bed_surface__medianSize'] = D50
        self._grid['node']['bed_surface__geometricMeanSize'] = geoMean
        self._grid['node']['bed_surface__standardDevSize'] = geoStd
       
        # GSD properties are in nodes. Now we map this information into links     
        # Grain sizes frequency
        self._grid['link']['bed_surface__grainSizeDistribution'] = \
            0.5 * (f_DEq[self.grid.nodes_at_link[:,0]] \
                   + f_DEq[self.grid.nodes_at_link[:,1]])
        self._grid['link']['bed_surface__grainSizeDistribution_Original'] = \
            copy.deepcopy(self._grid['link']['bed_surface__grainSizeDistribution'])
        self._grid['link']['bed_subsurface__grainSizeDistribution'] = \
            copy.deepcopy(self._grid['link']['bed_surface__grainSizeDistribution'])    
        self._grid['link']['bed_surface__medianSize'] = \
            0.5 * (D50[self.grid.nodes_at_link[:,0]] \
                   + D50[self.grid.nodes_at_link[:,1]])     
        self._grid['link']['bed_surface__geometricMeanSize'] = \
            0.5 * (geoMean[self.grid.nodes_at_link[:,0]] \
                   + geoMean[self.grid.nodes_at_link[:,1]])
        self._grid['link']['bed_surface__standardDevSize'] = \
            0.5 * (geoStd[self.grid.nodes_at_link[:,0]] \
                   + geoStd[self.grid.nodes_at_link[:,1]])    
        self._grid['link']['bed_surface__sandFraction'] = \
            0.5 * (FS[self.grid.nodes_at_link[:,0]] \
                   + FS[self.grid.nodes_at_link[:,1]])

        self.calculateActiveLayerThickness()
        
        # Sets the initial active layer thickness as zero meter deep.
        self._grid['node']['bed_surface__activeLayerThickness_previousTimeStep'] = \
            copy.deepcopy(self._grid['node']['bed_surface__activeLayerThickness']*0)
        self._grid['link']['bed_surface__activeLayerThickness_previousTimeStep'] = \
            copy.deepcopy(self._grid['link']['bed_surface__activeLayerThickness']*0)               
        # All equations will work with equivalent grains sizes
        # We droped the "Eq" in this definition
        self._D = DEq     

    def link_at_border_cells(self):

        k = self._grid.links_at_node
        
        l_r = k[:,0]  # list of links on the right face of each node
        l_l = k[:,2]  # list of links on the left face of each node
        l_u = k[:,1]  # list of links on the upper face of each node
        l_d = k[:,3]  # list of links on the lower face of each node
        
        lbr = l_l[self._grid.nodes_at_right_edge]   # Gives the list of links next to a border cell on the right side
        lbu = l_d[self._grid.nodes_at_top_edge]     # Gives the list of links next to a border cell on the top side
        lbl = l_r[self._grid.nodes_at_left_edge]    # Gives the list of links next to a border cell on the bottom side
        lbd = l_u[self._grid.nodes_at_bottom_edge]  # Gives the list of links next to a border cell on the left side
        
        self._link_at_boder_cells = np.sort(np.hstack((lbr,lbu,lbl,lbd)))
        
    def run_one_step(self):        
        """
        The component can be divided in two parts. In the first part all bed
        load transport and GSD calculations, including shear stress estimates, 
        are conducted. In the second part bed GSD and bed elevation can evolve
        
        --- First part ---
        Calculates shear stress and bed load transport rates across the grid.

        For one time step, this generates the shear stress across a given
        grid by accounting for the local water depth, bed surface slope, and
        water velocity at links. Then, based on this shear stress and local
        bed surface grain size distribution the bed load transport rate is 
        calculated at each link and mapped onto each node.Bed load grain size
        distributions are calculated when using Parker 1991 or Wilcock and 
        Crowe 2003 equations. Meyer-Peter and Muller and Fernandez Luque and
        van Beek models will only calculate the total bed load transport.

        Outputs the following bed surface properties through time at every link
        in the input grid: geometric mean size, GSD, median size, sand fraction, 
        std deviation size
        
        Also outputs the bed load GSD,bed load Rate, and shear stress through 
        time at every link in the input grid. The net bed load is output 
        through time at every node.
        
        """                                 
        if self._firstIteration is True:
            # Identify nodes with fixed elevations and surface GSD
            self.fixedNodeInfo()
            # Identify links with fixed bedload GSD
            self.fixedLinkInfo()
            # Identify the node upstream of the outlet
            self.outletNodeInfo() 
        
        # Unsteady shear stress is calculated every time step
        self.shearStress()
                               
        # Selects one bedload transport model to conduct all calculations       
        if self._bedloadEq == 'Parker1990':           
            self.bedloadParker1990()
        elif self._bedloadEq == 'MPM':           
            self.bedloadMeyerPeterMuller()   
        elif self._bedloadEq == 'FLvB':
            self.bedloadFernandezLuqueVanBeek()
        elif self._bedloadEq == 'WilcockAndCrowe':
            self.bedloadWilcockAndCrowe2003()

        #Calculates the net bedload transport from links into nodes
        self.calculateNetBedload()  
        
        # Maps bed load vectors from links (m2/s) onto nodes (m3/s) preserving
        # vector components.
        if self._output_vector is True:           
            self._tau_vector,self._tau_magnitude = self.vector_mapper(self._tau)       
            self._u_vector,self._u_magnitude = self.vector_mapper(self._u)
            self._qb_vector,self._qb_magnitude = self.bedload_mapper()
              
        """Second Part ---
        
        Erode grid topography.

        For one time step, this erodes into the grid topography acording to
        Exner equation.
        
        (1-λp) ∂Z/∂t = - (∂qbx/∂x + ∂qby/∂y)
         simplifying
         ∂Z/∂t = - (1 / (1-λp)) * (∂Qb/∂A)
         Z_t+1 = -(Δt * ΔQb)/(1-λp) + Z_t

        The grid field 'topographic__elevation' is altered each time step.
        """
        if self._evolveBed is True:
        
        # Enters to this routine only after the second time step to avoid
        # re-calculation. 
        # When the bed has changed new bed properties are calculated. This 
        # updated information will be used when calculating bedload rates.        
                                                                                                    
            # Update the bed elevation
            self.updateBedElevation()
            
            # Update the bed grain size distribution
            if (self._bedloadEq == 'Parker1990') \
                or (self._bedloadEq == 'WilcockAndCrowe'):
                if self._updateBedSurfaceGSD is True:
                    self.updateBedGSD()
                    self.updateBedProperties()
                self.map_GSD_from_link_to_node()
                
        self._firstIteration = False
        self._stratigraphyCycle += 1
        
    def fixedNodeInfo(self):
        """ Search and identify nodes defined as having a fixed elevation """
                
        # Gives the ID of the outlet node
        (self._fixedNodes_id,) = np.where(self._grid['node']['bed_surface__fixedElevation']==1)  
        (self._fixedSurfaceGSDNodes_id,) = np.where(self._grid['node']['bed_surface__fixed_grainSizeDistribution']==1) 
        
        # All connecting links to these nodes will be set as fixed too
        self._fixedLinks = np.unique(self._grid.links_at_node[self._fixedSurfaceGSDNodes_id])
        self._fixedLinks = self._fixedLinks[np.where(self._fixedLinks>=0)]
       
    def fixedLinkInfo(self):
        """ Search and identify links defined as having a fixed bed load GSD """
                
        # Gives the ID of the links
        if np.max(self._grid['link']['sediment_transport__imposed_bedloadGSD'] > 0):
            (a0,a1) = np.where(self._grid['link']['sediment_transport__imposed_bedloadGSD'] > 0)  
            self._fixedSurfaceGSDNodes_id_row = a0
            self._fixedSurfaceGSDNodes_id_col = a1
           
            # We use this flag to avoid extra calculations when there is no fixed 
            # bed load GSD
            if a0.shape[0] > 0:
                self._fixedLinkCalculate = True
        
    def calculateSandFraction(self):
        Psi_D = np.flip(np.log2(self._gsd[:,0]),axis=0)
        f_D = np.flip(self._gsd[:,1:],axis=0)
        FS = np.zeros(f_D.shape[1]) 
        i = np.max(np.where(Psi_D<=1))
        for j in np.arange(0,f_D.shape[1]):       
            FS[j] = ((1 - Psi_D[i]) /( Psi_D[i+1]-Psi_D[i] )) \
                *(f_D[i+1,j]-f_D[i,j]) + f_D[i,j]  
            
        return FS
    
    def calculateActiveLayerThickness(self):
        """ First for nodes """
        
        D = self._DOrig                     # Grain sizes
        Psi_D = np.flip(np.log2(D),axis=0)  # Grain sizes in Psi scale
        
        f_DEq = self._grid['node']['bed_surface__grainSizeDistribution']
        f_list = np.arange(0,self._grid.number_of_nodes)

        f = np.hstack((f_DEq,np.zeros([f_DEq.shape[0],1]))) # Grain sizes freq
        
        f_D = np.cumsum(np.flip(f,axis=1),axis=1) # Cumulative GSD in each node
        # Finds the index of the grain size smaller (i0) and larger (i1) than 
        # the fraction requested (fX) in each node
        i0 = np.argmin(f_D<=0.9,axis=1) - 1  
        i1 = np.argmax(f_D>0.9,axis=1)      
                
        Psi_DX = Psi_D[i0] + \
            ( (Psi_D[i1]-Psi_D[i0]) / (f_D[f_list,i1]-f_D[f_list,i0]) ) \
                * (0.9 - f_D[f_list,i0])
        D90_surface = 2**Psi_DX
        self._grid['node']['bed_surface__activeLayerThickness'] = 2 * D90_surface/1000
        
        """ Now for links"""
        
        f_DEq = self._grid['link']['bed_surface__grainSizeDistribution']
        f_list = np.arange(0,self._grid.number_of_links)
        
        f = np.hstack((f_DEq,np.zeros([f_DEq.shape[0],1]))) # Grain sizes freq
        
        f_D = np.cumsum(np.flip(f,axis=1),axis=1) # Cumulative GSD in each node
        # Finds the index of the grain size smaller (i0) and larger (i1) than 
        # the fraction requested (fX) in each node
        i0 = np.argmin(f_D<=0.9,axis=1) - 1  
        i1 = np.argmax(f_D>0.9,axis=1)      
                
        Psi_DX = Psi_D[i0] + \
            ( (Psi_D[i1]-Psi_D[i0]) / (f_D[f_list,i1]-f_D[f_list,i0]) ) \
                * (0.9 - f_D[f_list,i0])
        D90_surface = 2**Psi_DX
        self._grid['link']['bed_surface__activeLayerThickness'] = 2 * D90_surface/1000
    
    def outletNodeInfo(self):
        """ Search and identify the node upstream the outlet to apply boundary 
        conditions"""
                
        # Gives the ID of the outlet node
        (self._out_id,) = np.where(self._grid.status_at_node==1) 
        
        #nCols = self._grid.number_of_node_columns 
        #nRows = self._grid.number_of_node_rows
        # The outlet can be located on the bottom, top, right, or left edge of 
        # the Raster. Depending on the edge, the upstream node is then identify
        # upstream1 is upstream of the upstream node
            
        if np.isin(self._out_id,self._grid.nodes_at_right_edge).any():
            self._upstream_out_id = self._out_id - 1
            up_upstream_out_id = self._upstream_out_id - 1 
            up_up_upstream_out_id = up_upstream_out_id - 1 
            out_direction = 0
             
        elif np.isin(self._out_id,self._grid.nodes_at_top_edge).any():
            self._upstream_out_id = self._out_id - self._grid._shape[1]
            up_upstream_out_id = self._upstream_out_id - self._grid._shape[1]
            up_up_upstream_out_id = up_upstream_out_id - self._grid._shape[1]
            out_direction = 1
             
        elif np.isin(self._out_id,self._grid.nodes_at_left_edge).any():
            self._upstream_out_id = self._out_id + 1 
            up_upstream_out_id = self._upstream_out_id + 1
            up_up_upstream_out_id = up_upstream_out_id + 1
            out_direction = 2
            
        elif np.isin(self._out_id,self._grid.nodes_at_bottom_edge).any():
            self._upstream_out_id = self._out_id + self._grid._shape[1]
            up_upstream_out_id = self._upstream_out_id + self._grid._shape[1]
            up_up_upstream_out_id = up_upstream_out_id + self._grid._shape[1]
            out_direction = 3
             
        outletNodes = np.array((self._out_id,self._upstream_out_id,up_upstream_out_id))
        outletNodes = np.sort(outletNodes.flatten())
        outletLinks = np.unique(self._grid.links_at_node[outletNodes])
        outletLinks = outletLinks[np.where(outletLinks>=0)]
        
        # 4 comes from the number of levels that are affected by the OverlandFlow boundary
        # conditios. 4 rows or columns of links.
        outletLinks_4 = np.unique(self._grid.links_at_node[up_up_upstream_out_id])
        outletLinks_4 = outletLinks_4[np.where(outletLinks_4>=0)]
        outletLinks_4 = outletLinks_4 [~np.in1d(outletLinks_4,outletLinks)]
        
        outletLinksH = outletLinks[np.in1d(outletLinks,self._grid.horizontal_links)]
        outletLinksV = outletLinks[np.in1d(outletLinks,self._grid.vertical_links)] 
        self._outletLinks_4H = outletLinks_4[np.in1d(outletLinks_4,self._grid.horizontal_links)]
        self._outletLinks_4V = outletLinks_4[np.in1d(outletLinks_4,self._grid.vertical_links)] 
        
        nOutlets = self._out_id.shape[0]
        # 3 comes from the number of levels that we are adjusting.
        if (out_direction == 0) or (out_direction == 2):     # Outlet is at the right or left edge
            self._outletLinksH = np.reshape(outletLinksH,[nOutlets,3]).T     
            self._outletLinksV = np.reshape(outletLinksV,[3,nOutlets+1]).T 
            
        elif (out_direction == 1) or (out_direction == 3):   # Outlet is at the top or bottom edge
            self._outletLinksH = np.reshape(outletLinksH,[nOutlets+1,3])      
            self._outletLinksV = np.reshape(outletLinksV,[3,nOutlets]) 
                   
    def map_GSD_from_link_to_node(self):
        """Maps the bed surface grain size distribution from links to nodes.
        
        Given that the all our calculations are conducted in links we implemented
        thios function to display results in a raster or in nodes.
        
        """
        linksId = self._grid.links_at_node
        nNodes = self._grid.number_of_nodes
        F = self._grid['link']['bed_surface__grainSizeDistribution']
        F_nodes = 0.25 * (F[linksId[:,0],:] + F[linksId[:,1],:] + F[linksId[:,2],:] + F[linksId[:,3],:])
        F_nodes = F_nodes/np.reshape(np.sum(F_nodes,axis=1),[nNodes,1])
        
        # Revert any changes to the fixed GSD nodes and fixed elevations nodes
        fixedNodes = np.unique(np.hstack((self._fixedSurfaceGSDNodes_id,self._fixedNodes_id)))
        
        F_nodes[fixedNodes] = self._grid['node']['bed_surface__grainSizeDistribution_Original'][fixedNodes]  
        self._grid['node']['bed_surface__grainSizeDistribution'] = copy.deepcopy(F_nodes)
        
    def updateBedProperties(self):

        " Calculates the updated GSD properties"
        nLinks = self._grid.number_of_links
        D = self._DOrig     # Grain sizes
        DEq = self._D       # Equivalent grain sizes
        # Equivalent grain sizes frequency
        f_DEq = self._grid['link']['bed_surface__grainSizeDistribution']    
        # Grain sizes frequency
        f = np.hstack((f_DEq,np.zeros([f_DEq.shape[0],1])))                         
        
        """ Calculates D50 at each node based on the updated GSD"""
        Psi_D = np.flip(np.log2(D),axis=0)        # Grain sizes in Psi scale
        f_D = np.cumsum(np.flip(f,axis=1),axis=1) # Cumulative GSD in each link
        
        # Finds the index of the grain size smaller than 50% in each link
        i0 = np.argmin(f_D<=0.5,axis=1) - 1 
        # Finds the index of the grain size larger than 50% in each link
        i1 = np.argmax(f_D>0.5,axis=1) 
        linkList = np.arange(0,nLinks)
        
        Psi_D50 = Psi_D[i0] + \
            ( (Psi_D[i1]-Psi_D[i0]) / (f_D[linkList,i1]-f_D[linkList,i0]) ) \
                * (0.5 - f_D[linkList,i0])
        D50 = 2**Psi_D50                # Median grain size in each link
        
        # Calculates the geometric mean and standard deviation
        # Equivalent grain sizes in Psi scale
        Psi_DEq = np.log2(DEq)                          
        Psi_DEq_mean = np.sum(f_DEq*Psi_DEq,axis=1)
        geoMean = 2**Psi_DEq_mean       # Geometric mean size in each link
        
        # Equivalent grain sizes in Psi scale at each link
        Psi_DEq = np.tile(Psi_DEq,(nLinks,1)) 
        Psi_DEq_mean = np.reshape(Psi_DEq_mean,[nLinks,1])
        # Standard deviation at each node
        geoStd = 2**np.sqrt(np.sum( ((Psi_DEq - Psi_DEq_mean)**2) * f_DEq, axis=1)) 
        
        self._grid['link']['bed_surface__medianSize'] = D50
        self._grid['link']['bed_surface__geometricMeanSize'] = geoMean
        self._grid['link']['bed_surface__standardDevSize'] = geoStd        
              
        self._grid['link']['bed_surface__activeLayerThickness_previousTimeStep'] = \
            copy.deepcopy(self._grid['link']['bed_surface__activeLayerThickness'])            
        self.calculateActiveLayerThickness()
                    
    def shearStress(self):
        
        """ Unsteady shear stress calculated at links according to
        τ = rho * g * h * Sf
         
        Where Sf is the unsteady friction slope and is calculated as
        Sf = S0 - dh/ds - U/g du/ds - 1/g du/dt
        
        Alternatively, tau can be calculated as τ = rho * g * Rh * Sf
        but need to be manually selected (for consistency in the code)        

        The term ds indicates a certain direction, X and Y in this case. 
        All derivatives are first or second order approximations in each 
        direction. Most of grid parameters are read again because they may have
        change in time in overland flow or when the bed evolved.
        """
        
        #Reads the current topographic elevation
        z = self._grid.at_node['topographic__elevation']
        
        #First term in the friction slope - gradient of bed elevation 
        #S0 = - dz/ds
        self._dz_ds = - self._grid.calc_grad_at_link(z)
        
        #Second term in the friction slope - gradient of water depth
        h = self._grid['node']['surface_water__depth']
        dh_ds = self._grid.calc_grad_at_link(h)
        
        h_links = self._grid.at_link['surface_water__depth']

        # Third term in the friction slope - gradient of flow velocity
        # Velocity at current time step
        self._u = self._grid['link']['surface_water__velocity']
        
        # Velocity gradients are calculated in each direction
        du_ds = np.zeros_like(self._dz_ds)
        
        # In the X direction only horizontal gradients are valid. Here we use 
        # the map_mean_of_horizontal_links_to_node method. However, this tool 
        # is valid only for the X direction
        u_nodes_horizontal = self._grid.map_mean_of_horizontal_links_to_node(self._u)
        du_ds[self._grid.horizontal_links] = \
            self._grid.calc_grad_at_link(u_nodes_horizontal)[self._grid.horizontal_links]
        
        # In the Y direction an extra step is required. 
        u_nodes_vertical = self._grid.map_mean_of_vertical_links_to_node(self._u)
        u_nodes_vertical = \
            np.flip(np.flip((np.reshape(u_nodes_vertical,(self._grid._shape[0],self._grid._shape[1])))),axis=1)
        
        du_dsV = np.zeros((u_nodes_vertical.shape[0]-1,u_nodes_vertical.shape[1]))
        
        for i in np.arange(u_nodes_vertical.shape[0]-2,-1,-1):
            du_dsV[i,:] = (u_nodes_vertical[i,:] - u_nodes_vertical[i+1,:]) / self._grid.dy
        
        # Three operations are in the right side. Here we recover landLab 
        # indexing order 
        du_ds[self._grid.vertical_links] = np.flip(du_dsV.T,axis=1).flatten(order='F')        
        
        # Fourth term in the friction slope - rate of change of velocity in time
        u_t0 = self._grid['link']['surface_water__velocity_previous_time_step']
        du_dt = (self._u - u_t0) / self._grid._dt
        
        # Friction slope calculation at links including unsteady effects
        Sf = self._dz_ds - dh_ds - (self._u/self._g) * du_ds - 1/self._g * du_dt
        
        # And finally, the shear stress at links including unsteady effects
        if self._useHydraulicsRadiusInShearStress is True:
            # uses hydraulics ratio, so -shear stress is calcualtes as
            # τ = Rho * g * Rh * Sf
            
            # Different grid sizes (dx~=dy) are possible 
            A = np.zeros_like(h_links)
            A[self._grid.horizontal_links] = h_links[self._grid.horizontal_links] * self._grid.dx
            A[self._grid.vertical_links] = h_links[self._grid.vertical_links] * self._grid.dy
            
            P = np.zeros_like(h_links)
            P[self._grid.horizontal_links] = self._grid.dx + 2 * h_links[self._grid.horizontal_links]
            P[self._grid.vertical_links] = self._grid.dy + 2 * h_links[self._grid.vertical_links]
            Rh = A / P #Rh = wetted area / wetted perimeter#Rh = wetted area / wetted perimeter
            self._tau = self._rho * self._g * Rh * Sf
        else:
            # Equation is τ = Rho * g * h * Sf """
            self._tau = self._rho * self._g * h_links * Sf 
        
        # Apply a boundary condition of zero flux at link next to borders
        self._tau[self._link_at_boder_cells] = 0
        
        # Direction of flux will be recovered at the end of the bedload 
        # transport routine
        self._tauT = np.abs(self._tau)
        
        # Now we write the shear stress field into the grid so other components
        # or postprocess functions can access it
        self._grid['link']['surface_water__shearStress'] = self._tauT
        
    def bedloadParker1990(self):
        
        """ Surface-based bedload transport equation of Parker 1990
       
        G. Parker (1990) Surface-based bedload transport relation for gravel rivers, 
        Journal of Hydraulic Research, 28:4, 417-436, DOI: 10.1080/00221689009499058
        """

        # Variables definition - All these variables are updated each time step. 
        # Therefore, we read them again to ensure they are updated
        
        f = self._grid['link']['bed_surface__grainSizeDistribution']
        geoMean = self._grid['link']['bed_surface__geometricMeanSize']
        geoStd = self._grid['link']['bed_surface__standardDevSize']
        
        tauStar_sg = self._tauT/(self._rho*self._R*self._g*(geoMean/1000)) 
        self._phi_sgo = tauStar_sg/self._tauStar_rsgo
        
        self.strainFunctions()
        self._omega = 1+(np.log2(geoStd)/self._sigma0)*(self._omega0-1)
        
        Di_Dsg = np.tile(self._D,(self._nLinks,1))/(np.reshape(geoMean,[self._nLinks,1]))
        phi_sgo = np.reshape(self._phi_sgo,[self._phi_sgo.shape[0],1])
        omega = np.reshape(self._omega,[self._omega.shape[0],1])
        phi_i = omega* phi_sgo*(Di_Dsg)**-self._beta
        
        G = np.zeros_like(phi_i)
        
        # There are three intervals where G is evaluated
        (id0,id1) = np.where(phi_i > 1.59)
        if id0.shape[0]>0:
            G[id0,id1] = 5474 * (1 - 0.853/phi_i[id0,id1]) ** 4.5
        
        (id0,id1) = np.where( (phi_i >= 1) & (phi_i <= 1.59))
        if id0.shape[0]>0:
            G[id0,id1] = np.exp(14.2*(phi_i[id0,id1]-1)-9.28*(phi_i[id0,id1]-1) ** 2)
        
        (id0,id1) = np.where(phi_i < 1)
        if id0.shape[0]>0:        
            G[id0,id1] = phi_i[id0,id1] ** 14.2
            
        Gf = f * G
        GfSum = np.sum(Gf,axis=1)
        
        # Total bedload transport rate
        Wstar_s = 0.00218 * GfSum
        self._grid['link']['sediment_transport__bedloadRate'] = \
            (((np.sqrt(self._tauT/self._rho))**3*Wstar_s)/(self._R*self._g)) \
                * np.sign(self._tau)
        
        # Now that bed load rate has been calculated in all links we replace those that
        # were imposed. If there there are no imposed bedload GSD this will do nothing.
        if self._fixedLinkCalculate is True:
            linkId_fixedLinks = np.unique(self._fixedSurfaceGSDNodes_id_row)
                            
            self._grid['link']['sediment_transport__bedloadRate'][linkId_fixedLinks] = \
                self._grid['link']['sediment_transport__imposed_sediment_supply'][linkId_fixedLinks]
        
        # When calculating fractional sediment transport some fractions will 
        # be zero and a warning should be displayed because it will divide by 
        # zero. We know that it is not a problem because in those locations 
        # there is zero transport, so we deactivate the warning
        np.seterr(divide='ignore', invalid='ignore')
        
        # Frational bedload transport rate
        p = np.zeros_like(Gf)
        GfSum = np.transpose(np.tile(GfSum,(f.shape[1],1)))
        p = Gf / GfSum
        id0 = np.isnan(p)
        p[id0] = 0
        self._grid['link']['sediment_transport__bedloadGSD'] = p
        
        # Now that bedload GSD has been calculated in all links we replace those that
        # were imposed. If there there are no imposed bedload GSD this will do nothing.
        if self._fixedLinkCalculate is True:
            linkId_0 = self._fixedSurfaceGSDNodes_id_row
            linkId_1 = self._fixedSurfaceGSDNodes_id_col
                            
            self._grid['link']['sediment_transport__bedloadGSD'][linkId_0,linkId_1] = \
                self._grid['link']['sediment_transport__imposed_bedloadGSD'][linkId_0,linkId_1]
            
    def strainFunctions(self):
        
        omega0 = np.zeros_like(self._phi_sgo)
        sigma0 = np.zeros_like(self._phi_sgo)
            
        # There are three intervals where we can interpolate Omega and Sigma
        (I,) = np.where(self._phi_sgo <= 0.7639)
        if I.shape[0]>0:
            omega0[I] = self._oo[0]
            sigma0[I] = self._so[0]
        
        (I,) = np.where(self._phi_sgo > 231.2)
        if I.shape[0]>0:
            omega0[I] = np.interp(self._phi_sgo, self._po, self._oo)[I]
            sigma0[I] = np.interp(self._phi_sgo, self._po, self._oo)[I]
        
        (I,) = np.where( (self._phi_sgo > 0.7639) & (self._phi_sgo < 231.2))
        if I.shape[0]>0:        
            foo = interp1d(self._po, self._oo, kind='cubic')
            fso = interp1d(self._po, self._so, kind='cubic')
            omega0[I] = foo(self._phi_sgo[I])
            sigma0[I] = fso(self._phi_sgo[I])
        
        self._omega0 = omega0
        self._sigma0 = sigma0

    def bedloadWilcockAndCrowe2003(self): 
        
        """ Surface-based bedload transport equation of Wilcock and Crowe 2003
        
        Wilcock, P. R., & Crowe, J. C. (2003). Surface-based transport model 
        for mixed-size sediment. Journal of hydraulic engineering, 129(2), 
        120-128.
        """

        # Variables definition - All these variables are updated each time step. 
        # Therefore, we read them again to ensure they are updated
        f = self._grid['link']['bed_surface__grainSizeDistribution']
        geoMean = self._grid['link']['bed_surface__geometricMeanSize']
        FS = self._grid['link']['bed_surface__sandFraction']
        
        tauStar_sg = self._tauT/(self._rho*self._R*self._g*(geoMean/1000))                
        tauStar_rsg0 = 0.021 + 0.015 * np.exp(-20 * FS)                        
        phi_sg0 = tauStar_sg/tauStar_rsg0
        
        phi_i = np.zeros((phi_sg0.shape[0],self._D.shape[0]))
        
        G = np.zeros((phi_sg0.shape[0],self._D.shape[0]))
        p = np.zeros((phi_sg0.shape[0],self._D.shape[0]))
        
        for i in np.arange(0,self._D.shape[0]):
            b = 0.67 / (1 + np. exp(1.5 - self._D[i]/geoMean))
            phi_i[:,i] = phi_sg0 * (self._D[i]/(geoMean))**(-b)
            
        # There are two intervals where G is evaluated
        (id0,id1) = np.where(phi_i >= 1.35)
        if id0.shape[0]>0:
            G[id0,id1] = 14 * (1 - 0.894/ (phi_i[id0,id1]**0.5)) ** 4.5
        
        (id0,id1) = np.where(phi_i < 1.35)
        if id0.shape[0]>0:        
            G[id0,id1] = 0.002* phi_i[id0,id1] ** 7.5
            
        Wstar_i = f * G
        Wstar = np.sum(Wstar_i,axis=1)
        
        # Total bedload transport rate is calculated at each link
        # Notice that Wstar includes fi (Eq 2 in the paper)
        self._grid['link']['sediment_transport__bedloadRate'] = \
            (((np.sqrt(self._tauT/self._rho))**3*Wstar)/(self._R*self._g))  \
                * np.sign(self._tau)
        
        # Frational bedload transport rate is calculated as a fraction of the
        # total bedload transport rate at each link When calculating fractional 
        # sediment transport some fractions will be zero and a warning should 
        # be displayed because it will divide by zero. We know that it is not 
        # a problem becuase in those locations there is zero transport, so we 
        # deactivate the warning
        np.seterr(divide='ignore', invalid='ignore')
        p = Wstar_i/ np.reshape(Wstar,[Wstar.shape[0],1])
        id0 = np.isnan(p)
        p[id0] = 0
        self._grid['link']['sediment_transport__bedloadGSD'] = p    
        
        # Now that bedload GSD has been calculated in all links we replace those that
        # were imposed. If there there are no imposed bedload GSD this will do nothing.
        if self._fixedLinkCalculate is True:
            linkId_0 = self._fixedSurfaceGSDNodes_id_row
            linkId_1 = self._fixedSurfaceGSDNodes_id_col
                            
            self._grid['link']['sediment_transport__bedloadGSD'][linkId_0,linkId_1] = \
                self._grid['link']['sediment_transport__imposed_bedloadGSD'][linkId_0,linkId_1]
                
    def bedloadMeyerPeterMuller(self):   
        
        """ Surface-based bedload transport equation of Meyer-Peter and Müller
        
        Meyer-Peter, E. and Müller, R., 1948, Formulas for Bed-Load Transport, 
        Proceedings, 2nd Congress, International Association of Hydraulic 
        Research, Stockholm: 39-64.
        """        
        D50 = self._grid['link']['bed_surface__medianSize']
        tauStar = self._tauT/(self._rho*self._R*self._g*(D50/1000))
        tauCrStar = 0.047
        if self._variableCriticalShearStress == True:
            tauCrStar = self.criticalShearStress(tauCrStar)
        qbStar = np.where(tauStar-tauCrStar > 0 , 8 * np.abs(tauStar-tauCrStar)**(3/2) , 0) 
       
        self._grid['link']['sediment_transport__bedloadRate'] = \
            ( qbStar * (np.sqrt(self._R*self._g*(D50/1000)) * (D50/1000)) ) \
                * np.sign(self._tau)      
        
    def bedloadFernandezLuqueVanBeek(self): 
        
        """ Surface-based bedload transport equation of Fernandez Luque and 
        van Beek
        
        Fernandez Luque, R. and R. van Beek, 1976, Erosion and transport of 
        bedload sediment, Journal of Hydraulic Research, 14(2): 127-144.
        """
        D50 = self._grid['link']['bed_surface__medianSize']                
        tauStar = self._tauT/(self._rho*self._R*self._g*(D50/1000))
        tauCrStar = 0.045
        if self._variableCriticalShearStress == True:
            tauCrStar = self.criticalShearStress(tauCrStar)
        qbStar = np.where(tauStar-tauCrStar > 0, 5.7 * np.abs(tauStar-tauCrStar)**(3/2) ,0)
       
        # Total bedload transport rate is calculated at each link
        self._grid['link']['sediment_transport__bedloadRate'] = \
            ( qbStar * (np.sqrt(self._R*self._g*(D50/1000)) * (D50/1000)) ) \
                * np.sign(self._tau)      
    
    def criticalShearStress(self,tauCrStar):
        """ We calculate the current bed slope at links in case only when is used
        # to calculate a slope dependent critical shear stress """
        bedSlope = np.abs(self._dz_ds) # Direction is not important
        # Mueller et al. (2005) equation is used in steep slopes
        return np.where(bedSlope>0.03,2.18 * bedSlope + 0.021,tauCrStar) 
                         
    def calculateNetBedload(self):        
        """Calculates the net volumetric bedload coming from all links (m2/s) 
        onto nodes (m3/s).

        This method takes the volumetric bedload entering and exiting through a 
        face and determines the net volumetric bedload on a given node.
        
        """
        # Reads and modify the field sediment_transport__bedloadRate to account
        # for links where sediment supply is imposed.
        (self._Id_upSedSupp,)= np.nonzero(self._grid['link']['sediment_transport__imposed_sediment_supply'])
        self._grid['link']['sediment_transport__bedloadRate'][self._Id_upSedSupp] = \
            self._grid['link']['sediment_transport__imposed_sediment_supply'][self._Id_upSedSupp]

        qbx = np.sum( self._grid['link']['sediment_transport__bedloadRate'][self._grid.links_at_node[:,[0,2]]] * self._normal[:,[0,2]], axis = 1 ) * self._grid.dy
        qby = np.sum( self._grid['link']['sediment_transport__bedloadRate'][self._grid.links_at_node[:,[1,3]]] * self._normal[:,[1,3]], axis = 1 ) * self._grid.dx
                   
        net_qb_nodes = qbx + qby
        
        # At the boundary we don't have a link exiting, therefore we assume 
        # zero flux exiting. This is override in the Exner equation with a 
        # zero gradient boundary condition
        net_qb_nodes[self._grid.boundary_nodes] = 0
                
        self._grid['node']['sediment_transport__netBedload'] = net_qb_nodes

    def updateBedElevation(self):
        """ Applies the Exner equation and boundary conditions to predict
        the change in bed surface elevation"""
        
        # Bed elevation is updated using Exner equation       
        A = self._grid.dx * self._grid.dy
        DQb = self._grid['node']['sediment_transport__netBedload']
        z0 = self._grid['node']['topographic__elevation']
        Dz = -self._grid._dt / ((1 - self._lambda_p) * A) * DQb       
        
        self._grid['node']['topographic__elevation'] += Dz   
        
        # The outlet node has been modifed, but not according to the specified
        # boundary condition. Here, we return the outlet to the previous state
        self._grid['node']['topographic__elevation'][self._out_id] = z0[self._out_id]
        
        # Fixed nodes may have been modifed, but not according to the specified
        # condition. Here, we return the nodes to the orignal state
        self._grid['node']['topographic__elevation'][self._fixedNodes_id] = \
            self._grid['node']['topographic__elevation_original'][self._fixedNodes_id]       
        
        # Now we can apply the boundary conditions
        if self._outletBoundaryCondition == 'zeroGradient':
            Dz_outlet = Dz[self._upstream_out_id]
        elif self._outletBoundaryCondition == 'fixedValue':
            Dz_outlet = 0
            
        self._grid['node']['topographic__elevation'][self._out_id] = \
            z0[self._out_id] + Dz_outlet
        
        # Now we map data into links to update Bed GSD
        self._grid['link']['topographic__elevation'] = \
            self._grid.map_mean_of_link_nodes_to_link(self._grid['node']['topographic__elevation'])
                
        if self._trackStratigraphy is True:# Here we register how deep the deposited/eroded layer is
            self._grid['link']['bed_surface__surface_thickness_new_layer'] = \
                self._grid['link']['topographic__elevation'] - \
                    self._grid['link']['topographic__elevation_subsurface']
            
            # Checks if deposited material needs to be updated
            (self._deepLinksId,) = \
                np.where(self._grid['link']['bed_surface__surface_thickness_new_layer'] > self._newSurfaceLayerThickness)
            self._deepLinksId = self._deepLinksId[np.in1d(self._deepLinksId,self._grid.active_links)]
            self._deepLinksId = self._deepLinksId[~np.in1d(self._deepLinksId,self._fixedLinks)]
            
            if self._deepLinksId.shape[0]>0:
                if np.min(self._grid['link']['bed_surface__activeLayerThickness'][self._deepLinksId]) > self._newSurfaceLayerThickness:
                    print('Warning - New surface layer thickness is too thin compared to the active layer thickness')
                self._computeStratigraphy = True
                self._updateSubSurface = True
            
            # Checks if eroded material needs to be updated
            (self._erodedLinksId,) = \
                np.where(self._grid['link']['bed_surface__surface_thickness_new_layer'] < -self._newSurfaceLayerThickness)
            self._erodedLinksId = self._erodedLinksId[np.in1d(self._erodedLinksId,self._grid.active_links)]
            self._erodedLinksId = self._erodedLinksId[~np.in1d(self._erodedLinksId,self._fixedLinks)]
            
            if self._erodedLinksId.shape[0]>0:
                print('Found :',self._erodedLinksId,'\n')
                self._computeStratigraphy = True
                self._updateErodedSubSurface = True            
   
    def updateBedGSD(self):
        """ Uses the fractional Exner equation to update the bed GSD """
        
        # Here we create a number of variables that will be used in the 
        # following definitions to make the code a little bit cleaner
        # No deep copies are done to avoid unnecesary repetition.
        
        nLinks = self._grid.number_of_links
        nCols = self._grid.number_of_node_columns
        F = self._grid['link']['bed_surface__grainSizeDistribution']
        Fs = self._grid['link']['bed_subsurface__grainSizeDistribution'] 
        pl = self._grid['link']['sediment_transport__bedloadGSD']
        qbT = self._grid['link']['sediment_transport__bedloadRate'] # total bed load transport rate at each link
        La = np.reshape(self._grid['link']['bed_surface__activeLayerThickness'],[nLinks,1])
        Laold = np.reshape(self._grid['link']['bed_surface__activeLayerThickness_previousTimeStep'],[nLinks,1])
        
        lps = self._lambda_p
        dx = self._grid.dx
        dy = self._grid.dy
        alpha = self._alpha
        dt = self._grid._dt 
        
        dv = 2*nCols-1
        
        qbT = np.reshape(qbT,[nLinks,1])        
        qbTdev = np.zeros([nLinks,1])
        qbTdevNeg = np.zeros([nLinks,1])
        qjj1dev = np.zeros_like(pl)
        qjj1devNeg = np.zeros_like(pl)
        
        # Horizontal Links
        hlL = np.arange(0,nLinks-nCols+2,2*nCols-1)   # Horizontal Links at left edge
        hlR = np.arange(nCols-2,nLinks,2*nCols-1)     # Horizontal Links at right edge
        hl_Id = np.in1d(self._grid.horizontal_links,np.hstack((hlL,hlR)))    # Links at the top row
        hl = self._grid.horizontal_links[~hl_Id]
        
        # Vertical Links
        vl = self._grid.vertical_links[nCols:-nCols]   # Links within middle region
        vlB = self._grid.vertical_links[0:nCols]       # Links at the bottom row
        vlT = self._grid.vertical_links[-nCols:None]    # Links at the top row
        
        # First we assume that everywhere flow direction is in the positive direction
        qbTdev[hl] = alpha * (qbT[hl] - qbT[hl - 1]) / dy + \
                        (1 - alpha) * (qbT[hl + 1] - qbT[hl]) / dy
        qbTdev[hlL] = (qbT[hlL + 1] - qbT[hlL]) / dy
        qbTdev[hlR] = (qbT[hlR] - qbT[hlR - 1]) / dy
        
        qjj1dev[hl,:] = alpha * (qbT[hl] * pl[hl,:] - qbT[hl-1] * pl[hl-1,:]) / dy + \
                        (1 - alpha) * (qbT[hl + 1] * pl[hl + 1] - qbT[hl] * pl[hl]) / dy 
        qjj1dev[hlL,:] = (qbT[hlL + 1] * pl[hlL + 1,:] - qbT[hlL] * pl[hlL,:]) / dy 
        qjj1dev[hlR,:] = (qbT[hlR] * pl[hlR,:] - qbT[hlR - 1] * pl[hlR - 1,:]) / dy  
        
        qbTdev[vl] = alpha * (qbT[vl] - qbT[vl-dv]) / dx + \
                        (1 - alpha) * (qbT[vl+dv] - qbT[vl]) / dx
        qbTdev[vlB] = (qbT[vlB+dv] - qbT[vlB]) / dx
        qbTdev[vlT] = (qbT[vlT] - qbT[vlT-dv]) / dx
        
        qjj1dev[vl,:] = alpha * (qbT[vl] * pl[vl,:] - qbT[vl-dv] * pl[vl-dv,:]) / dx + \
                        (1 - alpha) * (qbT[vl + dv] * pl[vl + dv] - qbT[vl] * pl[vl]) / dx
        qjj1dev[vlB,:] = (qbT[vlB+dv] * pl[vlB+dv,:] - qbT[vlB] * pl[vlB,:]) / dx 
        qjj1dev[vlT,:] = (qbT[vlT] * pl[vlT,:] - qbT[vlT-dv] * pl[vlT-dv,:]) / dx 
        
        # Now we correct for flow at locations where it is flowing towards the negative direction
        (hl_neg_id,) = np.where(qbT[hl][:,0]<0) 
        (hlL_neg_id,) = np.where(qbT[hlL][:,0]<0) 
        (hlR_neg_id,) = np.where(qbT[hlR][:,0]<0) 
        
        qbTdevNeg[hl] = alpha * (qbT[hl] - qbT[hl + 1]) / dy + \
                        (1 - alpha) * (qbT[hl - 1] - qbT[hl]) / dy
        qbTdevNeg[hlL] = (qbT[hlL] - qbT[hlL - 1]) / dy
        qbTdevNeg[hlR] = (qbT[hlR - 1] - qbT[hlR]) / dy
        
        qbTdev[hl[hl_neg_id]] = -qbTdevNeg[hl[hl_neg_id]]
        qbTdev[hlL[hlL_neg_id]] = -qbTdevNeg[hlL[hlL_neg_id]]
        qbTdev[hlR[hlR_neg_id]] = -qbTdevNeg[hlR[hlR_neg_id]]
        
        qjj1devNeg[hl,:] = alpha * (qbT[hl] * pl[hl,:] - qbT[hl+1] * pl[hl+1,:]) / dy + \
                        (1 - alpha) * (qbT[hl - 1] * pl[hl - 1] - qbT[hl] * pl[hl]) / dy 
        qjj1devNeg[hlL,:] = (qbT[hlL] * pl[hlL,:] - qbT[hlL-1] * pl[hlL-1,:]) / dy 
        qjj1devNeg[hlR,:] = (qbT[hlR-1] * pl[hlR-1,:] - qbT[hlR] * pl[hlR,:]) / dy 
        
        qjj1dev[hl[hl_neg_id]] = -qjj1devNeg[hl[hl_neg_id]]
        qjj1dev[hlL[hlL_neg_id]] = -qjj1devNeg[hlL[hlL_neg_id]]
        qjj1dev[hlR[hlR_neg_id]] = -qjj1devNeg[hlR[hlR_neg_id]]
        
        (vl_neg_id,) = np.where(qbT[vl][:,0]<0) 
        (vlB_neg_id,) = np.where(qbT[vlB][:,0]<0) 
        (vlT_neg_id,) = np.where(qbT[vlT][:,0]<0) 
        
        qbTdevNeg[vl] = alpha * (qbT[vl] - qbT[vl+dv]) / dx + \
                        (1 - alpha) * (qbT[vl-dv] - qbT[vl]) / dx
        qbTdevNeg[vlB] = (qbT[vlB] - qbT[vlB + dv]) / dx
        qbTdevNeg[vlT] = (qbT[vlT - dv] - qbT[vlT]) / dx
        
        qbTdev[vl[vl_neg_id]] = -qbTdevNeg[vl[vl_neg_id]]
        qbTdev[vlB[vlB_neg_id]] = -qbTdevNeg[vlB[vlB_neg_id]]
        qbTdev[vlT[vlT_neg_id]] = -qbTdevNeg[vlT[vlT_neg_id]]
        
        qjj1devNeg[vl,:] = alpha * (qbT[vl] * pl[vl,:] - qbT[vl+dv] * pl[vl+dv,:]) / dx + \
                        (1 - alpha) * (qbT[vl - dv] * pl[vl - dv] - qbT[vl] * pl[vl]) / dx
        qjj1devNeg[vlB,:] = (qbT[vlB] * pl[vlB,:] - qbT[vlB+dv] * pl[vlB+dv,:]) / dx 
        qjj1devNeg[vlT,:] = (qbT[vlT-dv] * pl[vlT-dv,:] - qbT[vlT] * pl[vlT,:]) / dx 
        
        qjj1dev[vl[vl_neg_id]] = -qjj1devNeg[vl[vl_neg_id]]
        qjj1dev[vlB[vlB_neg_id]] = -qjj1devNeg[vlB[vlB_neg_id]]
        qjj1dev[vlT[vlT_neg_id]] = -qjj1devNeg[vlT[vlT_neg_id]]
        # Correction done 
        
        # Errors are introduced by roundig numbers. These lines can correct those errors
        qbTdev = np.where(np.abs(qbTdev)<1e-7,0,qbTdev) 
        qjj1dev = np.where(np.abs(qjj1dev)<1e-7,0,qjj1dev)  
        
        FIexc = copy.deepcopy(Fs)
        (id0,) = np.where(qbTdev[:,0] <= 0)
        FIexc[id0,:] = 0.7 * F[id0,:] + 0.3 * pl[id0,:]
          
        qjj2dev = FIexc * np.reshape(qbTdev,[nLinks,1])  
        
        Fnew = F + dt * (-qjj1dev+ qjj2dev) / (1 - lps) / La
        if (self._firstIteration is False) and (self._subsurfaceChanged is True):
            Fnew = Fnew+ (FIexc- F) / La * (La - Laold)
            self._subsurfaceChanged = False
        
        (id0,id1) = np.where(Fnew <= 0)    
        Fnew[id0,id1] = 0
        Fnew = Fnew/np.reshape(np.sum(Fnew,axis=1),[nLinks,1])
        Fnew = np.nan_to_num(Fnew)
        
        # Given the way in which OverlandFLow calculates flow near the outlets
        # we need to correct changes to GSD that are caused by the sudden drop in 
        # water depth. This affects always the 3 rows or columns of links near the 
        # outlet. Here we identify those 3 rows or columns
        # Corrects horizontal links
        for i in np.arange(0,self._outletLinks_4H.shape[0]):
            Fnew[self._outletLinksH[:,i],:] = Fnew[self._outletLinks_4H[i],:]
        
        # Corrects vertical links
        for i in np.arange(0,self._outletLinks_4V.shape[0]):
            Fnew[self._outletLinksV[:,i],:] = Fnew[self._outletLinks_4V[i],:] 
            
        # Now we update the bed surface GSD
        self._grid['link']['bed_surface__grainSizeDistribution'] = copy.deepcopy(Fnew)
           
        # If the bed is eroded below the original elevation it restores this 
        # initial GSD. First looks if there is any such link
        (erodedLinks_id,) = np.where(self._grid['link']['topographic__elevation'] < \
                                     self._grid['link']['topographic__elevation_original']) 
        
        if erodedLinks_id.shape[0]>0:
            self._grid['link']['bed_surface__grainSizeDistribution'][erodedLinks_id] = \
                self._grid['link']['bed_surface__grainSizeDistribution_Original'][erodedLinks_id]
        
        # Now, an eroded node cannot return to the original GSD if starts depositing again.
        # It will only use the original GSD if erodes deeper than the maximum that has been eroded
            self._grid['link']['topographic__elevation_original'][erodedLinks_id] = \
                copy.deepcopy(self._grid['link']['topographic__elevation'][erodedLinks_id])
        
        # Revert any changes to the fixed GSD nodes             
        self._grid['link']['bed_surface__grainSizeDistribution'][self._fixedLinks ] = \
            self._grid['link']['bed_surface__grainSizeDistribution_Original'][self._fixedLinks] 
        
        if self._trackStratigraphy is True:
            if (self._firstIteration is True) or \
                (self._stratigraphyCycle >= self._nCyclesToProcessStratigraphy) or \
                    (self._computeStratigraphy is True):
                self.stratigraphy()
 
    def stratigraphy(self):
        z = self._grid['link']['topographic__elevation']
        D = self._DOrig
        F = self._grid['link']['bed_surface__grainSizeDistribution']
        Fs = self._grid['link']['bed_subsurface__grainSizeDistribution'] 
        La = self._grid['link']['bed_surface__activeLayerThickness']
        dzl = self._grid['link']['bed_surface__surface_thickness_new_layer']
        
        activeLinks = self._grid.active_links
                
        if self._firstIteration is True:
            print('Creating files to store the  stratigraphy at time zero s','\n')
            # Creates a folder to store results and a file for each active link and node
            if os.path.exists(self._StratigraphyTempFilesPath):
                print('The folder') 
                print(self._StratigraphyTempFilesPath)
                print('Exists and it will be removed to store temporal files to process the stratigraphy');print(" ")
                shutil.rmtree(self._StratigraphyTempFilesPath)     
            
            os.mkdir(self._StratigraphyTempFilesPath)
            
            if os.path.exists(self._StratigraphyRawDataPath):
                print('The folder') 
                print(self._StratigraphyRawDataPath)
                print('Exists and it will be removed to store the stratigraphy at each active link');print(" ")
                shutil.rmtree(self._StratigraphyRawDataPath)
            
            os.mkdir(self._StratigraphyRawDataPath)
            
            # Now goes to the raw data folder to store the properties at time zero
            # The format is unfriendly but after the simulation is done a function will take
            # care of making it more friendly.
            os.chdir(self._StratigraphyRawDataPath)
                        
            for i in activeLinks:
                filename = 'link_' + str(i) + '.txt'
                data = np.hstack((self.t,z[i],F[i,:],Fs[i,:]))
                data = np.reshape(data,[1,data.shape[0]])
                with open(filename, 'ab') as f:
                    np.savetxt(f, data,'%.3f')                  
        
        # Here we store data
        os.chdir(self._cwd)
        os.chdir(self._StratigraphyTempFilesPath)
        for i in activeLinks:
            filename = 'link_' + str(i) + '.txt'
            data = np.hstack((self.t,z[i],dzl[i], La[i],F[i,:],Fs[i,:]))
            data = np.reshape(data,[1,data.shape[0]])
            with open(filename, 'ab') as f:
                np.savetxt(f, data,'%.3f') 
        
        # Here we update stratigraphy in case of deposition - only for new layers
        os.chdir(self._cwd)
        os.chdir(self._StratigraphyTempFilesPath)
        if self._updateSubSurface is True:
            print('Updating stratigraphy at time ',self.t, 's. New layer(s) created at links ',self._deepLinksId)
            for i in self._deepLinksId:
                linkData = np.loadtxt('link_' + str(i) + '.txt') 
                meanSubsurfaceGSD = np.mean(linkData[:,4:4+D.shape[0]-1],axis = 0) # 4 is the number of elements before F GSD
                self._grid['link']['bed_subsurface__grainSizeDistribution'][i,:] = \
                    meanSubsurfaceGSD/np.sum(meanSubsurfaceGSD) 
                os.remove('link_' + str(i) + '.txt')
                
            # Now updates the raw data for the links with deposition
            os.chdir(self._cwd)
            os.chdir(self._StratigraphyRawDataPath)
                    
            for i in self._deepLinksId:
                filename = 'link_' + str(i) + '.txt'
                data = np.hstack((self.t,z[i],F[i,:],self._grid['link']['bed_subsurface__grainSizeDistribution'][i,:]))
                data = np.reshape(data,[1,data.shape[0]])
                with open(filename, 'ab') as f:
                    np.savetxt(f, data,'%.3f') 
                       
            self._grid['link']['bed_surface__surface_thickness_new_layer'][self._deepLinksId] = \
                self._grid['link']['bed_surface__surface_thickness_new_layer'][self._deepLinksId] - self._newSurfaceLayerThickness
            self._grid['link']['topographic__elevation_subsurface'][self._deepLinksId] =\
                self._grid['link']['topographic__elevation_subsurface'][self._deepLinksId] + self._newSurfaceLayerThickness
            
            self._subsurfaceChanged = True
        
        # Here we update the stratigraphy in case of erosion - only for new layers
        os.chdir(self._cwd)
        
        if self._updateErodedSubSurface is True:
            print('Updating stratigraphy at time ',self.t, 's. Layer(s) eroded')
            os.chdir(self._StratigraphyTempFilesPath)
            for i in self._erodedLinksId:
                os.remove('link_' + str(i) + '.txt')            
            
            os.chdir(self._cwd)
            os.chdir(self._StratigraphyRawDataPath)
            
            for i in self._erodedLinksId:
                linkData = np.loadtxt('link_' + str(i) + '.txt')
                if linkData.shape[0] > 1:
                    self._grid['link']['bed_subsurface__grainSizeDistribution'][i,:] = \
                        linkData[linkData.shape[0]-2,9:None] # 9 is the number of elements before Fs GSD
                    os.remove('link_' + str(i) + '.txt')
                    linkData[-2,0:2] = np.array((self.t,z[i])) # Removes last layer
                    with open('link_' + str(i) + '.txt', 'ab') as f:
                        np.savetxt(f, linkData[0:-1,:],'%.3f') 
                    self._grid['link']['bed_surface__grainSizeDistribution'][i,:] = \
                        self._grid['link']['bed_subsurface__grainSizeDistribution'][i,:]            
        
            self._grid['link']['bed_surface__surface_thickness_new_layer'][self._erodedLinksId] = \
                self._grid['link']['bed_surface__surface_thickness_new_layer'][self._erodedLinksId] + self._newSurfaceLayerThickness
            self._grid['link']['topographic__elevation_subsurface'][self._erodedLinksId] =\
                self._grid['link']['topographic__elevation_subsurface'][self._erodedLinksId] - self._newSurfaceLayerThickness
            
            self._subsurfaceChanged = True
            
        os.chdir(self._cwd)
        self._stratigraphyCycle = 0
        self._computeStratigraphy = False
        self._updateSubSurface = False
        self._updateErodedSubSurface = False
        
    def bedload_mapper(self):
        
        """Maps the magnitude of the volumetric bedload per unit width values 
        from links (m2/s) onto nodes (m3/s).

        This method takes the bedload transport rates from links and calculates 
        the magnitude at a given node. 
             
        """
        qb_x_r = self._grid['link']['sediment_transport__bedloadRate'][self._grid.links_at_node[:,0]] * self._grid.dy
        (I,) = np.where(self._grid.links_at_node[:,0] < 0)
        qb_x_r[I] = 0  
        
        qb_x_l = self._grid['link']['sediment_transport__bedloadRate'][self._grid.links_at_node[:,2]] * self._grid.dy
        (I,) = np.where(self._grid.links_at_node[:,2] < 0)
        qb_x_l[I] = 0  
        
        qb_y_u = self._grid['link']['sediment_transport__bedloadRate'][self._grid.links_at_node[:,1]] * self._grid.dx
        (I,) = np.where(self._grid.links_at_node[:,1] < 0)
        qb_y_u[I] = 0  
        
        qb_y_l = self._grid['link']['sediment_transport__bedloadRate'][self._grid.links_at_node[:,3]] * self._grid.dx
        (I,) = np.where(self._grid.links_at_node[:,3] < 0)
        qb_y_l[I] = 0  
        
        qb_x = 0.5*(qb_x_r + qb_x_l)
        qb_y = 0.5*(qb_y_u + qb_y_l)
        
        self._bedloadRate_vector = np.transpose(np.vstack((qb_x,qb_y)))
        self._bedloadRate_magnitude = np.sqrt(qb_x**2 + qb_y**2)      
            
    def vector_mapper(self,vector):
        
        """Maps a vector, in this case shear stress or velocity, values from 
        links onto nodes preserving the components.

        This method takes the vectors values on links and determines the 
        vectors components in nodes.       
        
        Examples
        --------
        
        If we want to plot the velocity vector on top of the surface water 
        depth we can do (after running the main example)
        >>> from matplotlib import pyplot as plt 
        >>> (velocityVector,velocityMagnitude) = RBD.vector_mapper(RBD._u)
        >>> velocityVector_x = velocityVector[:,0]
        >>> velocityVector_x             
        array([ 0.125,  0.25 ,  0.25 ,  0.25 ,  0.125,  0.125,  0.25 ,  0.25 ,
        0.25 ,  0.125,  0.125,  0.25 ,  0.25 ,  0.25 ,  0.125,  0.125,
        0.25 ,  0.25 ,  0.25 ,  0.125,  0.125,  0.25 ,  0.25 ,  0.25 ,
        0.125])        
        >>> velocityVector_y = velocityVector[:,1]
        >>> velocityVector_y
        array([ 0.125,  0.125,  0.125,  0.125,  0.125,  0.25 ,  0.25 ,  0.25 ,
        0.25 ,  0.25 ,  0.25 ,  0.25 ,  0.25 ,  0.25 ,  0.25 ,  0.25 ,
        0.25 ,  0.25 ,  0.25 ,  0.25 ,  0.125,  0.125,  0.125,  0.125,
        0.125])
        >>> imshow_grid(grid, 'surface_water__depth',cmap='Blues',vmin=0,vmax=0.5)
        >>> plt.quiver(grid.x_of_node,grid.y_of_node,velocityVector_x,velocityVector_y,scale = 10)
              
        """        
        vector_x_r = vector[self._grid.links_at_node[:,0]]
        (I,) = np.where( self._grid.links_at_node[:,0] < 0)
        vector_x_r[I] = 0  
        
        vector_x_l = vector[self._grid.links_at_node[:,2]]
        (I,) = np.where( self._grid.links_at_node[:,2] < 0)
        vector_x_l[I] = 0  
        
        vector_y_u = vector[self._grid.links_at_node[:,1]]
        (I,) = np.where( self._grid.links_at_node[:,1] < 0)
        vector_y_u[I] = 0  
        
        vector_y_l = vector[self._grid.links_at_node[:,3]]
        (I,) = np.where( self._grid.links_at_node[:,3] < 0)
        vector_y_l[I] = 0  
        
        vector_x = 0.5*(vector_x_r + vector_x_l)
        vector_y = 0.5*(vector_y_u + vector_y_l)
        
        return np.transpose(np.vstack((vector_x,vector_y))), \
            np.sqrt(vector_x**2 + vector_y**2)
            
    def calculateDX(self, fX, mapped_in = 'node'):
        """ Calculates the bed surface and bed load grain size corresponding to
        any fraction. For example, 50%, which is the D50 or median grain size. 
        In that case use fX = 0.5
        
        This method takes the user specified fraction, from 0 to 1, and outputs
        the corresponding grain size in nodes or links. By default the node 
        option is used. Use mapped_in = 'link' to calculate in links

        Examples
        --------
        
        After running the main example, run
        >>> RBD.calculateDX(0.9)
        (array([ 29.17511963,  27.85761803,  27.85761803,  27.85761803,
         29.17511963,  29.17511963,  27.85761803,  27.85761803,
         27.85761803,  29.17511963,  29.17511963,  27.85761803,
         27.85761803,  27.85761803,  29.17511963,  29.17511963,
         27.85761803,  27.85761803,  27.85761803,  29.17511963,
         29.17511963,  27.85761803,  27.85761803,  27.85761803,  29.17511963]),
         array([ 27.04869462,  26.27655305,  97.00586026,  26.27655305,
         27.04869462,  27.59403443,  26.07884242,  97.00586026,
         26.07884242,  27.59403443,  27.59403443,  26.07884242,
         97.00586026,  26.07884242,  27.59403443,  27.59403443,
         26.07884242,  97.00586026,  26.07884242,  27.59403443,
         27.04869462,  25.55545386,  26.36216338,  25.55545386,  27.04869462]))
        
        >>> RBD.calculateDX(0.9,mapped_in = 'link')
        (array([ 28.64080227,  27.85761803,  27.85761803,  28.64080227,
         29.17511963,  27.85761803,  27.85761803,  27.85761803,
         29.17511963,  28.64080227,  27.85761803,  27.85761803,
         28.64080227,  29.17511963,  27.85761803,  27.85761803,
         ...        ,  ...        ,  ...        ,  ...        ,  
         28.64080227,  27.85761803,  27.85761803,  28.64080227]),
         array([ 26.85294086,   0.        ,   0.        ,  26.85294086,
         28.14885923,  25.03215792,  25.03215792,  25.03215792,
         28.14885923,  26.85294086,   0.        ,   0.        ,
         ...        ,  ...        ,  ...        ,  ...        ,  
          0.        ,   0.        ,  26.85294086,  28.14885923,
         25.03215792,  27.2670109 ,  25.03215792,  28.14885923,
         26.85294086,  25.03215792,  25.03215792,  26.85294086]))
        """
        nNodes = self._grid.number_of_nodes
        nodesList = np.arange(0,nNodes)
        nLinks = self._grid.number_of_links
        linksList = np.arange(0,nLinks)
        
        D = self._DOrig                     # Grain sizes
        Psi_D = np.flip(np.log2(D),axis=0)  # Grain sizes in Psi scale
        
        # First for the bed surface
        
        # Equivalent grain sizes frequency
        if mapped_in == 'link':
            f_DEq = self._grid['link']['bed_surface__grainSizeDistribution']
            f_list = linksList
        else:
            f_DEq = self._grid['node']['bed_surface__grainSizeDistribution']
            f_list = nodesList

        f = np.hstack((f_DEq,np.zeros([f_DEq.shape[0],1]))) # Grain sizes freq
        
        f_D = np.cumsum(np.flip(f,axis=1),axis=1) # Cumulative GSD in each node
        # Finds the index of the grain size smaller (i0) and larger (i1) than 
        # the fraction requested (fX) in each node
        i0 = np.argmin(f_D<=fX,axis=1) - 1  
        i1 = np.argmax(f_D>fX,axis=1)      
                
        Psi_DX = Psi_D[i0] + \
            ( (Psi_D[i1]-Psi_D[i0]) / (f_D[f_list,i1]-f_D[f_list,i0]) ) \
                * (fX - f_D[f_list,i0])
        DX_surface = 2**Psi_DX

        # Now for the bed load
        
        # Equivalent grain sizes frequency
        
        if mapped_in == 'link':
            f_DEq = self._grid['link']['sediment_transport__bedloadGSD']   
        else:
            f_DEq = self._grid['node']['sediment_transport__bedloadGSD'] 
        
        f = np.hstack((f_DEq,np.zeros([f_DEq.shape[0],1]))) # Grain sizes freq
        
        f_D = np.cumsum(np.flip(f,axis=1),axis=1) # Cumulative GSD in each node
        # Finds the index of the grain size smaller (i0) and larger (i1) than 
        # the fraction requested (fX) in each node
        i0 = np.argmin(f_D<=fX,axis=1) - 1  
        i1 = np.argmax(f_D>fX,axis=1)      
        
        Psi_DX = Psi_D[i0] + \
            ( (Psi_D[i1]-Psi_D[i0]) / (f_D[f_list,i1]-f_D[f_list,i0]) ) \
                * (fX - f_D[f_list,i0])
        DX_bedload = 2**Psi_DX
          
        return DX_surface, DX_bedload
        
    def formatGSDOutputs(self,bedloadGSD):
        """ Gives a more friendly format for the bed surface or bed load GSD. 
        Reads a bed load GSD, from links or nodes, and returns the GSD in
        cumulative percetage
        
        Examples
        --------
        
        After running the main example, run
        >>> RBD.formatGSDOutputs(grid['link']['sediment_transport__bedloadGSD'])
                   8         16     32
        Link_0   0.0  60.473028  100.0
        Link_1   0.0   0.000000    0.0
        Link_2   0.0   0.000000    0.0
        Link_3   0.0  60.473028  100.0
        Link_4   0.0  45.944616  100.0
        ...      ...  ...        ...
        Link_35  0.0  45.944616  100.0
        Link_36  0.0  60.473028  100.0
        Link_37  0.0  71.774474  100.0
        Link_38  0.0  71.774474  100.0
        Link_39  0.0  60.473028  100.0

        It also works for nodes
        >>> RBD.formatGSDOutputs(grid['node']['sediment_transport__bedloadGSD'])
                   8         16     32
        Node_0   0.0  58.765086  100.0
        Node_1   0.0  64.825377  100.0
        Node_2   0.0  34.588863   50.0
        Node_3   0.0  64.825377  100.0
        Node_4   0.0  58.765086  100.0
        ...      ...  ...        ...
        Node_20  0.0  58.765086  100.0
        Node_21  0.0  69.177726  100.0
        Node_22  0.0  64.235022  100.0
        Node_23  0.0  69.177726  100.0
        Node_24  0.0  58.765086  100.0

        """
        
        if bedloadGSD.shape[0] == self._grid.number_of_links:
            indexText = 'Link_'
        else:
            indexText = 'Node_'
        
        bedloadGSD = np.hstack((bedloadGSD,np.zeros([bedloadGSD.shape[0],1])))
        bedloadGSD = np.cumsum(np.fliplr(bedloadGSD),axis=1)*100
        
        D_ascending = np.sort(self._DOrig)
        columns = ["".join(item) for item in D_ascending.astype(str)]
        
        df = pd.DataFrame(bedloadGSD,columns=columns,
                           index = [indexText + str(i) 
                           for i in range(bedloadGSD.shape[0])])
        return df