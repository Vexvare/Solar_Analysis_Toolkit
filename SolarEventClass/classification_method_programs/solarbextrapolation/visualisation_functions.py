# -*- coding: utf-8 -*-
"""
Created on Mon Jul 27 14:39:00 2015

@author: alex_
"""


# Universal Imports
import numpy as np
import sunpy.map as mp
from astropy import units as u

# Visulisation Imports
from mayavi import mlab
from .mayavi_seed_streamlines import SeedStreamline, Streamline
from mayavi.tools.sources import vector_field
import matplotlib
matplotlib.use('Qt5Agg') 
import matplotlib.pyplot as plt
# Module Imports
#from classes import *
#from solarbextrapolation.map3dclasses import Map3D
from classification_method_programs.solarbextrapolation.utilities import decompose_ang_len

def visualise(aMap3D, **kwargs):
    """
    Basic function for visualising a vector field from an extrapolator.
    General usage involves passing boundary map and volume vector field and
    these are then aligned and plotted in mayavi.
    The vector field will be represented by streamlines generated from the
    given (or otherwise default) seed points.
    The boundary data should be rendered in approbriate colours for the given
    map data.

    Parameters
    ----------

    aMap3D : Map3D
        The 3D vector field from the extrapolator.

    boo_debug : boolean, optional
        If set, turns on logging functionality.

    seeds : numpy.array, optional
        If set, provides a list of manual seed points in the 3D vector field.

    boundary : sunpy.map, optional
        If set, provides the 2D map to place in the visulisation at the base of
    the volume.

    unit_length : `astropy.units.quantity.Quantity`, optional
        If set, provides the length of one unit in MayaVi for scaling maps.

    boundary_unit : `astropy.units.quantity.Quantity`, optional
        If set, provides a single unit for the x/y-axes of the boundary map.

    boundary_units : list, optional
        If set, provides a list of units for the x/y-axes of the boundary map.

    volume_unit : `astropy.units.quantity.Quantity`, optional
        If set, provides a single unit for the x/y/z-axes of the 3D vector field.

    volume_units : list, optional
        If set, provides a list of units for the x/y/z-axes of the 3D vector field.

    show_boundary_axes : boolean, optional
        If set, enables the display of the boundary map axes.

    show_volume_axes : boolean, optional
        If set, enables the display of the 3D vector field axes.

    """
    print('\n\nIn visualisation_functions... attempting to visualize the data... ')
    # Optional parameters
    boo_debug          = kwargs.get('debug', False)
    np_seeds           = kwargs.get('seeds', None)
    boundary           = kwargs.get('boundary', None)
    mayavi_unit_length = kwargs.get('unit_length', 1.0 * u.Mm) * 1.0
    boundary_unit      = kwargs.get('boundary_unit', mayavi_unit_length) * 1.0
    boundary_units     = kwargs.get('boundary_units', [ boundary_unit, boundary_unit, boundary_unit ])
    volume_unit        = kwargs.get('volume_unit', mayavi_unit_length) * 1.0
    volume_units       = kwargs.get('volume_units', [ volume_unit, volume_unit, volume_unit ])
    show_boundary_axes = kwargs.get('show_boundary_axes', True)
    show_volume_axes   = kwargs.get('show_volume_axes', True)

    # Setup the arc to length equivilence
    obs_distance = aMap3D.dsun - aMap3D.rsun_meters
    radian_length = [ (u.radian, u.meter, lambda x: obs_distance * x, lambda x: x / obs_distance) ]

    # Slice (scale) the fields to make the vectors usable in mayavi.
    int_slice_scale = 1

    npm_3d_sliced   = aMap3D.data[::int_slice_scale,::int_slice_scale,::int_slice_scale,:]

    # Plot the main vector field (volume).
    fig = mlab.figure(bgcolor = (1,1,1), fgcolor = (0,0,0)) #Just creates the scene, we will add the information to the scene next
    
    # Make 3D coords for ever point in the 3D grid.
    x_range = u.Quantity([ decompose_ang_len(aMap3D.xobsrange[0], equivalencies=radian_length),
                           decompose_ang_len(aMap3D.xobsrange[1], equivalencies=radian_length) ])
    y_range = u.Quantity([ decompose_ang_len(aMap3D.yobsrange[0], equivalencies=radian_length),
                           decompose_ang_len(aMap3D.yobsrange[1], equivalencies=radian_length) ])
    z_range = u.Quantity([ decompose_ang_len(aMap3D.zrange[0], equivalencies=radian_length),
                           decompose_ang_len(aMap3D.zrange[1], equivalencies=radian_length) ])
    
    x_range_scaled = (x_range/mayavi_unit_length).decompose().value
    y_range_scaled = (y_range/mayavi_unit_length).decompose().value
    z_range_scaled = (z_range/mayavi_unit_length).decompose().value

    X, Y, Z = np.mgrid[y_range_scaled[0]:y_range_scaled[1]:npm_3d_sliced.shape[0]*1j,    #Basically creates a empty box, or grid, of all of the spatial points that we will
                       x_range_scaled[0]:x_range_scaled[1]:npm_3d_sliced.shape[1]*1j,    # use to put all the magnetic field tube data into
                       z_range_scaled[0]:z_range_scaled[1]:npm_3d_sliced.shape[2]*1j]    # along with the base of the HMI map

    vec_field = vector_field(X, Y, Z, npm_3d_sliced[:,:,:,0], npm_3d_sliced[:,:,:,1], npm_3d_sliced[:,:,:,2],
                             name='Magnetic Vector Field', figure=fig)

    vec_field_mag = mlab.pipeline.extract_vector_norm(vec_field, name="Magnetic Field Magnitude")
    
    if show_volume_axes:
        # Label axes
        axes = mlab.axes()
        x_range_axis = decompose_ang_len((x_range/volume_units[0]).decompose(), equivalencies=radian_length, working_units=volume_units[0])#(x_range/volume_units[0]).decompose()
        y_range_axis = decompose_ang_len((y_range/volume_units[1]).decompose(), equivalencies=radian_length, working_units=volume_units[1])#(y_range/volume_units[1]).decompose()
        z_range_axis = decompose_ang_len((z_range/volume_units[2]).decompose(), equivalencies=radian_length, working_units=volume_units[2])#(z_range/volume_units[2]).decompose()
        '''
        print('\n')
        print('x_range: ' + str(x_range))
        print('y_range: ' + str(y_range))
        print('z_range: ' + str(z_range))
        print('\n')
        print('x_range_axis: ' + str(x_range_axis))
        print('y_range_axis: ' + str(y_range_axis))
        print('z_range_axis: ' + str(z_range_axis))
        print('\n')
        print('x_range_axis[0]: ' + str(x_range_axis[0]))
        print('y_range_axis[0]: ' + str(y_range_axis[0]))
        print('z_range_axis[0]: ' + str(z_range_axis[0]))
        print('\n')
        print('\nx_range_axis[0].value: ' + str(x_range_axis[0].value))
        print('x_range_axis[0].unit: ' + str(x_range_axis[0].unit))
        print('type(x_range_axis[0].unit): ' + str(type(x_range_axis[0].unit)))
        print('\n')
        '''
        #ROTATED DATA : 
        axes.axes.ranges = np.array([ y_range_axis[1], y_range_axis[0], x_range_axis[0], x_range_axis[1],  z_range_axis[0],  z_range_axis[1]])
        #NON - ROTATED DATA :
        #axes.axes.ranges = np.array([ x_range_axis[0], x_range_axis[1], y_range_axis[0], y_range_axis[1], z_range_axis[0],  z_range_axis[1]])
        
        axes.axes.use_ranges = True
        axes.axes.x_label = 'Solar Y (' + unit_label(volume_units[1]) + ')'
        axes.axes.y_label = 'Solar X (' + unit_label(volume_units[1]) + ')'
        axes.axes.z_label = 'Z (' + unit_label(volume_units[2]) + ')'
        
        axes.axes.font_factor = 0.80
        axes.axes.number_of_labels = 6
        
        axes.property.display_location = 'background'        
        
        'Adding Label Information and settings'
        axes.label_text_property.font_family = 'courier'
        axes.label_text_property.bold = False
        axes.label_text_property.italic = False
        axes.label_text_property.shadow = True
        axes.label_text_property.vertical_justification = 'bottom'
        

        'Adding Title Information and settings'
        axes.title_text_property.font_family = 'courier'
        axes.title_text_property.shadow = True
        axes.title_text_property.justification = 'left'
        axes.title_text_property.vertical_justification = 'bottom'
        axes.title_text_property.font_size = 2
        axes.title_text_property.bold = False
        
    # Plot the seed points
    if np_seeds is None:
        # Generate a plane for the streamline seed points, this plane is used to show/find the field lines that you wish to see in the program. 
        # It is a movable object that will find the exact field lines that will run through it
                                                            
        streamline = Streamline(streamline_type = 'tube')
        vec_field_mag.add_child(streamline)
        streamline.stream_tracer.integration_direction = 'both'
        streamline.seed.widget = streamline.seed.widget_list[2]
        streamline.seed.widget.resolution = 10
        streamline.seed.widget.representation = 'outline'
        streamline.tube_filter.radius = 0.5
        streamline.seed.widget.handle_size = 0.01
        # Some necessary points within the volume
        z = (0.15 * (z_range_scaled[1] - z_range_scaled[0])) + z_range_scaled[0]
        x_mid = (x_range_scaled[0] + x_range_scaled[1])/2.0
        y_mid = (y_range_scaled[0] + y_range_scaled[1])/2.0
        #rotated points
        streamline.seed.widget.normal_to_z_axis = True
        streamline.seed.widget.center = np.array([ y_mid,  x_mid,  z])
        streamline.seed.widget.point1 = np.array([ y_range_scaled[0],x_range_scaled[1],  z])
        streamline.seed.widget.point2 = np.array([ y_range_scaled[1],x_range_scaled[0],  z])
        streamline.seed.widget.origin = np.array([ y_range_scaled[0],x_range_scaled[0],  z])
        # Update the render
        scene = fig.scene
        scene.render()
    else:
        points = mlab.points3d(np_seeds[:,0], np_seeds[:,1], np_seeds[:,2])
        # Make the points smaller
        points.glyph.glyph.scale_factor = 10.0 #mayavi_scale
        # Make the points blue
        points.actor.property.color = (0.2,0,1)
        # Create the custom streamline object
        streamline = SeedStreamline(seed_points=np_seeds)
        # Add the streamline object to the plot and make it use the magentic field data,
        # by adding it as a child of the field we created earlier.
        # We add it to the magnitude field (which is in itself a child of bfield)
        # so that it picks up the scalar values and colours the lines.
        vec_field_mag.add_child(streamline)
    # Adjust some of the streamline appearance parameters
    streamline.module_manager.scalar_lut_manager.lut_mode = 'brg'
    streamline.module_manager.scalar_lut_manager.reverse_lut = True
    streamline.stream_tracer.integration_direction = 'both'
    streamline.stream_tracer.maximum_propagation = 500.0
    streamline.module_manager.scalar_lut_manager.show_legend = True #This legend is needed because we wish to see the strength of the magnetic fields at specific points
    streamline.module_manager.scalar_lut_manager.scalar_bar_representation.position = np.array([ 0.1,  0.1 ])
    streamline.update_pipeline()

    # Add the boundary data 2D map
    if boundary:

        x_range = u.Quantity([ decompose_ang_len(aMap3D.xrange[0], equivalencies=radian_length),
                               decompose_ang_len(aMap3D.xrange[1], equivalencies=radian_length) ])
        y_range = u.Quantity([ decompose_ang_len(aMap3D.yrange[0], equivalencies=radian_length),
                               decompose_ang_len(aMap3D.yrange[1], equivalencies=radian_length) ])
        
        x_range_scaled = (x_range/mayavi_unit_length).decompose().value
        y_range_scaled = (y_range/mayavi_unit_length).decompose().value
        
        # Create explicit points in 3D space
        X, Y = np.mgrid[y_range_scaled[0]:y_range_scaled[1]:boundary.data.shape[0]*1j,
                        x_range_scaled[0]:x_range_scaled[1]:boundary.data.shape[1]*1j]

        #print(boundary.__init__)
        # ONLY y ROTATION : 
        'I need to flip Y data into a new array of shape x=boundary.data.shape[1], y=boundary.data.shape[0]'
        temp_2d_boundary = np.zeros((boundary.data.shape[0],boundary.data.shape[1]))
        x=boundary.data.shape[1]-1
        y=boundary.data.shape[0]-1
        while x > 0 :
            while y > 0 :
                temp_2d_boundary[boundary.data.shape[0]-y,x] = boundary.data[y,x]
                y=y-1
            y=boundary.data.shape[0]-1
            x=x-1
        # Plot the current figure and add the boundary
        img_boundary = mlab.pipeline.array2d_source(X, Y, 
                                                    temp_2d_boundary, 
                                                    figure=fig
                                                    )

        img_boundary = mlab.pipeline.image_actor(img_boundary, figure = fig)
        img_boundary.actor.force_opaque = True
        img_boundary.module_manager.scalar_lut_manager.lut_mode = 'gist_earth'
        
        # Legend details
        img_boundary.module_manager.scalar_lut_manager.show_legend = False # This legend isn't needed, we dont need to view the 'Intensity' of each point; especially if we're gonna map a boundary around the coronal hole
        #img_boundary.module_manager.scalar_lut_manager.scalar_bar_representation.position = np.array([ 0.1,  0.1 ])

        # Place a small outline around the data cube
        outline = mlab.outline(line_width = 3)
        outline.outline_mode = 'cornered'
        
        # Show the axes if selected
        if show_boundary_axes:
            axes = mlab.axes()

            # Get the ranges of the boundary and scale to the selected units
            x_range = boundary.xrange.to(boundary_units[0].unit, equivalencies=radian_length)
            y_range = boundary.yrange.to(boundary_units[1].unit, equivalencies=radian_length)
            x_range_scaled = (x_range/boundary_units[0]).decompose().value
            y_range_scaled = (y_range/boundary_units[1]).decompose().value

            # Update the ranges manually to use custom units for the boundary
            axes.axes.ranges = np.array([ x_range_scaled[0],  x_range_scaled[1],  y_range_scaled[0],  y_range_scaled[1],  0,  0])
            axes.axes.use_ranges = True
            axes.axes.x_label = 'Solar X (' + unit_label(boundary_units[0]) + ')'
            axes.axes.y_label = 'Solar Y (' + unit_label(boundary_units[1]) + ')'
    return fig


def unit_label(quantity):
    """
    Small function to return a string label that is empty if value is 1.0 or is
    the given number otherwise.
    """
    if quantity.value == 1.0:
        return str(quantity.unit)
    return str(quantity)

