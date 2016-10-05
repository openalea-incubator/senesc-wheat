# -*- coding: latin-1 -*-

from __future__ import division # use "//" to do integer division

"""
    senescwheat.converter
    ~~~~~~~~~~~~~~~~~~~~~~~

    The module :mod:`senescwheat.converter` defines functions to:

        * convert :class:`dataframes <pandas.DataFrame>` to/from SenescWheat inputs or outputs format.
        * convert a :class:`MTG <openalea.mtg.mtg.MTG>` to/from SenescWheat inputs or outputs format.

    :copyright: Copyright 2014-2015 INRA-ECOSYS, see AUTHORS.
    :license: TODO, see LICENSE for details.

    .. seealso:: Barillot et al. 2015.
"""

"""
    Information about this versioned file:
        $LastChangedBy$
        $LastChangedDate$
        $LastChangedRevision$
        $URL$
        $Id$
"""

import numpy as np
import pandas as pd

import simulation


#: the name of the photosynthetic organs modeled by SenescWheat
SENESCWHEAT_PHOTOSYNTHETIC_ORGANS_NAMES = set(['internode', 'blade', 'sheath', 'peduncle', 'ear'])

#: the inputs needed by SenescWheat at roots scale
SENESCWHEAT_ROOTS_INPUTS = ['sucrose', 'amino_acids', 'mstruct', 'Nstruct', 'cytokinins']

#: the inputs needed by SenescWheat at element scale
SENESCWHEAT_ELEMENTS_INPUTS = ['green_area', 'proteins', 'mstruct', 'max_proteins', 'Nstruct', 'nitrates', 'amino_acids', 'starch', 'fructan', 'cytokinins', 'sucrose', 'is_growing']

#: the inputs needed by SenescWheat
SENESCWHEAT_INPUTS = SENESCWHEAT_ROOTS_INPUTS + SENESCWHEAT_ELEMENTS_INPUTS

#: the outputs computed by SenescWheat at roots scale
SENESCWHEAT_ROOTS_OUTPUTS = ['mstruct_death', 'mstruct', 'Nstruct', 'cytokinins']

#: the outputs computed by SenescWheat at elements scale
SENESCWHEAT_ELEMENTS_OUTPUTS = ['green_area', 'mstruct', 'Nstruct', 'starch', 'sucrose', 'fructan', 'proteins', 'amino_acids', 'cytokinins']

#: the outputs computed by SenescWheat
SENESCWHEAT_OUTPUTS = SENESCWHEAT_ROOTS_OUTPUTS + SENESCWHEAT_ELEMENTS_OUTPUTS

#: the inputs and outputs of SenescWheat at roots scale
SENESCWHEAT_ROOTS_INPUTS_OUTPUTS = sorted(list(set(SENESCWHEAT_ROOTS_INPUTS + SENESCWHEAT_ROOTS_OUTPUTS)))

#: the inputs and outputs of SenescWheat at elements scale
SENESCWHEAT_ELEMENTS_INPUTS_OUTPUTS = sorted(list(set(SENESCWHEAT_ELEMENTS_INPUTS + SENESCWHEAT_ELEMENTS_OUTPUTS)))

#: the inputs and outputs of SenescWheat.
SENESCWHEAT_INPUTS_OUTPUTS = SENESCWHEAT_INPUTS + SENESCWHEAT_OUTPUTS

#: the columns which define the topology of a roots in the input/output dataframe
ROOTS_TOPOLOGY_COLUMNS = ['plant', 'axis']

#: the columns which define the topology of an element in the input/output dataframe
ELEMENTS_TOPOLOGY_COLUMNS = ['plant', 'axis', 'metamer', 'organ', 'element']


def from_dataframes(roots_inputs, elements_inputs):
    """
    Convert inputs/outputs from Pandas dataframes to Senesc-Wheat format.

    :Parameters:

        - `roots_inputs` (:class:`pandas.DataFrame`) - Roots inputs dataframe to convert, with one line by roots.

        - `elements_inputs` (:class:`pandas.DataFrame`) - Elements inputs dataframe to convert, with one line by element.

    :Returns:
        The inputs/outputs in a dictionary.

    :Returns Type:
        :class:`dict` of :class:`dict`

    .. seealso:: see :attr:`simulation.Simulation.inputs` and :attr:`simulation.Simulation.outputs`
       for the structure of Senesc-Wheat inputs/outputs.

    """
    all_roots_dict = {}
    all_elements_dict = {}
    for (all_current_dict, current_dataframe, current_topology_columns) in ((all_roots_dict, roots_inputs, ROOTS_TOPOLOGY_COLUMNS),
                                                                            (all_elements_dict, elements_inputs, ELEMENTS_TOPOLOGY_COLUMNS)):
        current_columns = current_dataframe.columns.difference(current_topology_columns)
        for current_id, current_group in current_dataframe.groupby(current_topology_columns):
            current_series = current_group.loc[current_group.first_valid_index()]
            current_dict = current_series[current_columns].to_dict()
            all_current_dict[current_id] = current_dict

    return {'roots': all_roots_dict, 'elements': all_elements_dict}


def to_dataframes(data_dict):
    """
    Convert inputs/outputs from Senesc-Wheat format to Pandas dataframe.

    :Parameters:

        - `data_dict` (:class:`dict`) - The inputs/outputs in Senesc-Wheat format.

    :Returns:
        One dataframe for roots inputs/outputs and one dataframe for elements inputs/outputs.

    :Returns Type:
        :class:`tuple` of :class:`pandas.DataFrame`

    .. seealso:: see :attr:`simulation.Simulation.inputs` and :attr:`simulation.Simulation.outputs`
       for the structure of Senesc-Wheat inputs/outputs.

    """
    dataframes_dict = {}
    for (current_key, current_topology_columns, current_inputs_outputs_names) in (('roots', ROOTS_TOPOLOGY_COLUMNS, SENESCWHEAT_ROOTS_INPUTS_OUTPUTS),
                                                                                  ('elements', ELEMENTS_TOPOLOGY_COLUMNS, SENESCWHEAT_ELEMENTS_INPUTS_OUTPUTS)):
        current_data_dict = data_dict[current_key]
        current_ids_df = pd.DataFrame(current_data_dict.keys(), columns=current_topology_columns)
        current_data_df = pd.DataFrame(current_data_dict.values())
        current_df = pd.concat([current_ids_df, current_data_df], axis=1)
        current_df.sort_values(by=current_topology_columns, inplace=True)
        current_columns_sorted = current_topology_columns + [input_output for input_output in current_inputs_outputs_names if input_output in current_df.columns]
        current_df = current_df.reindex_axis(current_columns_sorted, axis=1, copy=False)
        current_df.reset_index(drop=True, inplace=True)
        dataframes_dict[current_key] = current_df

    return dataframes_dict['roots'], dataframes_dict['elements']


def from_MTG(g, roots_inputs, elements_inputs):
    """
    Convert a MTG to Senesc-Wheat inputs.
    Use data in `roots_inputs` and `elements_inputs` if `g` is incomplete.

    :Parameters:

            - g (:class:`openalea.mtg.mtg.MTG`) - A MTG which contains the inputs
              of Senesc-Wheat. These inputs are: :mod:`SENESCWHEAT_INPUTS`.

            - `roots_inputs` (:class:`pandas.DataFrame`) - Roots dataframe, with one line by roots.

            - `elements_inputs` (:class:`pandas.DataFrame`) - Elements dataframe, with one line by element.

    :Returns:
        The inputs of Senesc-Wheat.

    :Returns Type:
        :class:`dict` of :class:`dict`

    .. seealso:: see :attr:`simulation.Simulation.inputs` for the structure of Senesc-Wheat inputs.

    """
    all_roots_inputs_dict = {}
    all_elements_inputs_dict = {}

    roots_inputs_grouped = roots_inputs.groupby(ROOTS_TOPOLOGY_COLUMNS)
    elements_inputs_grouped = elements_inputs.groupby(ELEMENTS_TOPOLOGY_COLUMNS)

    # traverse the MTG recursively from top ...
    for plant_vid in g.components_iter(g.root):
        plant_index = int(g.index(plant_vid))
        for axis_vid in g.components_iter(plant_vid):
            axis_label = g.label(axis_vid)
            roots_id = (plant_index, axis_label)
            if roots_id in roots_inputs_grouped.groups:
                roots_inputs_group = roots_inputs_grouped.get_group(roots_id)
                roots_inputs_group_series = roots_inputs_group.loc[roots_inputs_group.first_valid_index(), SENESCWHEAT_ROOTS_INPUTS]
            else:
                roots_inputs_group_series = pd.Series()
            axis_properties = g.get_vertex_property(axis_vid)
            if 'roots' in axis_properties:
                roots_properties = axis_properties['roots']
                roots_inputs_dict = {}
                is_valid_roots = True
                for roots_input_name in SENESCWHEAT_ROOTS_INPUTS:
                    if roots_input_name in roots_properties:
                        # use the input from the MTG
                        roots_inputs_dict[roots_input_name] = roots_properties[roots_input_name]
                    else:
                        # use the input from the dataframe
                        if roots_input_name in roots_inputs_group_series:
                            roots_inputs_dict[roots_input_name] = roots_inputs_group_series[roots_input_name]
                        else:
                            is_valid_roots = False
                            break
                if is_valid_roots:
                    all_roots_inputs_dict[roots_id] = roots_inputs_dict
            else:
                roots_inputs_group_dict = roots_inputs_group_series.to_dict()
                if set(roots_inputs_group_dict).issuperset(SENESCWHEAT_ROOTS_INPUTS):
                    all_roots_inputs_dict[roots_id] = roots_inputs_group_dict

            for metamer_vid in g.components_iter(axis_vid):
                metamer_index = int(g.index(metamer_vid))
                for organ_vid in g.components_iter(metamer_vid):
                    organ_label = g.label(organ_vid)
                    if organ_label not in SENESCWHEAT_PHOTOSYNTHETIC_ORGANS_NAMES: continue
                    for element_vid in g.components_iter(organ_vid):
                        vertex_properties = g.get_vertex_property(element_vid)
                        element_label = g.label(element_vid)
                        element_id = (plant_index, axis_label, metamer_index, organ_label, element_label)
                        if element_id in elements_inputs_grouped.groups:
                            elements_inputs_group = elements_inputs_grouped.get_group(element_id)
                            elements_inputs_group_series = elements_inputs_group.loc[elements_inputs_group.first_valid_index(), SENESCWHEAT_ELEMENTS_INPUTS]
                        else:
                            elements_inputs_group_series = pd.Series()
                        element_inputs = {}
                        is_valid_element = True
                        for element_input_name in SENESCWHEAT_ELEMENTS_INPUTS:
                            if element_input_name in vertex_properties:
                                # use the input from the MTG
                                element_inputs[element_input_name] = vertex_properties[element_input_name]
                            else:
                                # use the input from the dataframe
                                if element_input_name in elements_inputs_group_series:
                                    element_inputs[element_input_name] = elements_inputs_group_series[element_input_name]
                                else:
                                    is_valid_element = False
                                    break
                        if is_valid_element:
                            all_elements_inputs_dict[element_id] = element_inputs

    return {'roots': all_roots_inputs_dict, 'elements': all_elements_inputs_dict}


def update_MTG(inputs, outputs, g):
    """
    Update a MTG from Senesc-Wheat inputs and outputs.

    :Parameters:

            - inputs (:class:`dict` of :class:`dict`) - Senesc-Wheat inputs.
            These inputs are: :mod:`SENESCWHEAT_INPUTS`.

            - outputs (:class:`dict` of :class:`dict`) - Senesc-Wheat outputs.
            These outputs are: :mod:`SENESCWHEAT_OUTPUTS`.

            - `g` (:class:`openalea.mtg.mtg.MTG`) - The MTG to update from the inputs and outputs of SenescWheat.

    .. seealso:: see :attr:`simulation.Simulation.inputs` and :attr:`simulation.Simulation.outputs` for the structure of Senesc-Wheat inputs.

    """
    # add the properties if needed
    property_names = g.property_names()
    for senescwheat_data_name in SENESCWHEAT_INPUTS_OUTPUTS:
        if senescwheat_data_name not in property_names:
            g.add_property(senescwheat_data_name)
    if 'roots' not in property_names:
        g.add_property('roots')

    roots_inputs_dict = inputs['roots']
    elements_inputs_dict = inputs['elements']
    roots_outputs_dict = outputs['roots']
    elements_outputs_dict = outputs['elements']

    # traverse the MTG recursively from top ...
    for plant_vid in g.components_iter(g.root):
        plant_index = int(g.index(plant_vid))
        for axis_vid in g.components_iter(plant_vid):
            axis_label = g.label(axis_vid)
            roots_id = (plant_index, axis_label)
            if roots_id not in roots_outputs_dict: continue
            # update the roots in the MTG
            if 'roots' not in g.get_vertex_property(axis_vid):
                g.property('roots')[axis_vid] = {}
            roots_properties = g.get_vertex_property(axis_vid)['roots']
            roots_properties.update(roots_inputs_dict[roots_id])
            roots_properties.update(roots_outputs_dict[roots_id])
            for metamer_vid in g.components_iter(axis_vid):
                metamer_index = int(g.index(metamer_vid))
                for organ_vid in g.components_iter(metamer_vid):
                    organ_label = g.label(organ_vid)
                    if organ_label not in SENESCWHEAT_PHOTOSYNTHETIC_ORGANS_NAMES: continue
                    for element_vid in g.components_iter(organ_vid):
                        element_label = g.label(element_vid)
                        element_id = (plant_index, axis_label, metamer_index, organ_label, element_label)
                        if element_id not in elements_outputs_dict: continue
                        # update the element in the MTG
                        element_inputs = elements_inputs_dict[element_id]
                        for element_input_name, element_input_value in element_inputs.iteritems():
                            g.property(element_input_name)[element_vid] = element_input_value
                        element_outputs = elements_outputs_dict[element_id]
                        for element_output_name, element_output_value in element_outputs.iteritems():
                            g.property(element_output_name)[element_vid] = element_output_value