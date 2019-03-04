#!/usr/bin/env python
# -*- coding: latin-1 -*-

from __future__ import division  # use '//' to do integer division

"""
    senescwheat.model
    ~~~~~~~~~~~~~~~~~~~

    Model of senescence.

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


class SenescenceModel(object):

    CONVERSION_FACTOR_20_TO_12 = 0.45 # modified_Arrhenius_equation(12)/modified_Arrhenius_equation(20)

    N_MOLAR_MASS = 14             #: Molar mass of nitrogen (g mol-1)
    SENESCENCE_ROOTS = 0 #3.5E-7 * CONVERSION_FACTOR_20_TO_12    #: Rate of root turnover at 12�C (s-1). Value at 20�C coming from Johnson and Thornley (1985), see also Asseng et al. (1997). TODO:
    # should be ontogenic
    # for vegetative stages, 0 in Asseng 1997, but not null in Johnson and Thornley
    FRACTION_N_MAX = {'blade': 0.5, 'stem': 0.425}  # Threshold of ([proteins]/[proteins]max) below which tissue death is triggered
    SENESCENCE_MAX_RATE = 0.2E-8 * CONVERSION_FACTOR_20_TO_12  # maximal senescence m� s-1 at 12�C (Tref)
    SENESCENCE_LENGTH_MAX_RATE = SENESCENCE_MAX_RATE / 3.5e-3  # maximal senescence m s-1 at 12�C (Tref)
    RATIO_N_MSTRUCT = {1:0.02, 2:0.02, 3:0.02, 4:0.02, 5:0.0175, 6:0.015, 7:0.01, 8:0.005, 9:0.005, 10:0.005, 11:0.005 } # Residual Mass of N in 1 g of mstruct at full senescence of the blade (from
    #  experiment NEMA)

    @classmethod
    def calculate_N_content_total(cls, proteins, amino_acids, nitrates, Nstruct, max_mstruct, mstruct, Nresidual):
        """ N content in the whole element (both green and senesced tissues).

        : Parameters:
            - `proteins` (:class:`float`) -protein concentration (�mol N proteins g-1 mstruct)
            - `amino_acids` (:class:`float`) - amino acids concentration (�mol N amino acids g-1 mstruct)
            - `nitrates` (:class:`float`) - nitrates concentration (�mol N nitrates g-1 mstruct)
            - `Nstruct` (:class:`float`) - structural N mass (g) - Should be constant during leaf life.
            - `max_mstruct` (:class:`float`) - structural mass maximal of the element i.e. structural mass of the whole element before senescence (g)
            - `mstruct` (:class:`float`) - structural mass (g)
            - `Nresidual` (:class:`float`) - residual mass of N in the senescent tissu (g)

        : Returns:
            N_content_total (between 0 and 1)

        :Returns Type:
            :class:`float`
        """
        return ( (proteins + amino_acids + nitrates)* 1E-6 * cls.N_MOLAR_MASS + Nresidual) / max_mstruct  + Nstruct/mstruct

    @classmethod
    def calculate_forced_relative_delta_green_area(cls, green_area_df, group_id, prev_green_area):
        """relative green_area variation due to senescence

        : Parameters:
            - `green_area_df` (:class:`pandas.core.frame.DataFrame`) - a pandas DataFrame containing the green area values for each photosynthetic element at each time
            - `group_id` (:class:`tuple`) - the group id to be used to select data in the DataFrame
            - `prev_green_area` (:class:`float`) - previous value of an organ green area (m-2)

        : Returns:
            new_green_area (m-2), relative_delta_green_area (dimensionless)

        :Returns Type:
            :class:`float`
        """
        new_green_area = green_area_df.get_group(group_id).green_area.values[0]
        relative_delta_green_area = (prev_green_area - new_green_area) / prev_green_area
        return new_green_area, relative_delta_green_area

    @classmethod
    def calculate_relative_delta_green_area(cls, organ_name, prev_green_area, proteins, max_proteins, delta_t, update_max_protein):
        """relative green_area variation due to senescence

        : Parameters:
            - `organ_name` (:class:`string`) - name of the organ to which belongs the element (used to distinguish lamina from stem organs)
            - `prev_green_area` (:class:`float`) - previous value of an organ green area (m-2)
            - `proteins` (:class:`float`) - protein concentration (�mol N proteins g-1 mstruct)
            - `max_proteins` (:class:`dict`) - a dictionnary where the maximal protein concentrations are stored by organ id
            - `delta_t` (:class:`float`) - value of the timestep (s)
            - `update_max_protein` (:class:`bool`) - whether to update the max proteins or not.

        : Returns:
            new_green_area (m-2), relative_delta_green_area (dimensionless)

        :Returns Type:
            :class:`float`

        .. todo:: remove update_max_protein

        """

        if organ_name == 'blade':
            fraction_N_max = cls.FRACTION_N_MAX['blade']
        else:
            fraction_N_max = cls.FRACTION_N_MAX['stem']

        # Overwrite max proteins
        if max_proteins < proteins and update_max_protein:
            max_proteins = proteins
            new_green_area = prev_green_area
            relative_delta_green_area = 0
        # Senescence if (actual proteins/max_proteins) < fraction_N_max
        elif (proteins / max_proteins) < fraction_N_max:
            senesced_area = min(prev_green_area, cls.SENESCENCE_MAX_RATE * delta_t)
            new_green_area = max(0, prev_green_area - senesced_area)
            relative_delta_green_area = senesced_area / prev_green_area
        else:
            new_green_area = prev_green_area
            relative_delta_green_area = 0
        return new_green_area, relative_delta_green_area, max_proteins

    # Temporaire
    @classmethod
    def calculate_relative_delta_senesced_length(cls, organ_name, prev_senesced_length,length, proteins, max_proteins, delta_t, update_max_protein):
        """relative green_area variation due to senescence

        : Parameters:
            - `organ_name` (:class:`string`) - name of the organ to which belongs the element (used to distinguish lamina from stem organs)
            - `prev_green_area` (:class:`float`) - previous value of an organ green area (m-2)
            - `proteins` (:class:`float`) - protein concentration (�mol N proteins g-1 mstruct)
            - `max_proteins` (:class:`dict`) - a dictionnary where the maximal protein concentrations are stored by organ id
            - `delta_t` (:class:`float`) - value of the timestep (s)
            - `update_max_protein` (:class:`bool`) - whether to update the max proteins or not.

        : Returns:
            new_green_area (m-2), relative_delta_green_area (dimensionless)

        :Returns Type:
            :class:`float`

        .. todo:: remove update_max_protein

        """

        if organ_name == 'blade':
            fraction_N_max = cls.FRACTION_N_MAX['blade']
        else:
            fraction_N_max = cls.FRACTION_N_MAX['stem']

        # Overwrite max proteins
        if max_proteins < proteins and update_max_protein:
            max_proteins = min( proteins, 4000 ) #TODO : put the number as parameter
            new_senesced_length = prev_senesced_length
            relative_delta_senesced_length = 0
        # Senescence if (actual proteins/max_proteins) < fraction_N_max
        elif (proteins / max_proteins) < fraction_N_max:
            senesced_length = cls.SENESCENCE_LENGTH_MAX_RATE * delta_t
            new_senesced_length = min(length, prev_senesced_length + senesced_length)
            if length == new_senesced_length:
                relative_delta_senesced_length = 1
            else:
                relative_delta_senesced_length = 1 - (length - new_senesced_length) / (length - prev_senesced_length)
        else:
            new_senesced_length = prev_senesced_length
            relative_delta_senesced_length = 0
        return new_senesced_length, relative_delta_senesced_length, max_proteins

    @classmethod
    def calculate_delta_mstruct_shoot(cls, relative_delta_green_area, prev_mstruct, prev_Nstruct):
        """delta of structural mass due to senescence of photosynthetic elements

        : Parameters:
            - `relative_delta_green_area` (:class:`float`) - relative variation of a photosynthetic element green area
            - `prev_mstruct` (:class:`float`) - previous value of an organ structural mass (g)
            - `prev_Nstruct` (:class:`float`) - previous value of an organ structural N (g)

        : Returns:
            new_mstruct (g), new_Nstruct (g)

        :Returns Type:
            :class:`float`
        """
        new_mstruct = prev_mstruct - prev_mstruct*relative_delta_green_area
        new_Nstruct = prev_Nstruct - prev_Nstruct*relative_delta_green_area
        return new_mstruct, new_Nstruct

    @classmethod
    def calculate_remobilisation(cls, metabolite, relative_delta_structure):
        """Metabolite remobilisation due to senescence over DELTA_T (�mol).
        : Parameters:
            - `relative_delta_structure` (:class:`float`) - could be relative variation of a photosynthetic element green area or relative variation of mstruct
        """
        return metabolite * relative_delta_structure

    @classmethod
    def calculate_remobilisation_proteins(cls, organ, element_index, metabolite, relative_delta_structure, ratio_N_mstruct_max):
        """Protein remobilisation due to senescence over DELTA_T. Part is remobilized as amino_acids (�mol N), the rest is increasing Nresidual (g).
        : Parameters:
            - `relative_delta_structure` (:class:`float`) - could be relative variation of a photosynthetic element green area or relative variation of mstruct
        : Returns:
            - Quantity of proteins remobilized either in amino acids, either in residual N
            - Quantity of proteins converted into amino_acids (�mol N)
            - Increment of Nresidual (g)
        """

        if organ != 'blade':
            remob_proteins = delta_aa = metabolite * relative_delta_structure
            delta_Nresidual = 0
        else:
            if ratio_N_mstruct_max <= cls.RATIO_N_MSTRUCT[element_index]: # then all the proteins are converted into Nresidual
                remob_proteins = metabolite
                delta_Nresidual = remob_proteins * 1E-6 * cls.N_MOLAR_MASS
                delta_aa = 0
            else: # then part of the proteins are converted into amino_acids, the rest is converted into Nresidual
                remob_proteins = metabolite * relative_delta_structure
                delta_aa = remob_proteins * 2/3.
                delta_Nresidual = (remob_proteins - delta_aa) * 1E-6 * cls.N_MOLAR_MASS
        return remob_proteins, delta_aa, delta_Nresidual

    @classmethod
    def calculate_roots_senescence(cls, mstruct, Nstruct):
        """Root senescence

        : Parameters:
            - `mstruct` (:class:`float`) - structural mass (g)
            - `Nstruct` (:class:`float`) - structural N (g)

        : Returns:
            Rate of mstruct loss by root senescence (g mstruct s-1), rate of Nstruct loss by root senescence (g Nstruct s-1)

        :Returns Type:
            :class:`float`
        """
        return mstruct * cls.SENESCENCE_ROOTS, Nstruct * cls.SENESCENCE_ROOTS

    @classmethod
    def calculate_relative_delta_mstruct_roots(cls, rate_mstruct_death, root_mstruct, delta_t):
        """Relative delta of root structural dry matter (g) over delta_t

        : Parameters:
            - `rate_mstruct_death` (:class:`float`) - Rate of mstruct loss by root senescence (g mstruct s-1)
            - `root_mstruct` (:class:`float`) - actual mstruct of roots (g)

        : Returns:
            relative_delta_mstruct

        :Returns Type:
            :class:`float`
        """
        return (rate_mstruct_death * delta_t) / root_mstruct
