# -*- coding: latin-1 -*-

from __future__ import division  # use '//' to do integer division

"""
    senescwheat.model
    ~~~~~~~~~~~~~~~~~~~

    Model of senescence.

    :copyright: Copyright 2014-2015 INRA-ECOSYS, see AUTHORS.
    :license: see LICENSE for details.

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
    CONVERSION_FACTOR_20_TO_12 = 0.45  # modified_Arrhenius_equation(12)/modified_Arrhenius_equation(20)

    N_MOLAR_MASS = 14  #: Molar mass of nitrogen (g mol-1)
    SENESCENCE_ROOTS_POSTFLOWERING = 3.5E-7 * CONVERSION_FACTOR_20_TO_12  #: Rate of root turnover at 12°C (s-1). Value at 20°C coming from Johnson and Thornley (1985), see also Asseng et al. (1997).
    SENESCENCE_ROOTS_PREFLOWERING = 0
    # TODO: should be ontogenic for vegetative stages, 0 in Asseng 1997, but not null in Johnson and Thornley
    FRACTION_N_MAX = {'blade': 0.5, 'stem': 0.425}  # Threshold of ([proteins]/[proteins]max) below which tissue death is triggered
    SENESCENCE_MAX_RATE = 0.2E-8 * CONVERSION_FACTOR_20_TO_12  # maximal senescence m² s-1 at 12°C (Tref)
    SENESCENCE_LENGTH_MAX_RATE = SENESCENCE_MAX_RATE / 3.5e-3  # maximal senescence m s-1 at 12°C (Tref)
    RATIO_N_MSTRUCT = {1: 0.02, 2: 0.02, 3: 0.02, 4: 0.02, 5: 0.0175, 6: 0.015, 7: 0.01, 8: 0.005, 9: 0.005, 10: 0.005, 11: 0.005}  #: N content in total organ mass (senesced + green) according to phytomer rank
    DEFAULT_RATIO_N_MSTRUCT = 0.005  #: default N content in total organ mass (senesced + green) if phytomer rank not found above
    AGE_EFFECT_SENESCENCE = 400  #: Age-induced senescence (degree-day calculated from elong-wheat as equivalent at 12°C)

    # Residual Mass of N in 1 g of mstruct at full senescence of the blade (from experiment NEMA)

    @classmethod
    def calculate_N_content_total(cls, proteins, amino_acids, nitrates, Nstruct, max_mstruct, mstruct, Nresidual):
        """ N content in the whole element (both green and senesced tissues).

        :param float proteins: protein concentration (µmol N proteins g-1 mstruct)
        :param float amino_acids: amino acids concentration (µmol N amino acids g-1 mstruct)
        :param float nitrates: nitrates concentration (µmol N nitrates g-1 mstruct)
        :param float Nstruct: structural N mass (g). Should be constant during leaf life.
        :param float max_mstruct: structural mass maximal of the element i.e. structural mass of the whole element before senescence (g)
        :param float mstruct: structural mass (g)
        :param float Nresidual: residual mass of N in the senescent tissu (g)

        :return: N_content_total (between 0 and 1)
        :rtype: float
        """
        return ((proteins + amino_acids + nitrates) * 1E-6 * cls.N_MOLAR_MASS + Nresidual) / max_mstruct + Nstruct / mstruct

    @classmethod
    def calculate_forced_relative_delta_green_area(cls, green_area_df, group_id, prev_green_area):
        """relative green_area variation due to senescence

        :param pandas.DataFrame green_area_df: a pandas DataFrame containing the green area values for each photosynthetic element at each time
        :param tuple group_id: the group id to be used to select data in the DataFrame
        :param float prev_green_area: previous value of an organ green area (m-2)

        :return: new_green_area (m-2), relative_delta_green_area (dimensionless)
        :rtype: tuple [float, float]
        """
        new_green_area = green_area_df.get_group(group_id).green_area.values[0]
        relative_delta_green_area = (prev_green_area - new_green_area) / prev_green_area
        return new_green_area, relative_delta_green_area

    @classmethod
    def calculate_relative_delta_green_area(cls, organ_name, prev_green_area, proteins, max_proteins, delta_t, update_max_protein):
        """relative green_area variation due to senescence

        :param str organ_name: name of the organ to which belongs the element (used to distinguish lamina from stem organs)
        :param float prev_green_area: previous value of an organ green area (m-2)
        :param float proteins: protein concentration (µmol N proteins g-1 mstruct)
        :param float max_proteins: maximal protein concentrations experienced by the organ (µmol N proteins g-1 mstruct)
        :param float delta_t: value of the timestep (s)
        :param bool update_max_protein: whether to update the max proteins or not.

        :return: new_green_area (m-2), relative_delta_green_area (dimensionless)
        :rtype: tuple [float, float]

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
        elif max_proteins == 0 or (proteins / max_proteins) < fraction_N_max:
            senesced_area = min(prev_green_area, cls.SENESCENCE_MAX_RATE * delta_t)
            new_green_area = max(0., prev_green_area - senesced_area)
            relative_delta_green_area = senesced_area / prev_green_area
        else:
            new_green_area = prev_green_area
            relative_delta_green_area = 0
        return new_green_area, relative_delta_green_area, max_proteins

    # Temporaire
    @classmethod
    def calculate_relative_delta_senesced_length(cls, organ_name, prev_senesced_length, length, proteins, max_proteins, delta_t, update_max_protein):
        """relative senesced length variation

        :param str organ_name: name of the organ to which belongs the element (used to distinguish lamina from stem organs)
        :param float prev_senesced_length: previous senesced length of an organ (m-2)
        :param float length: organ length (m)
        :param float proteins: protein concentration (µmol N proteins g-1 mstruct)
        :param float max_proteins: maximal protein concentrations experienced by the organ (µmol N proteins g-1 mstruct)
        :param float delta_t: value of the timestep (s)
        :param bool update_max_protein: whether to update the max proteins or not.

        :return: new_senesced_length (m), relative_delta_senesced_length (dimensionless), max_proteins (µmol N proteins g-1 mstruct)
        :rtype: tuple [float, float, float]
        
        .. todo:: remove update_max_protein
        """

        if organ_name == 'blade':
            fraction_N_max = cls.FRACTION_N_MAX['blade']
        else:
            fraction_N_max = cls.FRACTION_N_MAX['stem']

        # Overwrite max proteins
        if max_proteins < proteins and update_max_protein:
            max_proteins = proteins
            new_senesced_length = prev_senesced_length
            relative_delta_senesced_length = 0
        # Senescence if (actual proteins/max_proteins) < fraction_N_max
        elif max_proteins == 0 or (proteins / max_proteins) < fraction_N_max:
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

        :param float relative_delta_green_area: relative variation of a photosynthetic element green area (dimensionless)
        :param float prev_mstruct: previous value of an organ structural mass (g)
        :param float prev_Nstruct: previous value of an organ structural N (g)

        :return: new_mstruct (g), new_Nstruct (g)
        :rtype: tuple [float, float]
        """
        new_mstruct = prev_mstruct - prev_mstruct * relative_delta_green_area
        new_Nstruct = prev_Nstruct - prev_Nstruct * relative_delta_green_area
        return new_mstruct, new_Nstruct

    @classmethod
    def calculate_remobilisation(cls, metabolite, relative_delta_structure):
        """Metabolite remobilisation due to senescence over DELTA_T (µmol).
        
        :param float metabolite: amount of any metabolite to be remobilised (µmol) 
        :param float relative_delta_structure: could be relative variation of a photosynthetic element green area or relative variation of mstruct
        
        :return: metabolite remobilisation (µmol)
        :rtype: float
        """
        return metabolite * relative_delta_structure

    @classmethod
    def calculate_remobilisation_proteins(cls, organ, element_index, proteins, relative_delta_green_area, ratio_N_mstruct_max, full_remob):
        """Protein remobilisation due to senescence over DELTA_T. Part is remobilized as amino_acids (µmol N), the rest is increasing Nresidual (g).
        
        :param str organ: name of the organ
        :param int element_index: phytomer rank
        :param float proteins: amount of proteins (µmol N)
        :param float relative_delta_green_area: relative variation of a photosynthetic element green area
        :param float ratio_N_mstruct_max: N content in the whole element (both green and senesced tissues).
        :param bool full_remob: whether all proteins should be remobilised
        
        :return: Quantity of proteins remobilised either in amino acids, either in residual N (µmol),
                 Quantity of proteins converted into amino_acids (µmol N), 
                 Increment of Nresidual (g)
        :rtype: tuple [float, float, float]
        """

        if full_remob or organ != 'blade':
            remob_proteins = delta_amino_acids = proteins * relative_delta_green_area
            delta_Nresidual = 0
        else:
            if ratio_N_mstruct_max <= cls.RATIO_N_MSTRUCT.get(element_index, cls.DEFAULT_RATIO_N_MSTRUCT):  # then all the proteins are converted into Nresidual
                remob_proteins = proteins
                delta_Nresidual = remob_proteins * 1E-6 * cls.N_MOLAR_MASS
                delta_amino_acids = 0
            else:  # then part of the proteins are converted into amino_acids, the rest is converted into Nresidual
                remob_proteins = proteins * relative_delta_green_area
                delta_amino_acids = remob_proteins * 2 / 3.
                delta_Nresidual = (remob_proteins - delta_amino_acids) * 1E-6 * cls.N_MOLAR_MASS
        return remob_proteins, delta_amino_acids, delta_Nresidual

    @classmethod
    def calculate_roots_senescence(cls, mstruct, Nstruct, postflowering_stages):
        """Root senescence

        :param float mstruct: structural mass (g)
        :param float Nstruct: structural N (g)
        :param bool postflowering_stages: Option : True to run a simulation with postflo parameter

        :return: Rate of mstruct loss by root senescence (g mstruct s-1), rate of Nstruct loss by root senescence (g Nstruct s-1)
        :rtype: tuple [float, float]
        """
        if postflowering_stages:
            rate_senescence = cls.SENESCENCE_ROOTS_POSTFLOWERING
        else:
            rate_senescence = cls.SENESCENCE_ROOTS_PREFLOWERING
        return mstruct * rate_senescence, Nstruct * rate_senescence

    @classmethod
    def calculate_relative_delta_mstruct_roots(cls, rate_mstruct_death, root_mstruct, delta_t):
        """Relative delta of root structural dry matter (g) over delta_t

        :param float rate_mstruct_death: Rate of mstruct loss by root senescence (g mstruct s-1)
        :param float root_mstruct: actual mstruct of roots (g)
        :param float delta_t: value of the timestep (s)


        :return: relative_delta_mstruct (dimensionless)
        :rtype: float
        """
        return (rate_mstruct_death * delta_t) / root_mstruct
