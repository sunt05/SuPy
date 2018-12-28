
.. option:: aerodynamicresistancemethod

    :Description:
        Internal use. Please DO NOT modify
    :SUEWS-related variables:
        None
    :Dimensionality:
        0
    :Dimensionality remarks:
        Scalar


.. option:: ah_min

    :Description:
        Minimum QF values.
    :SUEWS-related variables:
        `AHMin_WD <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-ahmin-wd>`_, `AHMin_WE <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-ahmin-we>`_
    :Dimensionality:
        (2,)
    :Dimensionality remarks:
        2: Weekday and Weekend


.. option:: ah_slope_cooling

    :Description:
        Cooling slope of QF calculation.
    :SUEWS-related variables:
        `AHSlope_Cooling_WD <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-ahslope-cooling-wd>`_, `AHSlope_Cooling_WE <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-ahslope-cooling-we>`_
    :Dimensionality:
        (2,)
    :Dimensionality remarks:
        2: Weekday and Weekend


.. option:: ah_slope_heating

    :Description:
        Heating slope of QF calculation.
    :SUEWS-related variables:
        `AHSlope_Heating_WD <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-ahslope-heating-wd>`_, `AHSlope_Heating_WE <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-ahslope-heating-we>`_
    :Dimensionality:
        (2,)
    :Dimensionality remarks:
        2: Weekday and Weekend


.. option:: ahprof_24hr

    :Description:
        Hourly profile values used in energy use calculation.
    :SUEWS-related variables:
        `EnergyUseProfWD <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-energyuseprofwd>`_, `EnergyUseProfWE <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-energyuseprofwe>`_
    :Dimensionality:
        (24, 2)
    :Dimensionality remarks:
        24: hours of a day; 2: Weekday and Weekend


.. option:: alb

    :Description:
        Effective surface albedo (middle of the day value) for summertime.
    :SUEWS-related variables:
        `AlbedoMax <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-albedomax>`_
    :Dimensionality:
        (7,)
    :Dimensionality remarks:
        7: Seven SUEWS land cover types in the order of [paved; buildings; evergreen trees/shrubs; deciduous trees/shrubs; grass; bare soil and water]: https://suews-docs.readthedocs.io/en/latest/introduction.html


.. option:: albdectr_id

    :Description:
        Albedo of deciduous surface `DecTr` on day 0 of run
    :SUEWS-related variables:
        `albDecTr0 <https://suews-docs.readthedocs.io/en/latest/input_files/Initial_Conditions/Vegetation_parameters.html#cmdoption-arg-albdectr0>`_
    :Dimensionality:
        0
    :Dimensionality remarks:
        Scalar


.. option:: albevetr_id

    :Description:
        Albedo of evergreen surface `EveTr` on day 0 of run
    :SUEWS-related variables:
        `albEveTr0 <https://suews-docs.readthedocs.io/en/latest/input_files/Initial_Conditions/Vegetation_parameters.html#cmdoption-arg-albevetr0>`_
    :Dimensionality:
        0
    :Dimensionality remarks:
        Scalar


.. option:: albgrass_id

    :Description:
        Albedo of grass surface `Grass` on day 0 of run
    :SUEWS-related variables:
        `albGrass0 <https://suews-docs.readthedocs.io/en/latest/input_files/Initial_Conditions/Vegetation_parameters.html#cmdoption-arg-albgrass0>`_
    :Dimensionality:
        0
    :Dimensionality remarks:
        Scalar


.. option:: albmax_dectr

    :Description:
        Effective surface albedo (middle of the day value) for summertime.
    :SUEWS-related variables:
        `AlbedoMax <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-albedomax>`_
    :Dimensionality:
        0
    :Dimensionality remarks:
        Scalar


.. option:: albmax_evetr

    :Description:
        Effective surface albedo (middle of the day value) for summertime.
    :SUEWS-related variables:
        `AlbedoMax <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-albedomax>`_
    :Dimensionality:
        0
    :Dimensionality remarks:
        Scalar


.. option:: albmax_grass

    :Description:
        Effective surface albedo (middle of the day value) for summertime.
    :SUEWS-related variables:
        `AlbedoMax <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-albedomax>`_
    :Dimensionality:
        0
    :Dimensionality remarks:
        Scalar


.. option:: albmin_dectr

    :Description:
        Effective surface albedo (middle of the day value) for wintertime (not including snow).
    :SUEWS-related variables:
        `AlbedoMin <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-albedomin>`_
    :Dimensionality:
        0
    :Dimensionality remarks:
        Scalar


.. option:: albmin_evetr

    :Description:
        Effective surface albedo (middle of the day value) for wintertime (not including snow).
    :SUEWS-related variables:
        `AlbedoMin <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-albedomin>`_
    :Dimensionality:
        0
    :Dimensionality remarks:
        Scalar


.. option:: albmin_grass

    :Description:
        Effective surface albedo (middle of the day value) for wintertime (not including snow).
    :SUEWS-related variables:
        `AlbedoMin <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-albedomin>`_
    :Dimensionality:
        0
    :Dimensionality remarks:
        Scalar


.. option:: alpha_bioco2

    :Description:
        The mean apparent ecosystem quantum. Represents the initial slope of the light-response curve.
    :SUEWS-related variables:
        `alpha <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-alpha>`_
    :Dimensionality:
        (3,)
    :Dimensionality remarks:
        3: See variable description for specifics


.. option:: alpha_enh_bioco2

    :Description:
        Part of the `alpha` coeﬃcient related to the fraction of vegetation.
    :SUEWS-related variables:
        `alpha_enh <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-alpha-enh>`_
    :Dimensionality:
        (3,)
    :Dimensionality remarks:
        3: See variable description for specifics


.. option:: alt

    :Description:
        Used for both the radiation and water flow between grids.
    :SUEWS-related variables:
        `Alt <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-alt>`_
    :Dimensionality:
        0
    :Dimensionality remarks:
        Scalar


.. option:: baset

    :Description:
        Base Temperature for initiating growing degree days (GDD) for leaf growth. [°C]
    :SUEWS-related variables:
        `BaseT <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-baset>`_
    :Dimensionality:
        (3,)
    :Dimensionality remarks:
        3: See variable description for specifics


.. option:: basete

    :Description:
        Base temperature for initiating sensesance degree days (SDD) for leaf off. [°C]
    :SUEWS-related variables:
        `BaseTe <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-basete>`_
    :Dimensionality:
        (3,)
    :Dimensionality remarks:
        3: See variable description for specifics


.. option:: basethdd

    :Description:
        Base temperature for heating degree days [°C]
    :SUEWS-related variables:
        `BaseTHDD <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-basethdd>`_
    :Dimensionality:
        0
    :Dimensionality remarks:
        Scalar


.. option:: beta_bioco2

    :Description:
        The light-saturated gross photosynthesis of the canopy. [umol |m^-2| |s^-1| ]
    :SUEWS-related variables:
        `beta <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-beta>`_
    :Dimensionality:
        (3,)
    :Dimensionality remarks:
        3: See variable description for specifics


.. option:: beta_enh_bioco2

    :Description:
        Part of the `beta` coeﬃcient related to the fraction of vegetation.
    :SUEWS-related variables:
        `beta_enh <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-beta-enh>`_
    :Dimensionality:
        (3,)
    :Dimensionality remarks:
        3: See variable description for specifics


.. option:: bldgh

    :Description:
        Mean building height [m]
    :SUEWS-related variables:
        `H_Bldgs <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-h-bldgs>`_
    :Dimensionality:
        0
    :Dimensionality remarks:
        Scalar


.. option:: capmax_dec

    :Description:
        Maximum water storage capacity for upper surfaces (i.e. canopy)
    :SUEWS-related variables:
        `StorageMax <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-storagemax>`_
    :Dimensionality:
        0
    :Dimensionality remarks:
        Scalar


.. option:: capmin_dec

    :Description:
        Minimum water storage capacity for upper surfaces (i.e. canopy).
    :SUEWS-related variables:
        `StorageMin <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-storagemin>`_
    :Dimensionality:
        0
    :Dimensionality remarks:
        Scalar


.. option:: chanohm

    :Description:
        Bulk transfer coefficient for this surface to use in AnOHM [-]
    :SUEWS-related variables:
        `AnOHM_Ch <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-anohm-ch>`_
    :Dimensionality:
        (7,)
    :Dimensionality remarks:
        7: Seven SUEWS land cover types in the order of [paved; buildings; evergreen trees/shrubs; deciduous trees/shrubs; grass; bare soil and water]: https://suews-docs.readthedocs.io/en/latest/introduction.html


.. option:: cpanohm

    :Description:
        Volumetric heat capacity for this surface to use in AnOHM [J |m^-3|]
    :SUEWS-related variables:
        `AnOHM_Cp <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-anohm-cp>`_
    :Dimensionality:
        (7,)
    :Dimensionality remarks:
        7: Seven SUEWS land cover types in the order of [paved; buildings; evergreen trees/shrubs; deciduous trees/shrubs; grass; bare soil and water]: https://suews-docs.readthedocs.io/en/latest/introduction.html


.. option:: crwmax

    :Description:
        Maximum water holding capacity of snow [mm]
    :SUEWS-related variables:
        `CRWMax <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-crwmax>`_
    :Dimensionality:
        0
    :Dimensionality remarks:
        Scalar


.. option:: crwmin

    :Description:
        Minimum water holding capacity of snow [mm]
    :SUEWS-related variables:
        `CRWMin <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-crwmin>`_
    :Dimensionality:
        0
    :Dimensionality remarks:
        Scalar


.. option:: daywat

    :Description:
        Irrigation flag: 1 for on and 0 for off.
    :SUEWS-related variables:
        `DayWat(1) <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-daywat-1>`_, `DayWat(2) <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-daywat-2>`_, `DayWat(3) <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-daywat-3>`_, `DayWat(4) <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-daywat-4>`_, `DayWat(5) <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-daywat-5>`_, `DayWat(6) <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-daywat-6>`_, `DayWat(7) <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-daywat-7>`_
    :Dimensionality:
        (7,)
    :Dimensionality remarks:
        7: Seven days of a week: from Sunday to Saturday


.. option:: daywatper

    :Description:
        Fraction of properties using irrigation for each day of a week.
    :SUEWS-related variables:
        `DayWatPer(1) <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-daywatper-1>`_, `DayWatPer(2) <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-daywatper-2>`_, `DayWatPer(3) <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-daywatper-3>`_, `DayWatPer(4) <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-daywatper-4>`_, `DayWatPer(5) <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-daywatper-5>`_, `DayWatPer(6) <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-daywatper-6>`_, `DayWatPer(7) <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-daywatper-7>`_
    :Dimensionality:
        (7,)
    :Dimensionality remarks:
        7: Seven days of a week: from Sunday to Saturday


.. option:: decidcap_id

    :Description:
        Storage capacity of deciduous surface `DecTr` on day 0 of run.
    :SUEWS-related variables:
        `decidCap0 <https://suews-docs.readthedocs.io/en/latest/input_files/Initial_Conditions/Vegetation_parameters.html#cmdoption-arg-decidcap0>`_
    :Dimensionality:
        0
    :Dimensionality remarks:
        Scalar


.. option:: dectreeh

    :Description:
        Mean height of deciduous trees [m]
    :SUEWS-related variables:
        `H_DecTr <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-h-dectr>`_
    :Dimensionality:
        0
    :Dimensionality remarks:
        Scalar


.. option:: diagnose

    :Description:
        Internal use. Please DO NOT modify
    :SUEWS-related variables:
        None
    :Dimensionality:
        0
    :Dimensionality remarks:
        Scalar


.. option:: diagqn

    :Description:
        Internal use. Please DO NOT modify
    :SUEWS-related variables:
        None
    :Dimensionality:
        0
    :Dimensionality remarks:
        Scalar


.. option:: diagqs

    :Description:
        Internal use. Please DO NOT modify
    :SUEWS-related variables:
        None
    :Dimensionality:
        0
    :Dimensionality remarks:
        Scalar


.. option:: drainrt

    :Description:
        Drainage rate of bucket for LUMPS [mm |h^-1|]
    :SUEWS-related variables:
        `LUMPS_DrRate <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-lumps-drrate>`_
    :Dimensionality:
        0
    :Dimensionality remarks:
        Scalar


.. option:: ef_umolco2perj

    :Description:
        Emission factor for fuels used for building heating.
    :SUEWS-related variables:
        `EF_umolCO2perJ <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-ef-umolco2perj>`_
    :Dimensionality:
        0
    :Dimensionality remarks:
        Scalar


.. option:: emis

    :Description:
        Effective surface emissivity.
    :SUEWS-related variables:
        `Emissivity <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-emissivity>`_
    :Dimensionality:
        (7,)
    :Dimensionality remarks:
        7: Seven SUEWS land cover types in the order of [paved; buildings; evergreen trees/shrubs; deciduous trees/shrubs; grass; bare soil and water]: https://suews-docs.readthedocs.io/en/latest/introduction.html


.. option:: emissionsmethod

    :Description:
        Determines method for QF calculation.
    :SUEWS-related variables:
        `EmissionsMethod <https://suews-docs.readthedocs.io/en/latest/input_files/RunControl/scheme_options.html#cmdoption-arg-emissionsmethod>`_
    :Dimensionality:
        0
    :Dimensionality remarks:
        Scalar


.. option:: enddls

    :Description:
        End of the day light savings [DOY]
    :SUEWS-related variables:
        `EndDLS <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-enddls>`_
    :Dimensionality:
        0
    :Dimensionality remarks:
        Scalar


.. option:: enef_v_jkm

    :Description:
        Emission factor for heat [J k|m^-1|].
    :SUEWS-related variables:
        `EnEF_v_Jkm <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-enef-v-jkm>`_
    :Dimensionality:
        0
    :Dimensionality remarks:
        Scalar


.. option:: evapmethod

    :Description:
        Internal use. Please DO NOT modify
    :SUEWS-related variables:
        None
    :Dimensionality:
        0
    :Dimensionality remarks:
        Scalar


.. option:: evetreeh

    :Description:
        Mean height of evergreen trees [m]
    :SUEWS-related variables:
        `H_EveTr <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-h-evetr>`_
    :Dimensionality:
        0
    :Dimensionality remarks:
        Scalar


.. option:: faibldg

    :Description:
        Frontal area index for buildings [-]
    :SUEWS-related variables:
        `FAI_Bldgs <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-fai-bldgs>`_
    :Dimensionality:
        0
    :Dimensionality remarks:
        Scalar


.. option:: faidectree

    :Description:
        Frontal area index for deciduous trees [-]
    :SUEWS-related variables:
        `FAI_DecTr <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-fai-dectr>`_
    :Dimensionality:
        0
    :Dimensionality remarks:
        Scalar


.. option:: faievetree

    :Description:
        Frontal area index for evergreen trees [-]
    :SUEWS-related variables:
        `FAI_EveTr <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-fai-evetr>`_
    :Dimensionality:
        0
    :Dimensionality remarks:
        Scalar


.. option:: faut

    :Description:
        Fraction of irrigated area that is irrigated using automated systems
    :SUEWS-related variables:
        `Faut <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-faut>`_
    :Dimensionality:
        0
    :Dimensionality remarks:
        Scalar


.. option:: fcef_v_kgkm

    :Description:
        CO2 emission factor [kg |km^-1|]
    :SUEWS-related variables:
        `FcEF_v_kgkm <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-fcef-v-kgkm>`_
    :Dimensionality:
        0
    :Dimensionality remarks:
        Scalar


.. option:: flowchange

    :Description:
        Difference in input and output flows for water surface [mm |h^-1|]
    :SUEWS-related variables:
        `FlowChange <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-flowchange>`_
    :Dimensionality:
        0
    :Dimensionality remarks:
        Scalar


.. option:: frfossilfuel_heat

    :Description:
        Fraction of fossil fuels used for building heating [-]
    :SUEWS-related variables:
        `FrFossilFuel_Heat <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-frfossilfuel-heat>`_
    :Dimensionality:
        0
    :Dimensionality remarks:
        Scalar


.. option:: frfossilfuel_nonheat

    :Description:
        Fraction of fossil fuels used for building energy use [-]
    :SUEWS-related variables:
        `FrFossilFuel_NonHeat <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-frfossilfuel-nonheat>`_
    :Dimensionality:
        0
    :Dimensionality remarks:
        Scalar


.. option:: g1

    :Description:
        Related to maximum surface conductance [mm |s^-1|]
    :SUEWS-related variables:
        `G1 <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-g1>`_
    :Dimensionality:
        0
    :Dimensionality remarks:
        Scalar


.. option:: g2

    :Description:
        Related to Kdown dependence [W |m^-2|]
    :SUEWS-related variables:
        `G2 <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-g2>`_
    :Dimensionality:
        0
    :Dimensionality remarks:
        Scalar


.. option:: g3

    :Description:
        Related to VPD dependence [units depend on `gsModel`]
    :SUEWS-related variables:
        `G3 <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-g3>`_
    :Dimensionality:
        0
    :Dimensionality remarks:
        Scalar


.. option:: g4

    :Description:
        Related to VPD dependence [units depend on `gsModel`]
    :SUEWS-related variables:
        `G4 <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-g4>`_
    :Dimensionality:
        0
    :Dimensionality remarks:
        Scalar


.. option:: g5

    :Description:
        Related to temperature dependence [°C]
    :SUEWS-related variables:
        `G5 <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-g5>`_
    :Dimensionality:
        0
    :Dimensionality remarks:
        Scalar


.. option:: g6

    :Description:
        Related to soil moisture dependence [|mm^-1|]
    :SUEWS-related variables:
        `G6 <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-g6>`_
    :Dimensionality:
        0
    :Dimensionality remarks:
        Scalar


.. option:: gddfull

    :Description:
        The growing degree days (GDD) needed for full capacity of the leaf area index (LAI) [°C].
    :SUEWS-related variables:
        `GDDFull <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-gddfull>`_
    :Dimensionality:
        (3,)
    :Dimensionality remarks:
        3: See variable description for specifics


.. option:: gsmodel

    :Description:
        Formulation choice for conductance calculation.
    :SUEWS-related variables:
        `gsModel <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-gsmodel>`_
    :Dimensionality:
        0
    :Dimensionality remarks:
        Scalar


.. option:: humactivity_24hr

    :Description:
        Hourly profile values used in human activity calculation.
    :SUEWS-related variables:
        `ActivityProfWD <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-activityprofwd>`_, `ActivityProfWE <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-activityprofwe>`_
    :Dimensionality:
        (24, 2)
    :Dimensionality remarks:
        24: hours of a day; 2: Weekday and Weekend


.. option:: ie_a

    :Description:
        Coefficient for automatic irrigation model.
    :SUEWS-related variables:
        `Ie_a1 <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-ie-a1>`_, `Ie_a2 <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-ie-a2>`_, `Ie_a3 <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-ie-a3>`_
    :Dimensionality:
        (3,)
    :Dimensionality remarks:
        3: See variable description for specifics


.. option:: ie_end

    :Description:
        Day when irrigation ends [DOY]
    :SUEWS-related variables:
        `Ie_end <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-ie-end>`_
    :Dimensionality:
        0
    :Dimensionality remarks:
        Scalar


.. option:: ie_m

    :Description:
        Coefficient for manual irrigation model.
    :SUEWS-related variables:
        `Ie_m1 <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-ie-m1>`_, `Ie_m2 <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-ie-m2>`_, `Ie_m3 <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-ie-m3>`_
    :Dimensionality:
        (3,)
    :Dimensionality remarks:
        3: See variable description for specifics


.. option:: ie_start

    :Description:
        Day when irrigation starts [DOY]
    :SUEWS-related variables:
        `Ie_start <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-ie-start>`_
    :Dimensionality:
        0
    :Dimensionality remarks:
        Scalar


.. option:: internalwateruse_h

    :Description:
        Internal water use [mm |h^-1|]
    :SUEWS-related variables:
        `InternalWaterUse <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-internalwateruse>`_
    :Dimensionality:
        0
    :Dimensionality remarks:
        Scalar


.. option:: irrfracconif

    :Description:
        Fraction of evergreen trees that are irrigated [-]
    :SUEWS-related variables:
        `IrrFr_EveTr <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-irrfr-evetr>`_
    :Dimensionality:
        0
    :Dimensionality remarks:
        Scalar


.. option:: irrfracdecid

    :Description:
        Fraction of deciduous trees that are irrigated [-]
    :SUEWS-related variables:
        `IrrFr_DecTr <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-irrfr-dectr>`_
    :Dimensionality:
        0
    :Dimensionality remarks:
        Scalar


.. option:: irrfracgrass

    :Description:
        Fraction of `Grass` that is irrigated [-]
    :SUEWS-related variables:
        `IrrFr_Grass <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-irrfr-grass>`_
    :Dimensionality:
        0
    :Dimensionality remarks:
        Scalar


.. option:: kkanohm

    :Description:
        Thermal conductivity for this surface to use in AnOHM [W m |K^-1|]
    :SUEWS-related variables:
        `AnOHM_Kk <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-anohm-kk>`_
    :Dimensionality:
        (7,)
    :Dimensionality remarks:
        7: Seven SUEWS land cover types in the order of [paved; buildings; evergreen trees/shrubs; deciduous trees/shrubs; grass; bare soil and water]: https://suews-docs.readthedocs.io/en/latest/introduction.html


.. option:: kmax

    :Description:
        Maximum incoming shortwave radiation [W |m^-2|]
    :SUEWS-related variables:
        `Kmax <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-kmax>`_
    :Dimensionality:
        0
    :Dimensionality remarks:
        Scalar


.. option:: lai_id

    :Description:
        Initial LAI values.
    :SUEWS-related variables:
        `LAIinitialDecTr <https://suews-docs.readthedocs.io/en/latest/input_files/Initial_Conditions/Vegetation_parameters.html#cmdoption-arg-laiinitialdectr>`_, `LAIinitialEveTr <https://suews-docs.readthedocs.io/en/latest/input_files/Initial_Conditions/Vegetation_parameters.html#cmdoption-arg-laiinitialevetr>`_, `LAIinitialGrass <https://suews-docs.readthedocs.io/en/latest/input_files/Initial_Conditions/Vegetation_parameters.html#cmdoption-arg-laiinitialgrass>`_
    :Dimensionality:
        (3,)
    :Dimensionality remarks:
        3: See variable description for specifics


.. option:: laicalcyes

    :Description:
        Internal use. Please DO NOT modify
    :SUEWS-related variables:
        None
    :Dimensionality:
        0
    :Dimensionality remarks:
        Scalar


.. option:: laimax

    :Description:
        full leaf-on summertime value
    :SUEWS-related variables:
        `LAIMax <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-laimax>`_
    :Dimensionality:
        (3,)
    :Dimensionality remarks:
        3: See variable description for specifics


.. option:: laimin

    :Description:
        leaf-off wintertime value
    :SUEWS-related variables:
        `LAIMin <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-laimin>`_
    :Dimensionality:
        (3,)
    :Dimensionality remarks:
        3: See variable description for specifics


.. option:: laipower

    :Description:
        parameters required by LAI calculation.
    :SUEWS-related variables:
        `LeafGrowthPower1 <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-leafgrowthpower1>`_, `LeafGrowthPower2 <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-leafgrowthpower2>`_, `LeafOffPower1 <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-leafoffpower1>`_, `LeafOffPower2 <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-leafoffpower2>`_
    :Dimensionality:
        (4, 3)
    :Dimensionality remarks:
        4: See variable description for specifics; 3: Three vegetated land cover types (`EveTr`, `DecTr`, `Grass`)


.. option:: laitype

    :Description:
        LAI calculation choice.
    :SUEWS-related variables:
        `LAIEq <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-laieq>`_
    :Dimensionality:
        (3,)
    :Dimensionality remarks:
        3: See variable description for specifics


.. option:: lat

    :Description:
        Latitude [deg].
    :SUEWS-related variables:
        `lat <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-lat>`_
    :Dimensionality:
        0
    :Dimensionality remarks:
        Scalar


.. option:: lng

    :Description:
        longitude [deg]
    :SUEWS-related variables:
        `lng <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-lng>`_
    :Dimensionality:
        0
    :Dimensionality remarks:
        Scalar


.. option:: maxconductance

    :Description:
        The maximum conductance of each vegetation or surface type. [mm |s^-1|]
    :SUEWS-related variables:
        `MaxConductance <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-maxconductance>`_
    :Dimensionality:
        (3,)
    :Dimensionality remarks:
        3: See variable description for specifics


.. option:: maxqfmetab

    :Description:
        Maximum value for human heat emission. [W |m^-2|]
    :SUEWS-related variables:
        `MaxQFMetab <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-maxqfmetab>`_
    :Dimensionality:
        0
    :Dimensionality remarks:
        Scalar


.. option:: min_res_bioco2

    :Description:
        Minimum soil respiration rate (for cold-temperature limit) [umol |m^-2| |s^-1|].
    :SUEWS-related variables:
        `min_respi <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-min-respi>`_
    :Dimensionality:
        (3,)
    :Dimensionality remarks:
        3: See variable description for specifics


.. option:: minqfmetab

    :Description:
        Minimum value for human heat emission. [W |m^-2|]
    :SUEWS-related variables:
        `MinQFMetab <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-minqfmetab>`_
    :Dimensionality:
        0
    :Dimensionality remarks:
        Scalar


.. option:: narp_emis_snow

    :Description:
        Effective surface emissivity.
    :SUEWS-related variables:
        `Emissivity <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-emissivity>`_
    :Dimensionality:
        0
    :Dimensionality remarks:
        Scalar


.. option:: narp_trans_site

    :Description:
        Atmospheric transmissivity for NARP [-]
    :SUEWS-related variables:
        `NARP_Trans <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-narp-trans>`_
    :Dimensionality:
        0
    :Dimensionality remarks:
        Scalar


.. option:: netradiationmethod

    :Description:
        Determines method for calculation of radiation fluxes.
    :SUEWS-related variables:
        `NetRadiationMethod <https://suews-docs.readthedocs.io/en/latest/input_files/RunControl/scheme_options.html#cmdoption-arg-netradiationmethod>`_
    :Dimensionality:
        0
    :Dimensionality remarks:
        Scalar


.. option:: ohm_coef

    :Description:
        Coefficients for OHM calculation.
    :SUEWS-related variables:
        `a1 <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-a1>`_, `a2 <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-a2>`_, `a3 <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-a3>`_
    :Dimensionality:
        (8, 4, 3)
    :Dimensionality remarks:
        8: Seven SUEWS land cover types and one extra land cover type (currently NOT used); 4: SummerWet, SummerDry, WinterWet, WinterDry; 3: a1, a2, a3


.. option:: ohm_threshsw

    :Description:
        Temperature threshold determining whether summer/winter OHM coefficients are applied [°C]
    :SUEWS-related variables:
        `OHMThresh_SW <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-ohmthresh-sw>`_
    :Dimensionality:
        (8,)
    :Dimensionality remarks:
        8: Seven SUEWS land cover types and one extra land cover type (currently NOT used)


.. option:: ohm_threshwd

    :Description:
        Soil moisture threshold determining whether wet/dry OHM coefficients are applied [-]
    :SUEWS-related variables:
        `OHMThresh_WD <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-ohmthresh-wd>`_
    :Dimensionality:
        (8,)
    :Dimensionality remarks:
        8: Seven SUEWS land cover types and one extra land cover type (currently NOT used)


.. option:: ohmincqf

    :Description:
        Determines whether the storage heat flux calculation uses |Qstar| or ( |Qstar| +QF).
    :SUEWS-related variables:
        `OHMIncQF <https://suews-docs.readthedocs.io/en/latest/input_files/RunControl/scheme_options.html#cmdoption-arg-ohmincqf>`_
    :Dimensionality:
        0
    :Dimensionality remarks:
        Scalar


.. option:: pipecapacity

    :Description:
        Storage capacity of pipes [mm]
    :SUEWS-related variables:
        `PipeCapacity <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-pipecapacity>`_
    :Dimensionality:
        0
    :Dimensionality remarks:
        Scalar


.. option:: popdensdaytime

    :Description:
        Daytime population density (i.e. workers, tourists) [people |ha^-1|]
    :SUEWS-related variables:
        `PopDensDay <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-popdensday>`_
    :Dimensionality:
        0
    :Dimensionality remarks:
        Scalar


.. option:: popdensnighttime

    :Description:
        Night-time population density (i.e. residents) [people |ha^-1|]
    :SUEWS-related variables:
        `PopDensNight <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-popdensnight>`_
    :Dimensionality:
        0
    :Dimensionality remarks:
        Scalar


.. option:: popprof_24hr

    :Description:
        Hourly profile values used in dynamic population estimation.
    :SUEWS-related variables:
        `PopProfWD <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-popprofwd>`_, `PopProfWE <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-popprofwe>`_
    :Dimensionality:
        (24, 2)
    :Dimensionality remarks:
        24: hours of a day; 2: Weekday and Weekend


.. option:: pormax_dec

    :Description:
        full leaf-on summertime value Used only for `DecTr` (can affect roughness calculation)
    :SUEWS-related variables:
        `PorosityMax <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-porositymax>`_
    :Dimensionality:
        0
    :Dimensionality remarks:
        Scalar


.. option:: pormin_dec

    :Description:
        leaf-off wintertime value Used only for `DecTr` (can affect roughness calculation)
    :SUEWS-related variables:
        `PorosityMin <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-porositymin>`_
    :Dimensionality:
        0
    :Dimensionality remarks:
        Scalar


.. option:: porosity_id

    :Description:
        Porosity of deciduous vegetation on day 0 of run.
    :SUEWS-related variables:
        `porosity0 <https://suews-docs.readthedocs.io/en/latest/input_files/Initial_Conditions/Vegetation_parameters.html#cmdoption-arg-porosity0>`_
    :Dimensionality:
        0
    :Dimensionality remarks:
        Scalar


.. option:: preciplimit

    :Description:
        Limit for hourly snowfall when the ground is fully covered with snow [mm]
    :SUEWS-related variables:
        `PrecipLimSnow <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-preciplimsnow>`_
    :Dimensionality:
        0
    :Dimensionality remarks:
        Scalar


.. option:: preciplimitalb

    :Description:
        Limit for hourly precipitation when the ground is fully covered with snow. Then snow albedo is reset to AlbedoMax [mm]
    :SUEWS-related variables:
        `PrecipLimAlb <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-preciplimalb>`_
    :Dimensionality:
        0
    :Dimensionality remarks:
        Scalar


.. option:: qf0_beu

    :Description:
        Building energy use [W |m^-2|]
    :SUEWS-related variables:
        `QF0_BEU_WD <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-qf0-beu-wd>`_, `QF0_BEU_WE <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-qf0-beu-we>`_
    :Dimensionality:
        (2,)
    :Dimensionality remarks:
        2: Weekday and Weekend


.. option:: qf_a

    :Description:
        Base value for QF calculation.
    :SUEWS-related variables:
        `QF_A_WD <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-qf-a-wd>`_, `QF_A_WE <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-qf-a-we>`_
    :Dimensionality:
        (2,)
    :Dimensionality remarks:
        2: Weekday and Weekend


.. option:: qf_b

    :Description:
        Parameter related to heating degree days.
    :SUEWS-related variables:
        `QF_B_WD <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-qf-b-wd>`_, `QF_B_WE <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-qf-b-we>`_
    :Dimensionality:
        (2,)
    :Dimensionality remarks:
        2: Weekday and Weekend


.. option:: qf_c

    :Description:
        Parameter related to heating degree days.
    :SUEWS-related variables:
        `QF_C_WD <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-qf-c-wd>`_, `QF_C_WE <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-qf-c-we>`_
    :Dimensionality:
        (2,)
    :Dimensionality remarks:
        2: Weekday and Weekend


.. option:: radmeltfact

    :Description:
        Hourly radiation melt factor of snow [mm |w^-1| |h^-1|]
    :SUEWS-related variables:
        `RadMeltFactor <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-radmeltfactor>`_
    :Dimensionality:
        0
    :Dimensionality remarks:
        Scalar


.. option:: raincover

    :Description:
        Limit when surface totally covered with water for LUMPS [mm]
    :SUEWS-related variables:
        `LUMPS_Cover <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-lumps-cover>`_
    :Dimensionality:
        0
    :Dimensionality remarks:
        Scalar


.. option:: rainmaxres

    :Description:
        Maximum water bucket reservoir [mm] Used for LUMPS surface wetness control.
    :SUEWS-related variables:
        `LUMPS_MaxRes <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-lumps-maxres>`_
    :Dimensionality:
        0
    :Dimensionality remarks:
        Scalar


.. option:: resp_a

    :Description:
        Respiration coeﬃcient a.
    :SUEWS-related variables:
        `resp_a <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-resp-a>`_
    :Dimensionality:
        (3,)
    :Dimensionality remarks:
        3: See variable description for specifics


.. option:: resp_b

    :Description:
        Respiration coeﬃcient b - related to air temperature dependency.
    :SUEWS-related variables:
        `resp_b <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-resp-b>`_
    :Dimensionality:
        (3,)
    :Dimensionality remarks:
        3: See variable description for specifics


.. option:: roughlenheatmethod

    :Description:
        Determines method for calculating roughness length for heat.
    :SUEWS-related variables:
        `RoughLenHeatMethod <https://suews-docs.readthedocs.io/en/latest/input_files/RunControl/scheme_options.html#cmdoption-arg-roughlenheatmethod>`_
    :Dimensionality:
        0
    :Dimensionality remarks:
        Scalar


.. option:: roughlenmommethod

    :Description:
        Determines how aerodynamic roughness length (z0m) and zero displacement height (zdm) are calculated.
    :SUEWS-related variables:
        `RoughLenMomMethod <https://suews-docs.readthedocs.io/en/latest/input_files/RunControl/scheme_options.html#cmdoption-arg-roughlenmommethod>`_
    :Dimensionality:
        0
    :Dimensionality remarks:
        Scalar


.. option:: runofftowater

    :Description:
        Fraction of above-ground runoff flowing to water surface during flooding [-]
    :SUEWS-related variables:
        `RunoffToWater <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-runofftowater>`_
    :Dimensionality:
        0
    :Dimensionality remarks:
        Scalar


.. option:: s1

    :Description:
        A parameter related to soil moisture dependence [-]
    :SUEWS-related variables:
        `S1 <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-s1>`_
    :Dimensionality:
        0
    :Dimensionality remarks:
        Scalar


.. option:: s2

    :Description:
        A parameter related to soil moisture dependence [mm]
    :SUEWS-related variables:
        `S2 <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-s2>`_
    :Dimensionality:
        0
    :Dimensionality remarks:
        Scalar


.. option:: sathydraulicconduct

    :Description:
        Hydraulic conductivity for saturated soil [mm |s^-1|]
    :SUEWS-related variables:
        `SatHydraulicCond <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-sathydrauliccond>`_
    :Dimensionality:
        (7,)
    :Dimensionality remarks:
        7: Seven SUEWS land cover types in the order of [paved; buildings; evergreen trees/shrubs; deciduous trees/shrubs; grass; bare soil and water]: https://suews-docs.readthedocs.io/en/latest/introduction.html


.. option:: sddfull

    :Description:
        The sensesence degree days (SDD) needed to initiate leaf off. [°C]
    :SUEWS-related variables:
        `SDDFull <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-sddfull>`_
    :Dimensionality:
        (3,)
    :Dimensionality remarks:
        3: See variable description for specifics


.. option:: sfr

    :Description:
        Surface cover fractions.
    :SUEWS-related variables:
        `Fr_Bldgs <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-fr-bldgs>`_, `Fr_Bsoil <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-fr-bsoil>`_, `Fr_DecTr <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-fr-dectr>`_, `Fr_EveTr <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-fr-evetr>`_, `Fr_Grass <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-fr-grass>`_, `Fr_Paved <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-fr-paved>`_, `Fr_Water <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-fr-water>`_
    :Dimensionality:
        (7,)
    :Dimensionality remarks:
        7: Seven SUEWS land cover types in the order of [paved; buildings; evergreen trees/shrubs; deciduous trees/shrubs; grass; bare soil and water]: https://suews-docs.readthedocs.io/en/latest/introduction.html


.. option:: smdmethod

    :Description:
        Determines method for calculating soil moisture deficit (SMD).
    :SUEWS-related variables:
        `SMDMethod <https://suews-docs.readthedocs.io/en/latest/input_files/RunControl/scheme_options.html#cmdoption-arg-smdmethod>`_
    :Dimensionality:
        0
    :Dimensionality remarks:
        Scalar


.. option:: snowalb

    :Description:
        Initial snow albedo
    :SUEWS-related variables:
        `SnowAlb0 <https://suews-docs.readthedocs.io/en/latest/input_files/Initial_Conditions/Snow_related_parameters.html#cmdoption-arg-snowalb0>`_
    :Dimensionality:
        0
    :Dimensionality remarks:
        Scalar


.. option:: snowalbmax

    :Description:
        Effective surface albedo (middle of the day value) for summertime.
    :SUEWS-related variables:
        `AlbedoMax <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-albedomax>`_
    :Dimensionality:
        0
    :Dimensionality remarks:
        Scalar


.. option:: snowalbmin

    :Description:
        Effective surface albedo (middle of the day value) for wintertime (not including snow).
    :SUEWS-related variables:
        `AlbedoMin <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-albedomin>`_
    :Dimensionality:
        0
    :Dimensionality remarks:
        Scalar


.. option:: snowd

    :Description:
        Limit for the snow water equivalent when snow cover starts to be patchy [mm]
    :SUEWS-related variables:
        `SnowLimPatch <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-snowlimpatch>`_
    :Dimensionality:
        (7,)
    :Dimensionality remarks:
        7: Seven SUEWS land cover types in the order of [paved; buildings; evergreen trees/shrubs; deciduous trees/shrubs; grass; bare soil and water]: https://suews-docs.readthedocs.io/en/latest/introduction.html


.. option:: snowdens

    :Description:
        Initial snow density of each land cover.
    :SUEWS-related variables:
        `SnowDensBldgs <https://suews-docs.readthedocs.io/en/latest/input_files/Initial_Conditions/Snow_related_parameters.html#cmdoption-arg-snowdensbldgs>`_, `SnowDensPaved <https://suews-docs.readthedocs.io/en/latest/input_files/Initial_Conditions/Snow_related_parameters.html#cmdoption-arg-snowdenspaved>`_, `SnowDensDecTr <https://suews-docs.readthedocs.io/en/latest/input_files/Initial_Conditions/Snow_related_parameters.html#cmdoption-arg-snowdensdectr>`_, `SnowDensEveTr <https://suews-docs.readthedocs.io/en/latest/input_files/Initial_Conditions/Snow_related_parameters.html#cmdoption-arg-snowdensevetr>`_, `SnowDensGrass <https://suews-docs.readthedocs.io/en/latest/input_files/Initial_Conditions/Snow_related_parameters.html#cmdoption-arg-snowdensgrass>`_, `SnowDensBSoil <https://suews-docs.readthedocs.io/en/latest/input_files/Initial_Conditions/Snow_related_parameters.html#cmdoption-arg-snowdensbsoil>`_, `SnowDensWater <https://suews-docs.readthedocs.io/en/latest/input_files/Initial_Conditions/Snow_related_parameters.html#cmdoption-arg-snowdenswater>`_
    :Dimensionality:
        (7,)
    :Dimensionality remarks:
        7: Seven SUEWS land cover types in the order of [paved; buildings; evergreen trees/shrubs; deciduous trees/shrubs; grass; bare soil and water]: https://suews-docs.readthedocs.io/en/latest/introduction.html


.. option:: snowdensmax

    :Description:
        Maximum snow density [kg |m^-3|]
    :SUEWS-related variables:
        `SnowDensMax <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-snowdensmax>`_
    :Dimensionality:
        0
    :Dimensionality remarks:
        Scalar


.. option:: snowdensmin

    :Description:
        Fresh snow density [kg |m^-3|]
    :SUEWS-related variables:
        `SnowDensMin <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-snowdensmin>`_
    :Dimensionality:
        0
    :Dimensionality remarks:
        Scalar


.. option:: snowfrac

    :Description:
        Initial plan area fraction of snow on each land cover`
    :SUEWS-related variables:
        `SnowFracBldgs <https://suews-docs.readthedocs.io/en/latest/input_files/Initial_Conditions/Snow_related_parameters.html#cmdoption-arg-snowfracbldgs>`_, `SnowFracPaved <https://suews-docs.readthedocs.io/en/latest/input_files/Initial_Conditions/Snow_related_parameters.html#cmdoption-arg-snowfracpaved>`_, `SnowFracDecTr <https://suews-docs.readthedocs.io/en/latest/input_files/Initial_Conditions/Snow_related_parameters.html#cmdoption-arg-snowfracdectr>`_, `SnowFracEveTr <https://suews-docs.readthedocs.io/en/latest/input_files/Initial_Conditions/Snow_related_parameters.html#cmdoption-arg-snowfracevetr>`_, `SnowFracGrass <https://suews-docs.readthedocs.io/en/latest/input_files/Initial_Conditions/Snow_related_parameters.html#cmdoption-arg-snowfracgrass>`_, `SnowFracBSoil <https://suews-docs.readthedocs.io/en/latest/input_files/Initial_Conditions/Snow_related_parameters.html#cmdoption-arg-snowfracbsoil>`_, `SnowFracWater <https://suews-docs.readthedocs.io/en/latest/input_files/Initial_Conditions/Snow_related_parameters.html#cmdoption-arg-snowfracwater>`_
    :Dimensionality:
        (7,)
    :Dimensionality remarks:
        7: Seven SUEWS land cover types in the order of [paved; buildings; evergreen trees/shrubs; deciduous trees/shrubs; grass; bare soil and water]: https://suews-docs.readthedocs.io/en/latest/introduction.html


.. option:: snowlimbuild

    :Description:
        Limit of the snow water equivalent for snow removal from roads and roofs [mm]
    :SUEWS-related variables:
        `SnowLimRemove <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-snowlimremove>`_
    :Dimensionality:
        0
    :Dimensionality remarks:
        Scalar


.. option:: snowlimpaved

    :Description:
        Limit of the snow water equivalent for snow removal from roads and roofs [mm]
    :SUEWS-related variables:
        `SnowLimRemove <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-snowlimremove>`_
    :Dimensionality:
        0
    :Dimensionality remarks:
        Scalar


.. option:: snowpack

    :Description:
        Initial snow water equivalent on each land cover
    :SUEWS-related variables:
        `SnowPackBldgs <https://suews-docs.readthedocs.io/en/latest/input_files/Initial_Conditions/Snow_related_parameters.html#cmdoption-arg-snowpackbldgs>`_, `SnowPackPaved <https://suews-docs.readthedocs.io/en/latest/input_files/Initial_Conditions/Snow_related_parameters.html#cmdoption-arg-snowpackpaved>`_, `SnowPackDecTr <https://suews-docs.readthedocs.io/en/latest/input_files/Initial_Conditions/Snow_related_parameters.html#cmdoption-arg-snowpackdectr>`_, `SnowPackEveTr <https://suews-docs.readthedocs.io/en/latest/input_files/Initial_Conditions/Snow_related_parameters.html#cmdoption-arg-snowpackevetr>`_, `SnowPackGrass <https://suews-docs.readthedocs.io/en/latest/input_files/Initial_Conditions/Snow_related_parameters.html#cmdoption-arg-snowpackgrass>`_, `SnowPackBSoil <https://suews-docs.readthedocs.io/en/latest/input_files/Initial_Conditions/Snow_related_parameters.html#cmdoption-arg-snowpackbsoil>`_, `SnowPackWater <https://suews-docs.readthedocs.io/en/latest/input_files/Initial_Conditions/Snow_related_parameters.html#cmdoption-arg-snowpackwater>`_
    :Dimensionality:
        (7,)
    :Dimensionality remarks:
        7: Seven SUEWS land cover types in the order of [paved; buildings; evergreen trees/shrubs; deciduous trees/shrubs; grass; bare soil and water]: https://suews-docs.readthedocs.io/en/latest/introduction.html


.. option:: snowprof_24hr

    :Description:
        Hourly profile values used in snow clearing.
    :SUEWS-related variables:
        `SnowClearingProfWD <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-snowclearingprofwd>`_, `SnowClearingProfWE <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-snowclearingprofwe>`_
    :Dimensionality:
        (24, 2)
    :Dimensionality remarks:
        24: hours of a day; 2: Weekday and Weekend


.. option:: snowuse

    :Description:
        Determines whether the snow part of the model runs.
    :SUEWS-related variables:
        `SnowUse <https://suews-docs.readthedocs.io/en/latest/input_files/RunControl/scheme_options.html#cmdoption-arg-snowuse>`_
    :Dimensionality:
        0
    :Dimensionality remarks:
        Scalar


.. option:: snowwater

    :Description:
        Initial amount of liquid water in the snow on each land cover
    :SUEWS-related variables:
        `SnowWaterBldgsState <https://suews-docs.readthedocs.io/en/latest/input_files/Initial_Conditions/Snow_related_parameters.html#cmdoption-arg-snowwaterbldgsstate>`_, `SnowWaterPavedState <https://suews-docs.readthedocs.io/en/latest/input_files/Initial_Conditions/Snow_related_parameters.html#cmdoption-arg-snowwaterpavedstate>`_, `SnowWaterDecTrState <https://suews-docs.readthedocs.io/en/latest/input_files/Initial_Conditions/Snow_related_parameters.html#cmdoption-arg-snowwaterdectrstate>`_, `SnowWaterEveTrState <https://suews-docs.readthedocs.io/en/latest/input_files/Initial_Conditions/Snow_related_parameters.html#cmdoption-arg-snowwaterevetrstate>`_, `SnowWaterGrassState <https://suews-docs.readthedocs.io/en/latest/input_files/Initial_Conditions/Snow_related_parameters.html#cmdoption-arg-snowwatergrassstate>`_, `SnowWaterBSoilState <https://suews-docs.readthedocs.io/en/latest/input_files/Initial_Conditions/Snow_related_parameters.html#cmdoption-arg-snowwaterbsoilstate>`_, `SnowWaterWaterState <https://suews-docs.readthedocs.io/en/latest/input_files/Initial_Conditions/Snow_related_parameters.html#cmdoption-arg-snowwaterwaterstate>`_
    :Dimensionality:
        (7,)
    :Dimensionality remarks:
        7: Seven SUEWS land cover types in the order of [paved; buildings; evergreen trees/shrubs; deciduous trees/shrubs; grass; bare soil and water]: https://suews-docs.readthedocs.io/en/latest/introduction.html


.. option:: soildepth

    :Description:
        Depth of soil beneath the surface [mm]
    :SUEWS-related variables:
        `SoilDepth <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-soildepth>`_
    :Dimensionality:
        (7,)
    :Dimensionality remarks:
        7: Seven SUEWS land cover types in the order of [paved; buildings; evergreen trees/shrubs; deciduous trees/shrubs; grass; bare soil and water]: https://suews-docs.readthedocs.io/en/latest/introduction.html


.. option:: soilstore_id

    :Description:
        Initial water stored in soil beneath each land cover
    :SUEWS-related variables:
        `SoilstoreBldgsState <https://suews-docs.readthedocs.io/en/latest/input_files/Initial_Conditions/Soil_moisture_states.html#cmdoption-arg-soilstorebldgsstate>`_, `SoilstorePavedState <https://suews-docs.readthedocs.io/en/latest/input_files/Initial_Conditions/Soil_moisture_states.html#cmdoption-arg-soilstorepavedstate>`_, `SoilstoreDecTrState <https://suews-docs.readthedocs.io/en/latest/input_files/Initial_Conditions/Soil_moisture_states.html#cmdoption-arg-soilstoredectrstate>`_, `SoilstoreEveTrState <https://suews-docs.readthedocs.io/en/latest/input_files/Initial_Conditions/Soil_moisture_states.html#cmdoption-arg-soilstoreevetrstate>`_, `SoilstoreGrassState <https://suews-docs.readthedocs.io/en/latest/input_files/Initial_Conditions/Soil_moisture_states.html#cmdoption-arg-soilstoregrassstate>`_, `SoilstoreBSoilState <https://suews-docs.readthedocs.io/en/latest/input_files/Initial_Conditions/Soil_moisture_states.html#cmdoption-arg-soilstorebsoilstate>`_
    :Dimensionality:
        (7,)
    :Dimensionality remarks:
        7: Seven SUEWS land cover types in the order of [paved; buildings; evergreen trees/shrubs; deciduous trees/shrubs; grass; bare soil and water]: https://suews-docs.readthedocs.io/en/latest/introduction.html


.. option:: soilstorecap

    :Description:
        Limit value for `SoilDepth` [mm]
    :SUEWS-related variables:
        `SoilStoreCap <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-soilstorecap>`_
    :Dimensionality:
        (7,)
    :Dimensionality remarks:
        7: Seven SUEWS land cover types in the order of [paved; buildings; evergreen trees/shrubs; deciduous trees/shrubs; grass; bare soil and water]: https://suews-docs.readthedocs.io/en/latest/introduction.html


.. option:: stabilitymethod

    :Description:
        Defines which atmospheric stability functions are used.
    :SUEWS-related variables:
        `StabilityMethod <https://suews-docs.readthedocs.io/en/latest/input_files/RunControl/scheme_options.html#cmdoption-arg-stabilitymethod>`_
    :Dimensionality:
        0
    :Dimensionality remarks:
        Scalar


.. option:: startdls

    :Description:
        Start of the day light savings [DOY]
    :SUEWS-related variables:
        `StartDLS <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-startdls>`_
    :Dimensionality:
        0
    :Dimensionality remarks:
        Scalar


.. option:: state_id

    :Description:
        Initial wetness condition on each land cover
    :SUEWS-related variables:
        `BldgsState <https://suews-docs.readthedocs.io/en/latest/input_files/Initial_Conditions/Above_ground_state.html#cmdoption-arg-bldgsstate>`_, `PavedState <https://suews-docs.readthedocs.io/en/latest/input_files/Initial_Conditions/Above_ground_state.html#cmdoption-arg-pavedstate>`_, `DecTrState <https://suews-docs.readthedocs.io/en/latest/input_files/Initial_Conditions/Above_ground_state.html#cmdoption-arg-dectrstate>`_, `EveTrState <https://suews-docs.readthedocs.io/en/latest/input_files/Initial_Conditions/Above_ground_state.html#cmdoption-arg-evetrstate>`_, `GrassState <https://suews-docs.readthedocs.io/en/latest/input_files/Initial_Conditions/Above_ground_state.html#cmdoption-arg-grassstate>`_, `BSoilState <https://suews-docs.readthedocs.io/en/latest/input_files/Initial_Conditions/Above_ground_state.html#cmdoption-arg-bsoilstate>`_, `WaterState <https://suews-docs.readthedocs.io/en/latest/input_files/Initial_Conditions/Above_ground_state.html#cmdoption-arg-waterstate>`_
    :Dimensionality:
        (7,)
    :Dimensionality remarks:
        7: Seven SUEWS land cover types in the order of [paved; buildings; evergreen trees/shrubs; deciduous trees/shrubs; grass; bare soil and water]: https://suews-docs.readthedocs.io/en/latest/introduction.html


.. option:: statelimit

    :Description:
        Upper limit to the surface state. [mm]
    :SUEWS-related variables:
        `StateLimit <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-statelimit>`_
    :Dimensionality:
        (7,)
    :Dimensionality remarks:
        7: Seven SUEWS land cover types in the order of [paved; buildings; evergreen trees/shrubs; deciduous trees/shrubs; grass; bare soil and water]: https://suews-docs.readthedocs.io/en/latest/introduction.html


.. option:: storageheatmethod

    :Description:
        Determines method for calculating storage heat flux ΔQS.
    :SUEWS-related variables:
        `StorageHeatMethod <https://suews-docs.readthedocs.io/en/latest/input_files/RunControl/scheme_options.html#cmdoption-arg-storageheatmethod>`_
    :Dimensionality:
        0
    :Dimensionality remarks:
        Scalar


.. option:: storedrainprm

    :Description:
        Coefficients used in drainage calculation.
    :SUEWS-related variables:
        `DrainageCoef1 <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-drainagecoef1>`_, `DrainageCoef2 <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-drainagecoef2>`_, `DrainageEq <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-drainageeq>`_, `StorageMax <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-storagemax>`_, `StorageMin <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-storagemin>`_
    :Dimensionality:
        (6, 7)
    :Dimensionality remarks:
        6: See variable description for specifics; 7: Seven SUEWS land cover types in the order of [paved; buildings; evergreen trees/shrubs; deciduous trees/shrubs; grass; bare soil and water]: https://suews-docs.readthedocs.io/en/latest/introduction.html


.. option:: surfacearea

    :Description:
        Area of the grid [ha].
    :SUEWS-related variables:
        `SurfaceArea <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-surfacearea>`_
    :Dimensionality:
        0
    :Dimensionality remarks:
        Scalar


.. option:: t_critic_cooling

    :Description:
        Critical cooling temperature.
    :SUEWS-related variables:
        `TCritic_Cooling_WD <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-tcritic-cooling-wd>`_, `TCritic_Cooling_WE <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-tcritic-cooling-we>`_
    :Dimensionality:
        (2,)
    :Dimensionality remarks:
        2: Weekday and Weekend


.. option:: t_critic_heating

    :Description:
        Critical heating temperature.
    :SUEWS-related variables:
        `TCritic_Heating_WD <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-tcritic-heating-wd>`_, `TCritic_Heating_WE <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-tcritic-heating-we>`_
    :Dimensionality:
        (2,)
    :Dimensionality remarks:
        2: Weekday and Weekend


.. option:: tau_a

    :Description:
        Time constant for snow albedo aging in cold snow [-]
    :SUEWS-related variables:
        `tau_a <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-tau-a>`_
    :Dimensionality:
        0
    :Dimensionality remarks:
        Scalar


.. option:: tau_f

    :Description:
        Time constant for snow albedo aging in melting snow [-]
    :SUEWS-related variables:
        `tau_f <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-tau-f>`_
    :Dimensionality:
        0
    :Dimensionality remarks:
        Scalar


.. option:: tau_r

    :Description:
        Time constant for snow density ageing [-]
    :SUEWS-related variables:
        `tau_r <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-tau-r>`_
    :Dimensionality:
        0
    :Dimensionality remarks:
        Scalar


.. option:: tempmeltfact

    :Description:
        Hourly temperature melt factor of snow [mm |K^-1| |h^-1|]
    :SUEWS-related variables:
        `TempMeltFactor <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-tempmeltfactor>`_
    :Dimensionality:
        0
    :Dimensionality remarks:
        Scalar


.. option:: th

    :Description:
        Upper air temperature limit [°C]
    :SUEWS-related variables:
        `TH <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-th>`_
    :Dimensionality:
        0
    :Dimensionality remarks:
        Scalar


.. option:: theta_bioco2

    :Description:
        The convexity of the curve at light saturation.
    :SUEWS-related variables:
        `theta <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-theta>`_
    :Dimensionality:
        (3,)
    :Dimensionality remarks:
        3: See variable description for specifics


.. option:: timezone

    :Description:
        Time zone [h] for site relative to UTC (east is positive). This should be set according to the times given in the meteorological forcing file(s).
    :SUEWS-related variables:
        `Timezone <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-timezone>`_
    :Dimensionality:
        0
    :Dimensionality remarks:
        Scalar


.. option:: tl

    :Description:
        Lower air temperature limit [°C]
    :SUEWS-related variables:
        `TL <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-tl>`_
    :Dimensionality:
        0
    :Dimensionality remarks:
        Scalar


.. option:: trafficrate

    :Description:
        Traffic rate used for CO2 flux calculation.
    :SUEWS-related variables:
        `TrafficRate_WD <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-trafficrate-wd>`_, `TrafficRate_WE <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-trafficrate-we>`_
    :Dimensionality:
        (2,)
    :Dimensionality remarks:
        2: Weekday and Weekend


.. option:: trafficunits

    :Description:
        Units for the traffic rate for the study area. Not used in v2018a.
    :SUEWS-related variables:
        `TrafficUnits <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-trafficunits>`_
    :Dimensionality:
        0
    :Dimensionality remarks:
        Scalar


.. option:: traffprof_24hr

    :Description:
        Hourly profile values used in traffic activity calculation.
    :SUEWS-related variables:
        `TraffProfWD <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-traffprofwd>`_, `TraffProfWE <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-traffprofwe>`_
    :Dimensionality:
        (24, 2)
    :Dimensionality remarks:
        24: hours of a day; 2: Weekday and Weekend


.. option:: tstep

    :Description:
        Specifies the model time step [s].
    :SUEWS-related variables:
        `Tstep <https://suews-docs.readthedocs.io/en/latest/input_files/RunControl/Time_related_options.html#cmdoption-arg-tstep>`_
    :Dimensionality:
        0
    :Dimensionality remarks:
        Scalar


.. option:: veg_type

    :Description:
        Internal use. Please DO NOT modify
    :SUEWS-related variables:
        None
    :Dimensionality:
        0
    :Dimensionality remarks:
        Scalar


.. option:: waterdist

    :Description:
        Fraction of water redistribution
    :SUEWS-related variables:
        `ToBSoil <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-tobsoil>`_, `ToBldgs <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-tobldgs>`_, `ToDecTr <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-todectr>`_, `ToEveTr <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-toevetr>`_, `ToGrass <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-tograss>`_, `ToPaved <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-topaved>`_, `ToRunoff <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-torunoff>`_, `ToSoilStore <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-tosoilstore>`_, `ToWater <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-towater>`_
    :Dimensionality:
        (8, 6)
    :Dimensionality remarks:
        8: Seven SUEWS land cover types and Runoff/SoilStore as water receiver; 6: SUEWS land cover types other than water as water contributors


.. option:: waterusemethod

    :Description:
        Defines how external water use is calculated.
    :SUEWS-related variables:
        `WaterUseMethod <https://suews-docs.readthedocs.io/en/latest/input_files/RunControl/scheme_options.html#cmdoption-arg-waterusemethod>`_
    :Dimensionality:
        0
    :Dimensionality remarks:
        Scalar


.. option:: wetthresh

    :Description:
        Depth of water which determines whether evaporation occurs from a partially wet or completely wet surface [mm].
    :SUEWS-related variables:
        `WetThreshold <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-wetthreshold>`_
    :Dimensionality:
        (7,)
    :Dimensionality remarks:
        7: Seven SUEWS land cover types in the order of [paved; buildings; evergreen trees/shrubs; deciduous trees/shrubs; grass; bare soil and water]: https://suews-docs.readthedocs.io/en/latest/introduction.html


.. option:: wuprofa_24hr

    :Description:
        Hourly profile values used in automatic irrigation.
    :SUEWS-related variables:
        `WaterUseProfAutoWD <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-wateruseprofautowd>`_, `WaterUseProfAutoWE <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-wateruseprofautowe>`_
    :Dimensionality:
        (24, 2)
    :Dimensionality remarks:
        24: hours of a day; 2: Weekday and Weekend


.. option:: wuprofm_24hr

    :Description:
        Hourly profile values used in manual irrigation.
    :SUEWS-related variables:
        `WaterUseProfManuWD <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-wateruseprofmanuwd>`_, `WaterUseProfManuWE <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-wateruseprofmanuwe>`_
    :Dimensionality:
        (24, 2)
    :Dimensionality remarks:
        24: hours of a day; 2: Weekday and Weekend


.. option:: z

    :Description:
        Measurement height [m].
    :SUEWS-related variables:
        `z <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-z>`_
    :Dimensionality:
        0
    :Dimensionality remarks:
        Scalar


.. option:: z0m_in

    :Description:
        Roughness length for momentum [m]
    :SUEWS-related variables:
        `z0 <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-z0>`_
    :Dimensionality:
        0
    :Dimensionality remarks:
        Scalar


.. option:: zdm_in

    :Description:
        Zero-plane displacement [m]
    :SUEWS-related variables:
        `zd <https://suews-docs.readthedocs.io/en/latest/input_files/SUEWS_SiteInfo/Input_Options.html#cmdoption-arg-zd>`_
    :Dimensionality:
        0
    :Dimensionality remarks:
        Scalar
