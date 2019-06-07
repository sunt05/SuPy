
.. _df_state_var:

``df_state`` variables
============================



.. note:: Data structure of ``df_state`` is explained :ref:`here </data-structure/supy-io.ipynb#df_state_init:-model-initial-states>`.

.. option:: aerodynamicresistancemethod

    :Description:
        Internal use. Please DO NOT modify
    :Dimensionality:
        0
    :Dimensionality Remarks:
        Scalar
    :SUEWS-related variables:
        None


.. option:: ah_min

    :Description:
        Minimum QF values.
    :Dimensionality:
        (2,)
    :Dimensionality Remarks:
        2: {Weekday, Weekend}
    :SUEWS-related variables:
        :option:`AHMin_WD <suews:AHMin_WD>`, :option:`AHMin_WE <suews:AHMin_WE>`


.. option:: ah_slope_cooling

    :Description:
        Cooling slope of QF calculation.
    :Dimensionality:
        (2,)
    :Dimensionality Remarks:
        2: {Weekday, Weekend}
    :SUEWS-related variables:
        :option:`AHSlope_Cooling_WD <suews:AHSlope_Cooling_WD>`, :option:`AHSlope_Cooling_WE <suews:AHSlope_Cooling_WE>`


.. option:: ah_slope_heating

    :Description:
        Heating slope of QF calculation.
    :Dimensionality:
        (2,)
    :Dimensionality Remarks:
        2: {Weekday, Weekend}
    :SUEWS-related variables:
        :option:`AHSlope_Heating_WD <suews:AHSlope_Heating_WD>`, :option:`AHSlope_Heating_WE <suews:AHSlope_Heating_WE>`


.. option:: ahprof_24hr

    :Description:
        Hourly profile values used in energy use calculation.
    :Dimensionality:
        (24, 2)
    :Dimensionality Remarks:
        24: hours of a day

        2: {Weekday, Weekend}
    :SUEWS-related variables:
        :option:`EnergyUseProfWD <suews:EnergyUseProfWD>`, :option:`EnergyUseProfWE <suews:EnergyUseProfWE>`


.. option:: alb

    :Description:
        Effective surface albedo (middle of the day value) for summertime.
    :Dimensionality:
        (7,)
    :Dimensionality Remarks:
        7: { :term:`Paved`, :term:`Bldgs`, :term:`EveTr`, :term:`DecTr`, :term:`Grass`, :term:`BSoil`, :term:`Water`}
    :SUEWS-related variables:
        :option:`AlbedoMax <suews:AlbedoMax>`


.. option:: albdectr_id

    :Description:
        Albedo of deciduous surface `DecTr` on day 0 of run
    :Dimensionality:
        0
    :Dimensionality Remarks:
        Scalar
    :SUEWS-related variables:
        :option:`albDecTr0 <suews:albDecTr0>`


.. option:: albevetr_id

    :Description:
        Albedo of evergreen surface `EveTr` on day 0 of run
    :Dimensionality:
        0
    :Dimensionality Remarks:
        Scalar
    :SUEWS-related variables:
        :option:`albEveTr0 <suews:albEveTr0>`


.. option:: albgrass_id

    :Description:
        Albedo of grass surface `Grass` on day 0 of run
    :Dimensionality:
        0
    :Dimensionality Remarks:
        Scalar
    :SUEWS-related variables:
        :option:`albGrass0 <suews:albGrass0>`


.. option:: albmax_dectr

    :Description:
        Effective surface albedo (middle of the day value) for summertime.
    :Dimensionality:
        0
    :Dimensionality Remarks:
        Scalar
    :SUEWS-related variables:
        :option:`AlbedoMax <suews:AlbedoMax>`


.. option:: albmax_evetr

    :Description:
        Effective surface albedo (middle of the day value) for summertime.
    :Dimensionality:
        0
    :Dimensionality Remarks:
        Scalar
    :SUEWS-related variables:
        :option:`AlbedoMax <suews:AlbedoMax>`


.. option:: albmax_grass

    :Description:
        Effective surface albedo (middle of the day value) for summertime.
    :Dimensionality:
        0
    :Dimensionality Remarks:
        Scalar
    :SUEWS-related variables:
        :option:`AlbedoMax <suews:AlbedoMax>`


.. option:: albmin_dectr

    :Description:
        Effective surface albedo (middle of the day value) for wintertime (not including snow).
    :Dimensionality:
        0
    :Dimensionality Remarks:
        Scalar
    :SUEWS-related variables:
        :option:`AlbedoMin <suews:AlbedoMin>`


.. option:: albmin_evetr

    :Description:
        Effective surface albedo (middle of the day value) for wintertime (not including snow).
    :Dimensionality:
        0
    :Dimensionality Remarks:
        Scalar
    :SUEWS-related variables:
        :option:`AlbedoMin <suews:AlbedoMin>`


.. option:: albmin_grass

    :Description:
        Effective surface albedo (middle of the day value) for wintertime (not including snow).
    :Dimensionality:
        0
    :Dimensionality Remarks:
        Scalar
    :SUEWS-related variables:
        :option:`AlbedoMin <suews:AlbedoMin>`


.. option:: alpha_bioco2

    :Description:
        The mean apparent ecosystem quantum. Represents the initial slope of the light-response curve.
    :Dimensionality:
        (3,)
    :Dimensionality Remarks:
        3: { :term:`EveTr`, :term:`DecTr`, :term:`Grass`}
    :SUEWS-related variables:
        :option:`alpha <suews:alpha>`


.. option:: alpha_enh_bioco2

    :Description:
        Part of the `alpha` coefficient related to the fraction of vegetation.
    :Dimensionality:
        (3,)
    :Dimensionality Remarks:
        3: { :term:`EveTr`, :term:`DecTr`, :term:`Grass`}
    :SUEWS-related variables:
        :option:`alpha_enh <suews:alpha_enh>`


.. option:: alt

    :Description:
        Used for both the radiation and water flow between grids.
    :Dimensionality:
        0
    :Dimensionality Remarks:
        Scalar
    :SUEWS-related variables:
        :option:`Alt <suews:Alt>`


.. option:: baset

    :Description:
        Base Temperature for initiating growing degree days (GDD) for leaf growth. [°C]
    :Dimensionality:
        (3,)
    :Dimensionality Remarks:
        3: { :term:`EveTr`, :term:`DecTr`, :term:`Grass`}
    :SUEWS-related variables:
        :option:`BaseT <suews:BaseT>`


.. option:: basete

    :Description:
        Base temperature for initiating sensesance degree days (SDD) for leaf off. [°C]
    :Dimensionality:
        (3,)
    :Dimensionality Remarks:
        3: { :term:`EveTr`, :term:`DecTr`, :term:`Grass`}
    :SUEWS-related variables:
        :option:`BaseTe <suews:BaseTe>`


.. option:: basethdd

    :Description:
        Base temperature for heating degree days [°C]
    :Dimensionality:
        0
    :Dimensionality Remarks:
        Scalar
    :SUEWS-related variables:
        :option:`BaseTHDD <suews:BaseTHDD>`


.. option:: beta_bioco2

    :Description:
        The light-saturated gross photosynthesis of the canopy. [umol |m^-2| |s^-1| ]
    :Dimensionality:
        (3,)
    :Dimensionality Remarks:
        3: { :term:`EveTr`, :term:`DecTr`, :term:`Grass`}
    :SUEWS-related variables:
        :option:`beta <suews:beta>`


.. option:: beta_enh_bioco2

    :Description:
        Part of the `beta` coefficient related to the fraction of vegetation.
    :Dimensionality:
        (3,)
    :Dimensionality Remarks:
        3: { :term:`EveTr`, :term:`DecTr`, :term:`Grass`}
    :SUEWS-related variables:
        :option:`beta_enh <suews:beta_enh>`


.. option:: bldgh

    :Description:
        Mean building height [m]
    :Dimensionality:
        0
    :Dimensionality Remarks:
        Scalar
    :SUEWS-related variables:
        :option:`H_Bldgs <suews:H_Bldgs>`


.. option:: capmax_dec

    :Description:
        Maximum water storage capacity for upper surfaces (i.e. canopy)
    :Dimensionality:
        0
    :Dimensionality Remarks:
        Scalar
    :SUEWS-related variables:
        :option:`StorageMax <suews:StorageMax>`


.. option:: capmin_dec

    :Description:
        Minimum water storage capacity for upper surfaces (i.e. canopy).
    :Dimensionality:
        0
    :Dimensionality Remarks:
        Scalar
    :SUEWS-related variables:
        :option:`StorageMin <suews:StorageMin>`


.. option:: chanohm

    :Description:
        Bulk transfer coefficient for this surface to use in AnOHM [-]
    :Dimensionality:
        (7,)
    :Dimensionality Remarks:
        7: { :term:`Paved`, :term:`Bldgs`, :term:`EveTr`, :term:`DecTr`, :term:`Grass`, :term:`BSoil`, :term:`Water`}
    :SUEWS-related variables:
        :option:`AnOHM_Ch <suews:AnOHM_Ch>`


.. option:: cpanohm

    :Description:
        Volumetric heat capacity for this surface to use in AnOHM [J |m^-3|]
    :Dimensionality:
        (7,)
    :Dimensionality Remarks:
        7: { :term:`Paved`, :term:`Bldgs`, :term:`EveTr`, :term:`DecTr`, :term:`Grass`, :term:`BSoil`, :term:`Water`}
    :SUEWS-related variables:
        :option:`AnOHM_Cp <suews:AnOHM_Cp>`


.. option:: crwmax

    :Description:
        Maximum water holding capacity of snow [mm]
    :Dimensionality:
        0
    :Dimensionality Remarks:
        Scalar
    :SUEWS-related variables:
        :option:`CRWMax <suews:CRWMax>`


.. option:: crwmin

    :Description:
        Minimum water holding capacity of snow [mm]
    :Dimensionality:
        0
    :Dimensionality Remarks:
        Scalar
    :SUEWS-related variables:
        :option:`CRWMin <suews:CRWMin>`


.. option:: daywat

    :Description:
        Irrigation flag: 1 for on and 0 for off.
    :Dimensionality:
        (7,)
    :Dimensionality Remarks:
        7: {Sunday, Monday, Tuesday, Wednesday, Thursday, Friday, Saturday}
    :SUEWS-related variables:
        :option:`DayWat(1) <suews:DayWat(1)>`, :option:`DayWat(2) <suews:DayWat(2)>`, :option:`DayWat(3) <suews:DayWat(3)>`, :option:`DayWat(4) <suews:DayWat(4)>`, :option:`DayWat(5) <suews:DayWat(5)>`, :option:`DayWat(6) <suews:DayWat(6)>`, :option:`DayWat(7) <suews:DayWat(7)>`


.. option:: daywatper

    :Description:
        Fraction of properties using irrigation for each day of a week.
    :Dimensionality:
        (7,)
    :Dimensionality Remarks:
        7: {Sunday, Monday, Tuesday, Wednesday, Thursday, Friday, Saturday}
    :SUEWS-related variables:
        :option:`DayWatPer(1) <suews:DayWatPer(1)>`, :option:`DayWatPer(2) <suews:DayWatPer(2)>`, :option:`DayWatPer(3) <suews:DayWatPer(3)>`, :option:`DayWatPer(4) <suews:DayWatPer(4)>`, :option:`DayWatPer(5) <suews:DayWatPer(5)>`, :option:`DayWatPer(6) <suews:DayWatPer(6)>`, :option:`DayWatPer(7) <suews:DayWatPer(7)>`


.. option:: decidcap_id

    :Description:
        Storage capacity of deciduous surface `DecTr` on day 0 of run.
    :Dimensionality:
        0
    :Dimensionality Remarks:
        Scalar
    :SUEWS-related variables:
        :option:`decidCap0 <suews:decidCap0>`


.. option:: dectreeh

    :Description:
        Mean height of deciduous trees [m]
    :Dimensionality:
        0
    :Dimensionality Remarks:
        Scalar
    :SUEWS-related variables:
        :option:`H_DecTr <suews:H_DecTr>`


.. option:: diagnose

    :Description:
        Internal use. Please DO NOT modify
    :Dimensionality:
        0
    :Dimensionality Remarks:
        Scalar
    :SUEWS-related variables:
        None


.. option:: diagqn

    :Description:
        Internal use. Please DO NOT modify
    :Dimensionality:
        0
    :Dimensionality Remarks:
        Scalar
    :SUEWS-related variables:
        None


.. option:: diagqs

    :Description:
        Internal use. Please DO NOT modify
    :Dimensionality:
        0
    :Dimensionality Remarks:
        Scalar
    :SUEWS-related variables:
        None


.. option:: drainrt

    :Description:
        Drainage rate of bucket for LUMPS [mm |h^-1|]
    :Dimensionality:
        0
    :Dimensionality Remarks:
        Scalar
    :SUEWS-related variables:
        :option:`LUMPS_DrRate <suews:LUMPS_DrRate>`


.. option:: ef_umolco2perj

    :Description:
        Emission factor for fuels used for building heating.
    :Dimensionality:
        0
    :Dimensionality Remarks:
        Scalar
    :SUEWS-related variables:
        :option:`EF_umolCO2perJ <suews:EF_umolCO2perJ>`


.. option:: emis

    :Description:
        Effective surface emissivity.
    :Dimensionality:
        (7,)
    :Dimensionality Remarks:
        7: { :term:`Paved`, :term:`Bldgs`, :term:`EveTr`, :term:`DecTr`, :term:`Grass`, :term:`BSoil`, :term:`Water`}
    :SUEWS-related variables:
        :option:`Emissivity <suews:Emissivity>`


.. option:: emissionsmethod

    :Description:
        Determines method for QF calculation.
    :Dimensionality:
        0
    :Dimensionality Remarks:
        Scalar
    :SUEWS-related variables:
        :option:`EmissionsMethod <suews:EmissionsMethod>`


.. option:: enddls

    :Description:
        End of the day light savings [DOY]
    :Dimensionality:
        0
    :Dimensionality Remarks:
        Scalar
    :SUEWS-related variables:
        :option:`EndDLS <suews:EndDLS>`


.. option:: enef_v_jkm

    :Description:
        Emission factor for heat [J k|m^-1|].
    :Dimensionality:
        0
    :Dimensionality Remarks:
        Scalar
    :SUEWS-related variables:
        :option:`EnEF_v_Jkm <suews:EnEF_v_Jkm>`


.. option:: evapmethod

    :Description:
        Internal use. Please DO NOT modify
    :Dimensionality:
        0
    :Dimensionality Remarks:
        Scalar
    :SUEWS-related variables:
        None


.. option:: evetreeh

    :Description:
        Mean height of evergreen trees [m]
    :Dimensionality:
        0
    :Dimensionality Remarks:
        Scalar
    :SUEWS-related variables:
        :option:`H_EveTr <suews:H_EveTr>`


.. option:: faibldg

    :Description:
        Frontal area index for buildings [-]
    :Dimensionality:
        0
    :Dimensionality Remarks:
        Scalar
    :SUEWS-related variables:
        :option:`FAI_Bldgs <suews:FAI_Bldgs>`


.. option:: faidectree

    :Description:
        Frontal area index for deciduous trees [-]
    :Dimensionality:
        0
    :Dimensionality Remarks:
        Scalar
    :SUEWS-related variables:
        :option:`FAI_DecTr <suews:FAI_DecTr>`


.. option:: faievetree

    :Description:
        Frontal area index for evergreen trees [-]
    :Dimensionality:
        0
    :Dimensionality Remarks:
        Scalar
    :SUEWS-related variables:
        :option:`FAI_EveTr <suews:FAI_EveTr>`


.. option:: faut

    :Description:
        Fraction of irrigated area that is irrigated using automated systems
    :Dimensionality:
        0
    :Dimensionality Remarks:
        Scalar
    :SUEWS-related variables:
        :option:`Faut <suews:Faut>`


.. option:: fcef_v_kgkm

    :Description:
        CO2 emission factor [kg |km^-1|]
    :Dimensionality:
        0
    :Dimensionality Remarks:
        Scalar
    :SUEWS-related variables:
        :option:`FcEF_v_kgkm <suews:FcEF_v_kgkm>`


.. option:: flowchange

    :Description:
        Difference in input and output flows for water surface [mm |h^-1|]
    :Dimensionality:
        0
    :Dimensionality Remarks:
        Scalar
    :SUEWS-related variables:
        :option:`FlowChange <suews:FlowChange>`


.. option:: frfossilfuel_heat

    :Description:
        Fraction of fossil fuels used for building heating [-]
    :Dimensionality:
        0
    :Dimensionality Remarks:
        Scalar
    :SUEWS-related variables:
        :option:`FrFossilFuel_Heat <suews:FrFossilFuel_Heat>`


.. option:: frfossilfuel_nonheat

    :Description:
        Fraction of fossil fuels used for building energy use [-]
    :Dimensionality:
        0
    :Dimensionality Remarks:
        Scalar
    :SUEWS-related variables:
        :option:`FrFossilFuel_NonHeat <suews:FrFossilFuel_NonHeat>`


.. option:: g1

    :Description:
        Related to maximum surface conductance [mm |s^-1|]
    :Dimensionality:
        0
    :Dimensionality Remarks:
        Scalar
    :SUEWS-related variables:
        :option:`G1 <suews:G1>`


.. option:: g2

    :Description:
        Related to Kdown dependence [W |m^-2|]
    :Dimensionality:
        0
    :Dimensionality Remarks:
        Scalar
    :SUEWS-related variables:
        :option:`G2 <suews:G2>`


.. option:: g3

    :Description:
        Related to VPD dependence [units depend on `gsModel`]
    :Dimensionality:
        0
    :Dimensionality Remarks:
        Scalar
    :SUEWS-related variables:
        :option:`G3 <suews:G3>`


.. option:: g4

    :Description:
        Related to VPD dependence [units depend on `gsModel`]
    :Dimensionality:
        0
    :Dimensionality Remarks:
        Scalar
    :SUEWS-related variables:
        :option:`G4 <suews:G4>`


.. option:: g5

    :Description:
        Related to temperature dependence [°C]
    :Dimensionality:
        0
    :Dimensionality Remarks:
        Scalar
    :SUEWS-related variables:
        :option:`G5 <suews:G5>`


.. option:: g6

    :Description:
        Related to soil moisture dependence [|mm^-1|]
    :Dimensionality:
        0
    :Dimensionality Remarks:
        Scalar
    :SUEWS-related variables:
        :option:`G6 <suews:G6>`


.. option:: gddfull

    :Description:
        The growing degree days (GDD) needed for full capacity of the leaf area index (LAI) [°C].
    :Dimensionality:
        (3,)
    :Dimensionality Remarks:
        3: { :term:`EveTr`, :term:`DecTr`, :term:`Grass`}
    :SUEWS-related variables:
        :option:`GDDFull <suews:GDDFull>`


.. option:: gsmodel

    :Description:
        Formulation choice for conductance calculation.
    :Dimensionality:
        0
    :Dimensionality Remarks:
        Scalar
    :SUEWS-related variables:
        :option:`gsModel <suews:gsModel>`


.. option:: humactivity_24hr

    :Description:
        Hourly profile values used in human activity calculation.
    :Dimensionality:
        (24, 2)
    :Dimensionality Remarks:
        24: hours of a day

        2: {Weekday, Weekend}
    :SUEWS-related variables:
        :option:`ActivityProfWD <suews:ActivityProfWD>`, :option:`ActivityProfWE <suews:ActivityProfWE>`


.. option:: ie_a

    :Description:
        Coefficient for automatic irrigation model.
    :Dimensionality:
        (3,)
    :Dimensionality Remarks:
        3: { :term:`EveTr`, :term:`DecTr`, :term:`Grass`}
    :SUEWS-related variables:
        :option:`Ie_a1 <suews:Ie_a1>`, :option:`Ie_a2 <suews:Ie_a2>`, :option:`Ie_a3 <suews:Ie_a3>`


.. option:: ie_end

    :Description:
        Day when irrigation ends [DOY]
    :Dimensionality:
        0
    :Dimensionality Remarks:
        Scalar
    :SUEWS-related variables:
        :option:`Ie_end <suews:Ie_end>`


.. option:: ie_m

    :Description:
        Coefficient for manual irrigation model.
    :Dimensionality:
        (3,)
    :Dimensionality Remarks:
        3: { :term:`EveTr`, :term:`DecTr`, :term:`Grass`}
    :SUEWS-related variables:
        :option:`Ie_m1 <suews:Ie_m1>`, :option:`Ie_m2 <suews:Ie_m2>`, :option:`Ie_m3 <suews:Ie_m3>`


.. option:: ie_start

    :Description:
        Day when irrigation starts [DOY]
    :Dimensionality:
        0
    :Dimensionality Remarks:
        Scalar
    :SUEWS-related variables:
        :option:`Ie_start <suews:Ie_start>`


.. option:: internalwateruse_h

    :Description:
        Internal water use [mm |h^-1|]
    :Dimensionality:
        0
    :Dimensionality Remarks:
        Scalar
    :SUEWS-related variables:
        :option:`InternalWaterUse <suews:InternalWaterUse>`


.. option:: irrfracconif

    :Description:
        Fraction of evergreen trees that are irrigated [-]
    :Dimensionality:
        0
    :Dimensionality Remarks:
        Scalar
    :SUEWS-related variables:
        :option:`IrrFr_EveTr <suews:IrrFr_EveTr>`


.. option:: irrfracdecid

    :Description:
        Fraction of deciduous trees that are irrigated [-]
    :Dimensionality:
        0
    :Dimensionality Remarks:
        Scalar
    :SUEWS-related variables:
        :option:`IrrFr_DecTr <suews:IrrFr_DecTr>`


.. option:: irrfracgrass

    :Description:
        Fraction of `Grass` that is irrigated [-]
    :Dimensionality:
        0
    :Dimensionality Remarks:
        Scalar
    :SUEWS-related variables:
        :option:`IrrFr_Grass <suews:IrrFr_Grass>`


.. option:: kkanohm

    :Description:
        Thermal conductivity for this surface to use in AnOHM [W m |K^-1|]
    :Dimensionality:
        (7,)
    :Dimensionality Remarks:
        7: { :term:`Paved`, :term:`Bldgs`, :term:`EveTr`, :term:`DecTr`, :term:`Grass`, :term:`BSoil`, :term:`Water`}
    :SUEWS-related variables:
        :option:`AnOHM_Kk <suews:AnOHM_Kk>`


.. option:: kmax

    :Description:
        Maximum incoming shortwave radiation [W |m^-2|]
    :Dimensionality:
        0
    :Dimensionality Remarks:
        Scalar
    :SUEWS-related variables:
        :option:`Kmax <suews:Kmax>`


.. option:: lai_id

    :Description:
        Initial LAI values.
    :Dimensionality:
        (3,)
    :Dimensionality Remarks:
        3: { :term:`EveTr`, :term:`DecTr`, :term:`Grass`}
    :SUEWS-related variables:
        :option:`LAIinitialDecTr <suews:LAIinitialDecTr>`, :option:`LAIinitialEveTr <suews:LAIinitialEveTr>`, :option:`LAIinitialGrass <suews:LAIinitialGrass>`


.. option:: laicalcyes

    :Description:
        Internal use. Please DO NOT modify
    :Dimensionality:
        0
    :Dimensionality Remarks:
        Scalar
    :SUEWS-related variables:
        None


.. option:: laimax

    :Description:
        full leaf-on summertime value
    :Dimensionality:
        (3,)
    :Dimensionality Remarks:
        3: { :term:`EveTr`, :term:`DecTr`, :term:`Grass`}
    :SUEWS-related variables:
        :option:`LAIMax <suews:LAIMax>`


.. option:: laimin

    :Description:
        leaf-off wintertime value
    :Dimensionality:
        (3,)
    :Dimensionality Remarks:
        3: { :term:`EveTr`, :term:`DecTr`, :term:`Grass`}
    :SUEWS-related variables:
        :option:`LAIMin <suews:LAIMin>`


.. option:: laipower

    :Description:
        parameters required by LAI calculation.
    :Dimensionality:
        (4, 3)
    :Dimensionality Remarks:
        4: {`LeafGrowthPower1`, `LeafGrowthPower2`, `LeafOffPower1`, `LeafOffPower2`}

        3: { :term:`EveTr`, :term:`DecTr`, :term:`Grass`}
    :SUEWS-related variables:
        :option:`LeafGrowthPower1 <suews:LeafGrowthPower1>`, :option:`LeafGrowthPower2 <suews:LeafGrowthPower2>`, :option:`LeafOffPower1 <suews:LeafOffPower1>`, :option:`LeafOffPower2 <suews:LeafOffPower2>`


.. option:: laitype

    :Description:
        LAI calculation choice.
    :Dimensionality:
        (3,)
    :Dimensionality Remarks:
        3: { :term:`EveTr`, :term:`DecTr`, :term:`Grass`}
    :SUEWS-related variables:
        :option:`LAIEq <suews:LAIEq>`


.. option:: lat

    :Description:
        Latitude [deg].
    :Dimensionality:
        0
    :Dimensionality Remarks:
        Scalar
    :SUEWS-related variables:
        :option:`lat <suews:lat>`


.. option:: lng

    :Description:
        longitude [deg]
    :Dimensionality:
        0
    :Dimensionality Remarks:
        Scalar
    :SUEWS-related variables:
        :option:`lng <suews:lng>`


.. option:: maxconductance

    :Description:
        The maximum conductance of each vegetation or surface type. [mm |s^-1|]
    :Dimensionality:
        (3,)
    :Dimensionality Remarks:
        3: { :term:`EveTr`, :term:`DecTr`, :term:`Grass`}
    :SUEWS-related variables:
        :option:`MaxConductance <suews:MaxConductance>`


.. option:: maxqfmetab

    :Description:
        Maximum value for human heat emission. [W |m^-2|]
    :Dimensionality:
        0
    :Dimensionality Remarks:
        Scalar
    :SUEWS-related variables:
        :option:`MaxQFMetab <suews:MaxQFMetab>`


.. option:: min_res_bioco2

    :Description:
        Minimum soil respiration rate (for cold-temperature limit) [umol |m^-2| |s^-1|].
    :Dimensionality:
        (3,)
    :Dimensionality Remarks:
        3: { :term:`EveTr`, :term:`DecTr`, :term:`Grass`}
    :SUEWS-related variables:
        :option:`min_respi <suews:min_respi>`


.. option:: minqfmetab

    :Description:
        Minimum value for human heat emission. [W |m^-2|]
    :Dimensionality:
        0
    :Dimensionality Remarks:
        Scalar
    :SUEWS-related variables:
        :option:`MinQFMetab <suews:MinQFMetab>`


.. option:: narp_emis_snow

    :Description:
        Effective surface emissivity.
    :Dimensionality:
        0
    :Dimensionality Remarks:
        Scalar
    :SUEWS-related variables:
        :option:`Emissivity <suews:Emissivity>`


.. option:: narp_trans_site

    :Description:
        Atmospheric transmissivity for NARP [-]
    :Dimensionality:
        0
    :Dimensionality Remarks:
        Scalar
    :SUEWS-related variables:
        :option:`NARP_Trans <suews:NARP_Trans>`


.. option:: netradiationmethod

    :Description:
        Determines method for calculation of radiation fluxes.
    :Dimensionality:
        0
    :Dimensionality Remarks:
        Scalar
    :SUEWS-related variables:
        :option:`NetRadiationMethod <suews:NetRadiationMethod>`


.. option:: ohm_coef

    :Description:
        Coefficients for OHM calculation.
    :Dimensionality:
        (8, 4, 3)
    :Dimensionality Remarks:
        8: { :term:`Paved`, :term:`Bldgs`, :term:`EveTr`, :term:`DecTr`, :term:`Grass`, :term:`BSoil`, :term:`Water`, one extra land cover type (currently NOT used)}

        4: {SummerWet, SummerDry, WinterWet, WinterDry}

        3: {a1, a2, a3}
    :SUEWS-related variables:
        :option:`a1 <suews:a1>`, :option:`a2 <suews:a2>`, :option:`a3 <suews:a3>`


.. option:: ohm_threshsw

    :Description:
        Temperature threshold determining whether summer/winter OHM coefficients are applied [°C]
    :Dimensionality:
        (8,)
    :Dimensionality Remarks:
        8: { :term:`Paved`, :term:`Bldgs`, :term:`EveTr`, :term:`DecTr`, :term:`Grass`, :term:`BSoil`, :term:`Water`, one extra land cover type (currently NOT used)} 
    :SUEWS-related variables:
        :option:`OHMThresh_SW <suews:OHMThresh_SW>`


.. option:: ohm_threshwd

    :Description:
        Soil moisture threshold determining whether wet/dry OHM coefficients are applied [-]
    :Dimensionality:
        (8,)
    :Dimensionality Remarks:
        8: { :term:`Paved`, :term:`Bldgs`, :term:`EveTr`, :term:`DecTr`, :term:`Grass`, :term:`BSoil`, :term:`Water`, one extra land cover type (currently NOT used)} 
    :SUEWS-related variables:
        :option:`OHMThresh_WD <suews:OHMThresh_WD>`


.. option:: ohmincqf

    :Description:
        Determines whether the storage heat flux calculation uses |Qstar| or ( |Qstar| +QF).
    :Dimensionality:
        0
    :Dimensionality Remarks:
        Scalar
    :SUEWS-related variables:
        :option:`OHMIncQF <suews:OHMIncQF>`


.. option:: pipecapacity

    :Description:
        Storage capacity of pipes [mm]
    :Dimensionality:
        0
    :Dimensionality Remarks:
        Scalar
    :SUEWS-related variables:
        :option:`PipeCapacity <suews:PipeCapacity>`


.. option:: popdensdaytime

    :Description:
        Daytime population density (i.e. workers, tourists) [people |ha^-1|]
    :Dimensionality:
        0
    :Dimensionality Remarks:
        Scalar
    :SUEWS-related variables:
        :option:`PopDensDay <suews:PopDensDay>`


.. option:: popdensnighttime

    :Description:
        Night-time population density (i.e. residents) [people |ha^-1|]
    :Dimensionality:
        0
    :Dimensionality Remarks:
        Scalar
    :SUEWS-related variables:
        :option:`PopDensNight <suews:PopDensNight>`


.. option:: popprof_24hr

    :Description:
        Hourly profile values used in dynamic population estimation.
    :Dimensionality:
        (24, 2)
    :Dimensionality Remarks:
        24: hours of a day

        2: {Weekday, Weekend}
    :SUEWS-related variables:
        :option:`PopProfWD <suews:PopProfWD>`, :option:`PopProfWE <suews:PopProfWE>`


.. option:: pormax_dec

    :Description:
        full leaf-on summertime value Used only for `DecTr` (can affect roughness calculation)
    :Dimensionality:
        0
    :Dimensionality Remarks:
        Scalar
    :SUEWS-related variables:
        :option:`PorosityMax <suews:PorosityMax>`


.. option:: pormin_dec

    :Description:
        leaf-off wintertime value Used only for `DecTr` (can affect roughness calculation)
    :Dimensionality:
        0
    :Dimensionality Remarks:
        Scalar
    :SUEWS-related variables:
        :option:`PorosityMin <suews:PorosityMin>`


.. option:: porosity_id

    :Description:
        Porosity of deciduous vegetation on day 0 of run.
    :Dimensionality:
        0
    :Dimensionality Remarks:
        Scalar
    :SUEWS-related variables:
        :option:`porosity0 <suews:porosity0>`


.. option:: preciplimit

    :Description:
        Limit for hourly snowfall when the ground is fully covered with snow [mm]
    :Dimensionality:
        0
    :Dimensionality Remarks:
        Scalar
    :SUEWS-related variables:
        :option:`PrecipLimSnow <suews:PrecipLimSnow>`


.. option:: preciplimitalb

    :Description:
        Limit for hourly precipitation when the ground is fully covered with snow. Then snow albedo is reset to AlbedoMax [mm]
    :Dimensionality:
        0
    :Dimensionality Remarks:
        Scalar
    :SUEWS-related variables:
        :option:`PrecipLimAlb <suews:PrecipLimAlb>`


.. option:: qf0_beu

    :Description:
        Building energy use [W |m^-2|]
    :Dimensionality:
        (2,)
    :Dimensionality Remarks:
        2: {Weekday, Weekend}
    :SUEWS-related variables:
        :option:`QF0_BEU_WD <suews:QF0_BEU_WD>`, :option:`QF0_BEU_WE <suews:QF0_BEU_WE>`


.. option:: qf_a

    :Description:
        Base value for QF calculation.
    :Dimensionality:
        (2,)
    :Dimensionality Remarks:
        2: {Weekday, Weekend}
    :SUEWS-related variables:
        :option:`QF_A_WD <suews:QF_A_WD>`, :option:`QF_A_WE <suews:QF_A_WE>`


.. option:: qf_b

    :Description:
        Parameter related to heating degree days.
    :Dimensionality:
        (2,)
    :Dimensionality Remarks:
        2: {Weekday, Weekend}
    :SUEWS-related variables:
        :option:`QF_B_WD <suews:QF_B_WD>`, :option:`QF_B_WE <suews:QF_B_WE>`


.. option:: qf_c

    :Description:
        Parameter related to heating degree days.
    :Dimensionality:
        (2,)
    :Dimensionality Remarks:
        2: {Weekday, Weekend}
    :SUEWS-related variables:
        :option:`QF_C_WD <suews:QF_C_WD>`, :option:`QF_C_WE <suews:QF_C_WE>`


.. option:: radmeltfact

    :Description:
        Hourly radiation melt factor of snow [mm |w^-1| |h^-1|]
    :Dimensionality:
        0
    :Dimensionality Remarks:
        Scalar
    :SUEWS-related variables:
        :option:`RadMeltFactor <suews:RadMeltFactor>`


.. option:: raincover

    :Description:
        Limit when surface totally covered with water for LUMPS [mm]
    :Dimensionality:
        0
    :Dimensionality Remarks:
        Scalar
    :SUEWS-related variables:
        :option:`LUMPS_Cover <suews:LUMPS_Cover>`


.. option:: rainmaxres

    :Description:
        Maximum water bucket reservoir [mm] Used for LUMPS surface wetness control.
    :Dimensionality:
        0
    :Dimensionality Remarks:
        Scalar
    :SUEWS-related variables:
        :option:`LUMPS_MaxRes <suews:LUMPS_MaxRes>`


.. option:: resp_a

    :Description:
        Respiration coefficient a.
    :Dimensionality:
        (3,)
    :Dimensionality Remarks:
        3: { :term:`EveTr`, :term:`DecTr`, :term:`Grass`}
    :SUEWS-related variables:
        :option:`resp_a <suews:resp_a>`


.. option:: resp_b

    :Description:
        Respiration coefficient b - related to air temperature dependency.
    :Dimensionality:
        (3,)
    :Dimensionality Remarks:
        3: { :term:`EveTr`, :term:`DecTr`, :term:`Grass`}
    :SUEWS-related variables:
        :option:`resp_b <suews:resp_b>`


.. option:: roughlenheatmethod

    :Description:
        Determines method for calculating roughness length for heat.
    :Dimensionality:
        0
    :Dimensionality Remarks:
        Scalar
    :SUEWS-related variables:
        :option:`RoughLenHeatMethod <suews:RoughLenHeatMethod>`


.. option:: roughlenmommethod

    :Description:
        Determines how aerodynamic roughness length (z0m) and zero displacement height (zdm) are calculated.
    :Dimensionality:
        0
    :Dimensionality Remarks:
        Scalar
    :SUEWS-related variables:
        :option:`RoughLenMomMethod <suews:RoughLenMomMethod>`


.. option:: runofftowater

    :Description:
        Fraction of above-ground runoff flowing to water surface during flooding [-]
    :Dimensionality:
        0
    :Dimensionality Remarks:
        Scalar
    :SUEWS-related variables:
        :option:`RunoffToWater <suews:RunoffToWater>`


.. option:: s1

    :Description:
        A parameter related to soil moisture dependence [-]
    :Dimensionality:
        0
    :Dimensionality Remarks:
        Scalar
    :SUEWS-related variables:
        :option:`S1 <suews:S1>`


.. option:: s2

    :Description:
        A parameter related to soil moisture dependence [mm]
    :Dimensionality:
        0
    :Dimensionality Remarks:
        Scalar
    :SUEWS-related variables:
        :option:`S2 <suews:S2>`


.. option:: sathydraulicconduct

    :Description:
        Hydraulic conductivity for saturated soil [mm |s^-1|]
    :Dimensionality:
        (7,)
    :Dimensionality Remarks:
        7: { :term:`Paved`, :term:`Bldgs`, :term:`EveTr`, :term:`DecTr`, :term:`Grass`, :term:`BSoil`, :term:`Water`}
    :SUEWS-related variables:
        :option:`SatHydraulicCond <suews:SatHydraulicCond>`


.. option:: sddfull

    :Description:
        The sensesence degree days (SDD) needed to initiate leaf off. [°C]
    :Dimensionality:
        (3,)
    :Dimensionality Remarks:
        3: { :term:`EveTr`, :term:`DecTr`, :term:`Grass`}
    :SUEWS-related variables:
        :option:`SDDFull <suews:SDDFull>`


.. option:: sfr

    :Description:
        Surface cover fractions.
    :Dimensionality:
        (7,)
    :Dimensionality Remarks:
        7: { :term:`Paved`, :term:`Bldgs`, :term:`EveTr`, :term:`DecTr`, :term:`Grass`, :term:`BSoil`, :term:`Water`}
    :SUEWS-related variables:
        :option:`Fr_Bldgs <suews:Fr_Bldgs>`, :option:`Fr_Bsoil <suews:Fr_Bsoil>`, :option:`Fr_DecTr <suews:Fr_DecTr>`, :option:`Fr_EveTr <suews:Fr_EveTr>`, :option:`Fr_Grass <suews:Fr_Grass>`, :option:`Fr_Paved <suews:Fr_Paved>`, :option:`Fr_Water <suews:Fr_Water>`


.. option:: smdmethod

    :Description:
        Determines method for calculating soil moisture deficit (SMD).
    :Dimensionality:
        0
    :Dimensionality Remarks:
        Scalar
    :SUEWS-related variables:
        :option:`SMDMethod <suews:SMDMethod>`


.. option:: snowalb

    :Description:
        Initial snow albedo
    :Dimensionality:
        0
    :Dimensionality Remarks:
        Scalar
    :SUEWS-related variables:
        :option:`SnowAlb0 <suews:SnowAlb0>`


.. option:: snowalbmax

    :Description:
        Effective surface albedo (middle of the day value) for summertime.
    :Dimensionality:
        0
    :Dimensionality Remarks:
        Scalar
    :SUEWS-related variables:
        :option:`AlbedoMax <suews:AlbedoMax>`


.. option:: snowalbmin

    :Description:
        Effective surface albedo (middle of the day value) for wintertime (not including snow).
    :Dimensionality:
        0
    :Dimensionality Remarks:
        Scalar
    :SUEWS-related variables:
        :option:`AlbedoMin <suews:AlbedoMin>`


.. option:: snowdens

    :Description:
        Initial snow density of each land cover.
    :Dimensionality:
        (7,)
    :Dimensionality Remarks:
        7: { :term:`Paved`, :term:`Bldgs`, :term:`EveTr`, :term:`DecTr`, :term:`Grass`, :term:`BSoil`, :term:`Water`}
    :SUEWS-related variables:
        :option:`SnowDensBldgs <suews:SnowDensBldgs>`, :option:`SnowDensPaved <suews:SnowDensPaved>`, :option:`SnowDensDecTr <suews:SnowDensDecTr>`, :option:`SnowDensEveTr <suews:SnowDensEveTr>`, :option:`SnowDensGrass <suews:SnowDensGrass>`, :option:`SnowDensBSoil <suews:SnowDensBSoil>`, :option:`SnowDensWater <suews:SnowDensWater>`


.. option:: snowdensmax

    :Description:
        Maximum snow density [kg |m^-3|]
    :Dimensionality:
        0
    :Dimensionality Remarks:
        Scalar
    :SUEWS-related variables:
        :option:`SnowDensMax <suews:SnowDensMax>`


.. option:: snowdensmin

    :Description:
        Fresh snow density [kg |m^-3|]
    :Dimensionality:
        0
    :Dimensionality Remarks:
        Scalar
    :SUEWS-related variables:
        :option:`SnowDensMin <suews:SnowDensMin>`


.. option:: snowfrac

    :Description:
        Initial plan area fraction of snow on each land cover`
    :Dimensionality:
        (7,)
    :Dimensionality Remarks:
        7: { :term:`Paved`, :term:`Bldgs`, :term:`EveTr`, :term:`DecTr`, :term:`Grass`, :term:`BSoil`, :term:`Water`}
    :SUEWS-related variables:
        :option:`SnowFracBldgs <suews:SnowFracBldgs>`, :option:`SnowFracPaved <suews:SnowFracPaved>`, :option:`SnowFracDecTr <suews:SnowFracDecTr>`, :option:`SnowFracEveTr <suews:SnowFracEveTr>`, :option:`SnowFracGrass <suews:SnowFracGrass>`, :option:`SnowFracBSoil <suews:SnowFracBSoil>`, :option:`SnowFracWater <suews:SnowFracWater>`


.. option:: snowlimbldg

    :Description:
        Limit of the snow water equivalent for snow removal from roads and roofs [mm]
    :Dimensionality:
        0
    :Dimensionality Remarks:
        Scalar
    :SUEWS-related variables:
        :option:`SnowLimRemove <suews:SnowLimRemove>`


.. option:: snowlimpaved

    :Description:
        Limit of the snow water equivalent for snow removal from roads and roofs [mm]
    :Dimensionality:
        0
    :Dimensionality Remarks:
        Scalar
    :SUEWS-related variables:
        :option:`SnowLimRemove <suews:SnowLimRemove>`


.. option:: snowpack

    :Description:
        Initial snow water equivalent on each land cover
    :Dimensionality:
        (7,)
    :Dimensionality Remarks:
        7: { :term:`Paved`, :term:`Bldgs`, :term:`EveTr`, :term:`DecTr`, :term:`Grass`, :term:`BSoil`, :term:`Water`}
    :SUEWS-related variables:
        :option:`SnowPackBldgs <suews:SnowPackBldgs>`, :option:`SnowPackPaved <suews:SnowPackPaved>`, :option:`SnowPackDecTr <suews:SnowPackDecTr>`, :option:`SnowPackEveTr <suews:SnowPackEveTr>`, :option:`SnowPackGrass <suews:SnowPackGrass>`, :option:`SnowPackBSoil <suews:SnowPackBSoil>`, :option:`SnowPackWater <suews:SnowPackWater>`


.. option:: snowpacklimit

    :Description:
        Limit for the snow water equivalent when snow cover starts to be patchy [mm]
    :Dimensionality:
        (7,)
    :Dimensionality Remarks:
        7: { :term:`Paved`, :term:`Bldgs`, :term:`EveTr`, :term:`DecTr`, :term:`Grass`, :term:`BSoil`, :term:`Water`}
    :SUEWS-related variables:
        :option:`SnowLimPatch <suews:SnowLimPatch>`


.. option:: snowprof_24hr

    :Description:
        Hourly profile values used in snow clearing.
    :Dimensionality:
        (24, 2)
    :Dimensionality Remarks:
        24: hours of a day

        2: {Weekday, Weekend}
    :SUEWS-related variables:
        :option:`SnowClearingProfWD <suews:SnowClearingProfWD>`, :option:`SnowClearingProfWE <suews:SnowClearingProfWE>`


.. option:: snowuse

    :Description:
        Determines whether the snow part of the model runs.
    :Dimensionality:
        0
    :Dimensionality Remarks:
        Scalar
    :SUEWS-related variables:
        :option:`SnowUse <suews:SnowUse>`


.. option:: snowwater

    :Description:
        Initial amount of liquid water in the snow on each land cover
    :Dimensionality:
        (7,)
    :Dimensionality Remarks:
        7: { :term:`Paved`, :term:`Bldgs`, :term:`EveTr`, :term:`DecTr`, :term:`Grass`, :term:`BSoil`, :term:`Water`}
    :SUEWS-related variables:
        :option:`SnowWaterBldgsState <suews:SnowWaterBldgsState>`, :option:`SnowWaterPavedState <suews:SnowWaterPavedState>`, :option:`SnowWaterDecTrState <suews:SnowWaterDecTrState>`, :option:`SnowWaterEveTrState <suews:SnowWaterEveTrState>`, :option:`SnowWaterGrassState <suews:SnowWaterGrassState>`, :option:`SnowWaterBSoilState <suews:SnowWaterBSoilState>`, :option:`SnowWaterWaterState <suews:SnowWaterWaterState>`


.. option:: soildepth

    :Description:
        Depth of soil beneath the surface [mm]
    :Dimensionality:
        (7,)
    :Dimensionality Remarks:
        7: { :term:`Paved`, :term:`Bldgs`, :term:`EveTr`, :term:`DecTr`, :term:`Grass`, :term:`BSoil`, :term:`Water`}
    :SUEWS-related variables:
        :option:`SoilDepth <suews:SoilDepth>`


.. option:: soilstore_id

    :Description:
        Initial water stored in soil beneath each land cover
    :Dimensionality:
        (7,)
    :Dimensionality Remarks:
        7: { :term:`Paved`, :term:`Bldgs`, :term:`EveTr`, :term:`DecTr`, :term:`Grass`, :term:`BSoil`, :term:`Water`}
    :SUEWS-related variables:
        :option:`SoilstoreBldgsState <suews:SoilstoreBldgsState>`, :option:`SoilstorePavedState <suews:SoilstorePavedState>`, :option:`SoilstoreDecTrState <suews:SoilstoreDecTrState>`, :option:`SoilstoreEveTrState <suews:SoilstoreEveTrState>`, :option:`SoilstoreGrassState <suews:SoilstoreGrassState>`, :option:`SoilstoreBSoilState <suews:SoilstoreBSoilState>`


.. option:: soilstorecap

    :Description:
        Limit value for `SoilDepth` [mm]
    :Dimensionality:
        (7,)
    :Dimensionality Remarks:
        7: { :term:`Paved`, :term:`Bldgs`, :term:`EveTr`, :term:`DecTr`, :term:`Grass`, :term:`BSoil`, :term:`Water`}
    :SUEWS-related variables:
        :option:`SoilStoreCap <suews:SoilStoreCap>`


.. option:: stabilitymethod

    :Description:
        Defines which atmospheric stability functions are used.
    :Dimensionality:
        0
    :Dimensionality Remarks:
        Scalar
    :SUEWS-related variables:
        :option:`StabilityMethod <suews:StabilityMethod>`


.. option:: startdls

    :Description:
        Start of the day light savings [DOY]
    :Dimensionality:
        0
    :Dimensionality Remarks:
        Scalar
    :SUEWS-related variables:
        :option:`StartDLS <suews:StartDLS>`


.. option:: state_id

    :Description:
        Initial wetness condition on each land cover
    :Dimensionality:
        (7,)
    :Dimensionality Remarks:
        7: { :term:`Paved`, :term:`Bldgs`, :term:`EveTr`, :term:`DecTr`, :term:`Grass`, :term:`BSoil`, :term:`Water`}
    :SUEWS-related variables:
        :option:`BldgsState <suews:BldgsState>`, :option:`PavedState <suews:PavedState>`, :option:`DecTrState <suews:DecTrState>`, :option:`EveTrState <suews:EveTrState>`, :option:`GrassState <suews:GrassState>`, :option:`BSoilState <suews:BSoilState>`, :option:`WaterState <suews:WaterState>`


.. option:: statelimit

    :Description:
        Upper limit to the surface state. [mm]
    :Dimensionality:
        (7,)
    :Dimensionality Remarks:
        7: { :term:`Paved`, :term:`Bldgs`, :term:`EveTr`, :term:`DecTr`, :term:`Grass`, :term:`BSoil`, :term:`Water`}
    :SUEWS-related variables:
        :option:`StateLimit <suews:StateLimit>`


.. option:: storageheatmethod

    :Description:
        Determines method for calculating storage heat flux ΔQS.
    :Dimensionality:
        0
    :Dimensionality Remarks:
        Scalar
    :SUEWS-related variables:
        :option:`StorageHeatMethod <suews:StorageHeatMethod>`


.. option:: storedrainprm

    :Description:
        Coefficients used in drainage calculation.
    :Dimensionality:
        (6, 7)
    :Dimensionality Remarks:
        6: { `StorageMin`, `DrainageEq`, `DrainageCoef1`, `DrainageCoef2`, `StorageMax`, current storage}

        7: { :term:`Paved`, :term:`Bldgs`, :term:`EveTr`, :term:`DecTr`, :term:`Grass`, :term:`BSoil`, :term:`Water`}
    :SUEWS-related variables:
        :option:`DrainageCoef1 <suews:DrainageCoef1>`, :option:`DrainageCoef2 <suews:DrainageCoef2>`, :option:`DrainageEq <suews:DrainageEq>`, :option:`StorageMax <suews:StorageMax>`, :option:`StorageMin <suews:StorageMin>`


.. option:: surfacearea

    :Description:
        Area of the grid [ha].
    :Dimensionality:
        0
    :Dimensionality Remarks:
        Scalar
    :SUEWS-related variables:
        :option:`SurfaceArea <suews:SurfaceArea>`


.. option:: t_critic_cooling

    :Description:
        Critical cooling temperature.
    :Dimensionality:
        (2,)
    :Dimensionality Remarks:
        2: {Weekday, Weekend}
    :SUEWS-related variables:
        :option:`TCritic_Cooling_WD <suews:TCritic_Cooling_WD>`, :option:`TCritic_Cooling_WE <suews:TCritic_Cooling_WE>`


.. option:: t_critic_heating

    :Description:
        Critical heating temperature.
    :Dimensionality:
        (2,)
    :Dimensionality Remarks:
        2: {Weekday, Weekend}
    :SUEWS-related variables:
        :option:`TCritic_Heating_WD <suews:TCritic_Heating_WD>`, :option:`TCritic_Heating_WE <suews:TCritic_Heating_WE>`


.. option:: tau_a

    :Description:
        Time constant for snow albedo aging in cold snow [-]
    :Dimensionality:
        0
    :Dimensionality Remarks:
        Scalar
    :SUEWS-related variables:
        :option:`tau_a <suews:tau_a>`


.. option:: tau_f

    :Description:
        Time constant for snow albedo aging in melting snow [-]
    :Dimensionality:
        0
    :Dimensionality Remarks:
        Scalar
    :SUEWS-related variables:
        :option:`tau_f <suews:tau_f>`


.. option:: tau_r

    :Description:
        Time constant for snow density ageing [-]
    :Dimensionality:
        0
    :Dimensionality Remarks:
        Scalar
    :SUEWS-related variables:
        :option:`tau_r <suews:tau_r>`


.. option:: tempmeltfact

    :Description:
        Hourly temperature melt factor of snow [mm |K^-1| |h^-1|]
    :Dimensionality:
        0
    :Dimensionality Remarks:
        Scalar
    :SUEWS-related variables:
        :option:`TempMeltFactor <suews:TempMeltFactor>`


.. option:: th

    :Description:
        Upper air temperature limit [°C]
    :Dimensionality:
        0
    :Dimensionality Remarks:
        Scalar
    :SUEWS-related variables:
        :option:`TH <suews:TH>`


.. option:: theta_bioco2

    :Description:
        The convexity of the curve at light saturation.
    :Dimensionality:
        (3,)
    :Dimensionality Remarks:
        3: { :term:`EveTr`, :term:`DecTr`, :term:`Grass`}
    :SUEWS-related variables:
        :option:`theta <suews:theta>`


.. option:: timezone

    :Description:
        Time zone [h] for site relative to UTC (east is positive). This should be set according to the times given in the meteorological forcing file(s).
    :Dimensionality:
        0
    :Dimensionality Remarks:
        Scalar
    :SUEWS-related variables:
        :option:`Timezone <suews:Timezone>`


.. option:: tl

    :Description:
        Lower air temperature limit [°C]
    :Dimensionality:
        0
    :Dimensionality Remarks:
        Scalar
    :SUEWS-related variables:
        :option:`TL <suews:TL>`


.. option:: trafficrate

    :Description:
        Traffic rate used for CO2 flux calculation.
    :Dimensionality:
        (2,)
    :Dimensionality Remarks:
        2: {Weekday, Weekend}
    :SUEWS-related variables:
        :option:`TrafficRate_WD <suews:TrafficRate_WD>`, :option:`TrafficRate_WE <suews:TrafficRate_WE>`


.. option:: trafficunits

    :Description:
        Units for the traffic rate for the study area. Not used in v2018a.
    :Dimensionality:
        0
    :Dimensionality Remarks:
        Scalar
    :SUEWS-related variables:
        :option:`TrafficUnits <suews:TrafficUnits>`


.. option:: traffprof_24hr

    :Description:
        Hourly profile values used in traffic activity calculation.
    :Dimensionality:
        (24, 2)
    :Dimensionality Remarks:
        24: hours of a day

        2: {Weekday, Weekend}
    :SUEWS-related variables:
        :option:`TraffProfWD <suews:TraffProfWD>`, :option:`TraffProfWE <suews:TraffProfWE>`


.. option:: tstep

    :Description:
        Specifies the model time step [s].
    :Dimensionality:
        0
    :Dimensionality Remarks:
        Scalar
    :SUEWS-related variables:
        :option:`Tstep <suews:Tstep>`


.. option:: veg_type

    :Description:
        Internal use. Please DO NOT modify
    :Dimensionality:
        0
    :Dimensionality Remarks:
        Scalar
    :SUEWS-related variables:
        None


.. option:: waterdist

    :Description:
        Fraction of water redistribution
    :Dimensionality:
        (8, 6)
    :Dimensionality Remarks:
        8: { :term:`Paved`, :term:`Bldgs`, :term:`EveTr`, :term:`DecTr`, :term:`Grass`, :term:`BSoil`, :term:`Water`, one extra land cover type (currently NOT used)}

        6: { :term:`Paved`, :term:`Bldgs`, :term:`EveTr`, :term:`DecTr`, :term:`Grass`, :term:`BSoil`}
    :SUEWS-related variables:
        :option:`ToBSoil <suews:ToBSoil>`, :option:`ToBldgs <suews:ToBldgs>`, :option:`ToDecTr <suews:ToDecTr>`, :option:`ToEveTr <suews:ToEveTr>`, :option:`ToGrass <suews:ToGrass>`, :option:`ToPaved <suews:ToPaved>`, :option:`ToRunoff <suews:ToRunoff>`, :option:`ToSoilStore <suews:ToSoilStore>`, :option:`ToWater <suews:ToWater>`


.. option:: waterusemethod

    :Description:
        Defines how external water use is calculated.
    :Dimensionality:
        0
    :Dimensionality Remarks:
        Scalar
    :SUEWS-related variables:
        :option:`WaterUseMethod <suews:WaterUseMethod>`


.. option:: wetthresh

    :Description:
        Depth of water which determines whether evaporation occurs from a partially wet or completely wet surface [mm].
    :Dimensionality:
        (7,)
    :Dimensionality Remarks:
        7: { :term:`Paved`, :term:`Bldgs`, :term:`EveTr`, :term:`DecTr`, :term:`Grass`, :term:`BSoil`, :term:`Water`}
    :SUEWS-related variables:
        :option:`WetThreshold <suews:WetThreshold>`


.. option:: wuprofa_24hr

    :Description:
        Hourly profile values used in automatic irrigation.
    :Dimensionality:
        (24, 2)
    :Dimensionality Remarks:
        24: hours of a day

        2: {Weekday, Weekend}
    :SUEWS-related variables:
        :option:`WaterUseProfAutoWD <suews:WaterUseProfAutoWD>`, :option:`WaterUseProfAutoWE <suews:WaterUseProfAutoWE>`


.. option:: wuprofm_24hr

    :Description:
        Hourly profile values used in manual irrigation.
    :Dimensionality:
        (24, 2)
    :Dimensionality Remarks:
        24: hours of a day

        2: {Weekday, Weekend}
    :SUEWS-related variables:
        :option:`WaterUseProfManuWD <suews:WaterUseProfManuWD>`, :option:`WaterUseProfManuWE <suews:WaterUseProfManuWE>`


.. option:: z

    :Description:
        Measurement height [m].
    :Dimensionality:
        0
    :Dimensionality Remarks:
        Scalar
    :SUEWS-related variables:
        :option:`z <suews:z>`


.. option:: z0m_in

    :Description:
        Roughness length for momentum [m]
    :Dimensionality:
        0
    :Dimensionality Remarks:
        Scalar
    :SUEWS-related variables:
        :option:`z0 <suews:z0>`


.. option:: zdm_in

    :Description:
        Zero-plane displacement [m]
    :Dimensionality:
        0
    :Dimensionality Remarks:
        Scalar
    :SUEWS-related variables:
        :option:`zd <suews:zd>`

