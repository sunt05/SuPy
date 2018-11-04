
.. option:: aerodynamicresistancemethod

	:Dimensionality:
		0
	:Description:
		internal use, DO NOT change.
	:SUEWS-related variables:
		nan
    


.. option:: ah_min

	:Dimensionality:
		(2,)
	:Description:
		Minimum QF on weekdays [W |m^-2|]
Minimum QF on weekends [W |m^-2|]
	:SUEWS-related variables:
		AHMin_WD, AHMin_WE
    


.. option:: ah_slope_cooling

	:Dimensionality:
		(2,)
	:Description:
		Cooling slope of QF on weekdays [W |m^-2| |K^-1|]
Cooling slope of QF on weekends [W |m^-2| |K^-1|]
	:SUEWS-related variables:
		AHSlope_Cooling_WD, AHSlope_Cooling_WE
    


.. option:: ah_slope_heating

	:Dimensionality:
		(2,)
	:Description:
		Heating slope of QF on weekdays [W |m^-2| |K^-1|]
Heating slope of QF on weekends [W |m^-2| |K^-1|]
	:SUEWS-related variables:
		AHSlope_Heating_WD, AHSlope_Heating_WE
    


.. option:: ahprof_24hr

	:Dimensionality:
		(24, 2)
	:Description:
		Code linking to `EnergyUseProfWD` in `SUEWS_Profiles.txt`.
Code linking to `EnergyUseProfWE` in `SUEWS_Profiles.txt`.
	:SUEWS-related variables:
		EnergyUseProfWD, EnergyUseProfWE
    


.. option:: alb

	:Dimensionality:
		(7,)
	:Description:
		Effective surface albedo (middle of the day value) for summertime. View factors should be taken into account.
Example values [-] 0.85 Järvi et al. (2014) [Leena2014]_ 
Effective surface albedo (middle of the day value) for summertime, full leaf-on. View factors should be taken into account. Example values [-] 0.1 EveTr Oke (1987) [Ok87]_  0.18 DecTr Oke (1987) [Ok87]_  0.21 Grass Oke (1987) [Ok87]_ 
Effective albedo of the water surface. View factors should be taken into account. Example values [-] 0.1 Water Oke (1987) [Ok87]_ 
	:SUEWS-related variables:
		AlbedoMax
    


.. option:: albdectr_id

	:Dimensionality:
		0
	:Description:
		internal use, DO NOT change.
	:SUEWS-related variables:
		nan
    


.. option:: albevetr_id

	:Dimensionality:
		0
	:Description:
		internal use, DO NOT change.
	:SUEWS-related variables:
		nan
    


.. option:: albgrass_id

	:Dimensionality:
		0
	:Description:
		internal use, DO NOT change.
	:SUEWS-related variables:
		nan
    


.. option:: albmax_dectr

	:Dimensionality:
		0
	:Description:
		Effective surface albedo (middle of the day value) for summertime. View factors should be taken into account.
Example values [-] 0.85 Järvi et al. (2014) [Leena2014]_ 
Effective surface albedo (middle of the day value) for summertime, full leaf-on. View factors should be taken into account. Example values [-] 0.1 EveTr Oke (1987) [Ok87]_  0.18 DecTr Oke (1987) [Ok87]_  0.21 Grass Oke (1987) [Ok87]_ 
Effective albedo of the water surface. View factors should be taken into account. Example values [-] 0.1 Water Oke (1987) [Ok87]_ 
	:SUEWS-related variables:
		AlbedoMax
    


.. option:: albmax_evetr

	:Dimensionality:
		0
	:Description:
		Effective surface albedo (middle of the day value) for summertime. View factors should be taken into account.
Example values [-] 0.85 Järvi et al. (2014) [Leena2014]_ 
Effective surface albedo (middle of the day value) for summertime, full leaf-on. View factors should be taken into account. Example values [-] 0.1 EveTr Oke (1987) [Ok87]_  0.18 DecTr Oke (1987) [Ok87]_  0.21 Grass Oke (1987) [Ok87]_ 
Effective albedo of the water surface. View factors should be taken into account. Example values [-] 0.1 Water Oke (1987) [Ok87]_ 
	:SUEWS-related variables:
		AlbedoMax
    


.. option:: albmax_grass

	:Dimensionality:
		0
	:Description:
		Effective surface albedo (middle of the day value) for summertime. View factors should be taken into account.
Example values [-] 0.85 Järvi et al. (2014) [Leena2014]_ 
Effective surface albedo (middle of the day value) for summertime, full leaf-on. View factors should be taken into account. Example values [-] 0.1 EveTr Oke (1987) [Ok87]_  0.18 DecTr Oke (1987) [Ok87]_  0.21 Grass Oke (1987) [Ok87]_ 
Effective albedo of the water surface. View factors should be taken into account. Example values [-] 0.1 Water Oke (1987) [Ok87]_ 
	:SUEWS-related variables:
		AlbedoMax
    


.. option:: albmin_dectr

	:Dimensionality:
		0
	:Description:
		Effective surface albedo (middle of the day value) for wintertime (not including snow). View factors should be taken into account.

.. note::
    Not currently used for non-vegetated surfaces – set the same as AlbedoMax.
Example values [-] 0.18 Järvi et al. (2014) [Leena2014]_ 
Effective surface albedo (middle of the day value) for wintertime (not including snow), leaf-off. View factors should be taken into account. Example values [-] 0.1 EveTr Oke (1987) [Ok87]_  0.18 DecTr Oke (1987) [Ok87]_  0.21 Grass Oke (1987) [Ok87]_ 
View factors should be taken into account. Not currently used for water surface - set same as AlbedoMax.
	:SUEWS-related variables:
		AlbedoMin
    


.. option:: albmin_evetr

	:Dimensionality:
		0
	:Description:
		Effective surface albedo (middle of the day value) for wintertime (not including snow). View factors should be taken into account.

.. note::
    Not currently used for non-vegetated surfaces – set the same as AlbedoMax.
Example values [-] 0.18 Järvi et al. (2014) [Leena2014]_ 
Effective surface albedo (middle of the day value) for wintertime (not including snow), leaf-off. View factors should be taken into account. Example values [-] 0.1 EveTr Oke (1987) [Ok87]_  0.18 DecTr Oke (1987) [Ok87]_  0.21 Grass Oke (1987) [Ok87]_ 
View factors should be taken into account. Not currently used for water surface - set same as AlbedoMax.
	:SUEWS-related variables:
		AlbedoMin
    


.. option:: albmin_grass

	:Dimensionality:
		0
	:Description:
		Effective surface albedo (middle of the day value) for wintertime (not including snow). View factors should be taken into account.

.. note::
    Not currently used for non-vegetated surfaces – set the same as AlbedoMax.
Example values [-] 0.18 Järvi et al. (2014) [Leena2014]_ 
Effective surface albedo (middle of the day value) for wintertime (not including snow), leaf-off. View factors should be taken into account. Example values [-] 0.1 EveTr Oke (1987) [Ok87]_  0.18 DecTr Oke (1987) [Ok87]_  0.21 Grass Oke (1987) [Ok87]_ 
View factors should be taken into account. Not currently used for water surface - set same as AlbedoMax.
	:SUEWS-related variables:
		AlbedoMin
    


.. option:: alpha_bioco2

	:Dimensionality:
		(3,)
	:Description:
		The mean apparent ecosystem quantum. Represents the initial slope of the light-response curve.
	:SUEWS-related variables:
		alpha
    


.. option:: alpha_enh_bioco2

	:Dimensionality:
		(3,)
	:Description:
		Part of the `alpha` coeﬃcient related to the fraction of vegetation.
	:SUEWS-related variables:
		alpha_enh
    


.. option:: alt

	:Dimensionality:
		0
	:Description:
		Used for both the radiation and water flow between grids.

|NotAvail|
	:SUEWS-related variables:
		Alt
    


.. option:: baset

	:Dimensionality:
		(3,)
	:Description:
		See section 2.2 Järvi et al. (2011) [J11]_; Appendix A Järvi et al. (2014) [Leena2014]_. Example values [°C] 5 EveTr Järvi et al. (2011) [J11]_  5 DecTr Järvi et al. (2011) [J11]_  5 Grass Järvi et al. (2011) [J11]_ 
	:SUEWS-related variables:
		BaseT
    


.. option:: basete

	:Dimensionality:
		(3,)
	:Description:
		See section 2.2 Järvi et al. (2011) [J11]_ ; Appendix A Järvi et al. (2014) [Leena2014]_. Example values [°C] 10 EveTr Järvi et al. (2011) [J11]_  10 DecTr Järvi et al. (2011) [J11]_  10 Grass Järvi et al. (2011) [J11]_ 
	:SUEWS-related variables:
		BaseTe
    


.. option:: basethdd

	:Dimensionality:
		0
	:Description:
		Base temperature for heating degree days [°C]
	:SUEWS-related variables:
		BaseTHDD
    


.. option:: beta_bioco2

	:Dimensionality:
		(3,)
	:Description:
		The light-saturated gross photosynthesis of the canopy.
	:SUEWS-related variables:
		beta
    


.. option:: beta_enh_bioco2

	:Dimensionality:
		(3,)
	:Description:
		Part of the `beta` coeﬃcient related to the fraction of vegetation.
	:SUEWS-related variables:
		beta_enh
    


.. option:: bldgh

	:Dimensionality:
		0
	:Description:
		Mean building height [m]
	:SUEWS-related variables:
		H_Bldgs
    


.. option:: capmax_dec

	:Dimensionality:
		0
	:Description:
		Maximum water storage capacity for upper surfaces (i.e. canopy) Min and max values are to account for seasonal variation (e.g. leaf-on/leaf-off differences for vegetated surfaces).

.. note::
    Not currently used for non-vegetated surfaces – set the same as `StorageMin`.

Maximum water storage capacity for upper surfaces (i.e. canopy) Min/max values are to account for seasonal variation (e.g. leaf-off/leaf-on differences for vegetated surfaces) Only used for DecTr surfaces - set EveTr and Grass values the same as StorageMin. Example values [mm] 1.3 EveTr Breuer et al. (2003) [Br03]_  0.8 DecTr Breuer et al. (2003) [Br03]_  1.9 Grass Breuer et al. (2003) [Br03]_ 
Maximum water storage capacity for upper surfaces (i.e. canopy) Min and max values are to account for seasonal variation - not used for water surfaces so set same as StorageMin.
	:SUEWS-related variables:
		StorageMax
    


.. option:: capmin_dec

	:Dimensionality:
		0
	:Description:
		Minimum water storage capacity for upper surfaces (i.e. canopy). Min/max values are to account for seasonal variation (e.g. leaf-on/leaf-off differences for vegetated surfaces).

.. note::
    Not currently used for non-vegetated surfaces – set the same as `StorageMax`.
Minimum water storage capacity for upper surfaces (i.e. canopy). Min/max values are to account for seasonal variation (e.g. leaf-off/leaf-on differences for vegetated surfaces). Example values [mm] 1.3 EveTr Breuer et al. (2003) [Br03]_  0.3 DecTr Breuer et al. (2003) [Br03]_  1.9 Grass Breuer et al. (2003) [Br03]_ 
Minimum water storage capacity for upper surfaces (i.e. canopy). Min/max values are to account for seasonal variation - not used for water surfaces. Example values [mm] 0.5 Water
	:SUEWS-related variables:
		StorageMin
    


.. option:: chanohm

	:Dimensionality:
		(7,)
	:Description:
		Bulk transfer coefficient for this surface to use in AnOHM [-]
Bulk transfer coefficient for this surface to use in AnOHM [-]
Bulk transfer coefficient for this surface to use in AnOHM [-]
Bulk transfer coefficient for this surface to use in AnOHM [-]
	:SUEWS-related variables:
		AnOHM_Ch
    


.. option:: cpanohm

	:Dimensionality:
		(7,)
	:Description:
		Volumetric heat capacity for this surface to use in AnOHM [J |m^-3|]
Volumetric heat capacity for this surface to use in AnOHM [J |m^-3|]
Volumetric heat capacity for this surface to use in AnOHM [J |m^-3|]
Volumetric heat capacity for this surface to use in AnOHM [J |m^-3|]
	:SUEWS-related variables:
		AnOHM_Cp
    


.. option:: crwmax

	:Dimensionality:
		0
	:Description:
		Maximum water holding capacity of snow [mm]
	:SUEWS-related variables:
		CRWMax
    


.. option:: crwmin

	:Dimensionality:
		0
	:Description:
		Minimum water holding capacity of snow [mm]
	:SUEWS-related variables:
		CRWMin
    


.. option:: daywat

	:Dimensionality:
		(7,)
	:Description:
		Irrigation allowed on Sundays [1], if not [0]
Irrigation allowed on Mondays [1], if not [0]
Irrigation allowed on Tuesdays [1], if not [0]
Irrigation allowed on Wednesdays [1], if not [0]
Irrigation allowed on Thursdays [1], if not [0]
Irrigation allowed on Fridays [1], if not [0]
Irrigation allowed on Saturdays [1], if not [0]
	:SUEWS-related variables:
		DayWat(1), DayWat(2), DayWat(3), DayWat(4), DayWat(5), DayWat(6), DayWat(7)
    


.. option:: daywatper

	:Dimensionality:
		(7,)
	:Description:
		Fraction of properties using irrigation on Sundays [0-1]
Fraction of properties using irrigation on Mondays [0-1]
Fraction of properties using irrigation on Tuesdays [0-1]
Fraction of properties using irrigation on Wednesdays [0-1]
Fraction of properties using irrigation on Thursdays [0-1]
Fraction of properties using irrigation on Fridays [0-1]
Fraction of properties using irrigation on Saturdays [0-1]
	:SUEWS-related variables:
		DayWatPer(1), DayWatPer(2), DayWatPer(3), DayWatPer(4), DayWatPer(5), DayWatPer(6), DayWatPer(7)
    


.. option:: decidcap_id

	:Dimensionality:
		0
	:Description:
		internal use, DO NOT change.
	:SUEWS-related variables:
		nan
    


.. option:: dectreeh

	:Dimensionality:
		0
	:Description:
		Mean height of deciduous trees [m]
	:SUEWS-related variables:
		H_DecTr
    


.. option:: diagnose

	:Dimensionality:
		0
	:Description:
		internal use, DO NOT change.
	:SUEWS-related variables:
		nan
    


.. option:: diagqn

	:Dimensionality:
		0
	:Description:
		internal use, DO NOT change.
	:SUEWS-related variables:
		nan
    


.. option:: diagqs

	:Dimensionality:
		0
	:Description:
		internal use, DO NOT change.
	:SUEWS-related variables:
		nan
    


.. option:: dqndt

	:Dimensionality:
		0
	:Description:
		internal use, DO NOT change.
	:SUEWS-related variables:
		nan
    


.. option:: dqnsdt

	:Dimensionality:
		0
	:Description:
		internal use, DO NOT change.
	:SUEWS-related variables:
		nan
    


.. option:: drainrt

	:Dimensionality:
		0
	:Description:
		Drainage rate of bucket for LUMPS [mm |h^-1|] Used for LUMPS surface wetness control. Default recommended value of 0.25 mm |h^-1| from Loridan et al. (2011) [L2011]_ .
	:SUEWS-related variables:
		LUMPS_DrRate
    


.. option:: dt_since_start

	:Dimensionality:
		0
	:Description:
		internal use, DO NOT change.
	:SUEWS-related variables:
		nan
    


.. option:: ef_umolco2perj

	:Dimensionality:
		0
	:Description:
		Emission factor for fuels.
	:SUEWS-related variables:
		EF_umolCO2perJ
    


.. option:: emis

	:Dimensionality:
		(7,)
	:Description:
		Effective surface emissivity. View factors should be taken into account.
Effective surface emissivity. View factors should be taken into account Example values [-] 0.99 Järvi et al. (2014) [Leena2014]_ 
Effective surface emissivity. View factors should be taken into account. Example values [-] 0.98 EveTr Oke (1987) [Ok87]_  0.98 DecTr Oke (1987) [Ok87]_  0.93 Grass Oke (1987) [Ok87]_ 
Effective surface emissivity. View factors should be taken into account Example values [-] 0.95 Water Oke (1987) [Ok87]_ 
	:SUEWS-related variables:
		Emissivity
    


.. option:: emissionsmethod

	:Dimensionality:
		0
	:Description:
		Determines method for QF calculation.
	:SUEWS-related variables:
		EmissionsMethod
    


.. option:: enddls

	:Dimensionality:
		0
	:Description:
		End of the day light savings [DOY] See `Day_Light_Savings`.
	:SUEWS-related variables:
		EndDLS
    


.. option:: enef_v_jkm

	:Dimensionality:
		0
	:Description:
		Energy emission factor [J |km^-1| ]
	:SUEWS-related variables:
		EnEF_v_Jkm
    


.. option:: evapmethod

	:Dimensionality:
		0
	:Description:
		Determines method for evaporation calculation. Currently fixed. DO NOT change.
	:SUEWS-related variables:
		nan
    


.. option:: evetreeh

	:Dimensionality:
		0
	:Description:
		Mean height of evergreen trees [m]
	:SUEWS-related variables:
		H_EveTr
    


.. option:: faibldg

	:Dimensionality:
		0
	:Description:
		Frontal area index for buildings [-] Required if `RoughLenMomMethod` = 3 in `RunControl.nml` .
	:SUEWS-related variables:
		FAI_Bldgs
    


.. option:: faidectree

	:Dimensionality:
		0
	:Description:
		Frontal area index for deciduous trees [-] Required if `RoughLenMomMethod` = 3 in `RunControl.nml` .
	:SUEWS-related variables:
		FAI_DecTr
    


.. option:: faievetree

	:Dimensionality:
		0
	:Description:
		Frontal area index for evergreen trees [-] Required if `RoughLenMomMethod` = 3 in `RunControl.nml` .
	:SUEWS-related variables:
		FAI_EveTr
    


.. option:: faut

	:Dimensionality:
		0
	:Description:
		Fraction of irrigated area that is irrigated using automated systems (e.g. sprinklers).
	:SUEWS-related variables:
		Faut
    


.. option:: fcef_v_kgkm

	:Dimensionality:
		0
	:Description:
		CO2 emission factor [kg |km^-1| ]
	:SUEWS-related variables:
		FcEF_v_kgkm
    


.. option:: flowchange

	:Dimensionality:
		0
	:Description:
		Difference in input and output flows for water surface [mm |h^-1|] Used to indicate river or stream flow through the grid. Currently not fully tested!
	:SUEWS-related variables:
		FlowChange
    


.. option:: frfossilfuel_heat

	:Dimensionality:
		0
	:Description:
		Proportion of building energy use from fossil fuels rather than electricity.
	:SUEWS-related variables:
		FrFossilFuel_Heat
    


.. option:: frfossilfuel_nonheat

	:Dimensionality:
		0
	:Description:
		Fraction of Fossil Fuel for non heat.
	:SUEWS-related variables:
		FrFossilFuel_NonHeat
    


.. option:: g1

	:Dimensionality:
		0
	:Description:
		Related to maximum surface conductance [mm |s^-1|]
	:SUEWS-related variables:
		G1
    


.. option:: g2

	:Dimensionality:
		0
	:Description:
		Related to Kdown dependence [W |m^-2|]
	:SUEWS-related variables:
		G2
    


.. option:: g3

	:Dimensionality:
		0
	:Description:
		Related to VPD dependence [units depend on `gsModel` in `RunControl.nml`]
	:SUEWS-related variables:
		G3
    


.. option:: g4

	:Dimensionality:
		0
	:Description:
		Related to VPD dependence [units depend on `gsModel` in `RunControl.nml`]
	:SUEWS-related variables:
		G4
    


.. option:: g5

	:Dimensionality:
		0
	:Description:
		Related to temperature dependence [°C]
	:SUEWS-related variables:
		G5
    


.. option:: g6

	:Dimensionality:
		0
	:Description:
		Related to soil moisture dependence [m|m^-1|]
	:SUEWS-related variables:
		G6
    


.. option:: gdd_id

	:Dimensionality:
		(5,)
	:Description:
		internal use, DO NOT change.
	:SUEWS-related variables:
		nan
    


.. option:: gddfull

	:Dimensionality:
		(3,)
	:Description:
		This should be checked carefully for your study area using modelled LAI from the DailyState output file compared to known behaviour in the study area. See section 2.2 Järvi et al. (2011) [J11]_ ; Appendix A Järvi et al. (2014) [Leena2014]_ for more details. Example values [°C] 300 EveTr Järvi et al. (2011) [J11]_  300 DecTr Järvi et al. (2011) [J11]_  300 Grass Järvi et al. (2011) [J11]_ 
	:SUEWS-related variables:
		GDDFull
    


.. option:: gridiv

	:Dimensionality:
		0
	:Description:
		internal use, DO NOT change.
	:SUEWS-related variables:
		nan
    


.. option:: gsmodel

	:Dimensionality:
		0
	:Description:
		Formulation choice for conductance calculation.

- :code:`1` Järvi et al. (2011) [J11]_
- :code:`2` (**Recommended**) Ward et al. (2016)  [W16]_ 
	:SUEWS-related variables:
		gsModel
    


.. option:: hdd_id

	:Dimensionality:
		(12,)
	:Description:
		internal use, DO NOT change.
	:SUEWS-related variables:
		nan
    


.. option:: humactivity_24hr

	:Dimensionality:
		(24, 2)
	:Description:
		Code linking to `ActivityProfWD` in `SUEWS_Profiles.txt`.
Code linking to `ActivityProfWE` in `SUEWS_Profiles.txt`.
	:SUEWS-related variables:
		ActivityProfWD, ActivityProfWE
    


.. option:: icefrac

	:Dimensionality:
		(7,)
	:Description:
		internal use, DO NOT change.
	:SUEWS-related variables:
		nan
    


.. option:: ie_a

	:Dimensionality:
		(3,)
	:Description:
		Coefficient for automatic irrigation model [mm d -1 ]
Coefficient for automatic irrigation model [mm d -1 |K^-1|]
Coefficient for automatic irrigation model [mm d -2 ]
	:SUEWS-related variables:
		Ie_a1, Ie_a2, Ie_a3
    


.. option:: ie_end

	:Dimensionality:
		0
	:Description:
		Day when irrigation ends [DOY]
	:SUEWS-related variables:
		Ie_end
    


.. option:: ie_m

	:Dimensionality:
		(3,)
	:Description:
		Coefficient for manual irrigation model [mm d -1 ]
Coefficient for manual irrigation model [mm d -1 |K^-1|]
Coefficient for manual irrigation model [mm d -2 ]
	:SUEWS-related variables:
		Ie_m1, Ie_m2, Ie_m3
    


.. option:: ie_start

	:Dimensionality:
		0
	:Description:
		Day when irrigation starts [DOY]
	:SUEWS-related variables:
		Ie_start
    


.. option:: internalwateruse_h

	:Dimensionality:
		0
	:Description:
		Internal water use [mm |h^-1|]
	:SUEWS-related variables:
		InternalWaterUse
    


.. option:: irrfracconif

	:Dimensionality:
		0
	:Description:
		Fraction of evergreen trees that are irrigated [-] e.g. 50% of the evergreen trees/shrubs are irrigated
	:SUEWS-related variables:
		IrrFr_EveTr
    


.. option:: irrfracdecid

	:Dimensionality:
		0
	:Description:
		Fraction of deciduous trees that are irrigated [-]
	:SUEWS-related variables:
		IrrFr_DecTr
    


.. option:: irrfracgrass

	:Dimensionality:
		0
	:Description:
		Fraction of grass that is irrigated [-]
	:SUEWS-related variables:
		IrrFr_Grass
    


.. option:: kkanohm

	:Dimensionality:
		(7,)
	:Description:
		Thermal conductivity for this surface to use in AnOHM [W m |K^-1|]
Thermal conductivity for this surface to use in AnOHM [W m |K^-1|]
Thermal conductivity for this surface to use in AnOHM [W m |K^-1|]
Thermal conductivity for this surface to use in AnOHM [W m |K^-1|]
	:SUEWS-related variables:
		AnOHM_Kk
    


.. option:: kmax

	:Dimensionality:
		0
	:Description:
		Maximum incoming shortwave radiation [W |m^-2|]
	:SUEWS-related variables:
		Kmax
    


.. option:: lai_id

	:Dimensionality:
		(3,)
	:Description:
		internal use, DO NOT change.
	:SUEWS-related variables:
		nan
    


.. option:: laicalcyes

	:Dimensionality:
		0
	:Description:
		internal use, DO NOT change.
	:SUEWS-related variables:
		nan
    


.. option:: laimax

	:Dimensionality:
		(3,)
	:Description:
		full leaf-on summertime value Example values [|m^-2| |m^-2|] 5.1 EveTr Breuer et al. (2003) [Br03]_  5.5 DecTr Breuer et al. (2003) [Br03]_  5.9 Grass Breuer et al. (2003) [Br03]_ 
	:SUEWS-related variables:
		LAIMax
    


.. option:: laimin

	:Dimensionality:
		(3,)
	:Description:
		leaf-off wintertime value Example values [|m^-2| |m^-2|] 4. EveTr Järvi et al. (2011) [J11]_  1. DecTr Järvi et al. (2011) [J11]_  1.6 Grass Grimmond and Oke (1991) [3] and references therein
	:SUEWS-related variables:
		LAIMin
    


.. option:: laipower

	:Dimensionality:
		(4, 3)
	:Description:
		Example values LAIEq 0.03 Järvi et al. (2011) [J11]_ 0 0.04 Järvi et al. (2014) [Leena2014]_ 1
Example values [|K^-1|] LAIEq 0.0005 Järvi et al. (2011) [J11]_ 0 0.001 Järvi et al. (2014) [Leena2014]_ 1
Example values LAIEq 0.03 Järvi et al. (2011) [J11]_ 0 -1.5 Järvi et al. (2014) [Leena2014]_ 1
Example values [|K^-1|] LAIEq 0.0005 Järvi et al. (2011) [J11]_ 0 0.0015 Järvi et al. (2014) [Leena2014]_ 1
	:SUEWS-related variables:
		LeafGrowthPower1, LeafGrowthPower2, LeafOffPower1, LeafOffPower2
    


.. option:: laitype

	:Dimensionality:
		(3,)
	:Description:
		Options 0 Järvi et al. (2011) [J11]_  1 Järvi et al. (2014) [Leena2014]_  Coefficients are specified in the following four columns. N.B. North and South hemispheres are treated slightly differently.
	:SUEWS-related variables:
		LAIEq
    


.. option:: lat

	:Dimensionality:
		0
	:Description:
		Use coordinate system WGS84. Positive values are northern hemisphere (negative southern hemisphere). Used in radiation calculations.

.. note::
    If the total modelled area is small the latitude and longitude could be the same for each grid but small differences in radiation will not be determined.
    If you are defining the latitude and longitude differently between grids make certain that you provide enough decimal places.
	:SUEWS-related variables:
		lat
    


.. option:: lng

	:Dimensionality:
		0
	:Description:
		Use coordinate system WGS84. For compatibility with GIS, negative values are to the west, positive values are to the east (e.g. Vancouver = -123.12; Shanghai = 121.47)

.. note::
    this is a change of sign convention between v2016a and v2017a See `lat` for more details.
	:SUEWS-related variables:
		lng
    


.. option:: maxconductance

	:Dimensionality:
		(3,)
	:Description:
		Example values [mm |s^-1|] 7.4 EveTr Järvi et al. (2011) [J11]_  11.7 DecTr Järvi et al. (2011) [J11]_  33.1 Grass (unirrigated) Järvi et al. (2011) [J11]_  40. Grass (irrigated) Järvi et al. (2011) [J11]_ 
	:SUEWS-related variables:
		MaxConductance
    


.. option:: maxqfmetab

	:Dimensionality:
		0
	:Description:
		Maximum QF Metab value.
	:SUEWS-related variables:
		MaxQFMetab
    


.. option:: meltwaterstore

	:Dimensionality:
		(7,)
	:Description:
		Initial amount of liquid water in the snow
	:SUEWS-related variables:
		SnowWaterPavedState, SnowWaterBldgsState, SnowWaterEveTrState, SnowWaterDecTrState, SnowWaterGrassState, SnowWaterBSoilState, SnowWaterWaterState
    


.. option:: min_res_bioco2

	:Dimensionality:
		(3,)
	:Description:
		Minimum soil respiration rate (for cold-temperature limit)
	:SUEWS-related variables:
		min_respi
    


.. option:: minqfmetab

	:Dimensionality:
		0
	:Description:
		Minimum QF Metab value.
	:SUEWS-related variables:
		MinQFMetab
    


.. option:: narp_emis_snow

	:Dimensionality:
		0
	:Description:
		Effective surface emissivity. View factors should be taken into account.
Effective surface emissivity. View factors should be taken into account Example values [-] 0.99 Järvi et al. (2014) [Leena2014]_ 
Effective surface emissivity. View factors should be taken into account. Example values [-] 0.98 EveTr Oke (1987) [Ok87]_  0.98 DecTr Oke (1987) [Ok87]_  0.93 Grass Oke (1987) [Ok87]_ 
Effective surface emissivity. View factors should be taken into account Example values [-] 0.95 Water Oke (1987) [Ok87]_ 
	:SUEWS-related variables:
		Emissivity
    


.. option:: narp_trans_site

	:Dimensionality:
		0
	:Description:
		Atmospheric transmissivity for NARP [-] Value must in the range 0-1. Default recommended value of 1.
	:SUEWS-related variables:
		NARP_Trans
    


.. option:: netradiationmethod

	:Dimensionality:
		0
	:Description:
		Determines method for calculation of radiation fluxes.
	:SUEWS-related variables:
		NetRadiationMethod
    


.. option:: numcapita

	:Dimensionality:
		0
	:Description:
		population
	:SUEWS-related variables:
		PopDensDay, PopDensNight
    


.. option:: ohm_coef

	:Dimensionality:
		(8, 4, 3)
	:Description:
		Coefficient for Q* term [-]
Coefficient for dQ*/dt term [h]
Constant term [W |m^-2|]
	:SUEWS-related variables:
		a1, a2, a3
    


.. option:: ohm_threshsw

	:Dimensionality:
		(8,)
	:Description:
		Temperature threshold determining whether summer/winter OHM coefficients are applied [°C] If 5-day running mean air temperature is greater than or equal to this threshold, OHM coefficients for summertime are applied; otherwise coefficients for wintertime are applied. 
Temperature threshold determining whether summer/winter OHM coefficients are applied [°C] If 5-day running mean air temperature is greater than or equal to this threshold, OHM coefficients for summertime are applied; otherwise coefficients for wintertime are applied. Not actually used for Snow surface as winter wet conditions always assumed. 
Temperature threshold determining whether summer/winter OHM coefficients are applied [°C] If 5-day running mean air temperature is greater than or equal to this threshold, OHM coefficients for summertime are applied; otherwise coefficients for wintertime are applied. 
Temperature threshold determining whether summer/winter OHM coefficients are applied [°C] If 5-day running mean air temperature is greater than or equal to this threshold, OHM coefficients for summertime are applied; otherwise coefficients for wintertime are applied. 
	:SUEWS-related variables:
		OHMThresh_SW
    


.. option:: ohm_threshwd

	:Dimensionality:
		(8,)
	:Description:
		Soil moisture threshold determining whether wet/dry OHM coefficients are applied [-] If soil moisture (as a proportion of maximum soil moisture capacity) exceeds this threshold for bare soil and vegetated surfaces, OHM coefficients for wet conditions are applied; otherwise coefficients for dry coefficients are applied. Note that OHM coefficients for wet conditions are applied if the surface is wet. Not actually used for building and paved surfaces (as impervious). 
Soil moisture threshold determining whether wet/dry OHM coefficients are applied [-] If soil moisture (as a proportion of maximum soil moisture capacity) exceeds this threshold for bare soil and vegetated surfaces, OHM coefficients for wet conditions are applied; otherwise coefficients for dry coefficients are applied. Note that OHM coefficients for wet conditions are applied if the surface is wet. Not actually used for Snow surface as winter wet conditions always assumed. 
Soil moisture threshold determining whether wet/dry OHM coefficients are applied [-] If soil moisture (as a proportion of maximum soil moisture capacity) exceeds this threshold for bare soil and vegetated surfaces, OHM coefficients for wet conditions are applied; otherwise coefficients for dry coefficients are applied. Note that OHM coefficients for wet conditions are applied if the surface is wet. 
Soil moisture threshold determining whether wet/dry OHM coefficients are applied [-] If soil moisture (as a proportion of maximum soil moisture capacity) exceeds this threshold for bare soil and vegetated surfaces, OHM coefficients for wet conditions are applied; otherwise coefficients for dry coefficients are applied. Note that OHM coefficients for wet conditions are applied if the surface is wet. Not actually used for water surface (as no soil surface beneath). 
	:SUEWS-related variables:
		OHMThresh_WD
    


.. option:: ohmincqf

	:Dimensionality:
		0
	:Description:
		Determines whether the storage heat flux calculation uses |Qstar| or ( |Qstar| +QF).
	:SUEWS-related variables:
		OHMIncQF
    


.. option:: pipecapacity

	:Dimensionality:
		0
	:Description:
		Storage capacity of pipes [mm] Runoff amounting to less than the value specified here is assumed to be removed by pipes.
	:SUEWS-related variables:
		PipeCapacity
    


.. option:: popdensdaytime

	:Dimensionality:
		0
	:Description:
		Daytime population density (i.e. workers, tourists) [people ha -1 ] Population density is required if EmissionsMethod = 2 in `RunControl.nml` . The model will use the average of daytime and night-time population densities, unless only one is provided. If daytime population density is unknown, set to -999. 
	:SUEWS-related variables:
		PopDensDay
    


.. option:: popdensnighttime

	:Dimensionality:
		0
	:Description:
		Night-time population density (i.e. residents) [people ha -1 ] Population density is required if EmissionsMethod = 2 in `RunControl.nml` . The model will use the average of daytime and night-time population densities, unless only one is provided. If night-time population density is unknown, set to -999. 
	:SUEWS-related variables:
		PopDensNight
    


.. option:: popprof_24hr

	:Dimensionality:
		(24, 2)
	:Description:
		Code linking to `PopProfWD` in `SUEWS_Profiles.txt`.
Code linking to `PopProfWE` in `SUEWS_Profiles.txt`.
	:SUEWS-related variables:
		PopProfWD, PopProfWE
    


.. option:: pormax_dec

	:Dimensionality:
		0
	:Description:
		full leaf-on summertime value Used only for DecTr (can affect roughness calculation)
	:SUEWS-related variables:
		PorosityMax
    


.. option:: pormin_dec

	:Dimensionality:
		0
	:Description:
		leaf-off wintertime value Used only for DecTr (can affect roughness calculation)
	:SUEWS-related variables:
		PorosityMin
    


.. option:: porosity_id

	:Dimensionality:
		0
	:Description:
		internal use, DO NOT change.
	:SUEWS-related variables:
		nan
    


.. option:: preciplimit

	:Dimensionality:
		0
	:Description:
		Auer (1974) [Au74]_
	:SUEWS-related variables:
		PrecipLimSnow
    


.. option:: preciplimitalb

	:Dimensionality:
		0
	:Description:
		Limit for hourly precipitation when the ground is fully covered with snow. Then snow albedo is reset to AlbedoMax [mm]
	:SUEWS-related variables:
		PrecipLimAlb
    


.. option:: qf0_beu

	:Dimensionality:
		(2,)
	:Description:
		Weekday building energy use [W |m^-2|] Can be used for CO2 flux calculation.
Weekend building energy use [W |m^-2|] Can be used for CO2 flux calculation.
	:SUEWS-related variables:
		QF0_BEU_WD, QF0_BEU_WE
    


.. option:: qf_a

	:Dimensionality:
		(2,)
	:Description:
		Base value for QF on weekdays [W |m^-2| (Cap |ha^-1| |)^-1| ]
Base value for QF on weekends [W |m^-2| (Cap |ha^-1| |)^-1|]
	:SUEWS-related variables:
		QF_A_WD, QF_A_WE
    


.. option:: qf_b

	:Dimensionality:
		(2,)
	:Description:
		Parameter related to cooling degree days on weekdays [W |m^-2| |K^-1| (Cap |ha^-1| |)^-1|]
Parameter related to cooling degree days on weekends [W |m^-2| |K^-1| (Cap |ha^-1| |)^-1|]
	:SUEWS-related variables:
		QF_B_WD, QF_B_WE
    


.. option:: qf_c

	:Dimensionality:
		(2,)
	:Description:
		Parameter related to heating degree days on weekdays [W |m^-2| |K^-1| (Cap |ha^-1| |)^-1|]
Parameter related to heating degree days on weekends [W |m^-2| |K^-1| (Cap |ha^-1| |)^-1|]
	:SUEWS-related variables:
		QF_C_WD, QF_C_WE
    


.. option:: qn1_av

	:Dimensionality:
		0
	:Description:
		internal use, DO NOT change.
	:SUEWS-related variables:
		nan
    


.. option:: qn1_s_av

	:Dimensionality:
		0
	:Description:
		internal use, DO NOT change.
	:SUEWS-related variables:
		nan
    


.. option:: radmeltfact

	:Dimensionality:
		0
	:Description:
		Hourly radiation melt factor of snow [mm |w^-1| |h^-1|]
	:SUEWS-related variables:
		RadMeltFactor
    


.. option:: raincover

	:Dimensionality:
		0
	:Description:
		Limit when surface totally covered with water [mm] Used for LUMPS surface wetness control. Default recommended value of 1 mm from Loridan et al. (2011) [L2011]_ .
	:SUEWS-related variables:
		LUMPS_Cover
    


.. option:: rainmaxres

	:Dimensionality:
		0
	:Description:
		Maximum water bucket reservoir [mm] Used for LUMPS surface wetness control. Default recommended value of 10 mm from Loridan et al. (2011) [L2011]_ .
	:SUEWS-related variables:
		LUMPS_MaxRes
    


.. option:: resp_a

	:Dimensionality:
		(3,)
	:Description:
		Respiration coeﬃcient a
	:SUEWS-related variables:
		resp_a
    


.. option:: resp_b

	:Dimensionality:
		(3,)
	:Description:
		Respiration coeﬃcient b - related to air temperature dependency
	:SUEWS-related variables:
		resp_b
    


.. option:: roughlenheatmethod

	:Dimensionality:
		0
	:Description:
		Determines method for calculating roughness length for heat.
	:SUEWS-related variables:
		RoughLenHeatMethod
    


.. option:: roughlenmommethod

	:Dimensionality:
		0
	:Description:
		Determines how aerodynamic roughness length (z0m) and zero displacement height (zdm) are calculated.
	:SUEWS-related variables:
		RoughLenMomMethod
    


.. option:: runofftowater

	:Dimensionality:
		0
	:Description:
		Fraction of above-ground runoff flowing to water surface during flooding [-] Value must be in the range 0-1. Fraction of above-ground runoff that can flow to the water surface in the case of flooding.
	:SUEWS-related variables:
		RunoffToWater
    


.. option:: s1

	:Dimensionality:
		0
	:Description:
		Related to soil moisture dependence [-].

.. note::
    These will change in the future to ensure consistency with soil behaviour 
	:SUEWS-related variables:
		S1
    


.. option:: s2

	:Dimensionality:
		0
	:Description:
		Related to soil moisture dependence [mm].

.. note::
    These will change in the future to ensure consistency with soil behaviour
	:SUEWS-related variables:
		S2
    


.. option:: sathydraulicconduct

	:Dimensionality:
		(7,)
	:Description:
		Hydraulic conductivity for saturated soil [mm |s^-1|]
	:SUEWS-related variables:
		SatHydraulicCond
    


.. option:: sddfull

	:Dimensionality:
		(3,)
	:Description:
		This should be checked carefully for your study area using modelled LAI from the DailyState output file compared to known behaviour in the study area. See section 2.2 Järvi et al. (2011) [J11]_ ; Appendix A Järvi et al. (2014) [Leena2014]_ for more details. Example values [°C] -450 EveTr Järvi et al. (2011) [J11]_  -450 DecTr Järvi et al. (2011) [J11]_  -450 Grass Järvi et al. (2011) [J11]_ 
	:SUEWS-related variables:
		SDDFull
    


.. option:: sfr

	:Dimensionality:
		(7,)
	:Description:
		Surface cover fraction of buildings [-]
Surface cover fraction of bare soil or unmanaged land [-]
Surface cover fraction of deciduous trees and shrubs [-]
Surface cover fraction of evergreen trees and shrubs [-]
Surface cover fraction of grass [-]
Surface cover fraction of paved surfaces [-].
Surface cover fraction of open water [-] (e.g. river, lakes, ponds, swimming pools) 
	:SUEWS-related variables:
		Fr_Bldgs, Fr_Bsoil, Fr_DecTr, Fr_EveTr, Fr_Grass, Fr_Paved, Fr_Water
    


.. option:: smdmethod

	:Dimensionality:
		0
	:Description:
		Determines method for calculating soil moisture deficit (SMD).
	:SUEWS-related variables:
		SMDMethod
    


.. option:: snowalb

	:Dimensionality:
		0
	:Description:
		internal use, DO NOT change.
	:SUEWS-related variables:
		nan
    


.. option:: snowalbmax

	:Dimensionality:
		0
	:Description:
		Effective surface albedo (middle of the day value) for summertime. View factors should be taken into account.
Example values [-] 0.85 Järvi et al. (2014) [Leena2014]_ 
Effective surface albedo (middle of the day value) for summertime, full leaf-on. View factors should be taken into account. Example values [-] 0.1 EveTr Oke (1987) [Ok87]_  0.18 DecTr Oke (1987) [Ok87]_  0.21 Grass Oke (1987) [Ok87]_ 
Effective albedo of the water surface. View factors should be taken into account. Example values [-] 0.1 Water Oke (1987) [Ok87]_ 
	:SUEWS-related variables:
		AlbedoMax
    


.. option:: snowalbmin

	:Dimensionality:
		0
	:Description:
		Effective surface albedo (middle of the day value) for wintertime (not including snow). View factors should be taken into account.

.. note::
    Not currently used for non-vegetated surfaces – set the same as AlbedoMax.
Example values [-] 0.18 Järvi et al. (2014) [Leena2014]_ 
Effective surface albedo (middle of the day value) for wintertime (not including snow), leaf-off. View factors should be taken into account. Example values [-] 0.1 EveTr Oke (1987) [Ok87]_  0.18 DecTr Oke (1987) [Ok87]_  0.21 Grass Oke (1987) [Ok87]_ 
View factors should be taken into account. Not currently used for water surface - set same as AlbedoMax.
	:SUEWS-related variables:
		AlbedoMin
    


.. option:: snowd

	:Dimensionality:
		(7,)
	:Description:
		Meaning of this?
Limit of snow water equivalent when the surface surface is fully covered with snow. Not needed if `SnowUse` = 0 in `RunControl.nml` . Example values [mm] 190 EveTr Järvi et al. (2014) [Leena2014]_  190 DecTr Järvi et al. (2014) [Leena2014]_  190 Grass Järvi et al. (2014) [Leena2014]_ 
	:SUEWS-related variables:
		SnowLimPatch
    


.. option:: snowdens

	:Dimensionality:
		(7,)
	:Description:
		Initial snow density
	:SUEWS-related variables:
		nan
    


.. option:: snowdensmax

	:Dimensionality:
		0
	:Description:
		Maximum snow density [kg |m^-3|]
	:SUEWS-related variables:
		SnowDensMax
    


.. option:: snowdensmin

	:Dimensionality:
		0
	:Description:
		Fresh snow density [kg |m^-3|]
	:SUEWS-related variables:
		SnowDensMin
    


.. option:: snowfallcum

	:Dimensionality:
		0
	:Description:
		internal use, DO NOT change.
	:SUEWS-related variables:
		nan
    


.. option:: snowfrac

	:Dimensionality:
		(7,)
	:Description:
		Initial plan area fraction of snow
	:SUEWS-related variables:
		nan
    


.. option:: snowlimbuild

	:Dimensionality:
		0
	:Description:
		Meaning of this?
	:SUEWS-related variables:
		SnowLimRemove
    


.. option:: snowlimpaved

	:Dimensionality:
		0
	:Description:
		Meaning of this?
	:SUEWS-related variables:
		SnowLimRemove
    


.. option:: snowpack

	:Dimensionality:
		(7,)
	:Description:
		Initial snow water equivalent
	:SUEWS-related variables:
		nan
    


.. option:: snowprof_24hr

	:Dimensionality:
		(24, 2)
	:Description:
		Code for snow clearing profile (weekdays) Provides the link to column 1 of `SUEWS_Profiles.txt`. Value of integer is arbitrary but must match code specified in column 1 of `SUEWS_Profiles.txt`. e.g. 1 means use the characteristics specified in the row of input file `SUEWS_Profiles.txt` which has 1 in column 1 (Code).
Code for snow clearing profile (weekends) Provides the link to column 1 of `SUEWS_Profiles.txt`. Value of integer is arbitrary but must match code specified in column 1 of `SUEWS_Profiles.txt`. e.g. 1 means use the characteristics specified in the row of input file `SUEWS_Profiles.txt` which has 1 in column 1 (Code). Providing the same code for `SnowClearingProfWD` and `SnowClearingProfWE` would link to the same row in `SUEWS_Profiles.txt`, i.e. the same profile would be used for weekdays and weekends. 
	:SUEWS-related variables:
		SnowClearingProfWD, SnowClearingProfWE
    


.. option:: snowuse

	:Dimensionality:
		0
	:Description:
		Determines whether the snow part of the model runs.
	:SUEWS-related variables:
		SnowUse
    


.. option:: soildepth

	:Dimensionality:
		(7,)
	:Description:
		Depth of sub-surface soil store [mm] i.e. the depth of soil beneath the surface 
	:SUEWS-related variables:
		SoilDepth
    


.. option:: soilmoist_id

	:Dimensionality:
		(7,)
	:Description:
		internal use, DO NOT change.
	:SUEWS-related variables:
		nan
    


.. option:: soilstorecap

	:Dimensionality:
		(7,)
	:Description:
		SoilStoreCap must not be greater than SoilDepth.
	:SUEWS-related variables:
		SoilStoreCap
    


.. option:: stabilitymethod

	:Dimensionality:
		0
	:Description:
		Defines which atmospheric stability functions are used.
	:SUEWS-related variables:
		StabilityMethod
    


.. option:: startdls

	:Dimensionality:
		0
	:Description:
		Start of the day light savings [DOY] See `Day_Light_Savings`.
	:SUEWS-related variables:
		StartDLS
    


.. option:: state_id

	:Dimensionality:
		(7,)
	:Description:
		internal use, DO NOT change.
	:SUEWS-related variables:
		nan
    


.. option:: statelimit

	:Dimensionality:
		(7,)
	:Description:
		Currently only used for the water surface
Currently only used for the water surface
Surface state cannot exceed this value. Set to a large value (e.g. 20000 mm = 20 m) if the water body is substantial (lake, river, etc) or a small value (e.g. 10 mm) if water bodies are very shallow (e.g. fountains). WaterDepth (column 9) must not exceed this value.
	:SUEWS-related variables:
		StateLimit
    


.. option:: storageheatmethod

	:Dimensionality:
		0
	:Description:
		Determines method for calculating storage heat flux _QS.
	:SUEWS-related variables:
		StorageHeatMethod
    


.. option:: storedrainprm

	:Dimensionality:
		(6, 7)
	:Description:
		Coefficient used in `DrainageEq`
Example values DrainageEq 10 Coefficient D0 [mm |h^-1|] 3 Recommended [3] for Grass (irrigated) 0.013 Coefficient D0 [mm |h^-1|] 2 Recommended [3] for EveTr, DecTr, Grass (unirrigated)
Not currently used for water surface
Coefficient used in `DrainageEq`
Example values DrainageEq 3 Coefficient b [-] 3 Recommended [3] for Grass (irrigated) 1.71 Coefficient b [m|m^-1|] 2 Recommended [3] for EveTr, DecTr, Grass (unirrigated)
Not currently used for water surface
Formulation choice for drainage calculation.
Options 1 Falk and Niemczynowicz (1978) [FN78]_ 2 Halldin et al. (1979) [Ha79]_ (Rutter eqn corrected for c=0, see Calder & Wright (1986) [CW86]_ ) Recommended [3] for EveTr, DecTr, Grass (unirrigated) 3 Falk and Niemczynowicz (1978) [FN78]_ Recommended [3] for Grass (irrigated) Coefficients are specified in the following two columns.
Not currently used for water surface.
Maximum water storage capacity for upper surfaces (i.e. canopy) Min and max values are to account for seasonal variation (e.g. leaf-on/leaf-off differences for vegetated surfaces).

.. note::
    Not currently used for non-vegetated surfaces – set the same as `StorageMin`.

Maximum water storage capacity for upper surfaces (i.e. canopy) Min/max values are to account for seasonal variation (e.g. leaf-off/leaf-on differences for vegetated surfaces) Only used for DecTr surfaces - set EveTr and Grass values the same as StorageMin. Example values [mm] 1.3 EveTr Breuer et al. (2003) [Br03]_  0.8 DecTr Breuer et al. (2003) [Br03]_  1.9 Grass Breuer et al. (2003) [Br03]_ 
Maximum water storage capacity for upper surfaces (i.e. canopy) Min and max values are to account for seasonal variation - not used for water surfaces so set same as StorageMin.
Minimum water storage capacity for upper surfaces (i.e. canopy). Min/max values are to account for seasonal variation (e.g. leaf-on/leaf-off differences for vegetated surfaces).

.. note::
    Not currently used for non-vegetated surfaces – set the same as `StorageMax`.
Minimum water storage capacity for upper surfaces (i.e. canopy). Min/max values are to account for seasonal variation (e.g. leaf-off/leaf-on differences for vegetated surfaces). Example values [mm] 1.3 EveTr Breuer et al. (2003) [Br03]_  0.3 DecTr Breuer et al. (2003) [Br03]_  1.9 Grass Breuer et al. (2003) [Br03]_ 
Minimum water storage capacity for upper surfaces (i.e. canopy). Min/max values are to account for seasonal variation - not used for water surfaces. Example values [mm] 0.5 Water
	:SUEWS-related variables:
		DrainageCoef1, DrainageCoef2, DrainageEq, StorageMax, StorageMin
    


.. option:: surfacearea

	:Dimensionality:
		0
	:Description:
		Area of the grid [ha].
	:SUEWS-related variables:
		SurfaceArea
    


.. option:: t_critic_cooling

	:Dimensionality:
		(2,)
	:Description:
		Critical cooling temperature on weekdays [°C]
Critical cooling temperature on weekends [°C]
	:SUEWS-related variables:
		TCritic_Cooling_WD, TCritic_Cooling_WE
    


.. option:: t_critic_heating

	:Dimensionality:
		(2,)
	:Description:
		Critical heating temperature on weekdays [°C]
Critical heating temperature on weekends [°C]
	:SUEWS-related variables:
		TCritic_Heating_WD, TCritic_Heating_WE
    


.. option:: tair24hr

	:Dimensionality:
		(288,)
	:Description:
		internal use, DO NOT change.
	:SUEWS-related variables:
		nan
    


.. option:: tau_a

	:Dimensionality:
		0
	:Description:
		Time constant for snow albedo aging in cold snow [-]
	:SUEWS-related variables:
		tau_a
    


.. option:: tau_f

	:Dimensionality:
		0
	:Description:
		Time constant for snow albedo aging in melting snow [-]
	:SUEWS-related variables:
		tau_f
    


.. option:: tau_r

	:Dimensionality:
		0
	:Description:
		Time constant for snow density ageing [-]
	:SUEWS-related variables:
		tau_r
    


.. option:: tempmeltfact

	:Dimensionality:
		0
	:Description:
		Hourly temperature melt factor of snow [mm |K^-1| |h^-1|] (In previous model version, this parameter was 0.12)
	:SUEWS-related variables:
		TempMeltFactor
    


.. option:: th

	:Dimensionality:
		0
	:Description:
		Upper air temperature limit [°C]
	:SUEWS-related variables:
		TH
    


.. option:: theta_bioco2

	:Dimensionality:
		(3,)
	:Description:
		The convexity of the curve at light saturation.
	:SUEWS-related variables:
		theta
    


.. option:: timezone

	:Dimensionality:
		0
	:Description:
		Time zone [h] for site relative to UTC (east is positive). This should be set according to the times given in the meteorological forcing file(s).
	:SUEWS-related variables:
		Timezone
    


.. option:: tl

	:Dimensionality:
		0
	:Description:
		Lower air temperature limit [°C]
	:SUEWS-related variables:
		TL
    


.. option:: trafficrate

	:Dimensionality:
		(2,)
	:Description:
		Weekday traffic rate [veh km |m^-2| |s^-1|]. Can be used for CO2 flux calculation.
Weekend traffic rate [veh km |m^-2| |s^-1|]. Can be used for CO2 flux calculation.
	:SUEWS-related variables:
		TrafficRate_WD, TrafficRate_WE
    


.. option:: trafficunits

	:Dimensionality:
		0
	:Description:
		Traﬃc units choice.
	:SUEWS-related variables:
		TrafficUnits
    


.. option:: traffprof_24hr

	:Dimensionality:
		(24, 2)
	:Description:
		Code linking to `TraffProfWD` in `SUEWS_Profiles.txt`.
Code linking to `TraffProfWE` in `SUEWS_Profiles.txt`.
	:SUEWS-related variables:
		TraffProfWD, TraffProfWE
    


.. option:: tstep

	:Dimensionality:
		0
	:Description:
		internal use, DO NOT change.
	:SUEWS-related variables:
		nan
    


.. option:: tstep_prev

	:Dimensionality:
		0
	:Description:
		internal use, DO NOT change.
	:SUEWS-related variables:
		nan
    


.. option:: veg_type

	:Dimensionality:
		0
	:Description:
		internal use, DO NOT change.
	:SUEWS-related variables:
		nan
    


.. option:: waterdist

	:Dimensionality:
		(8, 6)
	:Description:
		Fraction of water going to `BSoil`
Fraction of water going to `Bldgs`
Fraction of water going to `DecTr`
Fraction of water going to `EveTr`
Fraction of water going to `Grass`
Fraction of water going to `Paved`
Fraction of water going to `Runoff`
Fraction of water going to `SoilStore`
Fraction of water going to `Water`
	:SUEWS-related variables:
		ToBSoil, ToBldgs, ToDecTr, ToEveTr, ToGrass, ToPaved, ToRunoff, ToSoilStore, ToWater
    


.. option:: waterusemethod

	:Dimensionality:
		0
	:Description:
		Defines how external water use is calculated.
	:SUEWS-related variables:
		WaterUseMethod
    


.. option:: wetthresh

	:Dimensionality:
		(7,)
	:Description:
		Depth of water which determines whether evaporation occurs from a partially wet or completely wet surface.
Depth of water which determines whether evaporation occurs from a partially wet or completely wet surface. Example values [mm] 1.8 EveTr 1. DecTr 2. Grass
Depth of water which determines whether evaporation occurs from a partially wet or completely wet surface. Example values [mm] 0.5 Water
	:SUEWS-related variables:
		WetThreshold
    


.. option:: wuday_id

	:Dimensionality:
		(9,)
	:Description:
		internal use, DO NOT change.
	:SUEWS-related variables:
		nan
    


.. option:: wuprofa_24hr

	:Dimensionality:
		(24, 2)
	:Description:
		Code for water use profile (automatic irrigation, weekdays) Provides the link to column 1 of `SUEWS_Profiles.txt`. Value of integer is arbitrary but must match code specified in column 1 of `SUEWS_Profiles.txt`. 
Code for water use profile (automatic irrigation, weekends) Provides the link to column 1 of `SUEWS_Profiles.txt`. Value of integer is arbitrary but must match code specified in column 1 of `SUEWS_Profiles.txt`. 
	:SUEWS-related variables:
		WaterUseProfAutoWD, WaterUseProfAutoWE
    


.. option:: wuprofm_24hr

	:Dimensionality:
		(24, 2)
	:Description:
		Code for water use profile (manual irrigation, weekdays) Provides the link to column 1 of `SUEWS_Profiles.txt`. Value of integer is arbitrary but must match code specified in column 1 of `SUEWS_Profiles.txt`. 
Code for water use profile (manual irrigation, weekends) Provides the link to column 1 of `SUEWS_Profiles.txt`. Value of integer is arbitrary but must match code specified in column 1 of `SUEWS_Profiles.txt`. 
	:SUEWS-related variables:
		WaterUseProfManuWD, WaterUseProfManuWE
    


.. option:: z

	:Dimensionality:
		0
	:Description:
		z must be greater than the displacement height. Forcing data should be representative of the local-scale, i.e. above the height of the roughness elements.
	:SUEWS-related variables:
		z
    


.. option:: z0m_in

	:Dimensionality:
		0
	:Description:
		Roughness length for momentum [m] Value supplied here is used if `RoughLenMomMethod` = 1 in `RunControl.nml` ; otherwise set to '-999' and a value will be calculated by the model (`RoughLenMomMethod` = 2, 3). 
	:SUEWS-related variables:
		z0
    


.. option:: zdm_in

	:Dimensionality:
		0
	:Description:
		Zero-plane displacement [m] Value supplied here is used if `RoughLenMomMethod` = 1 in `RunControl.nml` ; otherwise set to '-999' and a value will be calculated by the model (`RoughLenMomMethod` = 2, 3). 
	:SUEWS-related variables:
		zd
    
