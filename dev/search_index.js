var documenterSearchIndex = {"docs":
[{"location":"lib/library/#Library","page":"Library","title":"Library","text":"","category":"section"},{"location":"lib/library/","page":"Library","title":"Library","text":"Documentation for AstroForceModels.jl.","category":"page"},{"location":"lib/library/","page":"Library","title":"Library","text":"Modules = [AstroForceModels]","category":"page"},{"location":"lib/library/#AstroForceModels.AbstractSatelliteDragModel","page":"Library","title":"AstroForceModels.AbstractSatelliteDragModel","text":"Abstract satellite drag model to help compute drag forces.\n\n\n\n\n\n","category":"type"},{"location":"lib/library/#AstroForceModels.AbstractSatelliteSRPModel","page":"Library","title":"AstroForceModels.AbstractSatelliteSRPModel","text":"Abstract satellite drag model to help compute drag forces.\n\n\n\n\n\n","category":"type"},{"location":"lib/library/#AstroForceModels.CannonballFixedDrag","page":"Library","title":"AstroForceModels.CannonballFixedDrag","text":"Cannonball Fixed Drag Model struct Contains information to compute the ballistic coefficient of a cannonball drag model with a  fixed ballistic coefficient.\n\nFields\n\nradius::Number: The radius of the spacecraft.\nmass::Number: The mass of the spacecraft.\ndrag_coeff::Number: The drag coefficient of the spacecraft.\nballistic_coeff::Number: The fixed ballistic coefficient to use.\n\n\n\n\n\n","category":"type"},{"location":"lib/library/#AstroForceModels.CannonballFixedDrag-Tuple{Number, Number, Number}","page":"Library","title":"AstroForceModels.CannonballFixedDrag","text":"CannonballFixedDrag(radius::Number, mass::Number, drag_coeff::Number)\n\nConstructor for a fixed cannonball ballistic coefficient drag model.\n\nThe ballistic coefficient is computed with the following equation:\n\n            BC = CD * area/mass\n\nwhere area is the 2D projection of a sphere\n\n            area = π * r^2\n\nArguments\n\nradius::Number: The radius of the spacecraft.\nmass::Number: The mass of the spacecraft.\ndrag_coeff::Number: The drag coefficient of the spacecraft.\n\nReturns\n\ndrag_model::CannonballFixedDrag: A fixed ballistic coefficient drag model.\n\n\n\n\n\n","category":"method"},{"location":"lib/library/#AstroForceModels.CannonballFixedDrag-Tuple{Number}","page":"Library","title":"AstroForceModels.CannonballFixedDrag","text":"CannonballFixedDrag(ballistic_coeff::Number)\n\nConstructor for a fixed ballistic coefficient drag model.\n\nArguments\n\nballistic_coeff::Number: The fixed ballistic coefficient to use.\n\nReturns\n\ndrag_model::CannonballFixedDrag: A fixed ballistic coefficient drag model.\n\n\n\n\n\n","category":"method"},{"location":"lib/library/#AstroForceModels.CannonballFixedSRP","page":"Library","title":"AstroForceModels.CannonballFixedSRP","text":"Cannonball Fixed SRP Model struct Contains information to compute the reflectivity coefficient of a cannonball drag model with a  fixed reflectivity coefficient.\n\nFields\n\nradius::Number: The radius of the spacecraft.\nmass::Number: The mass of the spacecraft.\nreflectivity_coeff::Number: The reflectivity coefficient of the spacecraft.\nreflectivity_ballistic_coeff::Number: The fixed ballistic coefficient to use.\n\n\n\n\n\n","category":"type"},{"location":"lib/library/#AstroForceModels.CannonballFixedSRP-Tuple{Number, Number, Number}","page":"Library","title":"AstroForceModels.CannonballFixedSRP","text":"CannonballFixedSRP(radius::Number, mass::Number, drag_coeff::Number)\n\nConstructor for a fixed cannonball ballistic coefficient SRP model.\n\nThe ballistic coefficient is computed with the following equation:\n\n            RC = CD * area/mass\n\nwhere area is the 2D projection of a sphere\n\n            area = π * r^2\n\nArguments\n\nradius::Number: The radius of the spacecraft.\nmass::Number: The mass of the spacecraft.\nreflectivity_coeff::Number: The reflectivity coefficient of the spacecraft.\n\nReturns\n\nsrp_model::CannonballFixedSRP`: A fixed ballistic coefficient SRP model.\n\n\n\n\n\n","category":"method"},{"location":"lib/library/#AstroForceModels.CannonballFixedSRP-Tuple{Number}","page":"Library","title":"AstroForceModels.CannonballFixedSRP","text":"CannonballFixedSRP(ballistic_coeff::Number)\n\nConstructor for a fixed ballistic coefficient SRP model.\n\nArguments\n\nballistic_coeff::Number: The fixed ballistic coefficient to use.\n\nReturns\n\nsrp_model::CannonballFixedDrag: A fixed ballistic coefficient SRP model.\n\n\n\n\n\n","category":"method"},{"location":"lib/library/#AstroForceModels.CelestialBody-Tuple{Val{:Sun}}","page":"Library","title":"AstroForceModels.CelestialBody","text":"Creates a CelestialBody Object from a name Available: {Sun, Mercury, Venus, Earth, Mars, Jupiter, Saturn, Uranus, Neptune, Pluto, }\n\nArguments:\n\nname::String: Name of the Planet Object (Casing Doesn't Matter)\n(Optional) T::DataType: Datatype of PlanetBody Fields\n\nReturns: -CelestialBody: PlanetBody Object with Fields:     - name::String: Name of Object     - central_body::String: Name of Central Body Being Orbited     - jpl_code::Int: NAIF ID Code     - μ::T: Graviational Parameter [km/s]     - Req::T: Equatorial Radius [km]\n\n\n\n\n\n","category":"method"},{"location":"lib/library/#AstroForceModels.DragAstroModel","page":"Library","title":"AstroForceModels.DragAstroModel","text":"Drag Astro Model struct Contains information to compute the acceleration of a drag force on a spacecraft.\n\nFields\n\nsatellite_drag_model::AbstractSatelliteDragModel: The satellite drag model for computing the ballistic coefficient.\natmosphere_model::Symbol: The atmospheric model for computing the density.\neop_data::EopIau1980: Earth orientation parameters to help compute the density with the atmospheric model.\n\n\n\n\n\n","category":"type"},{"location":"lib/library/#AstroForceModels.SRPAstroModel","page":"Library","title":"AstroForceModels.SRPAstroModel","text":"SRP Astro Model struct Contains information to compute the acceleration of a SRP a spacecraft.\n\nFields\n\nsatellite_srp_model::AbstractSatelliteDragModel: The satellite srp model for computing the ballistic coefficient.\nsun_data::ThirdBodyModel: The data to compute the Sun's position.\n\n\n\n\n\n","category":"type"},{"location":"lib/library/#AstroForceModels.StateDragModel","page":"Library","title":"AstroForceModels.StateDragModel","text":"State Drag Model struct Empty struct used for when the simulation state includes the ballistic coefficient.\n\n\n\n\n\n","category":"type"},{"location":"lib/library/#AstroForceModels.StateSRPModel","page":"Library","title":"AstroForceModels.StateSRPModel","text":"State SRP Model struct Empty struct used for when the simulation state includes the reflectivity ballistic coefficient.\n\n\n\n\n\n","category":"type"},{"location":"lib/library/#AstroForceModels.ThirdBodyModel","page":"Library","title":"AstroForceModels.ThirdBodyModel","text":"Third Body Model Astro Model struct Contains information to compute the acceleration of a third body force acting on a spacecraft.\n\nFields\n\nbody::CelestialBody: Celestial body acting on the craft.\nephem_type::Symbol: Ephemeris type used to compute body's position. Options are currently :Vallado.\n\n\n\n\n\n","category":"type"},{"location":"lib/library/#AstroForceModels.ThirdBodyModel-Tuple{Number}","page":"Library","title":"AstroForceModels.ThirdBodyModel","text":"Convenience to compute the ephemeris position of a CelestialBody in a ThirdBodyModel Wraps get_position().\n\nArguments\n\ntime::Number: Current time of the simulation in seconds.\n\nReturns\n\nbody_position: SVector{3}: The 3-dimensional third body position in the J2000 frame.\n\n\n\n\n\n","category":"method"},{"location":"lib/library/#AstroForceModels.acceleration-Tuple{AbstractArray, ComponentArrays.ComponentVector, Number, DragAstroModel}","page":"Library","title":"AstroForceModels.acceleration","text":"'     acceleration(u::AbstractArray, p::ComponentVector, t::Number, drag_model::DragAstroModel)\n\nComputes the drag acceleration acting on a spacecraft given a drag model and current state and  parameters of an object.\n\nArguments\n\nu::AbstractArray: Current State of the simulation.\np::ComponentVector: Current parameters of the simulation.\nt::Number: Current time of the simulation.\ndrag_model::DragAstroModel: Drag model struct containing the relevant information to compute the acceleration.\n\nReturns\n\nacceleration: SVector{3}: The 3-dimensional drag acceleration acting on the spacecraft.\n\n\n\n\n\n","category":"method"},{"location":"lib/library/#AstroForceModels.acceleration-Tuple{AbstractArray, ComponentArrays.ComponentVector, Number, SRPAstroModel}","page":"Library","title":"AstroForceModels.acceleration","text":"acceleration(u::AbstractArray, p::ComponentVector, t::Number, srp_model::SRPAstroModel)\n\nComputes the drag acceleration acting on a spacecraft given a drag model and current state and  parameters of an object.\n\nArguments\n\nu::AbstractArray: Current State of the simulation.\np::ComponentVector: Current parameters of the simulation.\nt::Number: Current time of the simulation.\nsrp_model::SRPAstroModel: SRP model struct containing the relevant information to compute the acceleration.\n\nReturns\n\nacceleration: SVector{3}: The 3-dimensional srp acceleration acting on the spacecraft.\n\n\n\n\n\n","category":"method"},{"location":"lib/library/#AstroForceModels.acceleration-Tuple{AbstractArray, ComponentArrays.ComponentVector, Number, ThirdBodyModel}","page":"Library","title":"AstroForceModels.acceleration","text":"acceleration(u::AbstractArray, p::ComponentVector, t::Number, third_body_model::ThirdBodyModel)\n\nComputes the drag acceleration acting on a spacecraft given a drag model and current state and  parameters of an object.\n\nArguments\n\nu::AbstractArray: Current State of the simulation.\np::ComponentVector: Current parameters of the simulation.\nt::Number: Current time of the simulation.\nthird_body_model::ThirdBodyModel: Third body model struct containing the relevant information to compute the acceleration.\n\nReturns\n\nacceleration: SVector{3}: The 3-dimensional srp acceleration acting on the spacecraft.\n\n\n\n\n\n","category":"method"},{"location":"lib/library/#AstroForceModels.ballistic_coefficient-Tuple{AbstractArray, AbstractVector, Number, CannonballFixedDrag}","page":"Library","title":"AstroForceModels.ballistic_coefficient","text":"ballistic_coefficient(     u::AbstractArray,      p::ComponentVector,      t::Number,      model::CannonballFixedDrag)\n\nReturns the ballistic coefficient for a drag model given the model and current state of the simulation.\n\nArguments\n\n- `u::AbstractArray`: The current state of the simulation.\n- `p::ComponentVector`: The parameters of the simulation.\n- `t::Number`: The current time of the simulation.\n- `model::CannonballFixedDrag`: Drag model for the spacecraft.\n\nReturns\n\n-`ballistic_coeff::Number`: The current ballistic coefficient of the spacecraft.\n\n\n\n\n\n","category":"method"},{"location":"lib/library/#AstroForceModels.ballistic_coefficient-Tuple{AbstractArray, AbstractVector, Number, StateDragModel}","page":"Library","title":"AstroForceModels.ballistic_coefficient","text":"ballistic_coefficient(     u::AbstractArray,      p::ComponentVector,      t::Number,      model::StateDragModel)\n\nReturns the ballistic coefficient for a drag model given the model and current state  of the simulation.\n\nArguments\n\n- `u::AbstractArray`: The current state of the simulation.\n- `p::ComponentVector`: The parameters of the simulation.\n- `t::Number`: The current time of the simulation.\n- `model::StateDragModel`: The drag model for the spacecraft.\n\nReturns\n\n-`ballistic_coeff::Number`: The current ballistic coefficient of the spacecraft.\n\n\n\n\n\n","category":"method"},{"location":"lib/library/#AstroForceModels.compute_density","page":"Library","title":"AstroForceModels.compute_density","text":"compute_density(JD::Number, u::AbstractArray, eop_data::EopIau1980, Val(AtmosphereType::Symbol))\n\nComputes the atmospheric density at a point given the date, position, eop_data, and atmosphere type\n\nArguments\n\nJD::Number: The current time of the simulation in Julian days.\nu::AbstractArray: The current state of the simulation.\neop_data::EopIau1980: The earth orientation parameters.\nAtmosphereType::Symbol: The type of atmospheric model used to compute the density. Available    options are Jacchia-Bowman 2008 (:JB2008), Jacchia-Roberts 1971 (:JR1971), NRL MSIS 2000 (:MSIS2000),   Exponential (:ExpAtmo), and None (:None)\n\nReturns\n\nrho::Number: The density of the atmosphere at the provided time and point [kg/m^3].\n\n\n\n\n\n","category":"function"},{"location":"lib/library/#AstroForceModels.drag_accel-Tuple{AbstractArray, Number, Number, AbstractArray, Number}","page":"Library","title":"AstroForceModels.drag_accel","text":"drag_accel(u::AbstractArray, rho::Number, BC::Number, ω_vec::AbstractArray, t::Number, [DragModel]) -> SVector{3}{Number}\n\nCompute the Acceleration Atmospheric Drag\n\nThe atmosphere is treated as a solid revolving with the Earth and the apparent velocity of the satellite is computed using the transport theorem\n\n            𝐯_app = 𝐯 - 𝛚 x 𝐫\n\nThe acceleration from drag is then computed with a cannonball model as\n\n            𝐚 = 1/2 * ρ * BC * |𝐯_app|₂^2 * v̂\n\nnote: Note\nCurrently only fixed cannonball state based ballistic coefficients are supported, custom models can be created for higher fidelity.\n\nArguments\n\nu::AbstractArray: The current state of the spacecraft in the central body's inertial frame.\nrho::Number: Atmospheric density at (t, u) [kg/m^3].\nBC::Number: The ballistic coefficient of the satellite – (area/mass) * drag coefficient [kg/m^2].\nω_vec::AbstractArray: The angular velocity vector of Earth. Typically approximated as [0.0; 0.0; ω_Earth]\nt::Number: Current time of the simulation.\n\nReturns\n\nSVector{3}{Number}: Inertial acceleration from drag\n\n\n\n\n\n","category":"method"},{"location":"lib/library/#AstroForceModels.get_position-Union{Tuple{T}, Tuple{Val{:Vallado}, CelestialBody, T, Number}} where T<:Union{Nothing, SatelliteToolboxTransformations.EopIau1980, SatelliteToolboxTransformations.EopIau2000A}","page":"Library","title":"AstroForceModels.get_position","text":"Computes the position of the celestial body using Vallado's ephemeris\n\nArguments\n\nephem_type::Symbol: Ephemeris type used to compute body's position.\nbody::CelestialBody: Celestial body acting on the craft.\ntime::Number: Current time of the simulation in seconds.\n\nReturns\n\nbody_position: SVector{3}: The 3-dimensional third body position in the J2000 frame.\n\n\n\n\n\n","category":"method"},{"location":"lib/library/#AstroForceModels.reflectivity_ballistic_coefficient-Tuple{AbstractArray, AbstractVector, Number, CannonballFixedSRP}","page":"Library","title":"AstroForceModels.reflectivity_ballistic_coefficient","text":"reflectivityballisticcoefficient(     u::AbstractArray,      p::ComponentVector,      t::Number,      model::CannonballFixedSRP)\n\nReturns the ballistic coefficient for a SRP model given the model and current state of the simulation.\n\nArguments\n\n- `u::AbstractArray`: The current state of the simulation.\n- `p::ComponentVector`: The parameters of the simulation.\n- `t::Number`: The current time of the simulation.\n- `model::CannonballFixedSRP`: SRP model for the spacecraft.\n\nReturns\n\n-`ballistic_coeff::Number`: The current ballistic coefficient of the spacecraft.\n\n\n\n\n\n","category":"method"},{"location":"lib/library/#AstroForceModels.reflectivity_ballistic_coefficient-Tuple{AbstractArray, AbstractVector, Number, StateSRPModel}","page":"Library","title":"AstroForceModels.reflectivity_ballistic_coefficient","text":"ballistic_coefficient(     u::AbstractArray,      p::ComponentVector,      t::Number,      model::StateSRPModel)\n\nReturns the ballistic coefficient for a SRP model given the model and current state  of the simulation.\n\nArguments\n\n- `u::AbstractArray`: The current state of the simulation.\n- `p::ComponentVector`: The parameters of the simulation.\n- `t::Number`: The current time of the simulation.\n- `model::StateSRPModel`: The SRP model for the spacecraft.\n\nReturns\n\n-`ballistic_coeff::Number`: The current ballistic coefficient of the spacecraft.\n\n\n\n\n\n","category":"method"},{"location":"lib/library/#AstroForceModels.relativity_accel-Tuple{AbstractArray, Number, Number, AbstractArray, AbstractArray, Number, AbstractArray}","page":"Library","title":"AstroForceModels.relativity_accel","text":"\n\n\n\n","category":"method"},{"location":"lib/library/#AstroForceModels.shadow_model-Tuple{AbstractArray, AbstractArray, Val{:Cylindrical}}","page":"Library","title":"AstroForceModels.shadow_model","text":"shadow_model(sat_pos::AbstractArray, sun_pos::AbstractArray, R_Sun::Number, R_Earth::Number, t::Number, ShadowModel::Symbol)\n\nComputes the Lighting Factor of the Sun occur from the Umbra and Prenumbra of Earth's Shadow\n\nArguments\n\nsat_pos::AbstractArray: The current satellite position.\nsun_pos::AbstractArray: The current Sun position.\nR_Sun::Number: The radius of the Sun.\nR_Earth::Number: The radius of the Earth.\nShadowModel::Symbol: The Earth shadow model to use. Current Options – :Cylindrical, :Conical, :Conical_Simplified, :None\n\nReturns\n\nSVector{3}{Number}: Inertial acceleration from the 3rd body\n\n\n\n\n\n","category":"method"},{"location":"lib/library/#AstroForceModels.srp_accel-Tuple{AbstractArray, AbstractArray, Number}","page":"Library","title":"AstroForceModels.srp_accel","text":"srp_accel(u::AbstractArray, sun_pos::AbstractArray, R_Sun::Number, R_Earth::Number, Ψ::Number, RC::Number, t::Number, [SRPModel::Symbol]; ShadowModel::Symbol)\n\nCompute the Acceleration from Solar Radiaiton Pressure\n\nRadiation from the Sun reflects off the satellite's surface and transfers momentum perturbing the satellite's trajectory. This force can be computed using the a Cannonball model with the following equation\n\n            𝐚 = F * RC * Ψ * (AU/(R_sc_Sun))^2 * R̂_sc_Sun\n\nnote: Note\nCurrently only Cannonball SRP is supported, to use a higher fidelity drag either use a state varying function or compute the ballistic coefficient further upstream\n\nArguments\n\nu::AbstractArray: The current state of the spacecraft in the central body's inertial frame.\nsun_pos::AbstractArray: The current position of the Sun.\nR_Sun::Number: The radius of the Sun.\nR_Earth::Number: The radius of the Earth.\nΨ::Number: Solar Constant at 1 Astronomical Unit.\nRC::Number: The solar ballistic coefficient of the satellite – (Area/mass) * Reflectivity Coefficient [kg/m^2].\nt::Number: The current time of the Simulation\n\nOptional Arguments\n\nSRPModel::Symbol: SRP Model to use. Current Options – :Cannonball, :None\nShadowModel::Symbol: SRP Earth Shadow Model to use. Current Options – :Conical, :Conical_Simplified, :Cylinderical\n\nReturns\n\nSVector{3}{Number}: Inertial acceleration from the 3rd body\n\n\n\n\n\n","category":"method"},{"location":"lib/library/#AstroForceModels.third_body_accel-Tuple{AbstractArray, Number, AbstractArray}","page":"Library","title":"AstroForceModels.third_body_accel","text":"third_body_accel(u::AbstractArray, μ_body::Number, body_pos::AbstractArray, h::Number) -> SVector{3}{Number}\n\nCompute the Acceleration from a 3rd Body Represented as a Point Mass\n\nSince the central body is also being acted upon by the third body, the acceleration of body 𝐁 acting on  spacecraft 𝐀 in the orbiting body's 𝐂 is part of the force not acting on the central body\n\n            a = ∇UB(rA) - ∇UB(rC)\n\nArguments\n\nu::AbstractArray: The current state of the spacecraft in the central body's inertial frame.\nμ_body: Gravitation Parameter of the 3rd body.\nbody_pos::AbstractArray: The current position of the 3rd body in the central body's inertial frame.\n\nReturns\n\nSVector{3}{Number}: Inertial acceleration from the 3rd body\n\n\n\n\n\n","category":"method"},{"location":"man/usage/","page":"Usage","title":"Usage","text":"#TODO","category":"page"},{"location":"#AstroForceModels.jl","page":"Home","title":"AstroForceModels.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"This package contains the dominant astrodynamics forces affecting the orbital trajectory of a satellite for the SatelliteToolbox.jl ecosystem. Currently this package implements:","category":"page"},{"location":"","page":"Home","title":"Home","text":"[] Zonal Harmonics\n[] Solar Radiation Pressure\n[] Drag\n[] Third Body Gravity\n[] Relativistic\n[] Albedo\n[] Solid Tides","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"julia> using Pkg\njulia> Pkg.add(\"AstroForceModels\")","category":"page"},{"location":"man/api/","page":"API","title":"API","text":"#TODO","category":"page"}]
}
