@testset "JET Testing" begin
    rep = JET.test_package(
        SatelliteToolboxAtmosphericModels;
        toplevel_logger=nothing,
        target_modules=(@__MODULE__,),
    )
end
