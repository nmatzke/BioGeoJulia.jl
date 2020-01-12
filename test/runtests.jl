using Test, BGJ_Example, BioGeoJulia

@test hello("Julia") == "Hello, Julia"
@test domath(2.0) ≈ 7.0


@test hello_BioGeoJulia("Julia") == "BioGeoJulia says, hi Julia"
@test add_one_BioGeoJulia(2.0) ≈ 3.0
