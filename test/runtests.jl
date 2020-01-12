using Test, BioGeoJulia
using BioGeoJulia.BGJ_Example

@testset "BGJ_Example" begin
	@test hello("Julia") == "Hello, Julia"
	@test domath(2.0) ≈ 7.0
end


@testset "BioGeoJulia" begin
	@test hello_BioGeoJulia("Julia") == "BioGeoJulia says, hi Julia"
	@test add_one_BioGeoJulia(2.0) ≈ 3.0
end