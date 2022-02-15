### A Pluto.jl notebook ###
# v0.18.0

using Markdown
using InteractiveUtils

# ╔═╡ 6b1ad54f-61e4-490d-9032-7a557e8dc82f
md"""
## CHEME 5440/7770: Structural Analysis of the Urea Cycle Network (PS2)
"""

# ╔═╡ 7057c8e4-9e94-4a28-a885-07f5c96ebe39
html"""
<p style="font-size:20px;">Student name, Student name, Student name ... Student name</br>
Smith School of Chemical and Biomolecular Engineering, Cornell University, Ithaca NY 14850</p>
"""

# ╔═╡ f1e6f146-2984-4bff-8209-68dbf9e56192
#Kathleen Clifford Problem Set 2

# ╔═╡ 87a183bc-3857-4189-8103-18c46ff3245d
md"""
#### Build the stoichiometric array
"""

# ╔═╡ 4b647cd7-28d7-4630-8d61-63c0ead5c639


# ╔═╡ 630c76c7-1a17-4fe2-8770-9fd3f58932b4


# ╔═╡ 6970dab5-16bd-4898-b88d-723cb1b3d89e
md"""
#### Convex analysis: extreme pathway
"""

# ╔═╡ 8689cec2-f26a-438d-80d1-765cf07d90b2
#There were 5 extreme pathways  (rows of P), and 2 produced Urea, as seen by the fact that there was one non zero value in column 3, which was reaction 3 which produced urea.

#The reaction frequency for each reaction was 
#v1-3/5, v2-3/5, v3-2/5, v4-2/5, Fv5- 3/5, Rv5- 1/5, b1-2/5, b2-3/5, b3-3/5, b4-2/5, b5-3/5, Fb6-2/5, Rb6-2/5, b7-2/5, b8-2/5, b9-2/5, b10-2/5, b11- 2/5, b12-3/5, b13- 3/5, b14-2/5


# ╔═╡ 8f56cf8b-084d-472f-b578-a8554ad5788e
pidx = 7

# ╔═╡ b473b17e-3bf5-4b6c-af24-fe57b5a7e7e9
md"""
#### Metabolite connectivity array (MCA)
"""

# ╔═╡ ea09c6d6-1485-4d86-9eed-8a3418199148
#The diagnol of the MCA gives the number of reactions a particular metabolite participates in. In the following lines, the connectivity of the metabolites is ranked from most to least.

#"H20","Larginine","Lcitrulline","H","NADP","NADPH","O2","nitricoxide","2(NomegaLarginino)succinate","AMP","ATP","Lornithine","carbamoylphosphate","diphosphate","fumarate","phosphate","urea"


# ╔═╡ b7e5d1a6-57ed-4d09-a039-a4bd12386367
md"""
#### Reaction connectivity array (RCA)
"""

# ╔═╡ a958fbba-c172-4858-80ab-9f98dbc42d6c


# ╔═╡ 3a25bcad-7529-41c3-bfd0-90db0f192b38
#The diagnol of the RCA gives the number of metabolites participating in a reaction. The connectivity of the reactions are listed below, in order from most to least.

#Fv5,Rv5,v1,v3,v4,v2,"b₁","b₂","b₃","b₄","b₅","Fb₆","Rb₆","b₇","b₈","b₉","b₁₀","b₁₁","b₁₂","b₁₃","b₁₄"


# ╔═╡ 111bebcc-6fe2-440b-860c-9d18e092a617
#The correlation between reaction connectivity and extreme pathway reaction frequency was not able to be determined, I did not see a correlation. 

# ╔═╡ 267865de-1b5c-4579-861b-c6c46beb4739
function ingredients(path::String)
	
	# this is from the Julia source code (evalfile in base/loading.jl)
	# but with the modification that it returns the module instead of the last object
	name = Symbol("lib")
	m = Module(name)
	Core.eval(m,
        Expr(:toplevel,
             :(eval(x) = $(Expr(:core, :eval))($name, x)),
             :(include(x) = $(Expr(:top, :include))($name, x)),
             :(include(mapexpr::Function, x) = $(Expr(:top, :include))(mapexpr, $name, x)),
             :(include($path))))
	m
end

# ╔═╡ 67f5db98-88d0-11ec-27ac-b57538a166f4
begin
	# import some packages -
	using PlutoUI
	using PrettyTables
	using LinearAlgebra
	
	# setup paths -
	const _PATH_TO_NOTEBOOK = pwd()
	const _PATH_TO_SRC = joinpath(_PATH_TO_NOTEBOOK,"src")

	# load the PS2 code lib -
	lib = ingredients(joinpath(_PATH_TO_SRC, "Include.jl"));

	# return -
	nothing
end

# ╔═╡ 5338451e-3c4b-4030-bbbb-42eaf4209a89
begin


	# Setup a collection of reaction strings -
	reaction_array = Array{String,1}()

	# Setup a collection of reaction strings -
	reaction_array = Array{String,1}()

	# encode the reactions -
	#internal reactions
	push!(reaction_array,"v₁,ATP+LAspartate+Lcitrulline,AMP+diphosphate+2(NomegaLarginino)succinate,false")
	push!(reaction_array,"v₂,2(NomegaLarginino)succinate,fumarate+Larginine,false")
	push!(reaction_array,"v₃,Larginine+H20,Lornithine+urea,false")
	push!(reaction_array,"v₄,carbamoylphosphate+Lornithine,phosphate+Lcitrulline,false")
	push!(reaction_array,"v₅,2*Larginine+3*NADPH+3*H+4*O2,2*Lcitrulline+2*nitricoxide+3*NADP+4*H20,true")

	# exchange reactions -
	push!(reaction_array,"b₁,∅,carbamoylphosphate,false")
	push!(reaction_array,"b₂,∅,LAspartate,false")
	push!(reaction_array,"b₃,fumarate,∅,false")
	push!(reaction_array,"b₄,urea,∅,false")
	push!(reaction_array,"b₅,∅,ATP,false")
	push!(reaction_array,"b₆,∅,H20,true")
	push!(reaction_array,"b₇,∅,NADPH,false")
	push!(reaction_array,"b₈,∅,H,false")
	push!(reaction_array,"b₉,∅,O2,false")
	push!(reaction_array,"b₁₀,nitricoxide,∅,false")
	push!(reaction_array,"b₁₁,NADP,∅,false")
	push!(reaction_array,"b₁₂,AMP,∅,false")
	push!(reaction_array,"b₁₃,diphosphate,∅,false")
	push!(reaction_array,"b₁₄,phosphate,∅,false")
	
	

	# compute the stoichiometric matrix -
	(S, species_array, reaction_name_array) = lib.build_stoichiometric_matrix(reaction_array; expand=true);


	# show -
	nothing

	
end

# ╔═╡ 6a71c140-fd0a-4506-92b0-fe802fdca18a
(ℳ,ℛ) = size(S)

# ╔═╡ 36b15aa8-ca13-41e8-9bb6-a102b7e7b0b4
species_array

# ╔═╡ b5abeda2-9f63-4aa9-943e-c7f4ad0bfef1
S

# ╔═╡ 00bdd1d7-655d-4500-bfdd-499069e56490
reaction_name_array

# ╔═╡ 670ef5ac-aaae-4d2e-9d31-135417a8b198
species_array

# ╔═╡ 6e56cc9c-9fa8-405c-bfdd-191a02bee4aa
reaction_name_array

# ╔═╡ 9b2032ee-232f-4de5-96c0-3bcc9b78832c
B = S |> lib.binary_stoichiometric_matrix

# ╔═╡ 97b0763d-dcab-4afa-b660-52e18b3d523f
begin

	# compute the extreme pathways Tableu -
	PM = lib.expa(S)
	
	# P constaints the extreme pathways (rows) and 𝒩 is the "balanced" array (should be all zeros) -
	P = PM[:,1:ℛ]
	𝒩 = PM[:,(ℛ+1):end]

	# show -
	nothing
	
end

# ╔═╡ 80ea58f8-a77c-415e-92cb-86b836f0566e
P

# ╔═╡ 00b7bb48-c476-4d9b-940d-ca5dbed3477d
P[:,21]

# ╔═╡ 21b5a75a-3cb7-4d05-a8b2-639f3974febc
size(P)

# ╔═╡ 33bcfdc7-5dd4-4984-9f23-858a14e8bc10
Print(P)

# ╔═╡ 999ae1fd-5341-4f66-9db2-dec53fa0cd49
MCA = B*transpose(B)

# ╔═╡ e40b7e33-fae3-40a1-9a09-542483476c8f
diag(MCA)

# ╔═╡ 4520fc6e-7305-487e-924d-af22406e6d45
begin

RCA = transpose(B)*B
end

# ╔═╡ 5ba2c405-11be-4914-93be-379ed45eb63d
diag(RCA)

# ╔═╡ ab2bcfd5-3ba7-4388-8a3c-2cb95fba989a
html"""
<style>
main {
    max-width: 900px;
    width: 75%;
    margin: auto;
    font-family: "Roboto, monospace";
}

a {
    color: blue;
    text-decoration: none;
}
</style>"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
PrettyTables = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"

[compat]
PlutoUI = "~0.7.34"
PrettyTables = "~1.3.1"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.1"
manifest_format = "2.0"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "8eaf9f1b4921132a4cff3f36a1d9ba923b14a481"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.4"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "024fe24d83e4a5bf5fc80501a314ce0d1aa35597"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.DataAPI]]
git-tree-sha1 = "cc70b17275652eb47bc9e5f81635981f13cea5c8"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.9.0"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[deps.HypertextLiteral]]
git-tree-sha1 = "2b078b5a615c6c0396c77810d92ee8c6f470d238"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.3"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "8076680b162ada2a031f707ac7b4953e30667a37"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.2"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"

[[deps.Parsers]]
deps = ["Dates"]
git-tree-sha1 = "0b5cfbb704034b5b4c1869e36634438a047df065"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.2.1"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "UUIDs"]
git-tree-sha1 = "8979e9802b4ac3d58c503a20f2824ad67f9074dd"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.34"

[[deps.PrettyTables]]
deps = ["Crayons", "Formatting", "Markdown", "Reexport", "Tables"]
git-tree-sha1 = "dfb54c4e414caa595a1f2ed759b160f5a3ddcba5"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "1.3.1"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "TableTraits", "Test"]
git-tree-sha1 = "bb1064c9a84c52e277f1096cf41434b675cd368b"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.6.1"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
"""

# ╔═╡ Cell order:
# ╟─6b1ad54f-61e4-490d-9032-7a557e8dc82f
# ╟─7057c8e4-9e94-4a28-a885-07f5c96ebe39
# ╠═f1e6f146-2984-4bff-8209-68dbf9e56192
# ╟─87a183bc-3857-4189-8103-18c46ff3245d
# ╟─4b647cd7-28d7-4630-8d61-63c0ead5c639
# ╟─630c76c7-1a17-4fe2-8770-9fd3f58932b4
# ╠═5338451e-3c4b-4030-bbbb-42eaf4209a89
# ╠═6a71c140-fd0a-4506-92b0-fe802fdca18a
# ╠═36b15aa8-ca13-41e8-9bb6-a102b7e7b0b4
# ╠═b5abeda2-9f63-4aa9-943e-c7f4ad0bfef1
# ╠═9b2032ee-232f-4de5-96c0-3bcc9b78832c
# ╠═6970dab5-16bd-4898-b88d-723cb1b3d89e
# ╠═97b0763d-dcab-4afa-b660-52e18b3d523f
# ╠═80ea58f8-a77c-415e-92cb-86b836f0566e
# ╠═00b7bb48-c476-4d9b-940d-ca5dbed3477d
# ╠═8689cec2-f26a-438d-80d1-765cf07d90b2
# ╠═00bdd1d7-655d-4500-bfdd-499069e56490
# ╠═21b5a75a-3cb7-4d05-a8b2-639f3974febc
# ╠═33bcfdc7-5dd4-4984-9f23-858a14e8bc10
# ╠═8f56cf8b-084d-472f-b578-a8554ad5788e
# ╟─b473b17e-3bf5-4b6c-af24-fe57b5a7e7e9
# ╠═999ae1fd-5341-4f66-9db2-dec53fa0cd49
# ╠═e40b7e33-fae3-40a1-9a09-542483476c8f
# ╠═670ef5ac-aaae-4d2e-9d31-135417a8b198
# ╠═ea09c6d6-1485-4d86-9eed-8a3418199148
# ╟─b7e5d1a6-57ed-4d09-a039-a4bd12386367
# ╠═4520fc6e-7305-487e-924d-af22406e6d45
# ╠═a958fbba-c172-4858-80ab-9f98dbc42d6c
# ╠═5ba2c405-11be-4914-93be-379ed45eb63d
# ╠═6e56cc9c-9fa8-405c-bfdd-191a02bee4aa
# ╠═3a25bcad-7529-41c3-bfd0-90db0f192b38
# ╠═111bebcc-6fe2-440b-860c-9d18e092a617
# ╠═67f5db98-88d0-11ec-27ac-b57538a166f4
# ╠═267865de-1b5c-4579-861b-c6c46beb4739
# ╟─ab2bcfd5-3ba7-4388-8a3c-2cb95fba989a
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
