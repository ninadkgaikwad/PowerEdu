### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# ╔═╡ e2787ca2-f3fa-41bc-bc73-3adc1de1652a
using DataFrames

# ╔═╡ c2d09de5-f6a6-4de8-8cff-b95932e38923
using LinearAlgebra

# ╔═╡ 3e99b6c8-a1c9-43d5-9b54-06d0aaa49147
using PlutoUI

# ╔═╡ 29a33d0a-c5b9-4c85-9267-c8362de58b7f
include("src/SparseTechniques_Functions.jl");

# ╔═╡ 8ccc68a2-302c-4607-9622-318c523897cc
md"# Sparse Matrix Solution Techniques for Power Systems"

# ╔═╡ 37dfb22b-a604-45ef-8b78-6bcb192e5883
md" Julia implementations from Chapter 04: Sparse Matrix Solution Techniques from [Mariesa L. Crow's](https://www.linkedin.com/in/mariesa-crow-17a72895/) book [Computational Methods for Electric Power Systems](https://www.amazon.com/Computational-Methods-Electric-Systems-Engineering/dp/1032098228#customerReviews)"

# ╔═╡ 6fdad9d2-50de-4b55-bbeb-3bd3b431dff5
md"## Using Packages"

# ╔═╡ 2d8e35ff-8880-4237-a25e-4c4d0c07bbaa
md" Test Case Compressed Format Matrix (Elements are in order, sorted by row, then column)"

# ╔═╡ 4e923321-2836-4518-b3e6-0c28d838c45a
N = 5;

# ╔═╡ a1ca6ce6-57ee-4994-aa12-4a437d0f064e
N′ = N+1;

# ╔═╡ 8618c265-de20-4385-8173-12c23faf7b20
nnz = 12;

# ╔═╡ b4f4f208-6bf4-4194-a3a7-f7712d5cb38a
nnz′ = N^2 + 1;

# ╔═╡ 538bdf8c-6057-46fe-a802-475a816472e3
begin
	values = vec([-1, -2, 2, 8, 1, 3, -2, -3, 2, 1, 2, -4]);
	rows = vec([1 1 2 2 2 3 3 4 4 5 5 5]);
	cols = vec([1 3 1 2 4 3 5 2 3 1 2 5]);
	mat1 = DataFrame(val = values, i = rows, j = cols)
end

# ╔═╡ cabfd17c-3c44-4bfe-8212-26452b09cdbc
compressed2Full(mat1)

# ╔═╡ bcb71067-a9fe-49e6-b6ab-88895b93fa33
md"### Instantiate the $N$ vector"

# ╔═╡ 249d66ac-772f-44fb-a2f5-a00e95eefde8
firs, fics = [nnz′ * ones(Int64, N) for _ = 1:2];

# ╔═╡ daef88a8-db44-4a7d-be70-a25920b18726
NVec = DataFrame(FIR = firs, FIC = fics)

# ╔═╡ 42db1323-6675-4119-a0c7-c477e7c0e676
indices = vec(1:1:nnz);

# ╔═╡ b25d3b4c-9b74-4e23-af1c-4c189624d375
nrows, ncols = [N′ * ones(Int64, nnz) for _ in 1:2];

# ╔═╡ a047a9a6-7c8b-460f-af9d-6ab534f54e0f
nirs, nics = [nnz′ * ones(Int64, nnz) for _ in 1:2];

# ╔═╡ 93e59698-3702-4447-bf75-819f606e1f47
md"### Instantiate the $nnz$ vector"

# ╔═╡ e6692a0a-e39d-431f-8b82-abd1f5568bfb
nnzVec = DataFrame(idx = indices, VALUE = values, NROW = nrows, NCOL = ncols, NIR = nirs, NIC = nics)

# ╔═╡ dac8ab02-4334-455e-8702-0db4e31a5f61


# ╔═╡ 3d75f29b-d927-4459-8497-3e8d13c4799a


# ╔═╡ 49924d2f-9dcd-404d-96a2-61aa86afe6f0


# ╔═╡ ed1b3b9c-a3ba-4e08-ae69-43215d21bf8e


# ╔═╡ 8c1a2445-a323-458b-873f-81e49ae20aa0


# ╔═╡ caf47f37-63d0-47af-854c-9ce4c708b210


# ╔═╡ 3646b871-c83a-410a-af4b-b1ff4dd86540


# ╔═╡ 7fabe75c-5c86-4315-9b5c-0c6856ed4a8d


# ╔═╡ 8c064aa3-ae32-4736-bb54-5e1b9acf891b


# ╔═╡ d5f1ab89-ebc0-42fa-bbf4-fe9d5ac9c960


# ╔═╡ bdeae12a-e105-46b8-b91b-8bd06da8ffd6


# ╔═╡ b360eaea-452d-4f85-b8fe-832d51d96ee5


# ╔═╡ c062d472-9099-424f-a9d7-d62aa3bcc555


# ╔═╡ 7c994315-27c4-4ce0-9e65-bb47f1f459dc


# ╔═╡ 573cb2be-8e50-42e7-aa49-49cdb89e59c1


# ╔═╡ 00ae0a0d-2b17-4990-9e5b-5a1eb4ed7e60


# ╔═╡ 2dd4646b-9a80-4598-89c8-f34826a6ebff


# ╔═╡ deeaeafe-7701-4bef-a58d-340aef92feb1


# ╔═╡ fd1ebdc3-412a-4c40-9ef8-f5bfadac66f4


# ╔═╡ 9989fcde-17b3-4ef8-b6d9-89de778b8bbe


# ╔═╡ f1b454f3-ca44-44ea-b15a-d2c7ce68a8f8


# ╔═╡ 7c0eeee3-52a0-45d3-a9b4-bb0b497ae310


# ╔═╡ fcc50a39-3f19-4b41-8eb1-bed4144fcc57


# ╔═╡ c2774e13-82b1-439f-bd1f-06be228749be


# ╔═╡ 74031088-02c9-4ade-a1e5-36bee10951c7


# ╔═╡ 3fbb08a2-03b4-4c5c-b29b-9d41a1a203d6


# ╔═╡ 4f12df3d-1330-4adc-a74e-79a184f51615


# ╔═╡ 8543da9e-cb7a-409c-b924-61774c30a341


# ╔═╡ 54238595-c782-45b9-bb2c-2e60bf4bde9a


# ╔═╡ 0b2b3062-2689-49bf-9170-556287912f9c


# ╔═╡ f6ff7694-f21f-498b-821b-bae38529f045


# ╔═╡ fc18d362-9705-44dd-b321-6a8218b709bc


# ╔═╡ b46244ef-7c81-4185-b975-1fe8bac53c99


# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
DataFrames = "~1.5.0"
PlutoUI = "~0.7.51"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.9.0"
manifest_format = "2.0"
project_hash = "89604439512714c5bba0c40b17bfcaa99d67b9fd"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "8eaf9f1b4921132a4cff3f36a1d9ba923b14a481"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.4"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.Compat]]
deps = ["UUIDs"]
git-tree-sha1 = "7a60c856b9fa189eb34f5f8a6f6b5529b7942957"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.6.1"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.2+0"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.DataAPI]]
git-tree-sha1 = "8da84edb865b0b5b0100c0666a9bc9a0b71c553c"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.15.0"

[[deps.DataFrames]]
deps = ["Compat", "DataAPI", "Future", "InlineStrings", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrettyTables", "Printf", "REPL", "Random", "Reexport", "SentinelArrays", "SnoopPrecompile", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "aa51303df86f8626a962fccb878430cdb0a97eee"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.5.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "d1fff3a548102f48987a52a2e0d114fa97d730f0"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.13"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

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

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "c47c5fa4c5308f27ccaac35504858d8914e102f9"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.4"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "d75853a0bdbfb1ac815478bacd89cd27b550ace6"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.3"

[[deps.InlineStrings]]
deps = ["Parsers"]
git-tree-sha1 = "9cc2baf75c6d09f9da536ddf58eb2f29dedaf461"
uuid = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48"
version = "1.4.0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.InvertedIndices]]
git-tree-sha1 = "0dc7b50b8d436461be01300fd8cd45aa0274b038"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.3.0"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.3"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "7.84.0+0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.10.2+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+0"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "f66bdc5de519e8f8ae43bdc598782d35a25b1272"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.1.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.10.11"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.21+4"

[[deps.OrderedCollections]]
git-tree-sha1 = "d321bf2de576bf25ec4d3e4360faca399afca282"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.0"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "a5aef8d4a6e8d81f171b2bd4be5265b01384c74c"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.5.10"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.9.0"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "b478a748be27bd2f2c73a7690da219d0844db305"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.51"

[[deps.PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "a6062fe4063cdafe78f4a0a81cfffb89721b30e7"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.4.2"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "259e206946c293698122f63e2b513a7c99a244e8"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.1.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "7eb1686b4f04b82f96ed7a4ea5890a4f0c7a09f1"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.0"

[[deps.PrettyTables]]
deps = ["Crayons", "Formatting", "LaTeXStrings", "Markdown", "Reexport", "StringManipulation", "Tables"]
git-tree-sha1 = "213579618ec1f42dea7dd637a42785a608b1ea9c"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "2.2.4"

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
version = "0.7.0"

[[deps.SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "04bdff0b09c65ff3e06a05e3eb7b120223da3d39"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.4.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.SnoopPrecompile]]
deps = ["Preferences"]
git-tree-sha1 = "e760a70afdcd461cf01a575947738d359234665c"
uuid = "66db9d55-30c0-4569-8b51-7e840670fc0c"
version = "1.0.3"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "a4ada03f999bd01b3a25dcaa30b2d929fe537e00"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.1.0"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.9.0"

[[deps.StringManipulation]]
git-tree-sha1 = "46da2434b41f41ac3594ee9816ce5541c6096123"
uuid = "892a3eda-7b42-436c-8928-eab12a02cf0e"
version = "0.3.0"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "Pkg", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "5.10.1+6"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits", "Test"]
git-tree-sha1 = "1544b926975372da01227b382066ab70e574a3ec"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.10.1"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.Tricks]]
git-tree-sha1 = "aadb748be58b492045b4f56166b5188aa63ce549"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.7"

[[deps.URIs]]
git-tree-sha1 = "074f993b0ca030848b897beff716d93aca60f06a"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.4.2"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.7.0+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.48.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+0"
"""

# ╔═╡ Cell order:
# ╟─8ccc68a2-302c-4607-9622-318c523897cc
# ╟─37dfb22b-a604-45ef-8b78-6bcb192e5883
# ╟─6fdad9d2-50de-4b55-bbeb-3bd3b431dff5
# ╠═e2787ca2-f3fa-41bc-bc73-3adc1de1652a
# ╠═c2d09de5-f6a6-4de8-8cff-b95932e38923
# ╠═3e99b6c8-a1c9-43d5-9b54-06d0aaa49147
# ╠═29a33d0a-c5b9-4c85-9267-c8362de58b7f
# ╟─2d8e35ff-8880-4237-a25e-4c4d0c07bbaa
# ╠═4e923321-2836-4518-b3e6-0c28d838c45a
# ╠═a1ca6ce6-57ee-4994-aa12-4a437d0f064e
# ╠═8618c265-de20-4385-8173-12c23faf7b20
# ╠═b4f4f208-6bf4-4194-a3a7-f7712d5cb38a
# ╠═538bdf8c-6057-46fe-a802-475a816472e3
# ╠═cabfd17c-3c44-4bfe-8212-26452b09cdbc
# ╠═bcb71067-a9fe-49e6-b6ab-88895b93fa33
# ╠═249d66ac-772f-44fb-a2f5-a00e95eefde8
# ╠═daef88a8-db44-4a7d-be70-a25920b18726
# ╠═42db1323-6675-4119-a0c7-c477e7c0e676
# ╠═b25d3b4c-9b74-4e23-af1c-4c189624d375
# ╠═a047a9a6-7c8b-460f-af9d-6ab534f54e0f
# ╠═93e59698-3702-4447-bf75-819f606e1f47
# ╠═e6692a0a-e39d-431f-8b82-abd1f5568bfb
# ╠═dac8ab02-4334-455e-8702-0db4e31a5f61
# ╠═3d75f29b-d927-4459-8497-3e8d13c4799a
# ╠═49924d2f-9dcd-404d-96a2-61aa86afe6f0
# ╠═ed1b3b9c-a3ba-4e08-ae69-43215d21bf8e
# ╠═8c1a2445-a323-458b-873f-81e49ae20aa0
# ╠═caf47f37-63d0-47af-854c-9ce4c708b210
# ╠═3646b871-c83a-410a-af4b-b1ff4dd86540
# ╠═7fabe75c-5c86-4315-9b5c-0c6856ed4a8d
# ╠═8c064aa3-ae32-4736-bb54-5e1b9acf891b
# ╠═d5f1ab89-ebc0-42fa-bbf4-fe9d5ac9c960
# ╠═bdeae12a-e105-46b8-b91b-8bd06da8ffd6
# ╠═b360eaea-452d-4f85-b8fe-832d51d96ee5
# ╠═c062d472-9099-424f-a9d7-d62aa3bcc555
# ╠═7c994315-27c4-4ce0-9e65-bb47f1f459dc
# ╠═573cb2be-8e50-42e7-aa49-49cdb89e59c1
# ╠═00ae0a0d-2b17-4990-9e5b-5a1eb4ed7e60
# ╠═2dd4646b-9a80-4598-89c8-f34826a6ebff
# ╠═deeaeafe-7701-4bef-a58d-340aef92feb1
# ╠═fd1ebdc3-412a-4c40-9ef8-f5bfadac66f4
# ╠═9989fcde-17b3-4ef8-b6d9-89de778b8bbe
# ╠═f1b454f3-ca44-44ea-b15a-d2c7ce68a8f8
# ╠═7c0eeee3-52a0-45d3-a9b4-bb0b497ae310
# ╠═fcc50a39-3f19-4b41-8eb1-bed4144fcc57
# ╠═c2774e13-82b1-439f-bd1f-06be228749be
# ╠═74031088-02c9-4ade-a1e5-36bee10951c7
# ╠═3fbb08a2-03b4-4c5c-b29b-9d41a1a203d6
# ╠═4f12df3d-1330-4adc-a74e-79a184f51615
# ╠═8543da9e-cb7a-409c-b924-61774c30a341
# ╠═54238595-c782-45b9-bb2c-2e60bf4bde9a
# ╠═0b2b3062-2689-49bf-9170-556287912f9c
# ╠═f6ff7694-f21f-498b-821b-bae38529f045
# ╠═fc18d362-9705-44dd-b321-6a8218b709bc
# ╠═b46244ef-7c81-4185-b975-1fe8bac53c99
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
