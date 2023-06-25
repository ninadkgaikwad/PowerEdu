### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# ╔═╡ ecd816ef-75dd-45f5-b255-89d22324066e
using Pkg

# ╔═╡ 8899d627-5781-42dd-a8d1-02a6061465f5
Pkg.activate("sparseTechniques"); 

# ╔═╡ 3667b924-a7ab-4d5f-b08b-f5d9c94a2a25
using Base

# ╔═╡ e2787ca2-f3fa-41bc-bc73-3adc1de1652a
using DataFrames;

# ╔═╡ 7ed6c943-816c-422d-8128-125775f3f606
using LinearAlgebra;

# ╔═╡ 3e99b6c8-a1c9-43d5-9b54-06d0aaa49147
using PlutoUI;

# ╔═╡ a52fc0f5-1797-4fa9-8674-306a96da363b
using SparseArrays;

# ╔═╡ 29a33d0a-c5b9-4c85-9267-c8362de58b7f
include("src/SparseTechniques_Functions.jl");

# ╔═╡ 8ccc68a2-302c-4607-9622-318c523897cc
md"# Sparse Matrix Solution Techniques for Power Systems"

# ╔═╡ 37dfb22b-a604-45ef-8b78-6bcb192e5883
md" Julia implementations from Chapter 04: Sparse Matrix Solution Techniques from [Mariesa L. Crow's](https://www.linkedin.com/in/mariesa-crow-17a72895/) book [Computational Methods for Electric Power Systems](https://www.amazon.com/Computational-Methods-Electric-Systems-Engineering/dp/1032098228#customerReviews)"

# ╔═╡ 6fdad9d2-50de-4b55-bbeb-3bd3b431dff5
md"## Using Packages"

# ╔═╡ 0d452297-9850-45d2-93f3-759b2f0f5fb9
#include("src/Helper_Functions.jl");

# ╔═╡ 2d8e35ff-8880-4237-a25e-4c4d0c07bbaa
md" Test Case Compressed Format Matrix (Elements are in order, sorted by row, then column)"

# ╔═╡ 4e923321-2836-4518-b3e6-0c28d838c45a
N = 5;

# ╔═╡ 8618c265-de20-4385-8173-12c23faf7b20
nnz = 12;

# ╔═╡ d4a4cbf5-8c21-4793-99c6-ddd833b5d04e
md"### Defining the sparse matrix from the reference book."

# ╔═╡ 538bdf8c-6057-46fe-a802-475a816472e3
begin
	values = vec([-1, -2, 2, 8, 1, 3, -2, -3, 2, 1, 2, -4]);
	rows = vec([1 1 2 2 2 3 3 4 4 5 5 5]);
	cols = vec([1 3 1 2 4 3 5 2 3 1 2 5]);
	mat1 = DataFrame(Val = values, i = rows, j = cols)
end

# ╔═╡ 93e59698-3702-4447-bf75-819f606e1f47
md"### Use `sparmat` to build the $[N, nnz]$ vectors"

# ╔═╡ dac8ab02-4334-455e-8702-0db4e31a5f61
NVec, nnzVec = sparmat(mat1)

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


# ╔═╡ Cell order:
# ╟─8ccc68a2-302c-4607-9622-318c523897cc
# ╟─37dfb22b-a604-45ef-8b78-6bcb192e5883
# ╟─6fdad9d2-50de-4b55-bbeb-3bd3b431dff5
# ╟─ecd816ef-75dd-45f5-b255-89d22324066e
# ╟─8899d627-5781-42dd-a8d1-02a6061465f5
# ╟─3667b924-a7ab-4d5f-b08b-f5d9c94a2a25
# ╟─e2787ca2-f3fa-41bc-bc73-3adc1de1652a
# ╟─7ed6c943-816c-422d-8128-125775f3f606
# ╟─3e99b6c8-a1c9-43d5-9b54-06d0aaa49147
# ╟─a52fc0f5-1797-4fa9-8674-306a96da363b
# ╟─29a33d0a-c5b9-4c85-9267-c8362de58b7f
# ╟─0d452297-9850-45d2-93f3-759b2f0f5fb9
# ╟─2d8e35ff-8880-4237-a25e-4c4d0c07bbaa
# ╟─4e923321-2836-4518-b3e6-0c28d838c45a
# ╟─8618c265-de20-4385-8173-12c23faf7b20
# ╟─d4a4cbf5-8c21-4793-99c6-ddd833b5d04e
# ╟─538bdf8c-6057-46fe-a802-475a816472e3
# ╟─93e59698-3702-4447-bf75-819f606e1f47
# ╟─dac8ab02-4334-455e-8702-0db4e31a5f61
# ╟─3d75f29b-d927-4459-8497-3e8d13c4799a
# ╟─49924d2f-9dcd-404d-96a2-61aa86afe6f0
# ╟─ed1b3b9c-a3ba-4e08-ae69-43215d21bf8e
# ╟─8c1a2445-a323-458b-873f-81e49ae20aa0
# ╟─caf47f37-63d0-47af-854c-9ce4c708b210
# ╟─3646b871-c83a-410a-af4b-b1ff4dd86540
# ╟─7fabe75c-5c86-4315-9b5c-0c6856ed4a8d
# ╟─8c064aa3-ae32-4736-bb54-5e1b9acf891b
# ╟─d5f1ab89-ebc0-42fa-bbf4-fe9d5ac9c960
# ╟─bdeae12a-e105-46b8-b91b-8bd06da8ffd6
# ╟─b360eaea-452d-4f85-b8fe-832d51d96ee5
# ╟─c062d472-9099-424f-a9d7-d62aa3bcc555
# ╟─7c994315-27c4-4ce0-9e65-bb47f1f459dc
# ╟─573cb2be-8e50-42e7-aa49-49cdb89e59c1
# ╟─00ae0a0d-2b17-4990-9e5b-5a1eb4ed7e60
# ╟─2dd4646b-9a80-4598-89c8-f34826a6ebff
# ╟─deeaeafe-7701-4bef-a58d-340aef92feb1
# ╟─fd1ebdc3-412a-4c40-9ef8-f5bfadac66f4
# ╟─9989fcde-17b3-4ef8-b6d9-89de778b8bbe
# ╟─f1b454f3-ca44-44ea-b15a-d2c7ce68a8f8
# ╟─7c0eeee3-52a0-45d3-a9b4-bb0b497ae310
# ╟─fcc50a39-3f19-4b41-8eb1-bed4144fcc57
# ╟─c2774e13-82b1-439f-bd1f-06be228749be
# ╟─74031088-02c9-4ade-a1e5-36bee10951c7
# ╟─3fbb08a2-03b4-4c5c-b29b-9d41a1a203d6
# ╟─4f12df3d-1330-4adc-a74e-79a184f51615
# ╟─8543da9e-cb7a-409c-b924-61774c30a341
# ╟─54238595-c782-45b9-bb2c-2e60bf4bde9a
# ╟─0b2b3062-2689-49bf-9170-556287912f9c
# ╟─f6ff7694-f21f-498b-821b-bae38529f045
# ╟─fc18d362-9705-44dd-b321-6a8218b709bc
# ╟─b46244ef-7c81-4185-b975-1fe8bac53c99
