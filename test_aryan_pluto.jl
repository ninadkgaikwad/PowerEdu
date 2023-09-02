### A Pluto.jl notebook ###
# v0.19.27

#> [frontmatter]
#> title = "Sparse Data Structures and Linear System Techniques for Power Systems"
#> date = "2023-07-29"
#> tags = ["julia", " pluto", "power systems", " sparse matrix", " sparse techniques", "sparse data structures", "power system analysis", "sparse power flow", "sparse ybus", "sparse jacobian", "sparse LU factorization", "sparse LU decomposition"]
#> description = "Julia implementation of Sparse Data Structures and Algorithms for Solving Power Flow"

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ f5f02417-00b9-4c2d-ae80-ea31b5abda5e
using Pkg;

# ╔═╡ 26ee4f30-e1f7-45fa-8c72-50b811e5c8bc
Pkg.activate(".")

# ╔═╡ 8899d627-5781-42dd-a8d1-02a6061465f5
# ╠═╡ show_logs = false
Pkg.activate("sparseTechniques"); 

# ╔═╡ 01f8bcc0-633f-46df-9bee-ed3961da4b2a
using Revise

# ╔═╡ 896971e8-70f3-411b-90b7-276a3149a286
using PlutoUI;

# ╔═╡ 59abea63-093a-44e2-a63d-673fb8c5c294
using LaTeXStrings

# ╔═╡ e2787ca2-f3fa-41bc-bc73-3adc1de1652a
using DataFrames;

# ╔═╡ e8281b42-2efa-4f5d-9c79-1459efa303ab
using CSV;

# ╔═╡ 7ed6c943-816c-422d-8128-125775f3f606
using LinearAlgebra;

# ╔═╡ 672d5107-496d-407c-9908-da9bf9149f81
# ╠═╡ show_logs = false
using Plots;

# ╔═╡ a52fc0f5-1797-4fa9-8674-306a96da363b
using SparseArrays;

# ╔═╡ 75d68fa7-7a13-4356-a619-30cd7f662dc0
using Printf;

# ╔═╡ 6fe55bc7-11e4-4ef4-8c3e-8ca822b72535
include("src/Helper_Functions.jl");

# ╔═╡ 933a3383-5023-4df2-9a46-cdfdcb7087af
include("src/IEEE_CDF_Parser.jl");

# ╔═╡ 533cca13-03ef-4ad6-9d0c-1ed1f9176483
# ╠═╡ show_logs = false
include("src/SparseTechniques_Functions.jl");

# ╔═╡ 8ccc68a2-302c-4607-9622-318c523897cc
md"# Sparse Matrix Solution Techniques for Power Systems"

# ╔═╡ fade459d-91b9-440b-907b-766a4baa4a1a
Pkg.status();

# ╔═╡ 37dfb22b-a604-45ef-8b78-6bcb192e5883
md" Julia implementations from Chapter 04: Sparse Matrix Solution Techniques from [Mariesa L. Crow's](https://www.linkedin.com/in/mariesa-crow-17a72895/) book [Computational Methods for Electric Power Systems](https://www.amazon.com/Computational-Methods-Electric-Systems-Engineering/dp/1032098228#customerReviews)"

# ╔═╡ 5bf19ff7-4547-4718-bbc8-db563f2abc12
md"## System Selection"

# ╔═╡ e42b3750-5a62-45bf-88e3-c72c8478a1c0
md"""
Please choose your desired Power System from among the supported systems using this drop down menu: $(@bind systemName Select(["IEEE_14" => "IEEE 14", "IEEE_30" => "IEEE 30", "IEEE_57" => "IEEE 57", "IEEE_118" => "IEEE 118", "KOTH_3" => "Kothari 3 Bus System (Currently NOT Supported)", "IEEE_300" => "IEEE 300 (Currently NOT Supported.)"]))
"""

# ╔═╡ d99a50b4-4670-4559-b9b6-1e281f46e51a
md"### Other User Controls"

# ╔═╡ b3bc7313-0a3c-4996-a71c-602e240bc61e
md"""Do you want to see the Bus Data and Branch Data for your chosen system? $(@bind displaySystemData CheckBox(default=false) )

Do you want to see the sparse $Y_{Bus}$? $(@bind displayYBus CheckBox(default=true))

Do you want to see the sparse Jacobian? $(@bind displayJacobian CheckBox(default=true))

How to highlight the sparsity of data structures?
$(@bind returnValue Select(["print only" => "Only print out the % Sparsity Value.", "return value and print" => "Return % Sparsity Value and print it.", "return value only" => "Only return % Sparsity Value.", "nothing" => "Do NOT return or print."]))
"""

# ╔═╡ 6fdad9d2-50de-4b55-bbeb-3bd3b431dff5
md"## Using Packages and Modules"

# ╔═╡ 938ab111-da12-42ae-8da0-31aae092fae8
md" (Not Displayed)"

# ╔═╡ ecd816ef-75dd-45f5-b255-89d22324066e
# ╠═╡ disabled = true
#=╠═╡
using Pkg
  ╠═╡ =#

# ╔═╡ da3cfa14-9b0a-4b27-b680-de8f76ae58f1
md"### Housekeeping"

# ╔═╡ 071e9b51-c4e4-464f-81c1-5c4aa34511a3
md" (Not displayed)"

# ╔═╡ 3d75f29b-d927-4459-8497-3e8d13c4799a
folderInput = "data/";

# ╔═╡ 49924d2f-9dcd-404d-96a2-61aa86afe6f0
folder_processedData = "processedData/";

# ╔═╡ 8c1a2445-a323-458b-873f-81e49ae20aa0
# ╠═╡ show_logs = false
createFolderIfNotExisting(systemName, folder_processedData);

# ╔═╡ 29d49ed2-2af4-49be-a418-fd3ec06924cd
fileType_CDFFile = ".txt";

# ╔═╡ b39cb657-ff03-4ede-8b63-91780928aea3
filename_CDFFile = folderInput*systemName*"/"*systemName*"_Data"*fileType_CDFFile;

# ╔═╡ 6b1a9017-c902-429a-bf6c-a173c3ad190b
md"### Invoke the CDF Parser"

# ╔═╡ 3a394c37-70f7-476f-a830-d4287835fe59
md"to read the IEEE Common Data Format File for the chosen system"

# ╔═╡ 4017c4f8-9caa-4332-8c8d-994d696d7f5f
md" and converting relevant data into $pu$ values"

# ╔═╡ 54c04bc9-78d3-4b85-ab7b-abe2a6206596
CDF_DF_List = CDF_Parser(filename_CDFFile);

# ╔═╡ 1647fb65-2aec-4ff3-b695-e8ee2d190c25
CDF_DF_List_pu = CDF_pu_Converter(CDF_DF_List);

# ╔═╡ c330d28d-4cd1-4608-b540-cd7f51d0b376
busData_pu = CDF_DF_List_pu[2];

# ╔═╡ 093431aa-f04e-4e47-90c9-b1175e0aab8b
branchData_pu = CDF_DF_List_pu[3];

# ╔═╡ 4d9dd2f8-8776-4925-927f-825e8250324b
if displaySystemData
	# display(busData)
	(busData_pu = busData_pu, branchData_pu = branchData_pu)
	# display(branchData)
end

# ╔═╡ 0f6544d3-186f-4093-a5cf-e61f8bb36c94
md" ## Construct Sparse $Y_{Bus}$"

# ╔═╡ 3646b871-c83a-410a-af4b-b1ff4dd86540
sparYBus = constructSparseYBus(CDF_DF_List_pu);

# ╔═╡ bdaebe7a-ac92-434a-a3e6-4ff8cdba68fa
sparsityYBus = computeSparsity(sparYBus, returnValue=returnValue);

# ╔═╡ 459956ea-be4b-4d03-91c8-21ff1c77f973
@info sparsityYBus

# ╔═╡ e193276a-fa6a-481a-9b47-7d68882e52f2
begin

	YBusSpar2FullMags = abs.(spar2Full(sparYBus))
	
	markerSize1 = 4 ;
	
	pYBus = spy(YBusSpar2FullMags, 
		marker=:square, 
		markersize=markerSize1, 
		title=L"Sparsity Pattern of $Y_{Bus}$", 
		titlefont=("sans-serif", 18),
		gap="thiredy",
		color=:green);
	
	if returnValue == "return value only" || returnValue == "return value and print"
		annotate!(pYBus, ( (0.75, 0.1),
		text("Sparsity = $(@sprintf("%.2f", sparsityYBus/100))", 
		   	family="serif", 
            12, 
			:right, 
			:green) ))
	end
	
	
	plot(pYBus, layout=(1, 1))

end

# ╔═╡ ee890880-88c9-480a-bc95-ab0bb268efcd
if displayYBus
	# display(sparYBus.NVec)
	# display(sparYBus.MVec)
	# display(sparYBus.nnzVec)
	sparYBus
end

# ╔═╡ 311abf2d-5249-4e80-a97b-b809a79d7107
PowSysData = initializeVectors_pu(CDF_DF_List_pu);

# ╔═╡ 4727e699-d15e-48e6-b112-5526af6b4d8e
begin
	PSpecified = PowSysData.PSpecified;
	QSpecified = PowSysData.QSpecified;
	V = PowSysData.V;
	delta = PowSysData.delta;
end;

# ╔═╡ 4021a781-b33c-43a3-8893-3620306f7635
md" ### Compute Mismatches"

# ╔═╡ cc18a9dc-c35c-481e-a030-3fe919834b38
md" (Not Displayed)"

# ╔═╡ 8fda426a-5b55-4daa-8f0b-404290c14520
deltaP, deltaQ = computeMismatchesViaSparseYBus(PSpecified, QSpecified, V, delta, sparYBus);

# ╔═╡ 9aeffbaf-3223-45da-9d9f-6f7e6f9745e2
begin
	P = PSpecified - deltaP
	Q = QSpecified - deltaQ
end;

# ╔═╡ 4beab547-1b67-4de9-b667-f959fcd5177d
md" ## Construct Sparse Jacobian"

# ╔═╡ 77a35fcf-b866-48b1-8dc3-c303c799e9f5
sparJ = constructSparseJacobian(CDF_DF_List_pu, P, Q, V, delta, sparYBus);

# ╔═╡ 8eb2435e-eaa4-4d4a-8da0-f125c3125237
sparsityJ = computeSparsity(sparJ, returnValue=returnValue);

# ╔═╡ f42085e4-0615-4dca-8c2a-2e2eeaf1dcd7
if displayJacobian
	# display(sparJ)
	sparJ
end

# ╔═╡ 7fabe75c-5c86-4315-9b5c-0c6856ed4a8d
md"### Convert Sparse Jacobian into Sparse LU Factors"

# ╔═╡ 8c064aa3-ae32-4736-bb54-5e1b9acf891b
qluJ = sparLU(sparJ);

# ╔═╡ faeaf70f-a09d-4dae-8931-8f71035e2207
QJ = qluJ.Q;

# ╔═╡ 160ef94f-2ce8-4466-a080-99a046e4183f
sparsityQ = computeSparsity(QJ, returnValue=returnValue);

# ╔═╡ e5b94da1-be60-4d49-9a36-a0099c37fca9
md" #### Vizualising the sparse matrices"

# ╔═╡ 8a301f36-d59e-4cd1-9301-7257e3258b92
JSpar2Full = spar2Full(sparJ);

# ╔═╡ 1492fa9f-2ed8-4c64-a507-b23155df5bbd
QSpar2Full = spar2Full(QJ);

# ╔═╡ bdeae12a-e105-46b8-b91b-8bd06da8ffd6
begin

	markerSize = 1 ;
	
	pJ = spy(JSpar2Full, 
		marker=:square, 
		markersize=markerSize, 
		title=L"Sparsity Pattern of $J$", 
		color=:green);

	pQ = spy(QSpar2Full, 
	marker=:square, 
	markersize=markerSize, 
	title=L"Sparsity Pattern of $LU$ Factors", 
	color=:green);

	if returnValue == "return value only" || returnValue == "return value and print"

		annotate!(pJ, ( (0.85, -0.25) ,
			text("Sparsity = $(@sprintf("%.2f", sparsityJ/100))", 
				family="serif", 
				12, 
				:right, 
				:green) ) )

	

		annotate!(pQ, ( (0.85, -0.15),
			text("Fill-ins = $(qluJ.α)", 
				family="serif", 
				12, 
				:right, 
				:red)))
	
		annotate!(pQ, ( (0.85, -0.25),
				text("Sparsity = $(@sprintf("%.2f", sparsityQ/100))", 
					family="serif", 
					12, 
					:right, 
					:green) ))

	end

	plot(pJ, pQ, layout=(1, 2))

end

# ╔═╡ Cell order:
# ╟─8ccc68a2-302c-4607-9622-318c523897cc
# ╠═f5f02417-00b9-4c2d-ae80-ea31b5abda5e
# ╠═26ee4f30-e1f7-45fa-8c72-50b811e5c8bc
# ╠═fade459d-91b9-440b-907b-766a4baa4a1a
# ╠═896971e8-70f3-411b-90b7-276a3149a286
# ╟─37dfb22b-a604-45ef-8b78-6bcb192e5883
# ╟─5bf19ff7-4547-4718-bbc8-db563f2abc12
# ╟─e42b3750-5a62-45bf-88e3-c72c8478a1c0
# ╟─d99a50b4-4670-4559-b9b6-1e281f46e51a
# ╟─b3bc7313-0a3c-4996-a71c-602e240bc61e
# ╟─6fdad9d2-50de-4b55-bbeb-3bd3b431dff5
# ╟─938ab111-da12-42ae-8da0-31aae092fae8
# ╟─ecd816ef-75dd-45f5-b255-89d22324066e
# ╟─01f8bcc0-633f-46df-9bee-ed3961da4b2a
# ╟─59abea63-093a-44e2-a63d-673fb8c5c294
# ╟─8899d627-5781-42dd-a8d1-02a6061465f5
# ╟─e2787ca2-f3fa-41bc-bc73-3adc1de1652a
# ╟─e8281b42-2efa-4f5d-9c79-1459efa303ab
# ╟─7ed6c943-816c-422d-8128-125775f3f606
# ╟─672d5107-496d-407c-9908-da9bf9149f81
# ╟─a52fc0f5-1797-4fa9-8674-306a96da363b
# ╟─75d68fa7-7a13-4356-a619-30cd7f662dc0
# ╠═6fe55bc7-11e4-4ef4-8c3e-8ca822b72535
# ╟─933a3383-5023-4df2-9a46-cdfdcb7087af
# ╟─533cca13-03ef-4ad6-9d0c-1ed1f9176483
# ╟─da3cfa14-9b0a-4b27-b680-de8f76ae58f1
# ╟─071e9b51-c4e4-464f-81c1-5c4aa34511a3
# ╟─3d75f29b-d927-4459-8497-3e8d13c4799a
# ╟─49924d2f-9dcd-404d-96a2-61aa86afe6f0
# ╟─8c1a2445-a323-458b-873f-81e49ae20aa0
# ╟─29d49ed2-2af4-49be-a418-fd3ec06924cd
# ╟─b39cb657-ff03-4ede-8b63-91780928aea3
# ╟─6b1a9017-c902-429a-bf6c-a173c3ad190b
# ╟─3a394c37-70f7-476f-a830-d4287835fe59
# ╟─4017c4f8-9caa-4332-8c8d-994d696d7f5f
# ╠═54c04bc9-78d3-4b85-ab7b-abe2a6206596
# ╠═1647fb65-2aec-4ff3-b695-e8ee2d190c25
# ╠═c330d28d-4cd1-4608-b540-cd7f51d0b376
# ╠═093431aa-f04e-4e47-90c9-b1175e0aab8b
# ╟─4d9dd2f8-8776-4925-927f-825e8250324b
# ╟─0f6544d3-186f-4093-a5cf-e61f8bb36c94
# ╠═3646b871-c83a-410a-af4b-b1ff4dd86540
# ╠═bdaebe7a-ac92-434a-a3e6-4ff8cdba68fa
# ╠═459956ea-be4b-4d03-91c8-21ff1c77f973
# ╠═e193276a-fa6a-481a-9b47-7d68882e52f2
# ╟─ee890880-88c9-480a-bc95-ab0bb268efcd
# ╟─311abf2d-5249-4e80-a97b-b809a79d7107
# ╟─4727e699-d15e-48e6-b112-5526af6b4d8e
# ╟─4021a781-b33c-43a3-8893-3620306f7635
# ╟─cc18a9dc-c35c-481e-a030-3fe919834b38
# ╟─8fda426a-5b55-4daa-8f0b-404290c14520
# ╟─9aeffbaf-3223-45da-9d9f-6f7e6f9745e2
# ╟─4beab547-1b67-4de9-b667-f959fcd5177d
# ╠═77a35fcf-b866-48b1-8dc3-c303c799e9f5
# ╠═8eb2435e-eaa4-4d4a-8da0-f125c3125237
# ╟─f42085e4-0615-4dca-8c2a-2e2eeaf1dcd7
# ╟─7fabe75c-5c86-4315-9b5c-0c6856ed4a8d
# ╠═8c064aa3-ae32-4736-bb54-5e1b9acf891b
# ╠═faeaf70f-a09d-4dae-8931-8f71035e2207
# ╠═160ef94f-2ce8-4466-a080-99a046e4183f
# ╟─e5b94da1-be60-4d49-9a36-a0099c37fca9
# ╠═8a301f36-d59e-4cd1-9301-7257e3258b92
# ╠═1492fa9f-2ed8-4c64-a507-b23155df5bbd
# ╠═bdeae12a-e105-46b8-b91b-8bd06da8ffd6
