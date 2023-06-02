include("Test_Module_NKG.jl")

using .Test_Module_NKG

Sum1, Diff1 = Test_Module_NKG.Add_Sub_Numbers(10,10)

Diff = Test_Module_NKG.Test_Sub_NKG(10,10)

Sum = Test_Module_NKG.Test_Add_NKG(10,10)