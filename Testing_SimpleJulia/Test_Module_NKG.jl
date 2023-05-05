module Test_Module_NKG  

    export Test_Add_NKG , Test_Sub_NKG , Add_Sub_Numbers

    #export Test_Add_NKG , Test_Sub_NKG

    # Functions
    function Add_Sub_Numbers(Num1,Num2)

        Sum = Test_Add_NKG(Num1,Num2)

        Diff = Test_Sub_NKG(Num1,Num2)

        return Sum, Diff

    end


    # Including component files
    include("Test_FunctionFile_1.jl")

    include("Test_FunctionFile_2.jl")

end