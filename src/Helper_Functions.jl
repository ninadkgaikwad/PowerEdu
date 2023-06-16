# Helper_Functions.jl

function myprintln(verbose::Bool, args...)
    if verbose
        println(args...)
    end
end
