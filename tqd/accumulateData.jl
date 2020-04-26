using JSON, DataStructures

include("../lib/tqdParams.jl")

# Try to load the specified parameters
params = try
        JSON.parsefile(ARGS[1],dicttype=DataStructures.OrderedDict)
    catch err
        @warn err
        exit()
    end
dir = try pop!(params, "dir") catch err "" end
dir = replace(dir, r"\\^" => "")

loopKey, loopVal, noloopParams = initparams(params)
vals = Iterators.product(loopVal...)

improved = "improved" ∈ keys(noloopParams)
table = Array{Union{Float64,Int64,Missing}}(
            missing,
            reduce(*,length.(loopVal)),
            length(noloopParams)+length(loopKey)+3-(improved ? 1 : 0)
        )

noloopVals = copy(noloopParams)
_ = stepstoarray(pop!(noloopVals,"M"))
_ = pop!(noloopVals,"K")
improved && (_ = pop!(noloopVals,"improved"))

copy(noloopVals)

tableKeys = append!(copy(noloopVals.keys),loopKey)
append!(tableKeys,["r","σr","r̄","μr̄","σr̄"])
table[:,1:length(noloopVals.vals)] = repeat(noloopVals.vals',outer=size(table)[1])

rArr  = Array{Union{Float64,Missing}}(missing,length.(loopVal)...)
σrArr = copy(rArr)
r̄Arr  = copy(rArr)

vals = collect(vals)

for i ∈ 1:length(vals)
    val = vals[i]
    params = copy(noloopParams)
    alert = ""
    for i ∈ 1:length(loopKey)
        push!(params, loopKey[i]=>val[i])
        alert != "" && (alert *= ", ")
        alert *= "$(loopKey[i]) = $(val[i])"
    end
    println("Accumulating $alert")
    data = try
             JSON.parsefile(paramfilename(params, dir))
           catch err
               @warn err
               false
           end
    if data != false
        p = data["p"]
        σp = data["σp"]
        p̄  = data["p̄"]
        arrI = [findfirst(loopVal[i] .== val[i]) for i ∈ 1:length(loopKey)]
        rArr[ arrI... ] = (1-p[1])/2
        σrArr[ arrI... ] = σp[1]/2
        r̄Arr[ arrI... ] = (1-p̄[2][1])/2
        table[i,length(noloopVals.vals)+1:end] = [
            val...,
            rArr[arrI...],
            σrArr[arrI...],
            r̄Arr[arrI...],
            data["μr̄"],
            data["σr̄"]
        ]
    end
end


filename = replace(ARGS[1], r".*\/([^\/]*)\.json$" => s"\1")
if dir != ""
    filename = "$dir/$filename"
end

open("$filename.accumulated.csv","w") do f
    write(f,join(tableKeys,",")*"\n")
    for i ∈ 1:length(vals)
        write(f,join(table[i,:],",")*"\n")
    end
end
println("Accumulated data saved at $filename.accumulated.csv")

#=
For output in JSON instead
open("$filename.accumulated.json","w") do f
    JSON.print(f, DataStructures.OrderedDict([
        "key" => tableKeys,
        "val" => table
    ]))
end
=#
