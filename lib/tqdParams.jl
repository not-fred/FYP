using DataStructures, JSON

function initparams(params)

    loopKey      = Array{String}(undef,0)
    loopVal      = Array{Array}(undef,0)
    noloopParams = DataStructures.OrderedDict{String,Any}()

    for (key,val) ∈ params
        if length(val) > 1 && !isa(val, String) && key != "M"
            push!(loopKey, key)
            val = stepstoarray(val)
            push!(loopVal, val)
        else
            push!(noloopParams, key => val)
        end
    end

    return loopKey, loopVal, noloopParams

end

function paramfilename(p, dir::String="", ext::String=".json"; nameless=["U", "V"])
    dir = replace(dir, r"\\^" => "")
    name = dir != "" ? "$dir/" : ""
    for (k,v) ∈ p
        if !(k ∈ nameless)
            name *= "_$k("
            if length(v) > 1 && !isa(v, String)
                if "stop" ∈ keys(v)
                    step = "step" ∈ keys(v) ? deround(v["step"]) : 1
                    name *= "$(deround(v["start"])):$(step):$(deround(v["stop"]))"
                else
                    name *= "$(deround(v[1])):$(deround(v[2]-v[1])):$(deround(v[end]))"
                end
            else
                name *= "$(deround(v))"
            end
            name *= ")"
        end
    end
    return name * ext
end

function stepstoarray(val)
    if "stop" ∈ keys(val)
        step = "step" ∈ keys(val) ? val["step"] : 1
        return val["start"]:step:val["stop"]
    else
        return val
    end
end

function loadparams(val,loopKey,noloopParams)
    params = copy(noloopParams)
    alert = ""
    for i in 1:length(loopKey)
        push!(params, loopKey[i]=>val[i])
        alert != "" && (alert *= ", ")
        alert *= "$(loopKey[i]) = $(val[i])"
    end
    data = try
             JSON.parsefile(paramfilename(params, dir))
           catch err
               @warn err
               false
           end
    return data,alert
end


function deround(val)
    val = string(val)
    return replace(val, r"(\.\d*[1-9])0{3}0*[1-9](?:e-?\d+)?$" => s"\1")
end
