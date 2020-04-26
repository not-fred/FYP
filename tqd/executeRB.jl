using LsqFit, JSON, Dates, DataStructures, GLM
import Statistics.mean, Statistics.std, DataFrames.DataFrame

include("tqdRB.jl")
include("../lib/tqdParams.jl")
include("../lib/RB.jl")

# Detect use of MPI, do initialisation
# process if it exists
hasMPI = try using _MPI; true catch; false end

if hasMPI
    MPI.Init()
    comm = MPI.COMM_WORLD
    commSize = MPI.Comm_size(comm)
    rank = MPI.Comm_rank(comm)
else
    commSize = 1
    rank = 0
end

# Load the parameters in the first (or only) thread
if rank == 0
    t1 = now()
    # Try to load the specified parameters
    params = try
            JSON.parsefile(ARGS[1],dicttype=DataStructures.OrderedDict)
        catch err
            @warn err
            false
        end
    if params != false
        # Check directory to store the data if
        # specified, create directory if it doesn't exist
        dir = try pop!(params, "dir") catch err "" end
        dir != "" && !isdir(dir) && mkpath(dir)

        # Turn sequence lengths into array, since
        # it's under special instructions not to automatically
        # do that when passed through initparams (so it
        # won't be included in loopKey/loopVal)
        params["M"] = stepstoarray(params["M"])
        loopKey, loopVal, noloopParams = initparams(params)
        vals = collect(Iterators.product(loopVal...))

        # Divide out the work into the threads equally,
        # any remainder being left to the last thread
        states = Array{Array}(undef,commSize)
        each = length(vals) ÷ commSize
        for r ∈ 1:commSize
            eachTo = r*each
            r == commSize && (eachTo += length(vals)%commSize)
            states[r] = ((r-1)*each+1):eachTo
        end
    end
else
    # Initialise the variables in other threads
    # to assign the values later on
    dir = nothing
    params = nothing
    loopKey, loopVal, noloopParams = nothing, nothing, nothing
    vals = nothing
    states = nothing
end

hasMPI && (params = MPI.bcast(params,0,comm))
params == false && exit()

if hasMPI
    # Copy the variables over from the first thread
    dir = MPI.bcast(dir,0,comm)
    loopKey = MPI.bcast(loopKey,0,comm)
    loopVal = MPI.bcast(loopVal,0,comm)
    noloopParams = MPI.bcast(noloopParams,0,comm)
    vals = MPI.bcast(vals,0,comm)
    states = MPI.bcast(states,0,comm)
end

# If the run terminated halfway, resume it (…skipstates.json
# stores the runs that have been completed)
stateFilename = paramfilename(params, dir, ".skipStates.json")
skipStates = try JSON.parsefile(stateFilename) catch; [] end

for state ∈ states[rank+1]
    global skipStates
    state ∈ skipStates && continue

    # Load the parameters for the current step
    val = vals[state]
    params = copy(noloopParams)
    alert = ""
    for i ∈ 1:length(loopKey)
        push!(params, loopKey[i]=>val[i])
        alert != "" && (alert *= ", ")
        alert *= "$(loopKey[i]) = $(val[i])"
    end

    # Do an RB run
    println("Rank $rank started $alert on $(now())")
    @time begin
        F,τ,p̄,μr̄,σr̄ = tqdRB(
            params["M"],
            params["K"];
            [Symbol(k) => v for (k,v) ∈ params]...
        )
        print("Rank $rank took")
    end

    # Find the averages and do the fit
    μF = vec(mean(F,dims=2))
    σF = vec(std(F,dims=2))/√(params["K"]+1)
    p,σp = fitRB(params["M"],μF,σF)

    # Save the data
    open(paramfilename(params, dir, ".json"),"w") do f
        JSON.print(f, DataStructures.OrderedDict([
            "iF" => 1 .- F,
            "iμF" => 1 .- μF,
            "σF" => σF,
            "τ" => τ,
            "p" => p,
            "σp" => σp,
            "p̄" => p̄,
            "μr̄" => μr̄,
            "σr̄" => σr̄
        ]))
    end

    # Update completed files
    skipStates = try JSON.parsefile(stateFilename) catch; skipStates end
    push!(skipStates,state)
    open(stateFilename,"w") do f
        JSON.print(f, skipStates)
    end

    println("Parameters $(length(skipStates)) of $(length(vals)) finished")
end

if rank == commSize-1
    if hasMPI
        for r ∈ 0:commSize-2
            println(MPI.recv(r,commSize-1,comm)[1])
        end
    end
    try
        rm(stateFilename)
    catch err
        @warn err
    end
    println("RB experiment $(ARGS[1]) is finished, and it took $(now() - t1).")
else
    MPI.send("Rank $rank is finished",commSize-1,commSize-1,comm)
end
