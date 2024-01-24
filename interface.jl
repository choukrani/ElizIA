"""
    Lit un maillage exporté par GMSH.
"""
function readmesh(nomfile::String)
    fname = nameof(var"#self#")
    println("[$fname] Reading mesh $nomfile... ")

    # open the file and prepare an intput/output reader over it
    io = open(nomfile,"r")
    ### Lecture des noeuds
    while readline(io) != "\$Nodes"
    end
    Nbpt = parse(Int,readline(io))
    Coorneu     = zeros(Float64,Nbpt,2)
    Refneu      = zeros(Int,Nbpt)
    RefneuBis   = zeros(Int,Nbpt)
    for i=1:Nbpt
        # le . pour une application à chaque élément
        Coorneu[i,:] = parse.(Float64,split(readline(io)))[2:3] 
    end
    ### Lecture des éléments
    while readline(io) != "\$Elements"
    end
    Nbtri = parse(Int,readline(io))
    line = parse.(Int,split(readline(io)))
    while line[2] != 15
        line = parse.(Int,split(readline(io)))
        Nbtri -= 1
    end
    # noeuds des coins 
    while line[2] != 1 # on passe tous les 15
        # println("line : ", line)
        Refneu[line[end]] = line[end-2]
        line = parse.(Int,split(readline(io)))
        Nbtri -= 1
    end
    # noeuds du bord
    # on a max Nbpt arêtes *sur le bord*. Sera tronqué après
    Numaretes   = zeros(Int,Nbpt,2)
    Refaretes   = zeros(Int,Nbpt)
    Nbaretes    = 0        
    while line[2] != 2 # on passe tous les 1
        Nbaretes += 1
        RefneuBis[line[6:7]] .= line[4]
        Numaretes[Nbaretes,:] = line[6:7]
        Refaretes[Nbaretes] = line[4]
        line = parse.(Int,split(readline(io)))
        Nbtri -= 1
    end
    Numaretes = Numaretes[1:Nbaretes,:]
    Refaretes = Refaretes[1:Nbaretes]
    Refneu[Refneu .== 0] .= RefneuBis[Refneu .== 0]
    # triangles 
    Numtri = zeros(Int,Nbtri,3)
    Reftri = zeros(Int,Nbtri)
    for i=1:Nbtri-1
        Numtri[i,:] = line[end-2:end]
        Reftri[i] = line[end-3]
        line = parse.(Int,split(readline(io)))
    end
    # je sais, c'est moche, mais ça évite un if
    Numtri[Nbtri,:] = line[end-2:end]
    Reftri[Nbtri] = line[end-3]
    close(io)

    println("[$fname] ... done.")
    return Mesh(Nbpt,Nbtri,Coorneu,Refneu,Numtri,Reftri,Nbaretes,Numaretes,Refaretes)
end

"""
    Plot the mesh Mesh. 
    Call either with the name of the mesh, or 
    with the mesh itself. 
"""
function plotmesh(nom_maillage::String, titre="")
    mesh = readmesh(nom_maillage)
    plotmesh(mesh, titre)
end

function plotmesh(mesh::Mesh, titre="")
    p = plot(title=titre,legend=false)
    for itri in ProgressBar(1:mesh.Nbtri)
        for iar=1:3
            index = 3*(itri-1)+iar
            theline = [mesh.Coorneu[mesh.Numtri[itri,iar],:],mesh.Coorneu[mesh.Numtri[itri,iar%3+1],:]]
            plot!(p, [theline[1][1],theline[2][1]], [theline[1][2],theline[2][2]])
        end
    end
    println("MAX : ", max_value)
    display(p) 
end

"""
    Read a file with eigen elements. (TP3)
"""
function readelem(elemname::String)
    fname = nameof(var"#self#")
    println("[$fname] Reading " * elemname * "...")
    elem = readdlm(elemname)
    println("[$fname] ... done.")
    return elem 
end


"""

Simple loglog display of series of errors.
- listeh : vector of steps 
- errors : series to plot. One colums = one series
- names : vector of tags

"""
function plot_errors(listeh::Vector{Float64}, errors::VecOrMat{Float64}, names::Vector{String}, title="")
    fname = nameof(var"#self#")
    println("[$fname] Creating error plot... ")

    # padd(a,b) = [(b+a)/2 - 1.2 * (b-a)/2, (b+a)/2 + 1.2 * (b - a)/2]
    p = plot(xaxis=:log, yaxis=:log, xlabel="step", ylabel="error", legend=:topleft, title=title) #, aspect_ratio=:equal)
    abserr = abs.(errors)
    col = ["#B8255F", "#FF9933", "#AFB83B", "#6ACCBC", "#4073FF", "#AF38EB", "#808080", "#CCAC93"]
    mk = [:circle, :diamond, :xcross, :star5, :rect, :ltriangle, :pentagon]
    abserr[abserr .<= 1e-14] .= 1e-14
    for i=1:length(errors[1,:])
        if (maximum(abserr[:,i]) < 1e-10)
            println("[$fname] Skipping series $i, identically null")
        else
            plot!(p, listeh, abserr[:,i], shape=mk[i], ls=:solid, lc=col[i], mc=col[i], lw=1.0, label=names[i])
            if length(listeh) > 2 # handcrafted affine approximation
                xii = log10.(listeh)
                yii = log10.(abserr[:,i])
                X = sum(xii)/length(xii); Y = sum(yii)/length(yii)
                XX = sum([x * (x - X) for x in xii]); YY = sum([x * (y - Y) for (x,y) in zip(xii, yii)])
                a = YY / XX; b = Y - a * X
                poly_approx = 10 .^(a .* xii .+ b)
                plot!(p, listeh, poly_approx, shape=mk[i], ls=:dash,  lc=col[i], mc=col[i], lw=1.0, label="a = $(sprintf1("%2.5f",a))")
            end
        end
    end

    println("[$fname] ... done.")
    return p
end

"""
    Affichage d'une matrice avec noms de colonne / noms de lignes et éventuellement titre.
"""
function prettyprint(Data::VecOrMat, rownames::Vector{String}, colnames::Vector{String}, title="")

    if typeof(Data[1,1]) <: Int
        fmt = "d"
    else
        fmt = "f"
    end

    width = maximum([15, maximum(length.(rownames)), maximum(length.(colnames))])
    res = sprintf1("%$(width)s","") * prod([sprintf1("%$(width)s",name) for name in colnames]) * "\n"
    lres = length(res)-1
    linesep = "="^lres * "\n"
    if title==""
        res = linesep * res 
    else
        padtitle = floor(Int,(lres-length(title))/2)
        res = linesep * " "^padtitle * title * " "^(lres-padtitle-length(title)) * "\n" * res 
    end
    for line = 1:length(Data[:,1])
        res *= sprintf1("%$(width)s", rownames[line]) * prod([sprintf1("%$(width)$fmt", x) for x in Data[line,:]]) * "\n"
    end
    res *= linesep
    print(res)
end
