using Graphs

"""
    read_kep_file

Contruct G from a `.wmd` file from PrefLib.

# Parameters
* `wmd_file::String` : path of the `.wmd` file.
"""
function read_wmd_file(wmd_file::String)
    isfile(wmd_file) || throw(ArgumentError("$(wmd_file): file not found."))

    # Extract the graph structure from the .wmd file using a MetaGraph
    wmd_io = open(wmd_file, "r")

    # skip the first nine informative lines
    for i in 1:9
        readline(wmd_io)
    end

    # get the number of vertices and edges from the following lines
    splitted_line = split(readline(wmd_io), ' ')
    nb_vertices = parse(Int, splitted_line[4])
    splitted_line = split(readline(wmd_io), ' ')
    nb_edges = parse(Int, splitted_line[4])
    G = SimpleDiGraph(nb_vertices, 0)
    edge_weight = zeros(nb_vertices, nb_vertices)


    # skip next nb_vertices lines, which are only informative
    for i in 1:nb_vertices
        readline(wmd_io)
    end

    # read the set of edges
    while !eof(wmd_io)
        splitted_line = split(readline(wmd_io), ',')
        # /!\ Pairs are numbered from 0 in the second part of the file
        source = parse(Int, splitted_line[1])
        destination = parse(Int, splitted_line[2])
        weight = parse(Float64, splitted_line[3])

        # do not add an edge that has a zero weight
        if weight > 0
            add_edge!(G, source, destination)
            edge_weight[source, destination] = weight
        end
    end

    return G, edge_weight
end


"""
    gs_compatible

Renvoie true si un groupe sanguin de donneur est compatible avec celui d'un malade 

# Paramètres
* `gsd::String` : groupe sanguin d'un donneur
* `gsm::String` : groupe sanguin d'un malade
"""
function gs_compatible(gsd::String, gsm::String)
    if (gsd == "O") || (gsm == "AB")
        # O peut donner à {O,A,B,AB}, AB peut recevoir de {O,A,B,AB}
        return true
    elseif (gsd == "A") && (gsm == "A")
        #  A donne à {A,AB}
        return true
    elseif (gsd == "B") && (gsm == "B")
        # B donne à {B,AB}
        return true
    else
        # O ne peut pas recevoir de {A,B,AB}, A ne peut pas recevoir {B,AB}, B ne peut pas recevoir de {A,AB}
        return false
    end
end

"""
    read_dat_file

Contruct Gtilde from a `.dat` file from PrefLib.

# Parameters
* `dat_file::String` : path of the `.dat` file.
"""
function read_dat_file(dat_file::String)
    isfile(dat_file) || throw(ArgumentError(".dat file not found."))

    # Extraction des données individuelles depuis le fichier .dat
    file = readdlm(dat_file, '\n')
    nb_vertices = length(file) - 1
    gsm = Dict{Int64,String}() # groupe sanguin du malade
    gsd = Dict{Int64,String}() # groupe sanguin du donneur
    pra = Dict{Int64,Float64}() # PRA du malade
    for line in file[2:end]
        splitted_line = split(line, ',')
        pair = parse(Int, splitted_line[1])
        gsm[pair] = String(splitted_line[2])
        gsd[pair] = String(splitted_line[3])
        pra[pair] = parse(Float64, splitted_line[5])
    end

    # Construction du graphe de compatibilité a priori ̃G ; on fixe le poids des arcs à 1 par défaut (à modifier dans une autre fonction si souhaité) 
    Gtilde = SimpleDiGraph(nb_vertices, 0)
    edge_weight = zeros(nb_vertices, nb_vertices)
    for u in 1:nb_vertices
        for v in u+1:nb_vertices
            if gs_compatible(gsd[u], gsm[v])
                add_edge!(Gtilde, u, v)
                edge_weight[u, v] = 1
            end
            if gs_compatible(gsd[v], gsm[u])
                add_edge!(Gtilde, v, u)
                edge_weight[v, u] = 1
            end
        end
    end

    return Gtilde, edge_weight, gsm, gsd, pra
end