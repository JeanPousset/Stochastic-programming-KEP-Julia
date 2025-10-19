"""
    solve_cycle_infini

Résout le problème sans contrainte sur la taille des cycles pour une instance du problème de dons en dominos.

# Arguments
* `G::SimpleDiGraph` : graphe de compatibilité de l'instance
* `edge_weight::Matrix{Float64}` : matrice des poids des arcs
"""
function solve_cycle_infini(G::SimpleDiGraph, edge_weight::Matrix{Float64})
    model = Model(HiGHS.Optimizer);

    # Définition des variables : on ne peut pas faire n'importe quoi au moment de la définition des ensembles dans lesquels les indices varient, ci-dessous la façon la plus compacte, et en commentaire une possibilité qui passe par la définition de l'ensemble des arcs
    @variable(model, x[i in vertices(G), j in outneighbors(G,i)], Bin);

    # A = []
    # for e in edges(G)
    #   push!(A, (e.src,e.dst))
    # end
    # dans cette option, les variables ont un seul indice qui correspond à une paire (et non deux indices qui correspondent à un sommet), donc il faudra mettre les parenthèses partout
    # @variable(model, x[(i,j) in A]);

    # Définition de la fonction objectif
    @objective(model, Max, sum(edge_weight[e.src,e.dst]*x[e.src,e.dst] for e in edges(G)));

    # Définition des contraintes
    ## Contraintes de flux
    @constraint(model, [i in vertices(G)], sum(x[i,j] for j in outneighbors(G,i)) == sum(x[j,i] for j in inneighbors(G,i)));
    ## Pas plus d'un transfert par paire patient/donneur
    @constraint(model, [i in vertices(G)], sum(x[j,i] for j in inneighbors(G,i)) <= 1);

    # Résolution
    set_silent(model); # je désactive le slog de HiGHS pour que les sorties soient courtes, mais il peut être intéressant de les afficher pour essayer de comprendre ce que le solveur dit, le cours doit vous permettre de comprendre au moins une partie des messages
    timer = @timed optimize!(model);

    # Affichage des résultats
    println("\n-------------------------------")
    println("Résolution du modèle sans contrainte de longueur de cycle\n")
    println("Statut de la résolution : ", termination_status(model));
    println("Durée de la résolution : $(timer.time)")
    if (termination_status(model) == MOI.OPTIMAL)
        println("Nombre de transferts réalisés : ", objective_value(model));
        println("Liste des transferts réalisés :")
        for e in edges(G)
            if value(x[e.src,e.dst]) > 0
                print("$(e.src) -> $(e.dst) ; ");
            end
        end
    end
    println("\n-------------------------------\n")
end


"""
    solve_cycle_2

Résout le problème avec une contrainte sur la taille des cycles (au plus 2) pour une instance du problème de dons en dominos.

# Arguments
* `G::SimpleDiGraph` : graphe de compatibilité de l'instance
* `edge_weight::Matrix{Float64}` : matrice des poids des arcs
"""
function solve_cycle_2(G::SimpleDiGraph, edge_weight::Matrix{Float64})
    model = Model(HiGHS.Optimizer);
    

    # Il faut d'abord créer l'ensemble des arêtes du graphe non-orienté G2 correspondant à tous les échanges pouvant être faits entre deux paires patient-donneur
    mat_adj = adjacency_matrix(G); # matrice d'adjacence du graphe orienté, contient un nombre non nul si l'arc existe
    E2 = []; # liste des arêtes du graphe non-orienté
    voisins =  [];
    for i in vertices(G)
        push!(voisins, []);
    end
    for e in edges(G)
        # on ajoute que les arêtes (i,j) avec i<j pour éviter les répétitions
        if e.src < e.dst && mat_adj[e.dst,e.src] != 0
            push!(E2, (e.src,e.dst));
            push!(voisins[e.src], e.dst);
            push!(voisins[e.dst], e.src);
        end
    end

    # Définition des variables
    @variable(model, x[i in vertices(G), j in voisins[i] ; j > i], Bin);

    # Définition de la fonction objectif
    @objective(model, Max, sum((edge_weight[e[1],e[2]]+edge_weight[e[2],e[1]]) *x[e[1],e[2]] for e in E2));

    # Définition des contraintes
    ## Pas plus d'un transfert par paire patient/donneur
    @constraint(model, [i in vertices(G)], sum(x[i,j] for j in voisins[i] if j > i) + sum(x[j,i] for j in voisins[i] if j < i) <= 1);

    # Résolution
    set_silent(model)
    timer = @timed optimize!(model)

    # Affichage des résultats
    println("\n-------------------------------")
    println("Résolution du modèle avec cycles de taille 2\n")
    println("Statut de la résolution : ", termination_status(model));
    println("Durée de la résolution : $(timer.time)")
    if (termination_status(model) == MOI.OPTIMAL)
        println("Nombre de transferts réalisés : ", objective_value(model));
        println("Liste des transferts réalisés :")
        for e in E2
            if value(x[e[1],e[2]]) > 0
                print("$(e[1]) <-> $(e[2]) ; ");
            end
        end
    end
    println("\n-------------------------------\n")
end


"""
    solve_cycle_K

Résout le problème avec une contrainte sur la taille des cycles (au plus K) pour une instance du problème de dons en dominos.

# Arguments
* `G::SimpleDiGraph` : graphe de compatibilité de l'instance
* `edge_weight::Matrix{Float64}` : matrice des poids des arcs
* `K::Int` : taille maximale des cycles
"""
function solve_cycle_K(G::SimpleDiGraph, edge_weight::Matrix{Float64}, K::Int)
    model = Model(HiGHS.Optimizer);

    # Définition des variables
    H = collect(1:nv(G)); # ensemble d'hôpitaux
    @variable(model, transfert[i in vertices(G), j in outneighbors(G,i), h in H ; h <= i && h <= j], Bin);
    @variable(model, utilise_hopital[h in H], Bin);

    # Définition de la fonction objectif
    ## On ajoute un terme qui n'influe pas sur le nombre optimal de transferts  tant que les utilités sont entières ; cela permet de rompre une symétrie de plus en favorisant une solution optimale qui emploie un maximum d'hôpitaux. Cela évite notamment les solutions avec plusieurs cycles sur le même hôpital. Notez que si l'on considère plutôt le terme opposé, ça  favorise les solutions avec un petit nombre d'hôpitaux et on trouve une autre solution à 4 transferts pour la plus petite instance.
    @objective(model, Max, sum(edge_weight[e.src,e.dst]*transfert[e.src,e.dst, h] for e in edges(G) for h in H if h <= e.src && h <= e.dst) + 1/nv(G) * sum(utilise_hopital[h] for h in H));

    # Définition des contraintes
    ## Contraintes de flux dans chaque hôpital
    @constraint(model, [i in vertices(G), h in H ; h <= i], sum(transfert[i,j,h] for j in outneighbors(G,i) if j >= h) == sum(transfert[j,i,h] for j in inneighbors(G,i) if j >= h));

    ## Pas plus d'un transfert par paire patient/donneur
    @constraint(model, [i in vertices(G)], sum(transfert[j,i,h] for j in inneighbors(G,i) for h in H if h <= i && h <= j) <= 1);

    # Pas plus de K transferts par hôpital
    @constraint(model, [h in H], sum(transfert[e.src,e.dst,h] for e in edges(G) if e.src >= h && e.dst >= h) <= K);

    # On impose que la paire h soit impliquée dans un transfert de l'hôpital h si l'hôpital h est utilisé
    @constraint(model, [h in H], utilise_hopital[h] == sum(transfert[h,j,h] for j in outneighbors(G,h) if j >= h));

    # Résolution
    # Résolution
    set_silent(model)
    set_time_limit_sec(model, 120.0)
    timer = @timed optimize!(model)

    # Affichage des résultats
    println("\n-------------------------------")
    println("Résolution du modèle avec cycle de taille au plus $K\n")
    println("Statut de la résolution : ", termination_status(model));
    println("Durée de la résolution : $(timer.time)")
    if (has_values(model))
        println("Nombre de transferts réalisés : ", floor(Int, objective_value(model)));
        println("Liste des transferts réalisés :")
        for e in edges(G)
            val = sum(value(transfert[e.src,e.dst,h]) for h in H if h <= e.src && h <= e.dst)
            if val > 1.0e-4
                print("$(e.src) -> $(e.dst) ; ");
            end
        end
        println("\n-------------------------------\n")
    end
end
