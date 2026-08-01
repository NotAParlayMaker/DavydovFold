function open_database(path::AbstractString="davydovfoldon.sqlite")
    import SQLite
    db=SQLite.DB(path)
    for sql in ("CREATE TABLE IF NOT EXISTS proteins (id INTEGER PRIMARY KEY, sequence TEXT UNIQUE)",
                "CREATE TABLE IF NOT EXISTS simulations (id INTEGER PRIMARY KEY, protein_id INTEGER, trajectory_path TEXT, created_at TEXT)",
                "CREATE TABLE IF NOT EXISTS spectra (id INTEGER PRIMARY KEY, simulation_id INTEGER, data_path TEXT)",
                "CREATE TABLE IF NOT EXISTS experiments (id INTEGER PRIMARY KEY, name TEXT, data_path TEXT)")
        SQLite.execute(db,sql)
    end
    db
end
function store_simulation!(db,sequence,path)
    import SQLite
    SQLite.execute(db,"INSERT OR IGNORE INTO proteins(sequence) VALUES (?)",(sequence,))
    id=first(SQLite.Query(db,"SELECT id FROM proteins WHERE sequence = ?",(sequence,))).id
    SQLite.execute(db,"INSERT INTO simulations(protein_id,trajectory_path,created_at) VALUES (?,?,datetime('now'))",(id,path))
    SQLite.last_insert_rowid(db)
end

"""Start the lightweight JSON REST service. `POST /predict` accepts `sequence`."""
function serve(;host="127.0.0.1",port=8080,dbpath="davydovfoldon.sqlite")
    import HTTP, JSON3
    db=open_database(dbpath)
    HTTP.serve(host,port) do req
        if req.method=="POST" && req.target=="/predict"
            body=JSON3.read(req.body); seq=String(body.sequence); traj,_,_=run_soliton(seq;tmax=1.0)
            path=tempname()*".h5"; save_trajectory(path,traj); id=store_simulation!(db,seq,path)
            HTTP.Response(200,["Content-Type"=>"application/json"],JSON3.write((trajectory_id=id,path=path)))
        elseif req.method=="GET" && startswith(req.target,"/spectrum/")
            HTTP.Response(501,"spectrum persistence is not configured")
        else
            HTTP.Response(404,"not found")
        end
    end
end
