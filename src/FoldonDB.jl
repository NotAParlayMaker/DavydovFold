"""Open (or create) the DavydovFoldon SQLite database and its schema.

The returned handle is owned by the caller and should be closed with
`SQLite.close(db)` when it is no longer needed.
"""
function open_database(path::AbstractString="davydovfoldon.sqlite")
    import SQLite
    db = SQLite.DB(path)
    SQLite.execute(db, "PRAGMA foreign_keys = ON")
    statements = (
        """CREATE TABLE IF NOT EXISTS proteins (
               id INTEGER PRIMARY KEY,
               sequence TEXT NOT NULL UNIQUE CHECK(length(sequence) > 0)
           )""",
        """CREATE TABLE IF NOT EXISTS simulations (
               id INTEGER PRIMARY KEY,
               protein_id INTEGER NOT NULL REFERENCES proteins(id),
               trajectory_path TEXT NOT NULL,
               metadata_json TEXT NOT NULL DEFAULT '{}',
               created_at TEXT NOT NULL DEFAULT CURRENT_TIMESTAMP
           )""",
        """CREATE TABLE IF NOT EXISTS spectra (
               id INTEGER PRIMARY KEY,
               simulation_id INTEGER NOT NULL REFERENCES simulations(id),
               data_path TEXT NOT NULL,
               created_at TEXT NOT NULL DEFAULT CURRENT_TIMESTAMP
           )""",
        """CREATE TABLE IF NOT EXISTS experiments (
               id INTEGER PRIMARY KEY,
               name TEXT NOT NULL,
               data_path TEXT NOT NULL,
               metadata_json TEXT NOT NULL DEFAULT '{}'
           )""",
        "CREATE INDEX IF NOT EXISTS simulations_protein_idx ON simulations(protein_id)",
    )
    foreach(sql -> SQLite.execute(db, sql), statements)
    return db
end

"""Persist a trajectory reference and metadata, returning its database id."""
function store_simulation!(db, sequence::AbstractString, path::AbstractString;
                           metadata=NamedTuple())
    import JSON3, SQLite
    seq = uppercase(String(sequence))
    isempty(seq) && throw(ArgumentError("sequence must not be empty"))
    SQLite.execute(db, "INSERT OR IGNORE INTO proteins(sequence) VALUES (?)", (seq,))
    protein = first(SQLite.Query(db, "SELECT id FROM proteins WHERE sequence = ?", (seq,)))
    SQLite.execute(
        db,
        "INSERT INTO simulations(protein_id, trajectory_path, metadata_json) VALUES (?, ?, ?)",
        (protein.id, String(path), JSON3.write(metadata)),
    )
    return Int(SQLite.last_insert_rowid(db))
end

"""Start the local JSON service; `POST /predict` accepts a `sequence` field."""
function serve(; host::AbstractString="127.0.0.1", port::Integer=8080,
               dbpath::AbstractString="davydovfoldon.sqlite")
    import HTTP, JSON3
    1 <= port <= 65535 || throw(ArgumentError("port must be in 1:65535"))
    db = open_database(dbpath)
    HTTP.serve(String(host), Int(port)) do request
        if request.method == "POST" && request.target == "/predict"
            body = JSON3.read(request.body)
            hasproperty(body, :sequence) || return HTTP.Response(400, "missing sequence")
            sequence = String(body.sequence)
            trajectory, _, _ = run_soliton(sequence; tmax=1.0)
            path = tempname() * ".h5"
            save_trajectory(path, trajectory)
            id = store_simulation!(db, sequence, path; metadata=(temperature_K=310.0,))
            return HTTP.Response(200, ["Content-Type" => "application/json"],
                                 JSON3.write((trajectory_id=id, path=path)))
        end
        return HTTP.Response(404, "not found")
    end
end
