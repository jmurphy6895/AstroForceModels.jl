using Documenter
using AstroForceModels

makedocs(
    modules = [AstroForceModels],
    format = Documenter.HTML(
        prettyurls = !("local" in ARGS),
        canonical = "https//jmurphy6895.github.io/AstroForceModels.jl/stable/",
    ),
    sitename = "AstroForceModels.jl",
    authors = "Jordan Murphy",
    pages = [
        "Home" => "index.md",
        "Usage" => "man/usage.md",
        "API" => "man/api.md",
        "Library" => "lib/library.md",
    ],
)

deploydocs(
    repo = "github.com/jmurphy6895/AstroForceModels.jl.git",
    target = "build",
)
